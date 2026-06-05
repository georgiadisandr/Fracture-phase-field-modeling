// Thin driver around mesh::generate + the PFM solver.
//
// Responsibilities:
//   * load run configuration from a TOML file (path is argv[1])
//   * apply optional CLI overrides on top
//   * own gmsh::initialize / gmsh::finalize
//   * call mesh::generate
//   * write .msh / .vtk
//   * read FemInput from the live gmsh model using the loaded FemSpec
//   * loop load_factor 0 -> 1 in N_steps, calling pfm::solveStep per step
//   * maintain the crack-driving history field (irreversibility)
//   * write a numbered VTK snapshot per step (ParaView reads them as a series)
//   * optionally launch the Gmsh GUI at the end
//
// Usage:
//     phasefield_1 config.toml [--W 80 --steps 30 --no-gui ...]
//
// See Config.h for the full list of recognised CLI overrides and config.toml
// for the TOML schema.

#include "Config.h"
#include "GmshReader.h"
#include "History.h"
#include "MeshGeneration.h"
#include "Output.h"
#include "Solver.h"

#include <gmsh.h>

#include <Eigen/Dense>

#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {

void print_usage(const char* prog)
{
    std::cout <<
        "Usage: " << prog << " CONFIG.toml [options]\n"
        "\n"
        "  CONFIG.toml         path to TOML run configuration (required)\n"
        "\n"
        "Common CLI overrides (any TOML value can be overridden):\n"
        "  --out NAME          output base filename (no extension)\n"
        "  --W val             specimen width\n"
        "  --H val             specimen height (also resets y_crack = H/2\n"
        "                      unless --ycrack is also given)\n"
        "  --ycrack val        crack y-coordinate\n"
        "  --a val             crack length\n"
        "  --hfine val         fine mesh size near the crack\n"
        "  --hfar val          coarse mesh size in the bulk\n"
        "  --steps N           number of load-stepping increments\n"
        "  --gui / --no-gui    force Gmsh GUI on or off after the run\n"
        "  --help, -h          show this message and exit\n";
}

}  // namespace

int main(int argc, char** argv) try
{
    // ---- CLI: TOML path is mandatory ---------------------------------------
    if (argc < 2 ||
        std::string_view(argv[1]) == "--help" ||
        std::string_view(argv[1]) == "-h")
    {
        print_usage(argv[0]);
        return argc < 2 ? 1 : 0;
    }

    // Wall-clock timers (steady_clock so they aren't affected by NTP / DST).
    // t_start  = whole-run baseline
    // t_setup  = stamped after FEM input has been built, just before the
    //            first load step. The difference gives the setup time.
    using clock = std::chrono::steady_clock;
    const auto t_start = clock::now();

    // ---- Load + validate configuration -------------------------------------
    appcfg::AppConfig cfg = appcfg::loadConfigFromToml(argv[1]);
    appcfg::applyCliOverrides(cfg, argc, argv);
    appcfg::validate(cfg);

    std::cout << "[config] loaded \"" << argv[1] << "\"\n"
              << "         base_name = " << cfg.mesh.base_name << "\n"
              << "         geometry  W=" << cfg.mesh.W
              << "  H=" << cfg.mesh.H
              << "  a=" << cfg.mesh.a
              << "  y_crack=" << cfg.mesh.y_crack << "\n"
              << "         mesh size h_fine=" << cfg.mesh.h_fine
              << "  h_far="  << cfg.mesh.h_far  << "\n"
              << "         materials=" << cfg.fem.materials.size()
              << "  bcs=" << cfg.fem.bcs.size()
              << "  N_steps=" << cfg.N_steps << "\n";

    // ---- Mesh generation ----------------------------------------------------
    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::model::add(cfg.mesh.base_name);

    mesh::generate(cfg.mesh);

    const std::string msh = cfg.mesh.base_name + ".msh";
    const std::string vtk = cfg.mesh.base_name + ".vtk";
    gmsh::write(msh);
    gmsh::write(vtk);
    std::cout << "Wrote " << msh << " and " << vtk << "\n";

    // ---- Build FEM problem from the live gmsh model + the loaded FemSpec ---
    const FemInput d = fem::readFemInputFromGmsh(cfg.fem);

    std::cout << "[mesh] npoin = " << d.npoin
              << ",  nelem = "      << d.nelem
              << " (tri3="          << d.count_by_nnode(3)
              << ", quad4="         << d.count_by_nnode(4)
              << ", quad8="         << d.count_by_nnode(8)
              << "),  nvfix = "     << d.nvfix << "\n";

    Eigen::VectorXd u   = Eigen::VectorXd::Zero(2 * d.npoin);
    Eigen::VectorXd phi = Eigen::VectorXd::Zero(d.npoin);

    // Apply any [[initial_phi]] overrides resolved by GmshReader: each
    // entry sets phi(node) = value at step 0. The solver evolves phi from
    // there; no Dirichlet pinning is applied, so phi may relax if the
    // surrounding history field cannot hold it.
    if (!d.initial_phi_nodes.empty()) {
        for (const auto& e : d.initial_phi_nodes)
            phi(e.node) = e.value;
        std::cout << "[init] applied initial phi to "
                  << d.initial_phi_nodes.size()
                  << " node(s) from [[initial_phi]]\n";
    }

    // Crack-driving history field, enforcing irreversibility. Starts at zero
    // and is raised after every converged step (see updateHistory).
    pfm::HistoryField history = pfm::makeHistoryField(d);

    const auto stepFilename = [&](int step) {
        std::ostringstream ss;
        ss << cfg.mesh.base_name << "_step_"
           << std::setw(3) << std::setfill('0') << step << ".vtk";
        return ss.str();
    };

    // ---- Quasi-static load stepping ----------------------------------------
    // Walk load_factor from 0 -> 1 in N_steps equal increments. Each call to
    // solveStep finds Newton equilibrium at one load level; (u, phi) carry
    // over as warm-starts. After a converged step the history field is raised
    // so the phase field cannot heal at the next increment. On the first
    // non-converged step we halt.
    const int N_steps = cfg.N_steps;

    pfm::io::writeVTK(stepFilename(0), d, u, phi);
    std::cout << "[step  0]  load_factor = 0  (initial state written)\n";

    // Setup complete; everything after this is the solve.
    const auto t_solve_start = clock::now();
    {
        using namespace std::chrono;
        const double setup_s = duration<double>(t_solve_start - t_start).count();
        std::cout << "[time] setup (config + mesh + FEM input) = "
                  << std::fixed << std::setprecision(2)
                  << setup_s << " s\n"
                  << std::defaultfloat;
    }

    // Cumulative solve time -- the value reported at the end excludes
    // mesh generation and any post-processing GUI work.
    double solve_seconds = 0.0;

    // ---- Adaptive load-step subdivision ------------------------------------
    // Phase-field fracture is fragile during crack-growth events: a too-large
    // load increment can let subsolveU fail, which feeds spurious psi^+ into
    // the phase field (via H = max(H_stored, psi^+)) and contaminates phi.
    // The standard remedy is to walk load_factor adaptively:
    //   * nominal increment = 1 / N_steps
    //   * on failure: roll u/phi back to the last converged state, halve the
    //                 increment, retry (up to max_subdivs times)
    //   * on success: accept, then -- after a few consecutive successes --
    //                 grow the increment back toward the nominal
    // history is only mutated by main() after a converged step, so it never
    // needs an explicit rollback.
    const double lf_nominal  = 1.0 / N_steps;
    const int    max_subdivs = 6;    // smallest allowed = lf_nominal / 2^6
    const double lf_inc_min  = lf_nominal / (1 << max_subdivs);

    double lf_done        = 0.0;
    double lf_inc         = lf_nominal;
    int    step_index     = 0;       // counts ACCEPTED steps (drives VTK naming)
    int    subdiv_streak  = 0;       // subdivisions used for current failed step
    int    success_streak = 0;       // consecutive successes (drives re-growth)

    Eigen::VectorXd u_saved   = u;
    Eigen::VectorXd phi_saved = phi;

    while (lf_done < 1.0 - 1e-12) {
        const double lf_try = std::min(lf_done + lf_inc, 1.0);
        std::cout << "\n=== trying load_factor = " << lf_try
                  << "  (inc = " << lf_inc
                  << ", step_index = " << (step_index + 1)
                  << ", subdivs = " << subdiv_streak << ") ===\n";

        const auto t_step_start = clock::now();

        const pfm::SolverResult res =
            pfm::solveStep(d, u, phi, history, lf_try, cfg.solver);

        std::cout << "  result: "
                  << (res.converged ? "converged" : "DIVERGED")
                  << " in " << res.iters_used
                  << (cfg.solver.scheme == pfm::SolverScheme::Staggered
                          ? " staggered sweeps," : " Newton iters,")
                  << "  |R_final| = " << res.final_residual
                  << ",  max|u| = "   << u.cwiseAbs().maxCoeff()
                  << ",  max phi = "  << phi.maxCoeff() << "\n";

        if (res.converged) {
            // Accept: commit history, snapshot, advance lf_done, write VTK.
            pfm::updateHistory(d, u, history);
            u_saved   = u;
            phi_saved = phi;
            lf_done   = lf_try;
            ++step_index;
            subdiv_streak = 0;
            ++success_streak;

            const std::string fname = stepFilename(step_index);
            pfm::io::writeVTK(fname, d, u, phi);
            std::cout << "  [output] wrote " << fname << "\n";

            {
                using namespace std::chrono;
                const double step_s =
                    duration<double>(clock::now() - t_step_start).count();
                solve_seconds += step_s;
                std::cout << "  [time] step " << step_index << " = "
                          << std::fixed << std::setprecision(2)
                          << step_s << " s   (cumulative solve = "
                          << solve_seconds << " s)\n"
                          << std::defaultfloat;
            }

            // Grow the increment back toward nominal after the loading is
            // clearly behaving (two consecutive successes is a conservative
            // trigger; never exceed the user's nominal step size).
            if (success_streak >= 2 && lf_inc < lf_nominal) {
                lf_inc = std::min(lf_inc * 2.0, lf_nominal);
                std::cout << "  [adaptive] increment grown to " << lf_inc << "\n";
            }
        }
        else {
            // Reject: roll u/phi back to the last converged state and try a
            // smaller increment. Give up if we've subdivided too many times.
            u   = u_saved;
            phi = phi_saved;
            success_streak = 0;
            ++subdiv_streak;

            if (lf_inc * 0.5 < lf_inc_min || subdiv_streak > max_subdivs) {
                std::cerr << "[adaptive] cannot subdivide further (lf_inc = "
                          << lf_inc << ", subdivs used = " << subdiv_streak
                          << "); halting load-stepping loop at lf_done = "
                          << lf_done << ".\n";
                break;
            }

            lf_inc *= 0.5;
            std::cerr << "  [adaptive] non-converged step rolled back; "
                      << "increment halved to " << lf_inc
                      << " and retrying.\n";
        }
    }

    // Final timing summary -- printed regardless of whether the load
    // stepping completed all N_steps or bailed early on non-convergence.
    {
        using namespace std::chrono;
        const auto   t_end = clock::now();
        const double total_s = duration<double>(t_end - t_start).count();
        std::cout << "\n[time] summary:\n"
                  << std::fixed << std::setprecision(2)
                  << "         solve   = " << solve_seconds << " s\n"
                  << "         total   = " << total_s        << " s\n"
                  << std::defaultfloat;
    }

    if (cfg.show_gui) gmsh::fltk::run();

    gmsh::finalize();
    return 0;
}
catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    try { gmsh::finalize(); } catch (...) {}
    return 1;
}
