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
#include "ResidualPFM.h"
#include "Solver.h"

#include <gmsh.h>

#include <Eigen/Dense>

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <string_view>
#include <vector>

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

// ---------------------------------------------------------------------------
// Terminal logging: a streambuf that writes every character to TWO underlying
// buffers, so std::cout / std::cerr reach the console AND a log file at once.
// ---------------------------------------------------------------------------
class TeeBuf : public std::streambuf {
public:
    TeeBuf(std::streambuf* primary, std::streambuf* secondary)
        : primary_(primary), secondary_(secondary) {}

protected:
    int overflow(int ch) override {
        if (ch == traits_type::eof()) return traits_type::not_eof(ch);
        const int r1 = primary_   ? primary_->sputc(static_cast<char>(ch)) : ch;
        const int r2 = secondary_ ? secondary_->sputc(static_cast<char>(ch)) : ch;
        return (r1 == traits_type::eof() || r2 == traits_type::eof())
                   ? traits_type::eof() : ch;
    }
    int sync() override {
        const int r1 = primary_   ? primary_->pubsync() : 0;
        const int r2 = secondary_ ? secondary_->pubsync() : 0;
        return (r1 == 0 && r2 == 0) ? 0 : -1;
    }

private:
    std::streambuf* primary_;
    std::streambuf* secondary_;
};

// RAII: on construction redirect std::cout / std::cerr through tee buffers that
// also write to `path`; on destruction restore the original buffers FIRST (so
// the streams are valid again before the log file is closed).
class TerminalLogger {
public:
    TerminalLogger(const std::string& path, bool enabled) : enabled_(enabled) {
        if (!enabled_) return;
        file_.open(path, std::ios::out | std::ios::trunc);
        if (!file_) {                       // couldn't open -> run without a log
            enabled_ = false;
            std::cerr << "Warning: could not open log file \"" << path
                      << "\"; continuing without a terminal log.\n";
            return;
        }
        cout_old_ = std::cout.rdbuf();
        cerr_old_ = std::cerr.rdbuf();
        tee_out_  = std::make_unique<TeeBuf>(cout_old_, file_.rdbuf());
        tee_err_  = std::make_unique<TeeBuf>(cerr_old_, file_.rdbuf());
        std::cout.rdbuf(tee_out_.get());
        std::cerr.rdbuf(tee_err_.get());
    }

    ~TerminalLogger() {
        if (!enabled_) return;
        std::cout.flush();
        std::cerr.flush();
        std::cout.rdbuf(cout_old_);
        std::cerr.rdbuf(cerr_old_);
    }

    bool active() const { return enabled_; }

private:
    bool                     enabled_;
    std::ofstream            file_;
    std::streambuf*          cout_old_ = nullptr;
    std::streambuf*          cerr_old_ = nullptr;
    std::unique_ptr<TeeBuf>  tee_out_;
    std::unique_ptr<TeeBuf>  tee_err_;
};

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

    // Mirror all terminal output to "<base_name>_run.txt" (config: write_log).
    // Installed before the first print so the whole run is captured. Restored
    // automatically at scope exit. NOTE: Gmsh's own library messages print at
    // the C level and bypass std::cout, so they appear on the console only.
    const std::string log_path = cfg.mesh.base_name + "_run.txt";
    TerminalLogger logger(log_path, cfg.write_log);
    if (logger.active())
        std::cout << "[log] mirroring terminal output to " << log_path << "\n";

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

    // ---- Force-displacement logging (diagnostic) ---------------------------
    // The reaction on the loaded boundary equals the internal-force residual
    // at the constrained DOFs -- exactly the quantity zeroReactions() discards
    // from the convergence norm. We sum the vertical residual over the "loaded"
    // nodes (Dirichlet entries with a nonzero prescribed uy) and pair it with
    // the applied displacement load_factor * uy_full, written per converged
    // step. If the force rises, peaks, then drops while displacement keeps
    // increasing (or turns back), that is the snap-back signature.
    std::vector<int> loaded_ydof;     // global y-dof of each loaded node
    double           uy_full = 0.0;   // prescribed uy at full load
    for (int i = 0; i < d.nvfix; ++i) {
        if (d.iffix[i][1] != 0 && d.fixed[i][1] != 0.0) {
            loaded_ydof.push_back(2 * d.nofix[i] + 1);
            uy_full = d.fixed[i][1];
        }
    }
    std::cout << "[force-disp] logging reaction over " << loaded_ydof.size()
              << " loaded node(s); uy_full = " << uy_full << "\n";

    std::ofstream fd_csv(cfg.mesh.base_name + "_force_disp.csv");
    fd_csv << "step,load_factor,applied_uy,reaction_Fy,max_phi\n";

    // Compute the total vertical reaction at (u, phi) and append one CSV row.
    const auto logForceDisp = [&](int step, double load_factor) {
        const Eigen::VectorXd R_full =
            pfm::assembleGlobalResidual(d, u, phi, history);
        double Fy = 0.0;
        for (int g : loaded_ydof) Fy += R_full(g);
        fd_csv << step << ',' << load_factor << ','
               << (load_factor * uy_full) << ',' << Fy << ','
               << phi.maxCoeff() << '\n';
        fd_csv.flush();
    };

    // ---- Quasi-static load stepping ----------------------------------------
    // Walk load_factor from 0 -> 1 in N_steps equal increments. Each call to
    // solveStep finds Newton equilibrium at one load level; (u, phi) carry
    // over as warm-starts. After a converged step the history field is raised
    // so the phase field cannot heal at the next increment. On the first
    // non-converged step we halt.
    const int N_steps = cfg.N_steps;

    pfm::io::writeVTK(stepFilename(0), d, u, phi);
    logForceDisp(0, 0.0);
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
    const int max_subdivs = 6;    // smallest allowed increment = nominal / 2^6

    // Two-stage displacement control (Ambati-style) vs. the uniform schedule.
    // In two_stage mode the applied displacement u = lf * u_ref is advanced by
    // du_coarse until it reaches u_switch, then by du_fine; the adaptive
    // halving below still applies on non-convergence.
    const bool   two_stage = (cfg.step.mode == appcfg::StepMode::TwoStage);
    const double u_ref     = std::abs(uy_full);   // |full-load applied uy|

    if (two_stage && u_ref == 0.0)
        throw std::runtime_error(
            "run.step_mode = two_stage needs a nonzero prescribed Dirichlet uy "
            "(no loaded edge was found)");

    // Stage-dependent NOMINAL load-factor increment for the current u = done*u_ref.
    const auto nominal_inc = [&](double done) -> double {
        if (!two_stage) return 1.0 / N_steps;
        const double u_now = done * u_ref;
        const double du = (u_now < cfg.step.u_switch - 1e-15)
                              ? cfg.step.du_coarse : cfg.step.du_fine;
        return du / u_ref;
    };

    if (two_stage)
        std::cout << "[step] two-stage displacement schedule: du_coarse = "
                  << cfg.step.du_coarse << " up to u = " << cfg.step.u_switch
                  << ", then du_fine = " << cfg.step.du_fine
                  << "  (u_ref = " << u_ref << ", N_steps ignored)\n";
    else
        std::cout << "[step] uniform schedule: " << N_steps
                  << " equal load_factor increments\n";

    double lf_done        = 0.0;
    double lf_inc         = nominal_inc(0.0);
    int    step_index     = 0;       // counts ACCEPTED steps (drives VTK naming)
    int    subdiv_streak  = 0;       // subdivisions used for current failed step
    int    success_streak = 0;       // consecutive successes (drives re-growth)

    const int vtk_every     = cfg.vtk_every;  // write VTK every Nth accepted step
    int       last_vtk_step = 0;              // last step whose VTK was written (0)

    Eigen::VectorXd u_saved   = u;
    Eigen::VectorXd phi_saved = phi;

    while (lf_done < 1.0 - 1e-12) {
        // Current stage's nominal increment (changes when u crosses u_switch).
        const double nominal = nominal_inc(lf_done);
        // Never exceed the stage nominal -- this shrinks lf_inc the moment the
        // schedule steps from the coarse stage into the fine stage.
        if (lf_inc > nominal) lf_inc = nominal;

        double lf_try = std::min(lf_done + lf_inc, 1.0);
        // Two-stage: don't overshoot the switch displacement -- land exactly on
        // it so the coarse/fine boundary is clean.
        if (two_stage) {
            const double u_done = lf_done * u_ref;
            const double u_try  = lf_try  * u_ref;
            if (u_done < cfg.step.u_switch - 1e-15 &&
                u_try  > cfg.step.u_switch + 1e-15)
                lf_try = std::min(cfg.step.u_switch / u_ref, 1.0);
        }

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

            // Force-displacement curve is logged EVERY accepted step (cheap).
            logForceDisp(step_index, lf_done);

            // VTK snapshot only on the stride (step 0 already written before
            // the loop; the final step is guaranteed after the loop).
            if (vtk_every <= 1 || step_index % vtk_every == 0) {
                const std::string fname = stepFilename(step_index);
                pfm::io::writeVTK(fname, d, u, phi);
                last_vtk_step = step_index;
                std::cout << "  [output] wrote " << fname << "\n";
            }

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
            if (success_streak >= 2 && lf_inc < nominal) {
                lf_inc = std::min(lf_inc * 2.0, nominal);
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

            const double lf_inc_min = nominal / (1 << max_subdivs);
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

    // Always persist the final converged state, even if it fell between VTK
    // strides. u_saved / phi_saved hold the last accepted (u, phi).
    if (step_index > 0 && step_index != last_vtk_step) {
        const std::string fname = stepFilename(step_index);
        pfm::io::writeVTK(fname, d, u_saved, phi_saved);
        std::cout << "  [output] wrote final " << fname << "\n";
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
