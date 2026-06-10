#pragma once
//
// Config.h -- single point of truth for everything the driver needs to run.
//
// AppConfig bundles the existing problem-definition structs:
//   * mesh::Config         (specimen geometry & mesh sizing)
//   * fem::FemSpec         (element type, materials, BCs, group->material map)
//   * pfm::SolverSettings  (Newton tolerances, scheme, verbosity)
// plus a few driver-level knobs (load stepping, GUI flag).
//
// It is populated by two layers, in order:
//   1. loadConfigFromToml(path)   -- parse a TOML file into a fresh AppConfig.
//   2. applyCliOverrides(cfg, argc, argv)
//                                  -- override individual fields from the CLI.
//
// Both layers throw std::runtime_error / std::invalid_argument with a clear
// message on bad input; main() just catches std::exception.
//

#include "MeshGeneration.h"
#include "GmshReader.h"
#include "Solver.h"

#include <string>

namespace appcfg {

// Load-stepping mode.
//   Uniform  : N_steps equal load_factor increments 0 -> 1 (original behaviour).
//   TwoStage : displacement-controlled two-rate schedule (Ambati-style). The
//              applied displacement is advanced by du_coarse until it reaches
//              u_switch, then by du_fine to the end. The schedule is keyed on
//              the magnitude of the largest prescribed Dirichlet uy (the loaded
//              edge), converted to a load_factor increment internally. N_steps
//              is ignored in this mode. The adaptive halving safety net still
//              applies on non-convergence.
enum class StepMode { Uniform, TwoStage };

struct StepControl {
    StepMode mode      = StepMode::Uniform;
    double   du_coarse = 1.0e-5;   // coarse displacement increment (mm)
    double   u_switch  = 5.0e-3;   // switch displacement (mm)
    double   du_fine   = 1.0e-6;   // fine displacement increment (mm)
};

struct AppConfig {
    // Mesh + geometry knobs (W, H, a, y_crack, h_fine, h_far, r_fine, r_far,
    // base_name). Defaults come from mesh::Config.
    mesh::Config mesh;

    // Element type, materials library, physical-group-to-material map, and
    // Dirichlet / Neumann BCs. Defaults to an empty spec; TOML must populate it.
    fem::FemSpec fem;

    // Newton solver settings (tolerances, scheme, max_iter, max_staggered).
    pfm::SolverSettings solver;

    // Number of equal load increments for load_factor: 0 -> 1.
    // Used only when step.mode == StepMode::Uniform.
    int N_steps = 15;

    // Load-stepping schedule (uniform vs. two-stage displacement control).
    StepControl step;

    // Write a VTK snapshot every `vtk_every` ACCEPTED steps (1 = every step).
    // Step 0 (initial state) and the final converged step are always written.
    int vtk_every = 1;

    // Mirror all terminal output (std::cout / std::cerr) to a per-run text
    // file named "<base_name>_run.txt".
    bool write_log = true;

    // Launch the Gmsh GUI after the run.
    bool show_gui = true;
};

// Parse a TOML file into a fresh AppConfig. Throws on parse or validation
// errors. The TOML schema is documented at the top of the example
// `config.toml` shipped with this project.
[[nodiscard]] AppConfig loadConfigFromToml(const std::string& path);

// Apply CLI overrides on top of an already-loaded AppConfig. Recognised
// flags (all optional, all override the TOML value):
//
//   --out NAME      output base filename (no extension)
//   --W val         specimen width
//   --H val         specimen height (also resets y_crack = H/2 unless --ycrack
//                   is also given)
//   --ycrack val    crack y-coordinate
//   --a val         crack length
//   --hfine val     fine mesh size near the crack
//   --hfar val      coarse mesh size in the bulk
//   --steps N       number of load increments
//   --no-gui        do not open Gmsh GUI at end
//   --gui           force Gmsh GUI on at end
//
// Unknown flags produce a warning on stderr and are ignored. The TOML path
// (argv[1]) is *skipped* by this routine -- main() consumes it separately.
void applyCliOverrides(AppConfig& cfg, int argc, char** argv);

// Sanity-check a fully-populated AppConfig. Throws std::invalid_argument with
// a descriptive message if any field is out of bounds.
void validate(const AppConfig& cfg);

}  // namespace appcfg
