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
    int N_steps = 15;

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
