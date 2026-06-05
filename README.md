# PhaseField_test1 — 2D Phase-Field Fracture Solver

A C++ finite-element solver for **quasi-static brittle fracture** using the
phase-field method. It meshes a pre-cracked specimen, ramps an applied load
from 0 → 1 in increments, and at each increment solves the coupled
displacement + damage (phase-field) problem, writing a VTK snapshot per step
for ParaView.

- **Language / standard:** C++17
- **Dependencies:** [Eigen](https://eigen.tuxfamily.org) (dense + sparse linear
  algebra), [Gmsh SDK](https://gmsh.info) (mesh generation + I/O), and a vendored
  single-header [toml++](https://marzer.github.io/tomlplusplus/) (`tomlplusplus.hpp`).
- **Units:** lengths in mm, stresses in MPa, fracture energy `Gc` in N/mm.

---

## Building

The project uses CMake. Eigen and the Gmsh SDK paths are set near the top of
`CMakeLists.txt` — edit them to match your machine.

```sh
cmake -S . -B build              # configures (defaults to an optimized Release build)
cmake --build build              # builds build/phasefield_1.exe
```

On Windows the Gmsh DLL is copied next to the executable automatically.

> **Note on build type:** the code is Eigen-heavy (dense block products, sparse
> assembly, `SparseLU`). A Debug build runs 10–50× slower, so `CMakeLists.txt`
> defaults to `Release` (`-O3 -DNDEBUG`) unless you ask for Debug.

## Running

```sh
./build/phasefield_1.exe config.toml [options]
```

The TOML file is mandatory and defines the full run (geometry, mesh sizes,
material, boundary conditions, solver, load steps). Any TOML value can be
overridden on the command line, e.g. `--steps 30 --W 80 --no-gui`. Run with
`--help` for the full list. See `config.toml` for the documented schema.

**Outputs:** a base `.msh` / `.vtk` mesh plus a numbered VTK series
(`<base>_step_000.vtk`, `_001`, …) that ParaView opens as a time series.

---

## How it works (pipeline)

1. **Load config** (`Config.cpp`) — parse TOML, apply CLI overrides, validate.
2. **Generate mesh** (`MeshGeneration.cpp`) — build the cracked geometry in Gmsh
   with refinement near the crack, write `.msh` / `.vtk`.
3. **Read FEM input** (`GmshReader.cpp`) — pull nodes, elements, materials, and
   boundary conditions back out of the live Gmsh model into a `FemInput`.
4. **Load-step loop** (`main.cpp`) — walk `load_factor` 0 → 1. Each increment is
   solved to equilibrium by the solver; on non-convergence the increment is
   halved and retried (adaptive sub-stepping). After each converged step the
   irreversible history field is raised and a VTK snapshot is written.
5. **Solve** (`Solver.cpp`) — one Newton sequence per increment, using either the
   **monolithic** or **staggered** scheme.

---

## Module / function reference

### `main.cpp` — driver
Owns the whole run: loads config, runs Gmsh init/finalize, drives the adaptive
quasi-static load-stepping loop, calls the solver per step, maintains the
history field, writes per-step VTK, and prints timing.

### `Config.{h,cpp}` — run configuration
Parses the TOML file and holds it in an `AppConfig` (mesh `Config`, `FemSpec`,
solver settings, `N_steps`, `show_gui`).
- `loadConfigFromToml(path)` — parse a TOML file into a fresh `AppConfig`.
- `applyCliOverrides(cfg, argc, argv)` — overlay command-line overrides.
- `validate(cfg)` — sanity-check geometry, mesh, material, and solver values.

### `MeshGeneration.{h,cpp}` — geometry + meshing (Gmsh)
Builds the rectangular specimen with an embedded horizontal crack and the
refinement field around it. `Config` holds the knobs (`W`, `H`, `a`, `y_crack`,
`h_fine`, `h_far`, `r_fine`, `r_far`); `Geometry` holds the resulting Gmsh tags.
- `generate(cfg)` — top-level: build geometry, size field, mesh, physical groups.
- `configure_size_field(cfg, g)` — set fine/coarse element sizes vs. distance to crack.
- `configure_meshing(g)` — meshing algorithm / element-order options.
- `add_physical_groups(g)` — name boundaries (`Top`, `Bottom`, `Domain`, `CrackTip`, …).
- `print_mesh_stats()` — report node/element counts.

### `GmshReader.{h,cpp}` — model → FEM data
Translates the Gmsh model + the `FemSpec` (materials, BCs, group→material map,
optional Neumann / point-load / initial-phi specs) into a `FemInput`.
- `readFemInputFromGmsh(spec)` — reads, in order: nodes & coordinates, 2D
  elements (tri3 / quad4 / quad8), material resolution per group, Dirichlet BCs,
  Neumann tractions, concentrated point loads, and initial-phi overrides.

### `FEMINPUT.h` — central data structure
Header-only `FemInput` struct shared by every element routine: mesh
connectivity (`coord`, `conn`, `offset`, `matno`), problem constants (`ntype`,
`ndofn`, `ngaus`, `nstre`), material property tables + precomputed `MatCache`es,
Dirichlet data, and Neumann/point-load/initial-phi entries. Helper
`count_by_nnode(n)` counts elements with `n` nodes.

### `ConstitutiveModel.{h,cpp}` — material law + energy splits
The single home for the constitutive law (isotropic linear elasticity, 2D).
- `degradation(phi, k)` — phase-field degradation `g(phi) = (1-phi)^2 + k`.
- `SplitMaterial::from(E, nu, ntype)` — 2D elastic moduli (`mu`, `K2D`).
- `MatCache::build(E, nu, ntype, split)` — precompute per-material matrices once.
- `energySplit(eps, mc, model)` — returns `sigma^+`, `sigma^-`, tangents `C^+`,
  `C^-`, and the crack-driving energy `psi^+` for the chosen split model
  (`none` / `lancioni` / `amor` / `spectral`).

### FEM primitives
- `Gauss.{h,cpp}` — `triGauss(ngaus)` / `quadGauss(ngaus)`: Gauss-point
  coordinates and weights for triangles and quads.
- `ShapeFunc.{h,cpp}` — `shapeFunc(s, t, nnode)`: shape functions and their
  natural-coordinate derivatives at a point.
- `Jacob2.{h,cpp}` — `jacob2(...)`: element Jacobian, its determinant, Cartesian
  shape-function derivatives, and global Gauss-point coordinates.
- `ElementUtils.{h,cpp}` — `gaussPoints(nnode, ngaus)` plus the strain-displacement
  matrices `buildBu(...)` (displacement) and `buildBphi(...)` (phase field).

### `StiffnessPFM.{h,cpp}` — element + global assembly (the core)
Element-level physics (residual **and** tangent) in one shared routine, plus the
split forms used by the staggered scheme. `MatParams` holds `E, nu, Gc, l0, k`.
- `elementSystem(...)` — element residual `re` and (optionally) tangent `Ke` for
  the full coupled `[u | phi]` element.
- `assembleGlobalSystem(d, ...)` — assemble the global sparse tangent `K` and
  residual `R` (a `GlobalSystem`) in one pass.
- `elementUSystem(...)` / `assembleU_System(...)` — displacement sub-problem (phi frozen).
- `elementPhiSystem(...)` / `assemblePhi_System(...)` — phase-field sub-problem (u frozen).

### `ResidualPFM.{h,cpp}` — residual-only assembly
- `assembleGlobalResidual(d, u, phi, history)` — global internal-force residual
  (length `3*npoin`, layout `[u | phi]`) when the tangent is not needed, e.g. the
  staggered convergence check and monolithic post-loop residual.

### `History.{h,cpp}` — irreversibility
- `HistoryField` — per-element, per-Gauss-point crack-driving history `H`.
- `makeHistoryField(d)` — allocate and zero it.
- `updateHistory(d, u, history)` — after a converged step, raise `H` to
  `max(H_stored, psi^+)` so the crack cannot heal.

### `BoundaryConditions.{h,cpp}` — Dirichlet
- `applyDirichlet(K, R, ..., load_factor)` — impose prescribed displacements on
  the global system (scaled by `load_factor`).

### `NeumannBC.{h,cpp}` — external forces
- `assembleExternalForce(d)` — global external-force vector from surface
  tractions and concentrated point loads (scaled per step by `load_factor`).

### `Solver.{h,cpp}` — Newton solver for one load step
`solveStep(d, u, phi, history, load_factor, settings)` runs one Newton sequence,
returning a `SolverResult` (iterations, residual, converged flag). Two schemes
(`SolverSettings::scheme`):
- **Monolithic** — solve the full coupled `[u | phi]` system together. Fast
  quadratic convergence, but the tangent can be indefinite during crack growth.
- **Staggered** — alternate minimisation: solve `u` (phi frozen), then `phi`
  (u frozen), repeat until the coupled residual is small. More robust during
  crack propagation; linear outer convergence.

### `Output.{h,cpp}` — results
- `writeVTK(filename, d, u, phi)` — write a legacy VTK file with the displacement
  and phase-field fields for ParaView.

---


