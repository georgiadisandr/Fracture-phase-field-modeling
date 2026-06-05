// Config.cpp -- TOML loader + CLI override layer.
//
// Uses the vendored toml++ single header (tomlplusplus.hpp). The schema is
// kept simple and flat; see the example config.toml for a worked instance.

#include "Config.h"

#include "tomlplusplus.hpp"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace appcfg {

namespace {

template <typename T>
bool readScalar(const toml::table& tbl, std::string_view key, T& out)
{
    if (auto v = tbl[key].value<T>()) {
        out = *v;
        return true;
    }
    return false;
}

bool readDouble(const toml::table& tbl, std::string_view key, double& out)
{
    auto node = tbl[key];
    if (auto v = node.value<double>())  { out = *v;                      return true; }
    if (auto i = node.value<int64_t>()) { out = static_cast<double>(*i); return true; }
    return false;
}

bool readString(const toml::table& tbl, std::string_view key, std::string& out)
{
    if (auto v = tbl[key].value<std::string>()) {
        out = *v;
        return true;
    }
    return false;
}

std::vector<double> readDoubleArray(const toml::node_view<const toml::node>& node,
                                    std::string_view ctx)
{
    std::vector<double> out;
    if (!node) return out;
    const auto* arr = node.as_array();
    if (!arr)
        throw std::runtime_error(std::string(ctx) + ": expected an array");
    out.reserve(arr->size());
    for (const auto& e : *arr) {
        if (auto v = e.value<double>())      out.push_back(*v);
        else if (auto i = e.value<int64_t>()) out.push_back(static_cast<double>(*i));
        else
            throw std::runtime_error(std::string(ctx) + ": array element is not a number");
    }
    return out;
}

std::vector<int> readIntArray(const toml::node_view<const toml::node>& node,
                              std::string_view ctx)
{
    std::vector<int> out;
    if (!node) return out;
    const auto* arr = node.as_array();
    if (!arr)
        throw std::runtime_error(std::string(ctx) + ": expected an array");
    out.reserve(arr->size());
    for (const auto& e : *arr) {
        if (auto v = e.value<int64_t>())
            out.push_back(static_cast<int>(*v));
        else
            throw std::runtime_error(std::string(ctx) + ": array element is not an integer");
    }
    return out;
}

void loadMeshSection(const toml::table& root, mesh::Config& m)
{
    const auto* tbl = root["mesh"].as_table();
    if (!tbl) return;
    readDouble(*tbl, "W",       m.W);
    readDouble(*tbl, "H",       m.H);
    readDouble(*tbl, "a",       m.a);
    readDouble(*tbl, "y_crack", m.y_crack);
    readDouble(*tbl, "h_fine",  m.h_fine);
    readDouble(*tbl, "h_far",   m.h_far);
    readDouble(*tbl, "r_fine",  m.r_fine);
    readDouble(*tbl, "r_far",   m.r_far);
    readString(*tbl, "base_name", m.base_name);
    if (tbl->contains("H") && !tbl->contains("y_crack"))
        m.y_crack = m.H / 2.0;
}

void loadFemSection(const toml::table& root, fem::FemSpec& spec)
{
    if (const auto* tbl = root["fem"].as_table()) {
        int64_t tmp = 0;
        if (readScalar<int64_t>(*tbl, "ntype", tmp)) spec.ntype = static_cast<int>(tmp);
        if (readScalar<int64_t>(*tbl, "ndofn", tmp)) spec.ndofn = static_cast<int>(tmp);
        if (readScalar<int64_t>(*tbl, "ngaus", tmp)) spec.ngaus = static_cast<int>(tmp);
        if (readScalar<int64_t>(*tbl, "nstre", tmp)) spec.nstre = static_cast<int>(tmp);

        // Strain-energy split model. Case-insensitive; default = Amor.
        std::string es;
        if (readString(*tbl, "energy_split", es)) {
            for (auto& ch : es) ch = static_cast<char>(std::tolower(ch));
            if      (es == "none")     spec.split = pfm::SplitModel::None;
            else if (es == "lancioni") spec.split = pfm::SplitModel::Lancioni;
            else if (es == "amor")     spec.split = pfm::SplitModel::Amor;
            else if (es == "spectral") spec.split = pfm::SplitModel::Spectral;
            else
                throw std::runtime_error(
                    "fem.energy_split must be one of: none, lancioni, amor, "
                    "spectral (got \"" + es + "\")");
        }
    }

    std::unordered_map<std::string, int> name_to_index;

    if (const auto* arr = root["materials"].as_array()) {
        for (std::size_t i = 0; i < arr->size(); ++i) {
            const auto* mt = arr->get(i)->as_table();
            if (!mt)
                throw std::runtime_error("materials[" + std::to_string(i) + "] is not a table");
            fem::Material mat;
            mat.props = readDoubleArray((*mt)["props"],
                                        "materials[" + std::to_string(i) + "].props");
            if (mat.props.empty())
                throw std::runtime_error("materials[" + std::to_string(i)
                                          + "].props is missing or empty");
            spec.materials.push_back(std::move(mat));

            std::string name;
            if (readString(*mt, "name", name))
                name_to_index[name] = static_cast<int>(i);
        }
    }
    if (spec.materials.empty())
        throw std::runtime_error("[materials] must define at least one material");

    if (const auto* arr = root["material_for_group"].as_array()) {
        for (std::size_t i = 0; i < arr->size(); ++i) {
            const auto* mt = arr->get(i)->as_table();
            if (!mt)
                throw std::runtime_error("material_for_group[" + std::to_string(i)
                                          + "] is not a table");
            std::string group;
            if (!readString(*mt, "group", group))
                throw std::runtime_error("material_for_group[" + std::to_string(i)
                                          + "].group is required");

            int matidx = 0;
            if (auto v = (*mt)["material"].value<int64_t>()) {
                matidx = static_cast<int>(*v);
            } else if (auto s = (*mt)["material"].value<std::string>()) {
                auto it = name_to_index.find(*s);
                if (it == name_to_index.end())
                    throw std::runtime_error(
                        "material_for_group[" + std::to_string(i)
                        + "]: unknown material name \"" + *s + "\"");
                matidx = it->second;
            } else {
                throw std::runtime_error(
                    "material_for_group[" + std::to_string(i)
                    + "].material must be an integer index or a string name");
            }
            if (matidx < 0 || matidx >= static_cast<int>(spec.materials.size()))
                throw std::runtime_error(
                    "material_for_group[" + std::to_string(i)
                    + "].material index out of range");

            spec.material_for_group.emplace_back(std::move(group), matidx);
        }
    }

    if (const auto* arr = root["bcs"].as_array()) {
        for (std::size_t i = 0; i < arr->size(); ++i) {
            const auto* mt = arr->get(i)->as_table();
            if (!mt)
                throw std::runtime_error("bcs[" + std::to_string(i) + "] is not a table");
            fem::BCSpec bc;
            if (!readString(*mt, "physical_name", bc.physical_name))
                throw std::runtime_error("bcs[" + std::to_string(i)
                                          + "].physical_name is required");
            bc.flags  = readIntArray((*mt)["flags"],
                                     "bcs[" + std::to_string(i) + "].flags");
            bc.values = readDoubleArray((*mt)["values"],
                                        "bcs[" + std::to_string(i) + "].values");
            if (static_cast<int>(bc.flags.size())  != spec.ndofn ||
                static_cast<int>(bc.values.size()) != spec.ndofn)
                throw std::runtime_error("bcs[" + std::to_string(i)
                                          + "]: flags/values length must equal ndofn ("
                                          + std::to_string(spec.ndofn) + ")");
            spec.bcs.push_back(std::move(bc));
        }
    }

    // Neumann BCs: optional [[neumann]] entries, traction [tx, ty].
    if (const auto* arr = root["neumann"].as_array()) {
        for (std::size_t i = 0; i < arr->size(); ++i) {
            const auto* mt = arr->get(i)->as_table();
            if (!mt)
                throw std::runtime_error("neumann[" + std::to_string(i)
                                          + "] is not a table");
            fem::NeumannSpec nb;
            if (!readString(*mt, "physical_name", nb.physical_name))
                throw std::runtime_error("neumann[" + std::to_string(i)
                                          + "].physical_name is required");
            nb.traction = readDoubleArray((*mt)["traction"],
                                          "neumann[" + std::to_string(i)
                                          + "].traction");
            if (static_cast<int>(nb.traction.size()) != spec.ndofn)
                throw std::runtime_error("neumann[" + std::to_string(i)
                                          + "]: traction length must equal ndofn ("
                                          + std::to_string(spec.ndofn) + ")");
            spec.neumann.push_back(std::move(nb));
        }
    }

    // Point loads: optional [[point_load]] entries -- a Cartesian force
    // [Fx, Fy] applied at every node of a 0D physical group.
    if (const auto* arr = root["point_load"].as_array()) {
        for (std::size_t i = 0; i < arr->size(); ++i) {
            const auto* mt = arr->get(i)->as_table();
            if (!mt)
                throw std::runtime_error("point_load[" + std::to_string(i)
                                          + "] is not a table");
            fem::PointLoadSpec pl;
            if (!readString(*mt, "physical_name", pl.physical_name))
                throw std::runtime_error("point_load[" + std::to_string(i)
                                          + "].physical_name is required");
            pl.force = readDoubleArray((*mt)["force"],
                                       "point_load[" + std::to_string(i)
                                       + "].force");
            if (static_cast<int>(pl.force.size()) != spec.ndofn)
                throw std::runtime_error("point_load[" + std::to_string(i)
                                          + "]: force length must equal ndofn ("
                                          + std::to_string(spec.ndofn) + ")");
            spec.point_loads.push_back(std::move(pl));
        }
    }

    // Initial phi overrides: optional [[initial_phi]] entries. Each entry
    // pairs a physical group (any dimension) with a scalar phi value.
    if (const auto* arr = root["initial_phi"].as_array()) {
        for (std::size_t i = 0; i < arr->size(); ++i) {
            const auto* mt = arr->get(i)->as_table();
            if (!mt)
                throw std::runtime_error("initial_phi[" + std::to_string(i)
                                          + "] is not a table");
            fem::InitialPhiSpec ip;
            if (!readString(*mt, "physical_name", ip.physical_name))
                throw std::runtime_error("initial_phi[" + std::to_string(i)
                                          + "].physical_name is required");
            if (!readDouble(*mt, "value", ip.value))
                throw std::runtime_error("initial_phi[" + std::to_string(i)
                                          + "].value is required");
            spec.initial_phi.push_back(std::move(ip));
        }
    }
}

void loadSolverSection(const toml::table& root, pfm::SolverSettings& s)
{
    const auto* tbl = root["solver"].as_table();
    if (!tbl) return;
    readDouble(*tbl, "tol_rel",  s.tol_rel);
    readDouble(*tbl, "tol_abs",  s.tol_abs);
    int64_t mi = 0;
    if (readScalar<int64_t>(*tbl, "max_iter", mi))
        s.max_iter = static_cast<int>(mi);
    readScalar(*tbl, "verbose",  s.verbose);

    // Coupling scheme: "monolithic" (default) or "staggered". Case-insensitive.
    std::string scheme;
    if (readString(*tbl, "scheme", scheme)) {
        for (auto& c : scheme) c = static_cast<char>(std::tolower(c));
        if (scheme == "monolithic")
            s.scheme = pfm::SolverScheme::Monolithic;
        else if (scheme == "staggered")
            s.scheme = pfm::SolverScheme::Staggered;
        else
            throw std::runtime_error(
                "solver.scheme must be \"monolithic\" or \"staggered\" "
                "(got \"" + scheme + "\")");
    }

    int64_t ms = 0;
    if (readScalar<int64_t>(*tbl, "max_staggered", ms))
        s.max_staggered = static_cast<int>(ms);
}

void loadRunSection(const toml::table& root, AppConfig& cfg)
{
    const auto* tbl = root["run"].as_table();
    if (!tbl) return;
    int64_t ns = 0;
    if (readScalar<int64_t>(*tbl, "N_steps", ns))
        cfg.N_steps = static_cast<int>(ns);
    readScalar(*tbl, "show_gui", cfg.show_gui);
}

}  // namespace

// ===========================================================================
//  Public API
// ===========================================================================

AppConfig loadConfigFromToml(const std::string& path)
{
    AppConfig cfg;
    toml::table root;
    try {
        root = toml::parse_file(path);
    } catch (const toml::parse_error& e) {
        throw std::runtime_error("TOML parse error in \"" + path + "\": "
                                  + std::string(e.description()));
    }

    loadMeshSection  (root, cfg.mesh);
    loadFemSection   (root, cfg.fem);
    loadSolverSection(root, cfg.solver);
    loadRunSection   (root, cfg);

    return cfg;
}

void applyCliOverrides(AppConfig& cfg, int argc, char** argv)
{
    auto need = [&](int i, std::string_view flag) {
        if (i + 1 >= argc)
            throw std::invalid_argument(std::string(flag) + " requires a value");
    };

    bool saw_H = false, saw_ycrack = false;

    for (int i = 2; i < argc; ++i) {
        const std::string_view arg(argv[i]);
        if      (arg == "--gui")     { cfg.show_gui = true; }
        else if (arg == "--no-gui")  { cfg.show_gui = false; }
        else if (arg == "--out")     { need(i, arg); cfg.mesh.base_name = argv[++i]; }
        else if (arg == "--W")       { need(i, arg); cfg.mesh.W       = std::stod(argv[++i]); }
        else if (arg == "--H")       { need(i, arg); cfg.mesh.H       = std::stod(argv[++i]); saw_H = true; }
        else if (arg == "--ycrack")  { need(i, arg); cfg.mesh.y_crack = std::stod(argv[++i]); saw_ycrack = true; }
        else if (arg == "--a")       { need(i, arg); cfg.mesh.a       = std::stod(argv[++i]); }
        else if (arg == "--hfine")   { need(i, arg); cfg.mesh.h_fine  = std::stod(argv[++i]); }
        else if (arg == "--hfar")    { need(i, arg); cfg.mesh.h_far   = std::stod(argv[++i]); }
        else if (arg == "--steps")   { need(i, arg); cfg.N_steps      = std::stoi(argv[++i]); }
        else {
            std::cerr << "Warning: ignoring unknown argument \"" << arg << "\"\n";
        }
    }

    if (saw_H && !saw_ycrack)
        cfg.mesh.y_crack = cfg.mesh.H / 2.0;
}

void validate(const AppConfig& cfg)
{
    const auto& m = cfg.mesh;
    if (m.W <= 0.0 || m.H <= 0.0)
        throw std::invalid_argument("mesh.W and mesh.H must be positive");
    if (m.a <= 0.0 || m.a >= m.W)
        throw std::invalid_argument("mesh.a (crack length) must be in (0, W)");
    if (m.y_crack <= 0.0 || m.y_crack >= m.H)
        throw std::invalid_argument("mesh.y_crack must be in (0, H)");
    if (m.h_fine <= 0.0 || m.h_far <= 0.0 || m.h_fine > m.h_far)
        throw std::invalid_argument("require 0 < mesh.h_fine <= mesh.h_far");
    if (m.r_fine <= 0.0 || m.r_far <= 0.0 || m.r_fine > m.r_far)
        throw std::invalid_argument("require 0 < mesh.r_fine <= mesh.r_far");
    if (m.base_name.empty())
        throw std::invalid_argument("mesh.base_name must be non-empty");

    if (cfg.fem.materials.empty())
        throw std::invalid_argument("[materials] must define at least one material");
    if (cfg.fem.bcs.empty())
        throw std::invalid_argument("[[bcs]] must define at least one Dirichlet BC");

    for (std::size_t i = 0; i < cfg.fem.materials.size(); ++i) {
        if (cfg.fem.materials[i].props.size() < 5) {
            throw std::invalid_argument(
                "materials[" + std::to_string(i)
                + "].props must contain at least 5 values: { E, nu, Gc, l0, k }");
        }
    }

    if (cfg.N_steps <= 0)
        throw std::invalid_argument("run.N_steps must be > 0");
    if (cfg.solver.max_iter <= 0)
        throw std::invalid_argument("solver.max_iter must be > 0");
    if (cfg.solver.tol_rel <= 0.0 || cfg.solver.tol_abs <= 0.0)
        throw std::invalid_argument("solver.tol_rel and solver.tol_abs must be > 0");
    if (cfg.solver.scheme == pfm::SolverScheme::Staggered &&
        cfg.solver.max_staggered <= 0)
        throw std::invalid_argument(
            "solver.max_staggered must be > 0 when solver.scheme = \"staggered\"");
}

}  // namespace appcfg
