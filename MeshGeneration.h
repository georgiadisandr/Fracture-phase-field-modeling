#pragma once

#include <string>
#include <unordered_map>
#include <vector>

namespace mesh {

//Specimen configuration
struct Config {
    // All quantities are in millimetres.
    // Geometry
    double W       = 100.0;  // specimen width
    double H       = 40.0;   // specimen height
    double a       = 50.0;   // crack length
    double y_crack = 20.0;   // crack y-coordinate 

    // Mesh sizes
    double h_fine = 0.10;    // size near the crack
    double h_far  = 0.4;     // size in the bulk

    // Refinement zone radii
    double r_fine = 1.0;     // fully-fine zone around the crack faces
    double r_far  = 10.0;    // transition radius back to h_far

    // Output
    std::string base_name = "cracked_specimen_refined_line";
};

// Tags of the gmsh entities created by build_geometry.
struct Geometry {
    //Surfaces
    int upperSurf = 0; 
    int lowerSurf = 0;

    //tip of the crack   
    int p_tip    = 0;                
    
    //edges
    int lBottom  = 0, lTop      = 0;
    int lRightLo = 0, lRightUp  = 0;
    int lLeftUp  = 0, lLeftLo   = 0;
    int lCrackUp = 0, lCrackLo  = 0;
    int lLig     = 0; // ligament (crack tip -> right edge)                
};

// Build all gmsh entities for the cracked specimen and synchronize the CAD.
[[nodiscard]] Geometry build_geometry(const Config& cfg);

// Configure the background size field (Distance + Threshold on crack faces & tip).
void configure_size_field(const Config& cfg, const Geometry& g);

// Set recombination + meshing options for quad-dominant linear elements.
void configure_meshing(const Geometry& g);

// Add physical groups for boundaries, crack faces, ligament, tip, and the domain.
void add_physical_groups(const Geometry& g);

// Print mesh statistics (nodes, quads, tris) to stdout.
void print_mesh_stats();

// Caller is responsible for gmsh::initialize / gmsh::model::add / gmsh::write
// / gmsh::finalize around this call.
Geometry generate(const Config& cfg);

}  // namespace mesh
