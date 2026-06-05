#pragma once
//
// Small element-level helpers shared by every assembly routine (stiffness,
// residual, post-processing). They live in namespace fem because they only
// depend on the mesh / reference element, not on the physics.
//
// Quick map of who builds what:
//
//   Gauss.{h,cpp}         raw quadrature rules on the reference element
//   ShapeFunc.{h,cpp}     reference shape functions and their (xi,eta) derivs
//   Jacob2.{h,cpp}        reference -> physical mapping (cartd, djacb)
//   ElementUtils.{h,cpp}  glue that the physics calls: a Gauss-point list
//                         tailored to the element family, and the standard
//                         strain-displacement / phase-field B matrices.
//

#include <Eigen/Dense>
#include <vector>

namespace fem {

// One Gauss point on the reference element. eta is unused for 1D quads / lines
// but kept zero by callers for those (we don't use 1D right now).
struct GaussPoint {
    double xi  = 0.0;
    double eta = 0.0;
    double w   = 0.0;
};

// Gauss-point list for the chosen element family.
//   nnode = 3        : 3-point Dunavant rule on the reference triangle
//                       (ngaus_per_dir is ignored)
//   nnode = 4 or 8   : tensor-product Gauss-Legendre on [-1, 1]^2,
//                       n1d = ngaus_per_dir (defaults: 2 for Q4, 3 for Q8)
std::vector<GaussPoint> gaussPoints(int nnode, int ngaus_per_dir);

// Strain-displacement matrix B_u (3 x 2*nnode) in Voigt, engineering-shear.
// Column block for node i:
//        | dN_i/dx     0       |
//        |   0       dN_i/dy   |
//        | dN_i/dy   dN_i/dx   |
//
// cartd[idime][inode] is the Cartesian-derivative table produced by jacob2.
Eigen::MatrixXd buildBu(const std::vector<std::vector<double>>& cartd,
                        int nnode);

// Phase-field gradient matrix B_phi (ndime x nnode). Column i is grad N_i.
// Thin wrapper around cartd, kept here for symmetry with buildBu so call
// sites read uniformly.
Eigen::MatrixXd buildBphi(const std::vector<std::vector<double>>& cartd,
                          int ndime,
                          int nnode);

}  // namespace fem
