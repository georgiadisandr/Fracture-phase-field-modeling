#include "Jacob2.h"
#include <vector>
#include <iostream>

namespace fem {

Jacobian jacob2(const std::vector<std::vector<double>>& elcod,
                const ShapeData& sh_f,
                int nnode, int ndime)
{
    // const int ndime=2;
    Jacobian J;
    J.cartd.assign(ndime, std::vector<double>(nnode, 0.0));// cartesian Derivatives [∂N/∂x ∂N/∂y]
    J.gcpod.assign(ndime, 0.0); //Gauss point global coordinates Σ Ν_i*x_i

    // Global coordinates at this Gauss point
    for (int i = 0; i < ndime; ++i) {
        double acc = 0.0;
        for (int n = 0; n < nnode; ++n) {
            acc += elcod[i][n] * sh_f.Shape[n]; 
        }
        J.gcpod[i] = acc;
    }

    // Jacobian matrix (2D)
    double xjacm[2][2] = {{0.0, 0.0}, {0.0, 0.0}};//initialize J = [0]
    for (int i = 0; i < ndime; ++i) {
        for (int j = 0; j < ndime; ++j) {
            double acc = 0.0;
            for (int n = 0; n < nnode; ++n) {
                acc += sh_f.Deriv[i][n] * elcod[j][n]; //Σ Νd_i*x_i
            }
            xjacm[i][j] = acc;
        }
    }

    // Determinant det(J)
    J.djacb = xjacm[0][0] * xjacm[1][1] - xjacm[0][1] * xjacm[1][0];
    if (J.djacb <= 0.0) {
        std::cerr << "Zero or negative area for element " << '\n';
        return J;            // bail out before dividing by it
    }

    // Inverse — one reciprocal, four multiplies
    const double invDet = 1.0 / J.djacb;
    const double xjaci[2][2] = {
        {  xjacm[1][1] * invDet, -xjacm[0][1] * invDet },
        { -xjacm[1][0] * invDet,  xjacm[0][0] * invDet }
    };

    // Cartesian derivatives — write straight into J.cartd
    for (int i = 0; i < ndime; ++i) {
        for (int n = 0; n < nnode; ++n) {
            double acc = 0.0;
            for (int j = 0; j < ndime; ++j) {
                acc += xjaci[i][j] * sh_f.Deriv[j][n];
            }
            J.cartd[i][n] = acc;
        }
    }

    return J;
}

}  // namespace fem