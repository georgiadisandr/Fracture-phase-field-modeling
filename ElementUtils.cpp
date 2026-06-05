#include "ElementUtils.h"

#include "Gauss.h"

#include <stdexcept>
#include <string>

namespace fem {
std::vector<GaussPoint> gaussPoints(int nnode, int ngaus_per_dir)
{
    std::vector<GaussPoint> pts;

    if (nnode == 3) {
        const auto tri = quadrature::triGauss(3);
        pts.reserve(static_cast<std::size_t>(tri.ngaus()));
        for (int g = 0; g < tri.ngaus(); ++g)
            pts.push_back({ tri.xi[g], tri.eta[g], tri.w[g] });
        return pts;
    }

    if (nnode == 4 || nnode == 8) {
        const int n1 = (ngaus_per_dir > 0) ? ngaus_per_dir
                                           : (nnode == 4 ? 2 : 3);
                                           
        //if (ngaus_per_dir > 0) { n1 = ngaus_per_dir;}
        // else { if (nnode == 4) n1 = 2;
        //        else n1 = 3;}              

        const auto q = quadrature::quadGauss(n1);
        pts.reserve(static_cast<std::size_t>(n1) * n1);
        for (int i = 0; i < n1; ++i)
            for (int j = 0; j < n1; ++j)
                pts.push_back({ q.xi[i], q.xi[j], q.w[i] * q.w[j] });
        return pts;
    }

    throw std::runtime_error(
        "fem::gaussPoints: unsupported nnode = " + std::to_string(nnode));
}

Eigen::MatrixXd buildBu(const std::vector<std::vector<double>>& cartd, int nnode)
{
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 2 * nnode);
    for (int i = 0; i < nnode; ++i) {
        const double dNdx = cartd[0][i];
        const double dNdy = cartd[1][i];
        B(0, 2 * i)     = dNdx;
        B(1, 2 * i + 1) = dNdy;
        B(2, 2 * i)     = dNdy;
        B(2, 2 * i + 1) = dNdx;
    }
    return B;
}

Eigen::MatrixXd buildBphi(const std::vector<std::vector<double>>& cartd,
                          int ndime, int nnode)
{
    Eigen::MatrixXd B(ndime, nnode);
    for (int d = 0; d < ndime; ++d)
        for (int i = 0; i < nnode; ++i)
            B(d, i) = cartd[d][i];
    return B;
}

}  // namespace fem
