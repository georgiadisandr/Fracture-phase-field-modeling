#include "History.h"

#include "ConstitutiveModel.h"
#include "ElementUtils.h"
#include "Jacob2.h"
#include "ShapeFunc.h"
#include "StiffnessPFM.h"      // for pfm::MatParams

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace pfm {

HistoryField makeHistoryField(const FemInput& d)
{
    HistoryField history(d.nelem);
    for (int e = 0; e < d.nelem; ++e) {
        const int nnode = d.nnode_of(e);
        const int ngp   = static_cast<int>(
                              fem::gaussPoints(nnode, d.ngaus).size());
        history[e].assign(ngp, 0.0);
    }
    return history;
}

void updateHistory(const FemInput&        d,
                   const Eigen::VectorXd& u,
                   HistoryField&          history)
{
    if (u.size() != 2 * d.npoin)
        throw std::runtime_error("updateHistory: u size != 2*npoin");
    if (static_cast<int>(history.size()) != d.nelem)
        throw std::runtime_error("updateHistory: history size != nelem");

    const int ndime = d.ndime;
    std::vector<int> nodes;

    for (int e = 0; e < d.nelem; ++e) {
        nodes.assign(d.conn.begin() + d.offset[e],
                     d.conn.begin() + d.offset[e + 1]);
        const int nnode = static_cast<int>(nodes.size());

        // History only needs the strain-driven psi^+; the per-material cache
        // supplies the constitutive matrices. MatParams (Gc, l0, k) are not
        // used here -- they appear only in the residual / tangent assembly.
        const MatCache& mc = d.mat_caches[d.matno[e]];

        // Nodal coordinates (jacob2 layout) and element displacement.
        std::vector<std::vector<double>> elcod(ndime,
                                               std::vector<double>(nnode, 0.0));
        Eigen::VectorXd u_elem(2 * nnode);
        for (int i = 0; i < nnode; ++i) {
            const int node = nodes[i];
            for (int dim = 0; dim < ndime; ++dim)
                elcod[dim][i] = d.coord(node, dim);
            u_elem(2 * i)     = u(2 * node);
            u_elem(2 * i + 1) = u(2 * node + 1);
        }

        const auto gps = fem::gaussPoints(nnode, d.ngaus);
        if (static_cast<int>(history[e].size()) != static_cast<int>(gps.size()))
            throw std::runtime_error(
                "updateHistory: history[e] size != number of Gauss points");

        for (int g = 0; g < static_cast<int>(gps.size()); ++g) {
            const fem::ShapeData sh = fem::shapeFunc(gps[g].xi, gps[g].eta, nnode);
            const fem::Jacobian  J  = fem::jacob2(elcod, sh, nnode, ndime);
            const Eigen::MatrixXd Bu = fem::buildBu(J.cartd, nnode);

            const Eigen::Vector3d eps = Bu * u_elem;
            const EnergySplit     sp  = energySplit(eps, mc, d.split);

            // Irreversibility: H can only grow.
            history[e][g] = std::max(history[e][g], sp.psi_plus);
        }
    }
}

}  // namespace pfm
