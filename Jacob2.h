
#ifndef JACOB2_H
#define JACOB2_H
#include <vector>
#include "ShapeFunc.h"

namespace fem {

struct Jacobian {
    std::vector<std::vector<double>> cartd; //Cartesian derivatives
    std::vector<double> gcpod; //Gauss point global coordinates
    double djacb =0.0; //det of Jacobian matrix
};

Jacobian jacob2(
    const std::vector<std::vector<double>>& elcod, //element coordinates
    const ShapeData& sh_f,//shape functions in natural coordinates
    int nnode,int ndime);//number of node per element

}  // namespace fem

#endif //JACOB2_H
