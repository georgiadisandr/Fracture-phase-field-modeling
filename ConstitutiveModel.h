#pragma once
//
// Constitutive helpers shared by the PFM element routines (stiffness, residual,
// history). Keep ALL knowledge of the constitutive law in this one module.
//
// Conventions (isotropic linear elasticity, 2D):
//   * Strain Voigt order:  eps   = { eps_xx, eps_yy, gamma_xy }
//     (engineering shear, so gamma_xy = du/dy + dv/dx).
//   * Stress Voigt order:  sigma = { sigma_xx, sigma_yy, tau_xy }.
//   * Degradation       :  g(phi) = (1 - phi)^2 + k,   k a small residual.
//
// The phase field degrades only the "positive" (crack-driving) part sigma^+
// of the stress; sigma^- is protected. Which part is which is set by the
// SplitModel (see below). The element routines call energySplit() to obtain
// sigma^+, sigma^-, their tangents C^+, C^-, and the crack-driving energy
// psi^+ (which feeds the irreversible history field).
//

#include <Eigen/Dense>

namespace pfm {

// Phase-field degradation function:  g(phi) = (1 - phi)^2 + k.
double degradation(double phi, double k);

// ---------------------------------------------------------------------------
// Strain-energy split model. Selected per run via the TOML key
// `fem.energy_split`.
// ---------------------------------------------------------------------------
enum class SplitModel {
    None,        // 1. no split
    Lancioni,    // 2. deviatoric-volumetric (Lancioni & Royer-Carfagni)
    Amor,        // 3. Amor, Marigo, Maurini
    Spectral     // 4. spectral decomposition of the strain (Miehe)
};

// 2D elastic moduli used by every split model.
struct SplitMaterial {
    double mu  = 0.0;
    double K2D = 0.0;

    static SplitMaterial from(double E, double nu, int ntype);
};

// ---------------------------------------------------------------------------
// Per-material precomputed cache.
//
// Everything in here is a function of the material parameters (E, nu, ntype)
// and the split model -- not of the current strain or iteration. Building it
// ONCE up front and indexing by material id at the Gauss-point level removes
// the per-element-per-iteration cost of rebuilding D_d, Pi and C (and, for
// Amor / Lancioni / None, the closed-form C_plus / C_minus).
//
// The cached C_plus / C_minus pairs depend only on the sign of the volumetric
// strain theta:
//   _tens  -- value of C_plus / C_minus when theta > 0  (tensile  case)
//   _comp  -- value of C_plus / C_minus when theta <= 0 (compressive case)
//
// For Lancioni and None both pairs coincide (theta-independent).
// For Spectral the C_plus / C_minus fields are unused (the tangent is
// finite-differenced from spectralSigmaPlus and depends on the strain).
// ---------------------------------------------------------------------------
struct MatCache {
    SplitMaterial   sm;
    Eigen::Matrix3d C            = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d D_d          = Eigen::Matrix3d::Zero();

    Eigen::Matrix3d C_plus_tens  = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d C_minus_tens = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d C_plus_comp  = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d C_minus_comp = Eigen::Matrix3d::Zero();

    static MatCache build(double E, double nu, int ntype, SplitModel split);
};

// Result of an energy split at one Gauss point.
struct EnergySplit {
    Eigen::Vector3d sigma_plus  = Eigen::Vector3d::Zero();
    Eigen::Vector3d sigma_minus = Eigen::Vector3d::Zero();
    Eigen::Matrix3d C_plus      = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d C_minus     = Eigen::Matrix3d::Zero();
    double          psi_plus    = 0.0;
};

EnergySplit energySplit(const Eigen::Vector3d& eps,
                        const MatCache&        mc,
                        SplitModel             model);

}  // namespace pfm
