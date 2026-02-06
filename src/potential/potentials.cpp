#include "gem/potential/potentials.hpp"
#include <cmath>

namespace gem {

// CoulombPotential implementation
// ALL UNITS IN GeV (masses, energies) and fm (distances)

CoulombPotential::CoulombPotential(Real alpha_s, Real color_factor)
    : alpha_s_(alpha_s), color_factor_(color_factor) {}

Real CoulombPotential::evaluate(Real r) const {
    if (r < 1e-10) return 0.0;  // Regularize at origin
    // V = α_s × color_factor × ℏc / r
    // With r in fm, ℏc in GeV·fm, result in GeV
    return alpha_s_ * color_factor_ * constants::HBAR_C / r;
}

Matrix CoulombPotential::buildMatrix(const GaussianBasis& basis) const {
    // V_nm = α_s × color_factor × ℏc × <φ_n|1/r|φ_m>
    // <1/r> has units of 1/fm (from Gaussian integrals with ν in GeV²)
    // Actually, need to be careful: the Gaussian integrals give dimensionless results
    // when ν is in appropriate units. Let me reconsider.
    //
    // The basis functions are φ(r) = N exp(-ν r²) with r in fm, ν in GeV²
    // For this to work: ν r² must be dimensionless
    // ν [GeV²] × r² [fm²] needs conversion: ν r² = ν [GeV²] × r² [fm²] / (ℏc)² [GeV²·fm²]
    //
    // The coulomb integral <1/r> in fm^-1 gives:
    // V = α_s × color × ℏc [GeV·fm] × <1/r> [fm^-1] = α_s × color × ℏc / <r> [GeV]
    Matrix V = basis.coulombMatrix();
    V *= alpha_s_ * color_factor_ * constants::HBAR_C;
    return V;
}

// ConfinementPotential implementation

ConfinementPotential::ConfinementPotential(Real b, Real color_factor, Real V0)
    : b_(b), color_factor_(color_factor), V0_(V0) {}

Real ConfinementPotential::evaluate(Real r) const {
    // V = (-3b/4 × r + V_c) × color_factor
    // b is in GeV², r in fm, V_c in GeV
    // b [GeV²] × r [fm] / ℏc [GeV·fm] = b/ℏc × r [GeV]
    Real b_GeV_per_fm = b_ / constants::HBAR_C;
    Real linear_coeff = -0.75 * b_GeV_per_fm;  // -3b/4 in GeV/fm
    return linear_coeff * color_factor_ * r + V0_ * color_factor_;
}

Matrix ConfinementPotential::buildMatrix(const GaussianBasis& basis) const {
    // V_nm = (-3b/4 × <r> + V_c) × color_factor
    Real b_GeV_per_fm = b_ / constants::HBAR_C;
    Real linear_coeff = -0.75 * b_GeV_per_fm;  // -3b/4 in GeV/fm

    Matrix V_linear = basis.linearMatrix();
    V_linear *= linear_coeff * color_factor_;

    if (std::abs(V0_) > 1e-10) {
        Matrix S = basis.overlapMatrix();
        V_linear += V0_ * color_factor_ * S;
    }

    return V_linear;
}

// HyperfinePotential implementation

HyperfinePotential::HyperfinePotential(Real alpha_s_hyp, Real m1, Real m2,
                                       Real hyp_a, Real hyp_b_exp,
                                       Real color_factor, Real spin_factor)
    : alpha_s_(alpha_s_hyp), m1_(m1), m2_(m2),
      hyp_a_(hyp_a), hyp_b_exp_(hyp_b_exp),
      color_factor_(color_factor), spin_factor_(spin_factor) {
    computeTau();
}

void HyperfinePotential::computeTau() {
    // tau = (2*mu)^b / a, where mu is reduced mass in GeV
    // m1_, m2_ are now in GeV
    Real mu = m1_ * m2_ / (m1_ + m2_);  // reduced mass in GeV
    tau_ = std::pow(2.0 * mu, hyp_b_exp_) / hyp_a_;  // tau in GeV
}

void HyperfinePotential::setMasses(Real m1, Real m2) {
    m1_ = m1;
    m2_ = m2;
    computeTau();
}

Real HyperfinePotential::evaluate(Real r) const {
    // V_hyp = -(8π α_s^hyp / 3m_1 m_2) × (τ³/π^(3/2)) × exp(-τ²r²) × color × spin
    // With masses in GeV, τ in GeV, r in fm, result in GeV

    // Convert tau from GeV to fm^-1: tau_fm = tau_GeV / ℏc [GeV·fm]
    Real tau_fm = tau_ / constants::HBAR_C;  // tau in fm^-1

    Real tau_cubed = tau_fm * tau_fm * tau_fm;
    // Prefactor: 8π α_s / (3 m1 m2) with masses in GeV -> units of GeV^-2
    Real prefactor = (8.0 * constants::PI * alpha_s_ / 3.0) / (m1_ * m2_);
    Real normalization = tau_cubed / std::pow(constants::PI, 1.5);
    Real gaussian = std::exp(-tau_fm * tau_fm * r * r);
    // (ℏc)³ in GeV³·fm³ to convert tau_fm³ [fm^-3] to proper units
    Real hbar_c_cubed = constants::HBAR_C * constants::HBAR_C * constants::HBAR_C;

    // Result: [GeV^-2] × [fm^-3] × [GeV³·fm³] = GeV
    return -prefactor * normalization * gaussian * spin_factor_ * color_factor_ * hbar_c_cubed;
}

Matrix HyperfinePotential::buildMatrix(const GaussianBasis& basis) const {
    // Convert tau from GeV to fm^-1
    Real tau_fm = tau_ / constants::HBAR_C;

    Real prefactor = (8.0 * constants::PI * alpha_s_ / 3.0) / (m1_ * m2_);
    Real tau_cubed = tau_fm * tau_fm * tau_fm;
    Real normalization = tau_cubed / std::pow(constants::PI, 1.5);
    Real hbar_c_cubed = constants::HBAR_C * constants::HBAR_C * constants::HBAR_C;
    Real total_prefactor = -prefactor * normalization * spin_factor_ * color_factor_ * hbar_c_cubed;

    int n = basis.size();
    Matrix V(n, n);
    Real alpha = tau_fm * tau_fm;  // tau² in fm^-2

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            Real val = total_prefactor * basis.gaussian_potential(i, j, alpha);
            V(i, j) = val;
            V(j, i) = val;
        }
    }

    return V;
}

// CornellPotential implementation

CornellPotential::CornellPotential(Real alpha_s, Real b, Real color_factor, Real V0)
    : coulomb_(alpha_s, color_factor), confinement_(b, color_factor, V0) {}

Real CornellPotential::evaluate(Real r) const {
    return coulomb_.evaluate(r) + confinement_.evaluate(r);
}

Matrix CornellPotential::buildMatrix(const GaussianBasis& basis) const {
    return coulomb_.buildMatrix(basis) + confinement_.buildMatrix(basis);
}

// QuarkPotential implementation

QuarkPotential::QuarkPotential(const PotentialParameters& params, Real m1, Real m2,
                               Real color_factor, Real spin_factor)
    : coulomb_(params.alpha_s_coul, color_factor),
      confinement_(params.b, color_factor, params.V0),  // V0 already in GeV
      hyperfine_(params.alpha_s_hyp, m1, m2, params.hyp_a, params.hyp_b_exp,
                 color_factor, spin_factor) {}

Real QuarkPotential::evaluate(Real r) const {
    return coulomb_.evaluate(r) + confinement_.evaluate(r) + hyperfine_.evaluate(r);
}

Matrix QuarkPotential::buildMatrix(const GaussianBasis& basis) const {
    return coulomb_.buildMatrix(basis) +
           confinement_.buildMatrix(basis) +
           hyperfine_.buildMatrix(basis);
}

} // namespace gem
