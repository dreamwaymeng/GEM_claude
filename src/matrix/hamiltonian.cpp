#include "gem/matrix/hamiltonian.hpp"
#include <cmath>
#include <functional>
#include <array>

namespace gem {

// Helper structures and functions for Cholesky decomposition method
namespace {

// 2x2 transformation matrix for Jacobi coordinate transformation
struct TransformMatrix {
    Real T11, T12, T21, T22;

    // Compute T^{-1}
    TransformMatrix inverse() const {
        Real det = T11 * T22 - T12 * T21;
        return {T22 / det, -T12 / det, -T21 / det, T11 / det};
    }
};

// Compute transformation matrix from set (12)3 to set (ij)k
// For V_13: transform to (13)2 where r_13 = rho_2
// For V_23: transform to (23)1 where r_23 = rho_1
//
// Derived directly from Jacobi coordinate definitions:
// Set 3 (12)3: ρ₃ = r₁ - r₂, λ₃ = CM(12) - r₃
// Set 2 (13)2: ρ₂ = r₁ - r₃, λ₂ = CM(13) - r₂
// Set 1 (23)1: ρ₁ = r₂ - r₃, λ₁ = CM(23) - r₁
TransformMatrix computeTransformMatrix(Real m1, Real m2, Real m3, int target_pair) {
    Real M = m1 + m2 + m3;
    Real M12 = m1 + m2;
    Real M13 = m1 + m3;
    Real M23 = m2 + m3;

    Real T11, T12, T21, T22;

    if (target_pair == 13) {
        // Transform from (12)3 to (13)2
        // ρ₂ = (m₂/M₁₂) ρ₃ + λ₃
        // λ₂ = (m₁ M)/(M₁₂ M₁₃) ρ₃ - (m₃/M₁₃) λ₃
        T11 = m2 / M12;
        T12 = 1.0;
        T21 = m1 * M / (M12 * M13);
        T22 = -m3 / M13;
    } else {
        // Transform from (12)3 to (23)1
        // ρ₁ = -(m₁/M₁₂) ρ₃ + λ₃
        // λ₁ = -(m₂ M)/(M₁₂ M₂₃) ρ₃ - (m₃/M₂₃) λ₃
        T11 = -m1 / M12;
        T12 = 1.0;
        T21 = -m2 * M / (M12 * M23);
        T22 = -m3 / M23;
    }

    return {T11, T12, T21, T22};
}

// Compute matrix element using Cholesky decomposition
// Returns the potential integrals I_coul, I_lin, I_const, I_gauss
struct PotentialIntegrals {
    Real I_coul;    // Coulomb: 1/r
    Real I_lin;     // Linear: r
    Real I_const;   // Constant: 1
    Real I_gauss;   // Gaussian: exp(-tau^2 r^2)
    Real prefactor; // Common prefactor: N_bra * N_ket * pi^(3/2) / (l11 * l22)^3
};

PotentialIntegrals computeCholeskyIntegrals(
    Real nu_rho_bra, Real nu_lambda_bra,
    Real nu_rho_ket, Real nu_lambda_ket,
    const TransformMatrix& T_inv,
    Real tau_fm)  // tau in fm^-1 for hyperfine
{
    // Combined width matrix D = diag(nu_rho_sum, nu_lambda_sum)
    Real nu_rho_sum = nu_rho_bra + nu_rho_ket;
    Real nu_lambda_sum = nu_lambda_bra + nu_lambda_ket;

    // A = (T^-1)^T * D * T^-1
    // A_11 = T_inv_11^2 * nu_rho + T_inv_21^2 * nu_lambda
    // A_12 = T_inv_11 * T_inv_12 * nu_rho + T_inv_21 * T_inv_22 * nu_lambda
    // A_22 = T_inv_12^2 * nu_rho + T_inv_22^2 * nu_lambda
    Real A11 = T_inv.T11 * T_inv.T11 * nu_rho_sum + T_inv.T21 * T_inv.T21 * nu_lambda_sum;
    Real A12 = T_inv.T11 * T_inv.T12 * nu_rho_sum + T_inv.T21 * T_inv.T22 * nu_lambda_sum;
    Real A22 = T_inv.T12 * T_inv.T12 * nu_rho_sum + T_inv.T22 * T_inv.T22 * nu_lambda_sum;

    // Cholesky decomposition: A = L * L^T
    Real l11 = std::sqrt(A11);
    Real l21 = A12 / l11;
    Real l22 = std::sqrt(A22 - l21 * l21);

    // Normalizations
    Real N_rho_bra = std::pow(2.0 * nu_rho_bra / constants::PI, 0.75);
    Real N_lambda_bra = std::pow(2.0 * nu_lambda_bra / constants::PI, 0.75);
    Real N_rho_ket = std::pow(2.0 * nu_rho_ket / constants::PI, 0.75);
    Real N_lambda_ket = std::pow(2.0 * nu_lambda_ket / constants::PI, 0.75);
    Real N_bra = N_rho_bra * N_lambda_bra;
    Real N_ket = N_rho_ket * N_lambda_ket;

    // Correct formulas derived from completing the square in the 6D integral:
    // det(A) = (l11 * l22)^2  (note: A22 = l21^2 + l22^2 from Cholesky)
    //
    // Coulomb: <1/r> = N × 2π^{5/2} / (√A₂₂ × det(A))
    // Linear:  <r>   = N × 2π^{5/2} × √A₂₂ / det(A)²
    // Constant: <1>  = N × π³ / det(A)^{3/2}

    Real sqrt_A22 = std::sqrt(A22);  // A22 already computed above
    Real det_A = l11 * l11 * l22 * l22;  // (l11 * l22)^2

    // Common normalization factor
    Real norm_factor = N_bra * N_ket;

    // Potential integrals with correct formulas
    Real I_coul = norm_factor * 2.0 * std::pow(constants::PI, 2.5) / (sqrt_A22 * det_A);
    Real I_lin = norm_factor * 2.0 * std::pow(constants::PI, 2.5) * sqrt_A22 / (det_A * det_A);
    Real I_const = norm_factor * std::pow(constants::PI, 3.0) / std::pow(det_A, 1.5);

    // For backward compatibility, set prefactor to 1 since I_V now includes everything
    Real prefactor = 1.0;

    // Gaussian: <exp(-τ²r²)> = N × π³ / (det(A) + τ² A₂₂)^{3/2}
    // This comes from completing the square: A_eff + τ² = det(A)/A₂₂ + τ²
    Real tau2 = tau_fm * tau_fm;
    Real I_gauss = norm_factor * std::pow(constants::PI, 3.0)
                   / std::pow(det_A + tau2 * A22, 1.5);

    return {I_coul, I_lin, I_const, I_gauss, prefactor};
}

} // anonymous namespace

// TwoBodyHamiltonian implementation

TwoBodyHamiltonian::TwoBodyHamiltonian(const GaussianBasis& basis, Real reduced_mass,
                                       const QuarkPotential& potential)
    : basis_(basis), reduced_mass_(reduced_mass), potential_(potential) {
    int n = basis_.size();
    H_ = Matrix::Zero(n, n);
    S_ = Matrix::Zero(n, n);
    T_ = Matrix::Zero(n, n);
    V_ = Matrix::Zero(n, n);
}

void TwoBodyHamiltonian::build(Real color_factor, Real spin_factor) {
    int n = basis_.size();

    // Build overlap matrix
    S_ = basis_.overlapMatrix();

    // Build kinetic energy matrix
    T_ = basis_.kineticMatrix(reduced_mass_);

    // Build potential matrix
    V_ = potential_.buildMatrix(basis_);

    // Total Hamiltonian
    H_ = T_ + V_;
}

void TwoBodyHamiltonian::rebuild(const PotentialParameters& params) {
    // Update potential with new parameters (all in GeV units)
    potential_.coulomb().setAlphaS(params.alpha_s_coul);
    potential_.confinement().setStringTension(params.b);
    potential_.confinement().setV0(params.V0);  // V0 already in GeV
    // Note: hyperfine tau is computed from masses, not directly settable

    // Rebuild matrices
    V_ = potential_.buildMatrix(basis_);
    H_ = T_ + V_;
}

// ThreeBodyHamiltonian implementation

ThreeBodyHamiltonian::ThreeBodyHamiltonian(const GaussianBasis& basis_rho,
                                           const GaussianBasis& basis_lambda,
                                           Real mu_rho, Real mu_lambda,
                                           const std::array<Real, 3>& masses,
                                           const PotentialParameters& params)
    : basis_rho_(basis_rho), basis_lambda_(basis_lambda),
      mu_rho_(mu_rho), mu_lambda_(mu_lambda), masses_(masses), params_(params) {

    // Total basis size is n_rho × n_lambda
    int n_rho = basis_rho_.size();
    int n_lambda = basis_lambda_.size();
    int n_total = n_rho * n_lambda;

    H_ = Matrix::Zero(n_total, n_total);
    S_ = Matrix::Zero(n_total, n_total);
    T_ = Matrix::Zero(n_total, n_total);
    V_ = Matrix::Zero(n_total, n_total);
}

void ThreeBodyHamiltonian::build(int twice_S, bool pair12_identical) {
    buildKineticMatrix();
    buildPotentialMatrix(twice_S, pair12_identical);
    H_ = T_ + V_;
}

void ThreeBodyHamiltonian::buildKineticMatrix() {
    int n_rho = basis_rho_.size();
    int n_lambda = basis_lambda_.size();

    // Kinetic and overlap matrices for each coordinate
    Matrix T_rho = basis_rho_.kineticMatrix(mu_rho_);
    Matrix T_lambda = basis_lambda_.kineticMatrix(mu_lambda_);
    Matrix S_rho = basis_rho_.overlapMatrix();
    Matrix S_lambda = basis_lambda_.overlapMatrix();

    // Build product space matrices
    // Index mapping: (i_rho, i_lambda) -> i_rho * n_lambda + i_lambda
    for (int i_rho = 0; i_rho < n_rho; ++i_rho) {
        for (int i_lambda = 0; i_lambda < n_lambda; ++i_lambda) {
            int i = i_rho * n_lambda + i_lambda;

            for (int j_rho = 0; j_rho < n_rho; ++j_rho) {
                for (int j_lambda = 0; j_lambda < n_lambda; ++j_lambda) {
                    int j = j_rho * n_lambda + j_lambda;

                    // Overlap: S_ij = S_rho(i_rho, j_rho) × S_lambda(i_lambda, j_lambda)
                    S_(i, j) = S_rho(i_rho, j_rho) * S_lambda(i_lambda, j_lambda);

                    // Kinetic: T_ij = T_rho × S_lambda + S_rho × T_lambda
                    T_(i, j) = T_rho(i_rho, j_rho) * S_lambda(i_lambda, j_lambda)
                             + S_rho(i_rho, j_rho) * T_lambda(i_lambda, j_lambda);
                }
            }
        }
    }
}

void ThreeBodyHamiltonian::buildPotentialMatrix(int twice_S, bool pair12_identical) {
    int n_rho = basis_rho_.size();
    int n_lambda = basis_lambda_.size();

    // Color factor for baryon (antisymmetric color singlet)
    Real color_factor = constants::color_factor::BARYON_QQ;

    // Spin factors for different pairs depend on total spin and quark content
    // Sum rule: Σ⟨sᵢ·sⱼ⟩ = [S(S+1) - 9/4]/2
    //   S=1/2: sum = -3/4
    //   S=3/2: sum = +3/4
    //
    // For S=1/2 baryons, there are two cases:
    //   Λ-type (all different flavors, e.g., uds, udc):
    //     Pair 12 in spin singlet: ⟨s₁·s₂⟩ = -3/4, ⟨s₁·s₃⟩ = ⟨s₂·s₃⟩ = 0
    //   Σ-type (two identical quarks in pair 12, e.g., uuc, ddc):
    //     Fermi statistics requires pair 12 in spin triplet (antisym color × sym spin × sym space)
    //     ⟨s₁·s₂⟩ = +1/4, ⟨s₁·s₃⟩ = ⟨s₂·s₃⟩ = -1/2 (from sum rule)
    //
    // For S=3/2 baryons: all pairs symmetric (s=1), ⟨sᵢ·sⱼ⟩ = +1/4

    // Get masses
    Real m1 = masses_[0], m2 = masses_[1], m3 = masses_[2];

    Real spin_factor_12, spin_factor_13, spin_factor_23;
    if (twice_S == 1) {
        // S = 1/2 baryon
        if (pair12_identical) {
            // Σ-type: identical quarks in pair 12 must be in spin triplet
            spin_factor_12 = constants::spin_factor::SPIN_1;  // +1/4
            spin_factor_13 = -0.5;  // From sum rule: -3/4 - 1/4 = -1, split equally
            spin_factor_23 = -0.5;
        } else {
            // Λ-type: pair 12 in spin singlet
            spin_factor_12 = constants::spin_factor::SPIN_0;  // -3/4
            spin_factor_13 = 0.0;
            spin_factor_23 = 0.0;
        }
    } else {
        // S = 3/2: all pairs symmetric
        spin_factor_12 = constants::spin_factor::SPIN_1;  // +1/4
        spin_factor_13 = constants::spin_factor::SPIN_1;
        spin_factor_23 = constants::spin_factor::SPIN_1;
    }

    // Initialize potential matrix
    V_ = Matrix::Zero(n_rho * n_lambda, n_rho * n_lambda);

    // ========== Pair 1-2: diagonal in ρ ==========
    {
        QuarkPotential pot_12(params_, m1, m2, color_factor, spin_factor_12);
        Matrix V_rho_12 = pot_12.buildMatrix(basis_rho_);
        Matrix S_lambda = basis_lambda_.overlapMatrix();

        for (int i_rho = 0; i_rho < n_rho; ++i_rho) {
            for (int i_lambda = 0; i_lambda < n_lambda; ++i_lambda) {
                int i = i_rho * n_lambda + i_lambda;
                for (int j_rho = 0; j_rho < n_rho; ++j_rho) {
                    for (int j_lambda = 0; j_lambda < n_lambda; ++j_lambda) {
                        int j = j_rho * n_lambda + j_lambda;
                        V_(i, j) += V_rho_12(i_rho, j_rho) * S_lambda(i_lambda, j_lambda);
                    }
                }
            }
        }
    }

    // ========== Pairs 1-3 and 2-3: use Cholesky decomposition method ==========

    // Compute transformation matrices (and their inverses)
    TransformMatrix T_13 = computeTransformMatrix(m1, m2, m3, 13);
    TransformMatrix T_23 = computeTransformMatrix(m1, m2, m3, 23);
    TransformMatrix T_13_inv = T_13.inverse();
    TransformMatrix T_23_inv = T_23.inverse();

    // Compute hyperfine tau for each pair (in fm^-1)
    Real mu_13 = m1 * m3 / (m1 + m3);
    Real mu_23 = m2 * m3 / (m2 + m3);
    Real tau_13 = std::pow(2.0 * mu_13, params_.hyp_b_exp) / params_.hyp_a;
    Real tau_23 = std::pow(2.0 * mu_23, params_.hyp_b_exp) / params_.hyp_a;
    Real tau_13_fm = tau_13 / constants::HBAR_C;
    Real tau_23_fm = tau_23 / constants::HBAR_C;

    // Hyperfine prefactors (note: I_gauss already includes the Gaussian integral)
    // V_hyp = prefactor * I_gauss where I_gauss = (pi * l11^2 / (l11^2 + tau^2))^(3/2)
    Real hyp_prefactor_13 = -(8.0 * constants::PI * params_.alpha_s_hyp / 3.0) / (m1 * m3)
                           * std::pow(tau_13_fm, 3) / std::pow(constants::PI, 1.5)
                           * spin_factor_13 * color_factor
                           * constants::HBAR_C * constants::HBAR_C * constants::HBAR_C;
    Real hyp_prefactor_23 = -(8.0 * constants::PI * params_.alpha_s_hyp / 3.0) / (m2 * m3)
                           * std::pow(tau_23_fm, 3) / std::pow(constants::PI, 1.5)
                           * spin_factor_23 * color_factor
                           * constants::HBAR_C * constants::HBAR_C * constants::HBAR_C;

    // Coulomb and confinement prefactors
    Real coul_prefactor = params_.alpha_s_coul * color_factor * constants::HBAR_C;
    Real lin_prefactor = -0.75 * params_.b / constants::HBAR_C * color_factor;
    Real const_prefactor = params_.V0 * color_factor;

    // Build V_13 and V_23 using Cholesky decomposition method
    for (int i_rho = 0; i_rho < n_rho; ++i_rho) {
        Real nu_i_rho = basis_rho_.nu(i_rho);

        for (int i_lambda = 0; i_lambda < n_lambda; ++i_lambda) {
            Real nu_i_lambda = basis_lambda_.nu(i_lambda);
            int i = i_rho * n_lambda + i_lambda;

            for (int j_rho = 0; j_rho < n_rho; ++j_rho) {
                Real nu_j_rho = basis_rho_.nu(j_rho);

                for (int j_lambda = 0; j_lambda < n_lambda; ++j_lambda) {
                    Real nu_j_lambda = basis_lambda_.nu(j_lambda);
                    int j = j_rho * n_lambda + j_lambda;

                    // ===== Pair 1-3 using Cholesky =====
                    PotentialIntegrals integrals_13 = computeCholeskyIntegrals(
                        nu_i_rho, nu_i_lambda, nu_j_rho, nu_j_lambda,
                        T_13_inv, tau_13_fm);

                    Real V_coul_13 = coul_prefactor * integrals_13.prefactor * integrals_13.I_coul;
                    Real V_lin_13 = lin_prefactor * integrals_13.prefactor * integrals_13.I_lin;
                    Real V_const_13 = const_prefactor * integrals_13.prefactor * integrals_13.I_const;
                    Real V_hyp_13 = hyp_prefactor_13 * integrals_13.prefactor * integrals_13.I_gauss;

                    V_(i, j) += V_coul_13 + V_lin_13 + V_const_13 + V_hyp_13;

                    // ===== Pair 2-3 using Cholesky =====
                    PotentialIntegrals integrals_23 = computeCholeskyIntegrals(
                        nu_i_rho, nu_i_lambda, nu_j_rho, nu_j_lambda,
                        T_23_inv, tau_23_fm);

                    Real V_coul_23 = coul_prefactor * integrals_23.prefactor * integrals_23.I_coul;
                    Real V_lin_23 = lin_prefactor * integrals_23.prefactor * integrals_23.I_lin;
                    Real V_const_23 = const_prefactor * integrals_23.prefactor * integrals_23.I_const;
                    Real V_hyp_23 = hyp_prefactor_23 * integrals_23.prefactor * integrals_23.I_gauss;

                    V_(i, j) += V_coul_23 + V_lin_23 + V_const_23 + V_hyp_23;
                }
            }
        }
    }

    // ========== Three-body force ==========
    // V_3body = -C_3 / (m1 * m2 * m3)
    // No color, spin, or isospin dependence
    // Matrix element: <Φ|V_3body|Φ'> = V_3body * S
    constexpr Real C_3 = 2.02e-3;  // GeV^4
    Real V_3body = -C_3 / (m1 * m2 * m3);
    V_ += V_3body * S_;
}

// MatrixElements implementation

Real MatrixElements::kinetic(const GaussianBasis& basis, int n, int m, Real mu) {
    return basis.kinetic(n, m, mu);
}

Real MatrixElements::coulomb(const GaussianBasis& basis, int n, int m) {
    return basis.coulomb(n, m);
}

Real MatrixElements::linear(const GaussianBasis& basis, int n, int m) {
    return basis.linear(n, m);
}

Real MatrixElements::gaussian(const GaussianBasis& basis, int n, int m, Real alpha) {
    return basis.gaussian_potential(n, m, alpha);
}

Real MatrixElements::twoBodyInThreeBody(const GaussianBasis& basis_rho,
                                         const GaussianBasis& basis_lambda,
                                         int n_rho, int n_lambda, int m_rho, int m_lambda,
                                         Real a, Real b,
                                         const std::function<Real(Real)>& potential) {
    // This is a complex integral that typically requires numerical methods
    // or specialized analytical formulas for Gaussian integrands

    // For Gaussian basis with Gaussian potential:
    // The integral can be done analytically
    // For other potentials, numerical integration may be needed

    // Placeholder - full implementation would compute:
    // <φ_n_rho(ρ) φ_n_lambda(λ) | V(sqrt(a²ρ² + b²λ²)) | φ_m_rho(ρ) φ_m_lambda(λ)>

    return 0.0;  // Placeholder
}

} // namespace gem
