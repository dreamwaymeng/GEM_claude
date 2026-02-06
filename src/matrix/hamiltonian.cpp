#include "gem/matrix/hamiltonian.hpp"
#include <cmath>

namespace gem {

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

void ThreeBodyHamiltonian::build(int twice_S) {
    buildKineticMatrix();
    buildPotentialMatrix(twice_S);
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

void ThreeBodyHamiltonian::buildPotentialMatrix(int twice_S) {
    int n_rho = basis_rho_.size();
    int n_lambda = basis_lambda_.size();

    // Color factor for baryon (antisymmetric 3̄)
    Real color_factor = constants::color_factor::BARYON_QQ;

    // Spin factor depends on total spin
    Real spin_factor;
    if (twice_S == 1) {
        // S = 1/2: mixed symmetry, average spin-spin interaction
        spin_factor = constants::spin_factor::SPIN_0;  // Simplified
    } else {
        // S = 3/2: symmetric, all pairs have s=1
        spin_factor = constants::spin_factor::SPIN_1;
    }

    // Get masses
    Real m1 = masses_[0], m2 = masses_[1], m3 = masses_[2];
    Real m12 = m1 + m2;

    // Pair coefficients for (ρ, λ)
    // r12² = ρ²
    // r13² = (m2/m12)² ρ² + λ² + 2(m2/m12) ρ·λ
    // r23² = (m1/m12)² ρ² + λ² - 2(m1/m12) ρ·λ

    // For l=0, the angular integration simplifies
    // <V(r_ij)> where r_ij² = a² ρ² + b² λ² + 2ab ρλ cos(θ)

    // Build potential for each pair
    V_ = Matrix::Zero(n_rho * n_lambda, n_rho * n_lambda);

    // This is a simplified implementation for l=0 only
    // Full implementation would include angular momentum coupling

    // Pair 1-2: pure ρ dependence
    Matrix V_12 = Matrix::Zero(n_rho * n_lambda, n_rho * n_lambda);
    {
        QuarkPotential pot_12(params_, m1, m2, color_factor, spin_factor);
        Matrix V_rho_12 = pot_12.buildMatrix(basis_rho_);
        Matrix S_lambda = basis_lambda_.overlapMatrix();

        for (int i_rho = 0; i_rho < n_rho; ++i_rho) {
            for (int i_lambda = 0; i_lambda < n_lambda; ++i_lambda) {
                int i = i_rho * n_lambda + i_lambda;
                for (int j_rho = 0; j_rho < n_rho; ++j_rho) {
                    for (int j_lambda = 0; j_lambda < n_lambda; ++j_lambda) {
                        int j = j_rho * n_lambda + j_lambda;
                        V_12(i, j) = V_rho_12(i_rho, j_rho) * S_lambda(i_lambda, j_lambda);
                    }
                }
            }
        }
    }

    // For pairs 1-3 and 2-3, the potential depends on both ρ and λ
    // This requires a more complex integration or approximation
    // Here we use a simplified approach treating them similarly

    V_ = V_12;  // Simplified: only 1-2 pair fully implemented

    // Full implementation would add V_13 and V_23 with proper transformations
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
