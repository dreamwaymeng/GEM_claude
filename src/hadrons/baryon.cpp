#include "gem/hadrons/baryon.hpp"
#include <cmath>
#include <stdexcept>

namespace gem {

Baryon::Baryon(QuarkFlavor flavor1, QuarkFlavor flavor2, QuarkFlavor flavor3,
               int twice_total_spin)
    : m1_(constants::getQuarkMass(flavor1)),
      m2_(constants::getQuarkMass(flavor2)),
      m3_(constants::getQuarkMass(flavor3)),
      twice_total_spin_(twice_total_spin),
      pair12_identical_(flavor1 == flavor2),  // Σ-type if same flavor
      jacobi_(m1_, m2_, m3_),
      binding_energy_(0.0),
      valid_(false) {
    if (twice_total_spin != 1 && twice_total_spin != 3) {
        throw std::invalid_argument("Baryon total spin must be 1/2 (twice_S=1) or 3/2 (twice_S=3)");
    }
    initialize();
}

Baryon::Baryon(Real m1, Real m2, Real m3, int twice_total_spin)
    : m1_(m1), m2_(m2), m3_(m3),
      twice_total_spin_(twice_total_spin),
      pair12_identical_(std::abs(m1 - m2) < 1e-6),  // Assume identical if same mass
      jacobi_(m1, m2, m3),
      binding_energy_(0.0),
      valid_(false) {
    if (twice_total_spin != 1 && twice_total_spin != 3) {
        throw std::invalid_argument("Baryon total spin must be 1/2 (twice_S=1) or 3/2 (twice_S=3)");
    }
    initialize();
}

void Baryon::initialize() {
    // Set default parameters for both coordinates - r_min/r_max in fm
    basis_params_rho_.n_basis = constants::default_basis::N_BASIS;
    basis_params_rho_.r_min = constants::default_basis::R_MIN;
    basis_params_rho_.r_max = constants::default_basis::R_MAX;

    basis_params_lambda_ = basis_params_rho_;

    pot_params_.alpha_s_coul = constants::default_potential::ALPHA_S_COULOMB;
    pot_params_.alpha_s_hyp = constants::default_potential::ALPHA_S_HYPERFINE;
    pot_params_.b = constants::default_potential::B_STRING;
    pot_params_.V0 = constants::default_potential::V0;
    pot_params_.hyp_a = constants::default_potential::HYPERFINE_A;
    pot_params_.hyp_b_exp = constants::default_potential::HYPERFINE_B_EXP;
}

void Baryon::setBasis(const BasisParameters& params_rho, const BasisParameters& params_lambda) {
    basis_params_rho_ = params_rho;
    basis_params_lambda_ = params_lambda;
}

void Baryon::setBasis(int n_basis, Real r_min, Real r_max) {
    basis_params_rho_.n_basis = n_basis;
    basis_params_rho_.r_min = r_min;
    basis_params_rho_.r_max = r_max;
    basis_params_lambda_ = basis_params_rho_;
}

void Baryon::setPotential(const PotentialParameters& params) {
    pot_params_ = params;
}

void Baryon::setConfig(const Config& config) {
    basis_params_rho_ = config.getBasisParameters();
    basis_params_lambda_ = config.getBasisParameters();
    pot_params_ = config.getPotentialParameters();
}

Real Baryon::calculateMass() {
    // Create bases for ρ and λ coordinates
    basis_rho_ = std::make_unique<GaussianBasis>(basis_params_rho_);
    basis_lambda_ = std::make_unique<GaussianBasis>(basis_params_lambda_);

    // Build three-body Hamiltonian
    std::array<Real, 3> masses = {m1_, m2_, m3_};
    ThreeBodyHamiltonian hamiltonian(*basis_rho_, *basis_lambda_,
                                     jacobi_.mu_rho(), jacobi_.mu_lambda(),
                                     masses, pot_params_);
    hamiltonian.build(twice_total_spin_, pair12_identical_);

    // Solve eigenvalue problem
    EigenvalueSolver solver(EigenvalueSolver::Method::CholeskyReduce);
    EigenResult result = solver.solve(hamiltonian);

    if (!result.converged) {
        valid_ = false;
        return 0.0;
    }

    valid_ = true;
    eigenvalues_ = result.eigenvalues;
    wave_function_ = result.eigenvectors.col(0);

    // Ground state binding energy
    binding_energy_ = eigenvalues_(0);

    // Total mass = constituent masses + binding energy
    return totalMass();
}

} // namespace gem
