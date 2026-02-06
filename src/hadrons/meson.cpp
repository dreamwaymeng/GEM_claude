#include "gem/hadrons/meson.hpp"
#include <cmath>
#include <stdexcept>

namespace gem {

Meson::Meson(QuarkFlavor flavor1, QuarkFlavor flavor2, int total_spin)
    : m1_(constants::getQuarkMass(flavor1)),
      m2_(constants::getQuarkMass(flavor2)),
      total_spin_(total_spin),
      jacobi_(m1_, m2_),
      binding_energy_(0.0),
      rms_radius_(0.0),
      valid_(false) {
    if (total_spin < 0 || total_spin > 1) {
        throw std::invalid_argument("Meson total spin must be 0 or 1");
    }
    initialize();
}

Meson::Meson(Real m1, Real m2, int total_spin)
    : m1_(m1), m2_(m2), total_spin_(total_spin),
      jacobi_(m1, m2),
      binding_energy_(0.0),
      rms_radius_(0.0),
      valid_(false) {
    if (total_spin < 0 || total_spin > 1) {
        throw std::invalid_argument("Meson total spin must be 0 or 1");
    }
    initialize();
}

void Meson::initialize() {
    // Set default parameters - r_min/r_max in fm
    basis_params_.n_basis = constants::default_basis::N_BASIS;
    basis_params_.r_min = constants::default_basis::R_MIN;
    basis_params_.r_max = constants::default_basis::R_MAX;

    pot_params_.alpha_s_coul = constants::default_potential::ALPHA_S_COULOMB;
    pot_params_.alpha_s_hyp = constants::default_potential::ALPHA_S_HYPERFINE;
    pot_params_.b = constants::default_potential::B_STRING;
    pot_params_.V0 = constants::default_potential::V0;
    pot_params_.hyp_a = constants::default_potential::HYPERFINE_A;
    pot_params_.hyp_b_exp = constants::default_potential::HYPERFINE_B_EXP;
}

void Meson::setBasis(const BasisParameters& params) {
    basis_params_ = params;
}

void Meson::setBasis(int n_basis, Real r_min, Real r_max) {
    basis_params_.n_basis = n_basis;
    basis_params_.r_min = r_min;
    basis_params_.r_max = r_max;
}

void Meson::setPotential(const PotentialParameters& params) {
    pot_params_ = params;
}

void Meson::setConfig(const Config& config) {
    basis_params_ = config.getBasisParameters();
    pot_params_ = config.getPotentialParameters();
}

Real Meson::spinFactor() const {
    // S=0 (pseudoscalar): <s1·s2> = -3/4
    // S=1 (vector): <s1·s2> = +1/4
    return (total_spin_ == 0) ? constants::spin_factor::SPIN_0
                              : constants::spin_factor::SPIN_1;
}

Real Meson::calculateMass() {
    // Create basis
    basis_ = std::make_unique<GaussianBasis>(basis_params_);

    // Create potential with color and spin factors
    Real color_factor = colorFactor();
    Real spin_factor = this->spinFactor();

    potential_ = std::make_unique<QuarkPotential>(
        pot_params_, m1_, m2_, color_factor, spin_factor
    );

    // Build Hamiltonian
    Real mu = jacobi_.reducedMass();
    TwoBodyHamiltonian hamiltonian(*basis_, mu, *potential_);
    hamiltonian.build(color_factor, spin_factor);

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

    // Ground state binding energy (relative to continuum)
    binding_energy_ = eigenvalues_(0);

    // Calculate RMS radius
    Matrix R_matrix = basis_->linearMatrix();
    Matrix S_matrix = basis_->overlapMatrix();
    Real r_expect = wave_function_.transpose() * R_matrix * wave_function_;
    Real norm = wave_function_.transpose() * S_matrix * wave_function_;
    Real r_mean = r_expect / norm;

    // <r²> calculation
    Real r2 = 0.0;
    for (int n = 0; n < basis_->size(); ++n) {
        for (int m = 0; m < basis_->size(); ++m) {
            Real nu_sum = basis_->nu(n) + basis_->nu(m);
            Real ovlp = basis_->overlap(n, m);
            r2 += wave_function_(n) * wave_function_(m) * (1.5 / nu_sum) * ovlp;
        }
    }
    r2 /= norm;
    rms_radius_ = std::sqrt(r2);

    // Total mass = constituent masses + binding energy
    return totalMass();
}

Real Meson::excitedStateEnergy(int n) const {
    if (!valid_ || n >= eigenvalues_.size()) {
        return 0.0;
    }
    return m1_ + m2_ + eigenvalues_(n);
}

} // namespace gem
