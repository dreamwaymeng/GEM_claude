#include "gem/basis/gaussian_basis.hpp"
#include <cmath>
#include <stdexcept>

namespace gem {

GaussianBasis::GaussianBasis(int n_basis, Real nu_min, Real nu_max)
    : n_basis_(n_basis), nu_(n_basis), norm_(n_basis) {
    if (n_basis < 1) {
        throw std::invalid_argument("n_basis must be at least 1");
    }
    if (nu_min <= 0 || nu_max <= 0) {
        throw std::invalid_argument("Width parameters must be positive");
    }
    if (nu_min >= nu_max && n_basis > 1) {
        throw std::invalid_argument("nu_min must be less than nu_max");
    }
    initializeGeometricSeries(nu_min, nu_max);
}

GaussianBasis::GaussianBasis(const BasisParameters& params)
    : GaussianBasis(params.n_basis, params.nu_min(), params.nu_max()) {}

void GaussianBasis::initializeGeometricSeries(Real nu_min, Real nu_max) {
    // Generate width parameters in geometric series: ν_n = ν_1 * r^(n-1)
    // where r = (ν_max/ν_min)^(1/(n-1)) for n > 1
    if (n_basis_ == 1) {
        nu_[0] = nu_min;
    } else {
        Real ratio = std::pow(nu_max / nu_min, 1.0 / (n_basis_ - 1));
        for (int n = 0; n < n_basis_; ++n) {
            nu_[n] = nu_min * std::pow(ratio, n);
        }
    }

    // Compute normalization constants for l=0: N_n = (2ν_n/π)^(3/4)
    for (int n = 0; n < n_basis_; ++n) {
        norm_[n] = std::pow(2.0 * nu_[n] / constants::PI, 0.75);
    }
}

Real GaussianBasis::evaluate(int n, Real r) const {
    return norm_[n] * std::exp(-nu_[n] * r * r);
}

Real GaussianBasis::overlap(int n, int m) const {
    // <φ_n|φ_m> = N_n * N_m * (π / (ν_n + ν_m))^(3/2)
    Real nu_sum = nu_[n] + nu_[m];
    return norm_[n] * norm_[m] * std::pow(constants::PI / nu_sum, 1.5);
}

Real GaussianBasis::kinetic(int n, int m, Real mu) const {
    // T = -ℏ²/(2μ) ∇²
    // ∇² exp(-νr²) = (4ν²r² - 6ν) exp(-νr²)
    // <φ_n|∇²|φ_m> = -6 ν_n ν_m / (ν_n + ν_m) × <φ_n|φ_m>
    // <φ_n|T|φ_m> = (ℏc)² / (2μ) * 6 * ν_n * ν_m / (ν_n + ν_m) * <φ_n|φ_m>
    // Note: μ in GeV, ν in fm^-2, (ℏc)² in (GeV·fm)², result in GeV
    Real nu_sum = nu_[n] + nu_[m];
    Real nu_prod = nu_[n] * nu_[m];
    Real ovlp = overlap(n, m);

    // ℏ²/(2μ) with μ in GeV: (ℏc)² = 0.0389 (GeV·fm)²
    Real factor = constants::HBAR_C_SQUARED / (2.0 * mu);

    return factor * 6.0 * nu_prod / nu_sum * ovlp;
}

Real GaussianBasis::coulomb(int n, int m) const {
    // <φ_n|1/r|φ_m> = N_n * N_m * 2π / (ν_n + ν_m)
    Real nu_sum = nu_[n] + nu_[m];
    return norm_[n] * norm_[m] * 2.0 * constants::PI / nu_sum;
}

Real GaussianBasis::linear(int n, int m) const {
    // <φ_n|r|φ_m> = N_n * N_m * 2π / (ν_n + ν_m)²
    // Derivation: ∫ r³ exp(-αr²) dr = 1/(2α²), so with 4π prefactor:
    // 4π × 1/(2(ν_n+ν_m)²) = 2π/(ν_n+ν_m)²
    Real nu_sum = nu_[n] + nu_[m];
    return norm_[n] * norm_[m] * 2.0 * constants::PI / (nu_sum * nu_sum);
}

Real GaussianBasis::gaussian_potential(int n, int m, Real alpha) const {
    // <φ_n|exp(-α*r²)|φ_m> = N_n * N_m * (π / (ν_n + ν_m + α))^(3/2)
    Real nu_sum_alpha = nu_[n] + nu_[m] + alpha;
    return norm_[n] * norm_[m] * std::pow(constants::PI / nu_sum_alpha, 1.5);
}

Real GaussianBasis::delta_smeared(int n, int m, Real sigma) const {
    // Smeared delta function: δ_σ(r) = (α/π)^(3/2) exp(-α*r²) where α = 1/σ²
    // <φ_n|δ_σ|φ_m> = (α/π)^(3/2) * <φ_n|exp(-α*r²)|φ_m>
    //              = (α/π)^(3/2) * N_n * N_m * (π / (ν_n + ν_m + α))^(3/2)
    Real alpha = 1.0 / (sigma * sigma);
    Real prefactor = std::pow(alpha / constants::PI, 1.5);
    return prefactor * gaussian_potential(n, m, alpha);
}

Matrix GaussianBasis::overlapMatrix() const {
    Matrix S(n_basis_, n_basis_);
    for (int n = 0; n < n_basis_; ++n) {
        for (int m = 0; m <= n; ++m) {
            Real val = overlap(n, m);
            S(n, m) = val;
            S(m, n) = val;
        }
    }
    return S;
}

Matrix GaussianBasis::kineticMatrix(Real mu) const {
    Matrix T(n_basis_, n_basis_);
    for (int n = 0; n < n_basis_; ++n) {
        for (int m = 0; m <= n; ++m) {
            Real val = kinetic(n, m, mu);
            T(n, m) = val;
            T(m, n) = val;
        }
    }
    return T;
}

Matrix GaussianBasis::coulombMatrix() const {
    Matrix V(n_basis_, n_basis_);
    for (int n = 0; n < n_basis_; ++n) {
        for (int m = 0; m <= n; ++m) {
            Real val = coulomb(n, m);
            V(n, m) = val;
            V(m, n) = val;
        }
    }
    return V;
}

Matrix GaussianBasis::linearMatrix() const {
    Matrix V(n_basis_, n_basis_);
    for (int n = 0; n < n_basis_; ++n) {
        for (int m = 0; m <= n; ++m) {
            Real val = linear(n, m);
            V(n, m) = val;
            V(m, n) = val;
        }
    }
    return V;
}

Matrix GaussianBasis::deltaMatrix(Real sigma) const {
    Matrix D(n_basis_, n_basis_);
    for (int n = 0; n < n_basis_; ++n) {
        for (int m = 0; m <= n; ++m) {
            Real val = delta_smeared(n, m, sigma);
            D(n, m) = val;
            D(m, n) = val;
        }
    }
    return D;
}

} // namespace gem
