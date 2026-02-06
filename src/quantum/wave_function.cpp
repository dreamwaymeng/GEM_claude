#include "gem/quantum/wave_function.hpp"
#include "gem/basis/gaussian_basis.hpp"
#include <cmath>
#include <stdexcept>

namespace gem {

// DiscreteWaveFunction implementation

DiscreteWaveFunction::DiscreteWaveFunction(const ComplexVector& color, const ComplexVector& spin)
    : color_(color), spin_(spin) {}

Real DiscreteWaveFunction::colorFactor() const {
    return GellMann::colorFactor(color_);
}

Real DiscreteWaveFunction::spinFactor() const {
    return SpinOperators::spinFactor(spin_);
}

bool DiscreteWaveFunction::isNormalized(Real tol) const {
    Complex color_norm = color_.squaredNorm();
    Complex spin_norm = spin_.squaredNorm();
    return std::abs(color_norm.real() - 1.0) < tol &&
           std::abs(spin_norm.real() - 1.0) < tol;
}

void DiscreteWaveFunction::normalize() {
    Real color_norm = std::sqrt(color_.squaredNorm());
    Real spin_norm = std::sqrt(spin_.squaredNorm());
    if (color_norm > 0) color_ /= color_norm;
    if (spin_norm > 0) spin_ /= spin_norm;
}

// WaveFunctionFactory implementation

DiscreteWaveFunction WaveFunctionFactory::meson(int total_spin, int Sz) {
    ComplexVector color = GellMann::mesonColorSinglet();
    ComplexVector spin;

    if (total_spin == 0) {
        spin = SpinOperators::spinSinglet();
    } else if (total_spin == 1) {
        spin = SpinOperators::spinTriplet(Sz);
    } else {
        throw std::invalid_argument("Meson total spin must be 0 or 1");
    }

    return DiscreteWaveFunction(color, spin);
}

DiscreteWaveFunction WaveFunctionFactory::baryon(int twice_total_spin, int twice_Sz) {
    ComplexVector color = GellMann::baryonColorSinglet();
    ComplexVector spin = mixedSymmetrySpinState(twice_total_spin, twice_Sz);
    return DiscreteWaveFunction(color, spin);
}

ComplexVector WaveFunctionFactory::colorSinglet(int n_quarks) {
    switch (n_quarks) {
        case 2:
            return GellMann::mesonColorSinglet();
        case 3:
            return GellMann::baryonColorSinglet();
        default:
            throw std::invalid_argument("Color singlet not implemented for this quark number");
    }
}

ComplexVector WaveFunctionFactory::symmetricSpinState(int n_quarks, int twice_S) {
    if (n_quarks == 2) {
        if (twice_S == 0) return SpinOperators::spinSinglet();
        if (twice_S == 2) return SpinOperators::spinTriplet(0);
    }
    throw std::invalid_argument("Symmetric spin state not implemented");
}

ComplexVector WaveFunctionFactory::mixedSymmetrySpinState(int twice_S, int twice_Sz) {
    // For 3 quarks, the spin wave function is in 2⊗2⊗2 = 8 dimensional space
    // Total spin can be S=1/2 (mixed symmetry) or S=3/2 (symmetric)
    ComplexVector wf = ComplexVector::Zero(8);

    // Indices: |s1,s2,s3> where si = 0 (up) or 1 (down)
    // |000> = 0, |001> = 1, |010> = 2, |011> = 3,
    // |100> = 4, |101> = 5, |110> = 6, |111> = 7

    if (twice_S == 1) {
        // S = 1/2 (mixed symmetry)
        if (twice_Sz == 1) {
            // Sz = +1/2
            // Mixed symmetry state: ρ-type
            // |1/2, 1/2>_ρ = (1/√6)(2|↑↑↓> - |↑↓↑> - |↓↑↑>)
            wf(1) = 2.0 / std::sqrt(6.0);   // |↑↑↓> = |001>
            wf(2) = -1.0 / std::sqrt(6.0);  // |↑↓↑> = |010>
            wf(4) = -1.0 / std::sqrt(6.0);  // |↓↑↑> = |100>
        } else if (twice_Sz == -1) {
            // Sz = -1/2
            wf(3) = -2.0 / std::sqrt(6.0);  // |↑↓↓> = |011>
            wf(5) = 1.0 / std::sqrt(6.0);   // |↓↑↓> = |101>
            wf(6) = 1.0 / std::sqrt(6.0);   // |↓↓↑> = |110>
        }
    } else if (twice_S == 3) {
        // S = 3/2 (symmetric)
        if (twice_Sz == 3) {
            wf(0) = 1.0;  // |↑↑↑>
        } else if (twice_Sz == 1) {
            Real norm = 1.0 / std::sqrt(3.0);
            wf(1) = norm;  // |↑↑↓>
            wf(2) = norm;  // |↑↓↑>
            wf(4) = norm;  // |↓↑↑>
        } else if (twice_Sz == -1) {
            Real norm = 1.0 / std::sqrt(3.0);
            wf(3) = norm;  // |↑↓↓>
            wf(5) = norm;  // |↓↑↓>
            wf(6) = norm;  // |↓↓↑>
        } else if (twice_Sz == -3) {
            wf(7) = 1.0;  // |↓↓↓>
        }
    }

    return wf;
}

// FullWaveFunction implementation

FullWaveFunction::FullWaveFunction(const GaussianBasis& basis, const Vector& coefficients,
                                   const DiscreteWaveFunction& discrete)
    : basis_(basis), coefficients_(coefficients), discrete_(discrete) {}

Real FullWaveFunction::evaluateSpatial(Real r) const {
    Real value = 0.0;
    for (int n = 0; n < basis_.size(); ++n) {
        value += coefficients_(n) * basis_.evaluate(n, r);
    }
    return value;
}

Real FullWaveFunction::meanRadius() const {
    // <r> = Σ_{nm} c_n c_m <φ_n|r|φ_m>
    Matrix R = basis_.linearMatrix();
    return coefficients_.transpose() * R * coefficients_;
}

Real FullWaveFunction::meanRadiusSquared() const {
    // <r²> = Σ_{nm} c_n c_m <φ_n|r²|φ_m>
    // For Gaussian basis, this can be computed analytically
    Real r2 = 0.0;
    for (int n = 0; n < basis_.size(); ++n) {
        for (int m = 0; m < basis_.size(); ++m) {
            // <φ_n|r²|φ_m> = (3/2) / (ν_n + ν_m) * <φ_n|φ_m>
            Real nu_sum = basis_.nu(n) + basis_.nu(m);
            Real ovlp = basis_.overlap(n, m);
            r2 += coefficients_(n) * coefficients_(m) * (1.5 / nu_sum) * ovlp;
        }
    }
    return r2;
}

Real FullWaveFunction::rmsRadius() const {
    return std::sqrt(meanRadiusSquared());
}

} // namespace gem
