#include "gem/quantum/spin_operators.hpp"
#include <cmath>
#include <stdexcept>

namespace gem {

const Complex SpinOperators::I_(0.0, 1.0);

Matrix2c SpinOperators::sigma(int i) {
    Matrix2c mat = Matrix2c::Zero();

    switch (i) {
        case 1:  // σ_x
            mat(0, 1) = 1.0;
            mat(1, 0) = 1.0;
            break;
        case 2:  // σ_y
            mat(0, 1) = -I_;
            mat(1, 0) = I_;
            break;
        case 3:  // σ_z
            mat(0, 0) = 1.0;
            mat(1, 1) = -1.0;
            break;
        default:
            throw std::invalid_argument("Pauli index must be 1, 2, or 3");
    }
    return mat;
}

Matrix2c SpinOperators::spin(int i) {
    return sigma(i) / 2.0;
}

std::array<Matrix2c, 3> SpinOperators::allSigma() {
    return {sigma(1), sigma(2), sigma(3)};
}

ComplexMatrix SpinOperators::spinDotSpin() {
    // Compute (σ_1 · σ_2) / 4 in the 2⊗2 = 4 dimensional space
    // s_1 · s_2 = (1/4) Σ_i (σ_i ⊗ σ_i)
    ComplexMatrix result = ComplexMatrix::Zero(4, 4);
    auto sigmas = allSigma();

    for (int i = 0; i < 3; ++i) {
        // Kronecker product σ_i ⊗ σ_i
        for (int a = 0; a < 2; ++a) {
            for (int b = 0; b < 2; ++b) {
                for (int c = 0; c < 2; ++c) {
                    for (int d = 0; d < 2; ++d) {
                        int row = a * 2 + c;
                        int col = b * 2 + d;
                        result(row, col) += sigmas[i](a, b) * sigmas[i](c, d);
                    }
                }
            }
        }
    }
    return result / 4.0;
}

Real SpinOperators::spinFactor(int twice_S) {
    // For two spin-1/2 particles: <s_1 · s_2> = [S(S+1) - 3/2] / 2
    Real S = static_cast<Real>(twice_S) / 2.0;
    return (S * (S + 1.0) - 1.5) / 2.0;
}

Real SpinOperators::spinFactor(int twice_s1, int twice_s2, int twice_S) {
    Real s1 = static_cast<Real>(twice_s1) / 2.0;
    Real s2 = static_cast<Real>(twice_s2) / 2.0;
    Real S = static_cast<Real>(twice_S) / 2.0;
    return (S * (S + 1.0) - s1 * (s1 + 1.0) - s2 * (s2 + 1.0)) / 2.0;
}

Real SpinOperators::spinFactor(const ComplexVector& spin_wf) {
    ComplexMatrix s1s2 = spinDotSpin();
    Complex factor = spin_wf.adjoint() * s1s2 * spin_wf;
    return factor.real();
}

ComplexVector SpinOperators::spinSinglet() {
    // |0,0> = (1/√2)(|↑↓> - |↓↑>)
    // |↑↑> = |0>, |↑↓> = |1>, |↓↑> = |2>, |↓↓> = |3>
    ComplexVector wf = ComplexVector::Zero(4);
    Real norm = 1.0 / std::sqrt(2.0);
    wf(1) = norm;   // |↑↓>
    wf(2) = -norm;  // -|↓↑>
    return wf;
}

ComplexVector SpinOperators::spinTriplet(int Sz) {
    ComplexVector wf = ComplexVector::Zero(4);

    switch (Sz) {
        case 1:  // |1,+1> = |↑↑>
            wf(0) = 1.0;
            break;
        case 0: {  // |1,0> = (1/√2)(|↑↓> + |↓↑>)
            Real norm = 1.0 / std::sqrt(2.0);
            wf(1) = norm;
            wf(2) = norm;
            break;
        }
        case -1:  // |1,-1> = |↓↓>
            wf(3) = 1.0;
            break;
        default:
            throw std::invalid_argument("Sz must be -1, 0, or 1 for triplet");
    }
    return wf;
}

ComplexVector SpinOperators::coupledSpinState(int twice_S, int twice_Sz) {
    if (twice_S == 0 && twice_Sz == 0) {
        return spinSinglet();
    } else if (twice_S == 2) {
        return spinTriplet(twice_Sz / 2);
    }
    throw std::invalid_argument("Invalid S, Sz for two spin-1/2 particles");
}

Real SpinOperators::clebschGordan_half_half(int twice_m1, int twice_m2, int twice_J, int twice_M) {
    // CG coefficients for j1=1/2, j2=1/2
    // Check selection rule: M = m1 + m2
    if (twice_M != twice_m1 + twice_m2) {
        return 0.0;
    }

    const Real sqrt2_inv = 1.0 / std::sqrt(2.0);

    if (twice_J == 0 && twice_M == 0) {
        // J=0, M=0: <1/2,m1,1/2,m2|0,0>
        if (twice_m1 == 1 && twice_m2 == -1) return sqrt2_inv;
        if (twice_m1 == -1 && twice_m2 == 1) return -sqrt2_inv;
    } else if (twice_J == 2) {
        // J=1
        if (twice_M == 2) {
            // M=1: <1/2,1/2,1/2,1/2|1,1>
            if (twice_m1 == 1 && twice_m2 == 1) return 1.0;
        } else if (twice_M == 0) {
            // M=0: <1/2,m1,1/2,m2|1,0>
            if (twice_m1 == 1 && twice_m2 == -1) return sqrt2_inv;
            if (twice_m1 == -1 && twice_m2 == 1) return sqrt2_inv;
        } else if (twice_M == -2) {
            // M=-1: <1/2,-1/2,1/2,-1/2|1,-1>
            if (twice_m1 == -1 && twice_m2 == -1) return 1.0;
        }
    }
    return 0.0;
}

} // namespace gem
