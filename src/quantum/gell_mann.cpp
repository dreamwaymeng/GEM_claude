#include "gem/quantum/gell_mann.hpp"
#include <cmath>

namespace gem {

bool GellMann::tables_initialized_ = false;
std::array<std::array<std::array<Real, 8>, 8>, 8> GellMann::f_abc_;
std::array<std::array<std::array<Real, 8>, 8>, 8> GellMann::d_abc_;

Matrix3c GellMann::lambda(int a) {
    Matrix3c mat = Matrix3c::Zero();
    const Complex I(0.0, 1.0);
    const Real sqrt3 = std::sqrt(3.0);

    switch (a) {
        case 1:  // λ_1
            mat(0, 1) = 1.0;
            mat(1, 0) = 1.0;
            break;
        case 2:  // λ_2
            mat(0, 1) = -I;
            mat(1, 0) = I;
            break;
        case 3:  // λ_3
            mat(0, 0) = 1.0;
            mat(1, 1) = -1.0;
            break;
        case 4:  // λ_4
            mat(0, 2) = 1.0;
            mat(2, 0) = 1.0;
            break;
        case 5:  // λ_5
            mat(0, 2) = -I;
            mat(2, 0) = I;
            break;
        case 6:  // λ_6
            mat(1, 2) = 1.0;
            mat(2, 1) = 1.0;
            break;
        case 7:  // λ_7
            mat(1, 2) = -I;
            mat(2, 1) = I;
            break;
        case 8:  // λ_8
            mat(0, 0) = 1.0 / sqrt3;
            mat(1, 1) = 1.0 / sqrt3;
            mat(2, 2) = -2.0 / sqrt3;
            break;
        default:
            throw std::invalid_argument("Gell-Mann index must be 1-8");
    }
    return mat;
}

std::array<Matrix3c, 8> GellMann::allLambda() {
    std::array<Matrix3c, 8> matrices;
    for (int a = 1; a <= 8; ++a) {
        matrices[a - 1] = lambda(a);
    }
    return matrices;
}

ComplexMatrix GellMann::lambdaDotLambda() {
    // Compute Σ_a (λ_a ⊗ λ_a) in the 3⊗3 = 9 dimensional space
    ComplexMatrix result = ComplexMatrix::Zero(9, 9);
    auto lambdas = allLambda();

    for (int a = 0; a < 8; ++a) {
        // Kronecker product λ_a ⊗ λ_a
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    for (int l = 0; l < 3; ++l) {
                        int row = i * 3 + k;
                        int col = j * 3 + l;
                        result(row, col) += lambdas[a](i, j) * lambdas[a](k, l);
                    }
                }
            }
        }
    }
    return result;
}

Real GellMann::colorFactor(const ComplexVector& color_wf) {
    ComplexMatrix lambdaDot = lambdaDotLambda();
    Complex factor = color_wf.adjoint() * lambdaDot * color_wf;
    return factor.real() / 4.0;  // Divide by 4 for standard normalization
}

ComplexVector GellMann::mesonColorSinglet() {
    // |1> = (1/√3)(|rr̄> + |gḡ> + |bb̄>)
    // In the 3⊗3̄ space, this is indexed as |i,j> where j is the antiquark color
    // For meson: antiquark carries anti-color, so |r,r̄> = |0,0>, |g,ḡ> = |1,1>, |b,b̄> = |2,2>
    ComplexVector wf = ComplexVector::Zero(9);
    Real norm = 1.0 / std::sqrt(3.0);
    wf(0) = norm;  // |r,r̄> -> index 0*3+0 = 0
    wf(4) = norm;  // |g,ḡ> -> index 1*3+1 = 4
    wf(8) = norm;  // |b,b̄> -> index 2*3+2 = 8
    return wf;
}

ComplexVector GellMann::baryonColorSinglet() {
    // |1> = (1/√6) ε_{ijk} |i,j,k>
    // This is a 27-dimensional vector for 3⊗3⊗3
    // Indices: |i,j,k> -> i*9 + j*3 + k
    ComplexVector wf = ComplexVector::Zero(27);
    Real norm = 1.0 / std::sqrt(6.0);

    // ε_{012} = +1: |r,g,b>
    wf(0 * 9 + 1 * 3 + 2) = norm;
    // ε_{120} = +1: |g,b,r>
    wf(1 * 9 + 2 * 3 + 0) = norm;
    // ε_{201} = +1: |b,r,g>
    wf(2 * 9 + 0 * 3 + 1) = norm;
    // ε_{021} = -1: |r,b,g>
    wf(0 * 9 + 2 * 3 + 1) = -norm;
    // ε_{102} = -1: |g,r,b>
    wf(1 * 9 + 0 * 3 + 2) = -norm;
    // ε_{210} = -1: |b,g,r>
    wf(2 * 9 + 1 * 3 + 0) = -norm;

    return wf;
}

void GellMann::initializeTables() {
    if (tables_initialized_) return;

    // Initialize to zero
    for (auto& plane : f_abc_) {
        for (auto& row : plane) {
            row.fill(0.0);
        }
    }
    for (auto& plane : d_abc_) {
        for (auto& row : plane) {
            row.fill(0.0);
        }
    }

    const Real sqrt3 = std::sqrt(3.0);

    // Non-zero f_abc (totally antisymmetric)
    // f_{123} = 1
    f_abc_[0][1][2] = 1.0;
    // f_{147} = 1/2
    f_abc_[0][3][6] = 0.5;
    // f_{156} = -1/2
    f_abc_[0][4][5] = -0.5;
    // f_{246} = 1/2
    f_abc_[1][3][5] = 0.5;
    // f_{257} = 1/2
    f_abc_[1][4][6] = 0.5;
    // f_{345} = 1/2
    f_abc_[2][3][4] = 0.5;
    // f_{367} = -1/2
    f_abc_[2][5][6] = -0.5;
    // f_{458} = √3/2
    f_abc_[3][4][7] = sqrt3 / 2.0;
    // f_{678} = √3/2
    f_abc_[5][6][7] = sqrt3 / 2.0;

    // Fill by antisymmetry
    for (int a = 0; a < 8; ++a) {
        for (int b = 0; b < 8; ++b) {
            for (int c = 0; c < 8; ++c) {
                if (a < b && b < c && f_abc_[a][b][c] != 0.0) {
                    Real val = f_abc_[a][b][c];
                    // Even permutations: same sign
                    f_abc_[b][c][a] = val;
                    f_abc_[c][a][b] = val;
                    // Odd permutations: opposite sign
                    f_abc_[a][c][b] = -val;
                    f_abc_[b][a][c] = -val;
                    f_abc_[c][b][a] = -val;
                }
            }
        }
    }

    // Non-zero d_abc (totally symmetric)
    // d_{118} = d_{228} = d_{338} = 1/√3
    d_abc_[0][0][7] = 1.0 / sqrt3;
    d_abc_[1][1][7] = 1.0 / sqrt3;
    d_abc_[2][2][7] = 1.0 / sqrt3;
    // d_{448} = d_{558} = d_{668} = d_{778} = -1/(2√3)
    d_abc_[3][3][7] = -1.0 / (2.0 * sqrt3);
    d_abc_[4][4][7] = -1.0 / (2.0 * sqrt3);
    d_abc_[5][5][7] = -1.0 / (2.0 * sqrt3);
    d_abc_[6][6][7] = -1.0 / (2.0 * sqrt3);
    // d_{888} = -1/√3
    d_abc_[7][7][7] = -1.0 / sqrt3;
    // d_{146} = d_{157} = 1/2
    d_abc_[0][3][5] = 0.5;
    d_abc_[0][4][6] = 0.5;
    // d_{247} = d_{256} = 1/2
    d_abc_[1][3][6] = 0.5;
    d_abc_[1][4][5] = 0.5;
    // d_{344} = d_{355} = 1/2
    d_abc_[2][3][3] = 0.5;
    d_abc_[2][4][4] = 0.5;
    // d_{366} = d_{377} = -1/2
    d_abc_[2][5][5] = -0.5;
    d_abc_[2][6][6] = -0.5;

    // Fill by symmetry
    for (int a = 0; a < 8; ++a) {
        for (int b = 0; b < 8; ++b) {
            for (int c = 0; c < 8; ++c) {
                if (d_abc_[a][b][c] != 0.0) {
                    Real val = d_abc_[a][b][c];
                    d_abc_[a][c][b] = val;
                    d_abc_[b][a][c] = val;
                    d_abc_[b][c][a] = val;
                    d_abc_[c][a][b] = val;
                    d_abc_[c][b][a] = val;
                }
            }
        }
    }

    tables_initialized_ = true;
}

Real GellMann::structureConstant_f(int a, int b, int c) {
    initializeTables();
    if (a < 1 || a > 8 || b < 1 || b > 8 || c < 1 || c > 8) {
        return 0.0;
    }
    return f_abc_[a - 1][b - 1][c - 1];
}

Real GellMann::structureConstant_d(int a, int b, int c) {
    initializeTables();
    if (a < 1 || a > 8 || b < 1 || b > 8 || c < 1 || c > 8) {
        return 0.0;
    }
    return d_abc_[a - 1][b - 1][c - 1];
}

} // namespace gem
