#ifndef GEM_QUANTUM_SPIN_OPERATORS_HPP
#define GEM_QUANTUM_SPIN_OPERATORS_HPP

#include "../core/types.hpp"
#include <array>
#include <cmath>

namespace gem {

/**
 * @brief SU(2) Pauli matrices and spin algebra
 *
 * Pauli matrices σ_i (i=1,2,3 or x,y,z):
 * σ_1 = |0 1|  σ_2 = |0 -i|  σ_3 = |1  0|
 *       |1 0|        |i  0|        |0 -1|
 *
 * Spin operator: s = σ/2
 * Commutation: [σ_i, σ_j] = 2i ε_{ijk} σ_k
 */
class SpinOperators {
public:
    /**
     * @brief Get Pauli matrix σ_i (i = 1, 2, 3)
     */
    static Matrix2c sigma(int i);

    /**
     * @brief Get spin-1/2 operator s_i = σ_i / 2
     */
    static Matrix2c spin(int i);

    /**
     * @brief Get all 3 Pauli matrices
     */
    static std::array<Matrix2c, 3> allSigma();

    /**
     * @brief Compute s_1 · s_2 = (σ_1 · σ_2) / 4 in 4-dimensional product space
     * Returns a 4x4 matrix for spin_1 ⊗ spin_2
     */
    static ComplexMatrix spinDotSpin();

    /**
     * @brief Spin factor for two spin-1/2 particles coupled to total spin S
     * <S| s_1 · s_2 |S> = [S(S+1) - 3/2] / 2
     * S=0: -3/4,  S=1: +1/4
     */
    static Real spinFactor(int twice_S);

    /**
     * @brief Spin factor for two particles with spins s1 and s2 coupled to S
     * <S| s_1 · s_2 |S> = [S(S+1) - s1(s1+1) - s2(s2+1)] / 2
     */
    static Real spinFactor(int twice_s1, int twice_s2, int twice_S);

    /**
     * @brief Compute the spin factor <s_1 · s_2> for given spin wave function
     */
    static Real spinFactor(const ComplexVector& spin_wf);

    /**
     * @brief Generate spin singlet state for two spin-1/2 particles (S=0)
     * |0,0> = (1/√2)(|↑↓> - |↓↑>)
     */
    static ComplexVector spinSinglet();

    /**
     * @brief Generate spin triplet states for two spin-1/2 particles (S=1)
     * |1,+1> = |↑↑>
     * |1,0>  = (1/√2)(|↑↓> + |↓↑>)
     * |1,-1> = |↓↓>
     * @param Sz z-component of spin (1, 0, or -1)
     */
    static ComplexVector spinTriplet(int Sz);

    /**
     * @brief Generate coupled spin state for two spin-1/2 particles
     * @param twice_S Total spin times 2 (0 or 2)
     * @param twice_Sz Sz times 2 (-S to +S in steps of 2)
     */
    static ComplexVector coupledSpinState(int twice_S, int twice_Sz);

    /**
     * @brief Clebsch-Gordan coefficient <j1,m1,j2,m2|J,M>
     * For spin-1/2 + spin-1/2 coupling
     */
    static Real clebschGordan_half_half(int twice_m1, int twice_m2, int twice_J, int twice_M);

private:
    static const Complex I_;  // Imaginary unit
};

} // namespace gem

#endif // GEM_QUANTUM_SPIN_OPERATORS_HPP
