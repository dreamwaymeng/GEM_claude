#ifndef GEM_QUANTUM_GELL_MANN_HPP
#define GEM_QUANTUM_GELL_MANN_HPP

#include "../core/types.hpp"
#include <array>
#include <cmath>

namespace gem {

/**
 * @brief SU(3) Gell-Mann matrices and color algebra
 *
 * The 8 Gell-Mann matrices λ_a (a=1,...,8) are the generators of SU(3).
 * They satisfy: [λ_a, λ_b] = 2i f_{abc} λ_c
 *              {λ_a, λ_b} = (4/3) δ_{ab} + 2 d_{abc} λ_c
 *              Tr(λ_a λ_b) = 2 δ_{ab}
 */
class GellMann {
public:
    /**
     * @brief Get Gell-Mann matrix λ_a (a = 1 to 8)
     */
    static Matrix3c lambda(int a);

    /**
     * @brief Get all 8 Gell-Mann matrices
     */
    static std::array<Matrix3c, 8> allLambda();

    /**
     * @brief Compute λ_i · λ_j = Σ_a λ_a^(i) λ_a^(j)
     * For quarks i and j in the fundamental representation
     * Returns a 9x9 matrix in the product space (color_i ⊗ color_j)
     */
    static ComplexMatrix lambdaDotLambda();

    /**
     * @brief Color factor for quark-antiquark in singlet state
     * <1| λ_i · λ_j / 4 |1> = -4/3 for meson
     */
    static constexpr Real mesonColorFactor() { return -4.0 / 3.0; }

    /**
     * @brief Color factor for two quarks in antisymmetric 3-bar state
     * <3-bar| λ_i · λ_j / 4 |3-bar> = -2/3 for baryon
     */
    static constexpr Real baryonColorFactor() { return -2.0 / 3.0; }

    /**
     * @brief Color factor for two quarks in symmetric 6 state
     * <6| λ_i · λ_j / 4 |6> = +1/3
     */
    static constexpr Real sextetColorFactor() { return 1.0 / 3.0; }

    /**
     * @brief Compute the color factor <λ_i · λ_j> / 4 for given color wave function
     * @param color_wf Color wave function coefficients (9 components for 2 quarks)
     */
    static Real colorFactor(const ComplexVector& color_wf);

    /**
     * @brief Generate the color singlet wave function for q-qbar (meson)
     * |1> = (1/√3) Σ_i |i, ī> = (1/√3)(|rr̄> + |gḡ> + |bb̄>)
     */
    static ComplexVector mesonColorSinglet();

    /**
     * @brief Generate the color singlet wave function for qqq (baryon)
     * |1> = (1/√6) ε_{ijk} |i,j,k> (antisymmetric)
     */
    static ComplexVector baryonColorSinglet();

    /**
     * @brief Structure constants f_abc (totally antisymmetric)
     */
    static Real structureConstant_f(int a, int b, int c);

    /**
     * @brief Symmetric structure constants d_abc
     */
    static Real structureConstant_d(int a, int b, int c);

private:
    static void initializeTables();
    static bool tables_initialized_;
    static std::array<std::array<std::array<Real, 8>, 8>, 8> f_abc_;
    static std::array<std::array<std::array<Real, 8>, 8>, 8> d_abc_;
};

} // namespace gem

#endif // GEM_QUANTUM_GELL_MANN_HPP
