#ifndef GEM_QUANTUM_WAVE_FUNCTION_HPP
#define GEM_QUANTUM_WAVE_FUNCTION_HPP

#include "../core/types.hpp"
#include "../basis/gaussian_basis.hpp"
#include "gell_mann.hpp"
#include "spin_operators.hpp"
#include <vector>

namespace gem {

/**
 * @brief Discrete wave function components (color × spin × flavor)
 *
 * The full hadron wave function can be written as:
 * Ψ = ψ_space(r) × χ_color × χ_spin × χ_flavor
 *
 * This class handles the discrete parts (color, spin, flavor).
 */
class DiscreteWaveFunction {
public:
    /**
     * @brief Create discrete wave function with given color and spin states
     */
    DiscreteWaveFunction(const ComplexVector& color, const ComplexVector& spin);

    /**
     * @brief Get color wave function
     */
    const ComplexVector& color() const { return color_; }

    /**
     * @brief Get spin wave function
     */
    const ComplexVector& spin() const { return spin_; }

    /**
     * @brief Compute color factor <λ_i · λ_j> / 4
     */
    Real colorFactor() const;

    /**
     * @brief Compute spin factor <s_i · s_j>
     */
    Real spinFactor() const;

    /**
     * @brief Compute combined color-spin factor
     */
    Real colorSpinFactor() const { return colorFactor() * spinFactor(); }

    /**
     * @brief Check if wave function is normalized
     */
    bool isNormalized(Real tol = 1e-10) const;

    /**
     * @brief Normalize the wave function
     */
    void normalize();

private:
    ComplexVector color_;
    ComplexVector spin_;
};

/**
 * @brief Factory class for creating standard hadron wave functions
 */
class WaveFunctionFactory {
public:
    /**
     * @brief Create meson (qq̄) wave function
     * @param total_spin Total spin (0 or 1)
     * @param Sz z-component of total spin
     */
    static DiscreteWaveFunction meson(int total_spin, int Sz = 0);

    /**
     * @brief Create baryon (qqq) wave function
     * @param total_spin Total spin (1 or 3, meaning 1/2 or 3/2)
     * @param Sz z-component (times 2)
     */
    static DiscreteWaveFunction baryon(int twice_total_spin, int twice_Sz = 1);

    /**
     * @brief Create color singlet wave function for given number of quarks
     */
    static ComplexVector colorSinglet(int n_quarks);

    /**
     * @brief Create symmetric spin state
     */
    static ComplexVector symmetricSpinState(int n_quarks, int twice_S);

    /**
     * @brief Create mixed symmetry spin state (for baryons)
     */
    static ComplexVector mixedSymmetrySpinState(int twice_S, int twice_Sz);
};

/**
 * @brief Full spatial-discrete wave function
 */
class FullWaveFunction {
public:
    FullWaveFunction(const GaussianBasis& basis, const Vector& coefficients,
                     const DiscreteWaveFunction& discrete);

    /**
     * @brief Get spatial part coefficients
     */
    const Vector& spatialCoefficients() const { return coefficients_; }

    /**
     * @brief Get discrete wave function
     */
    const DiscreteWaveFunction& discrete() const { return discrete_; }

    /**
     * @brief Get underlying basis
     */
    const GaussianBasis& basis() const { return basis_; }

    /**
     * @brief Evaluate spatial wave function at point r
     */
    Real evaluateSpatial(Real r) const;

    /**
     * @brief Compute <r> expectation value
     */
    Real meanRadius() const;

    /**
     * @brief Compute <r²> expectation value
     */
    Real meanRadiusSquared() const;

    /**
     * @brief Compute RMS radius √<r²>
     */
    Real rmsRadius() const;

private:
    const GaussianBasis& basis_;
    Vector coefficients_;
    DiscreteWaveFunction discrete_;
};

} // namespace gem

#endif // GEM_QUANTUM_WAVE_FUNCTION_HPP
