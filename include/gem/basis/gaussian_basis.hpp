#ifndef GEM_BASIS_GAUSSIAN_BASIS_HPP
#define GEM_BASIS_GAUSSIAN_BASIS_HPP

#include "../core/types.hpp"
#include "../core/constants.hpp"
#include <vector>
#include <cmath>

namespace gem {

/**
 * @brief Gaussian basis function for the Gaussian Expansion Method
 *
 * The basis function is: φ_n(r) = N_n * exp(-ν_n * r^2)
 * where N_n is the normalization constant and ν_n is the width parameter.
 *
 * For l=0 (s-wave): N_n = (2ν_n/π)^(3/4)
 */
class GaussianBasis {
public:
    /**
     * @brief Construct a Gaussian basis set using geometric series
     * @param n_basis Number of basis functions
     * @param nu_min Minimum width parameter (fm^-2)
     * @param nu_max Maximum width parameter (fm^-2)
     */
    GaussianBasis(int n_basis, Real nu_min, Real nu_max);

    /**
     * @brief Construct from BasisParameters
     */
    explicit GaussianBasis(const BasisParameters& params);

    /**
     * @brief Get the number of basis functions
     */
    int size() const { return n_basis_; }

    /**
     * @brief Get the width parameter for basis function n
     */
    Real nu(int n) const { return nu_[n]; }

    /**
     * @brief Get the normalization constant for basis function n
     * For l=0: N_n = (2ν_n/π)^(3/4)
     */
    Real normalization(int n) const { return norm_[n]; }

    /**
     * @brief Evaluate basis function n at distance r
     */
    Real evaluate(int n, Real r) const;

    /**
     * @brief Overlap integral <φ_n|φ_m>
     * = N_n * N_m * (π / (ν_n + ν_m))^(3/2)
     */
    Real overlap(int n, int m) const;

    /**
     * @brief Kinetic energy integral <φ_n|T|φ_m> where T = -ℏ²/(2μ) ∇²
     * = (ℏ²/2μ) * 6 * ν_n * ν_m / (ν_n + ν_m) * <φ_n|φ_m>
     * @param mu Reduced mass in GeV
     */
    Real kinetic(int n, int m, Real mu) const;

    /**
     * @brief Coulomb potential integral <φ_n|1/r|φ_m>
     * = N_n * N_m * 2π / (ν_n + ν_m)
     */
    Real coulomb(int n, int m) const;

    /**
     * @brief Linear potential integral <φ_n|r|φ_m>
     * = N_n * N_m * 2π / (ν_n + ν_m)²
     */
    Real linear(int n, int m) const;

    /**
     * @brief Gaussian potential integral <φ_n|exp(-α*r²)|φ_m>
     * = N_n * N_m * (π / (ν_n + ν_m + α))^(3/2)
     */
    Real gaussian_potential(int n, int m, Real alpha) const;

    /**
     * @brief Delta function (smeared) integral for hyperfine interaction
     * Using Gaussian smearing: δ_σ(r) = (α/π)^(3/2) exp(-α*r²) where α = 1/σ²
     */
    Real delta_smeared(int n, int m, Real sigma) const;

    /**
     * @brief Build the overlap matrix S_nm = <φ_n|φ_m>
     */
    Matrix overlapMatrix() const;

    /**
     * @brief Build the kinetic energy matrix T_nm
     */
    Matrix kineticMatrix(Real mu) const;

    /**
     * @brief Build the Coulomb potential matrix V_nm = <φ_n|1/r|φ_m>
     */
    Matrix coulombMatrix() const;

    /**
     * @brief Build the linear potential matrix V_nm = <φ_n|r|φ_m>
     */
    Matrix linearMatrix() const;

    /**
     * @brief Build the smeared delta function matrix
     */
    Matrix deltaMatrix(Real sigma) const;

private:
    int n_basis_;
    std::vector<Real> nu_;    // Width parameters
    std::vector<Real> norm_;  // Normalization constants

    void initializeGeometricSeries(Real nu_min, Real nu_max);
};

} // namespace gem

#endif // GEM_BASIS_GAUSSIAN_BASIS_HPP
