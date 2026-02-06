#ifndef GEM_SOLVER_EIGENVALUE_SOLVER_HPP
#define GEM_SOLVER_EIGENVALUE_SOLVER_HPP

#include "../core/types.hpp"
#include "../matrix/hamiltonian.hpp"

namespace gem {

/**
 * @brief Eigenvalue solver for the generalized eigenvalue problem
 *
 * Solves: H c = E S c
 * where H is the Hamiltonian matrix, S is the overlap matrix,
 * and E are the eigenvalues (energies).
 */
class EigenvalueSolver {
public:
    /**
     * @brief Solver methods
     */
    enum class Method {
        GeneralizedEigen,  // Eigen's generalized eigenvalue solver
        CholeskyReduce,    // Reduce to standard form via Cholesky decomposition
        SvdStabilize       // SVD stabilization for near-singular overlap
    };

    EigenvalueSolver(Method method = Method::CholeskyReduce);

    /**
     * @brief Solve the generalized eigenvalue problem H c = E S c
     * @param H Hamiltonian matrix
     * @param S Overlap matrix
     * @return EigenResult containing eigenvalues and eigenvectors
     */
    EigenResult solve(const Matrix& H, const Matrix& S);

    /**
     * @brief Solve using the Hamiltonian object
     */
    EigenResult solve(const Hamiltonian& hamiltonian);

    /**
     * @brief Get the ground state energy
     */
    Real groundStateEnergy() const { return eigenvalues_(0); }

    /**
     * @brief Get the ground state wave function coefficients
     */
    Vector groundStateCoeffs() const { return eigenvectors_.col(0); }

    /**
     * @brief Get eigenvalue at index i
     */
    Real eigenvalue(int i) const { return eigenvalues_(i); }

    /**
     * @brief Get eigenvector at index i
     */
    Vector eigenvector(int i) const { return eigenvectors_.col(i); }

    /**
     * @brief Get all eigenvalues
     */
    const Vector& eigenvalues() const { return eigenvalues_; }

    /**
     * @brief Get all eigenvectors (columns)
     */
    const Matrix& eigenvectors() const { return eigenvectors_; }

    /**
     * @brief Set numerical tolerance for overlap matrix conditioning
     */
    void setTolerance(Real tol) { tolerance_ = tol; }

    /**
     * @brief Check if solution converged properly
     */
    bool hasConverged() const { return converged_; }

    /**
     * @brief Get number of eigenvalues found
     */
    int numEigenvalues() const { return eigenvalues_.size(); }

private:
    Method method_;
    Real tolerance_;
    bool converged_;

    Vector eigenvalues_;
    Matrix eigenvectors_;

    /**
     * @brief Solve using generalized eigenvalue solver
     */
    EigenResult solveGeneralized(const Matrix& H, const Matrix& S);

    /**
     * @brief Solve using Cholesky decomposition to reduce to standard form
     * S = L L^T, then solve (L^-1 H L^-T) y = E y, c = L^-T y
     */
    EigenResult solveCholesky(const Matrix& H, const Matrix& S);

    /**
     * @brief Solve with SVD stabilization for ill-conditioned overlap
     */
    EigenResult solveSvdStabilized(const Matrix& H, const Matrix& S);
};

/**
 * @brief Result container with additional analysis
 */
struct SolverResult {
    EigenResult eigen;
    Real ground_state_energy;
    Vector ground_state_coeffs;

    // Analysis results
    Real normalization;          // Should be 1
    Real energy_expectation;     // <ψ|H|ψ>/<ψ|S|ψ> should equal eigenvalue
    Real overlap_condition;      // Condition number of S

    SolverResult() : ground_state_energy(0.0), normalization(0.0),
                     energy_expectation(0.0), overlap_condition(0.0) {}
};

/**
 * @brief Analyze the solution quality
 */
SolverResult analyzeResult(const EigenResult& result, const Matrix& H, const Matrix& S);

} // namespace gem

#endif // GEM_SOLVER_EIGENVALUE_SOLVER_HPP
