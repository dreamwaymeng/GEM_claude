#include "gem/solver/eigenvalue_solver.hpp"
#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>
#include <Eigen/SVD>
#include <algorithm>
#include <cmath>

namespace gem {

EigenvalueSolver::EigenvalueSolver(Method method)
    : method_(method), tolerance_(1e-12), converged_(false) {}

EigenResult EigenvalueSolver::solve(const Matrix& H, const Matrix& S) {
    switch (method_) {
        case Method::GeneralizedEigen:
            return solveGeneralized(H, S);
        case Method::CholeskyReduce:
            return solveCholesky(H, S);
        case Method::SvdStabilize:
            return solveSvdStabilized(H, S);
        default:
            return solveCholesky(H, S);
    }
}

EigenResult EigenvalueSolver::solve(const Hamiltonian& hamiltonian) {
    return solve(hamiltonian.matrix(), hamiltonian.overlapMatrix());
}

EigenResult EigenvalueSolver::solveGeneralized(const Matrix& H, const Matrix& S) {
    EigenResult result;
    int n = H.rows();

    // Use Eigen's generalized eigenvalue solver
    Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> solver(H, S);

    if (solver.info() != Eigen::Success) {
        result.converged = false;
        return result;
    }

    eigenvalues_ = solver.eigenvalues();
    eigenvectors_ = solver.eigenvectors();

    // Sort by eigenvalue (should already be sorted, but verify)
    std::vector<int> indices(n);
    for (int i = 0; i < n; ++i) indices[i] = i;
    std::sort(indices.begin(), indices.end(),
              [this](int a, int b) { return eigenvalues_(a) < eigenvalues_(b); });

    Vector sorted_evals(n);
    Matrix sorted_evecs(n, n);
    for (int i = 0; i < n; ++i) {
        sorted_evals(i) = eigenvalues_(indices[i]);
        sorted_evecs.col(i) = eigenvectors_.col(indices[i]);
    }
    eigenvalues_ = sorted_evals;
    eigenvectors_ = sorted_evecs;

    result.eigenvalues = eigenvalues_;
    result.eigenvectors = eigenvectors_;
    result.converged = true;
    converged_ = true;

    return result;
}

EigenResult EigenvalueSolver::solveCholesky(const Matrix& H, const Matrix& S) {
    EigenResult result;
    int n = H.rows();

    // Cholesky decomposition: S = L L^T
    Eigen::LLT<Matrix> llt(S);
    if (llt.info() != Eigen::Success) {
        // Fall back to SVD method if Cholesky fails
        return solveSvdStabilized(H, S);
    }

    Matrix L = llt.matrixL();
    Matrix L_inv = L.inverse();
    Matrix L_inv_T = L_inv.transpose();

    // Transform to standard eigenvalue problem: H' y = E y
    // where H' = L^-1 H L^-T
    Matrix H_prime = L_inv * H * L_inv_T;

    // Solve standard eigenvalue problem
    Eigen::SelfAdjointEigenSolver<Matrix> solver(H_prime);

    if (solver.info() != Eigen::Success) {
        result.converged = false;
        return result;
    }

    eigenvalues_ = solver.eigenvalues();
    // Transform eigenvectors back: c = L^-T y
    eigenvectors_ = L_inv_T * solver.eigenvectors();

    // Normalize eigenvectors with respect to S
    for (int i = 0; i < n; ++i) {
        Vector c = eigenvectors_.col(i);
        Real norm = std::sqrt(c.transpose() * S * c);
        if (norm > tolerance_) {
            eigenvectors_.col(i) /= norm;
        }
    }

    result.eigenvalues = eigenvalues_;
    result.eigenvectors = eigenvectors_;
    result.converged = true;
    converged_ = true;

    return result;
}

EigenResult EigenvalueSolver::solveSvdStabilized(const Matrix& H, const Matrix& S) {
    EigenResult result;
    int n = H.rows();

    // Eigenvalue decomposition of overlap matrix: S = U Λ U^T
    // (S is symmetric positive semi-definite)
    Eigen::SelfAdjointEigenSolver<Matrix> s_solver(S);
    if (s_solver.info() != Eigen::Success) {
        result.converged = false;
        return result;
    }

    Vector lambda = s_solver.eigenvalues();   // eigenvalues in ascending order
    Matrix U = s_solver.eigenvectors();

    // Find effective rank: count eigenvalues above threshold
    // Use relative tolerance based on largest eigenvalue
    Real max_lambda = lambda(n - 1);  // largest eigenvalue (ascending order)
    Real threshold = tolerance_ * max_lambda;

    int rank = 0;
    for (int i = 0; i < n; ++i) {
        if (lambda(i) > threshold) {
            ++rank;
        }
    }

    if (rank == 0) {
        result.converged = false;
        return result;
    }

    // Build transformation matrix X that projects into well-conditioned subspace
    // X = U_r Λ_r^(-1/2), where U_r contains only eigenvectors with λ > threshold
    // Then S^(-1/2) ≈ X X^T in the reduced space
    Matrix X(n, rank);
    int col = 0;
    for (int i = 0; i < n; ++i) {
        if (lambda(i) > threshold) {
            X.col(col) = U.col(i) / std::sqrt(lambda(i));
            ++col;
        }
    }

    // Transform Hamiltonian to reduced space: H_red = X^T H X (rank × rank matrix)
    Matrix H_reduced = X.transpose() * H * X;

    // Solve standard eigenvalue problem in reduced space
    Eigen::SelfAdjointEigenSolver<Matrix> h_solver(H_reduced);
    if (h_solver.info() != Eigen::Success) {
        result.converged = false;
        return result;
    }

    // Get eigenvalues (these are the physical energies)
    Vector evals_reduced = h_solver.eigenvalues();
    Matrix evecs_reduced = h_solver.eigenvectors();

    // Transform eigenvectors back to original space: c = X y
    Matrix evecs_full = X * evecs_reduced;

    // Store results (only rank eigenvalues are meaningful)
    eigenvalues_ = evals_reduced;
    eigenvectors_ = evecs_full;

    // Normalize eigenvectors with respect to S
    for (int i = 0; i < rank; ++i) {
        Vector c = eigenvectors_.col(i);
        Real norm = std::sqrt(c.transpose() * S * c);
        if (norm > 1e-15) {
            eigenvectors_.col(i) /= norm;
        }
    }

    result.eigenvalues = eigenvalues_;
    result.eigenvectors = eigenvectors_;
    result.converged = true;
    converged_ = true;

    return result;
}

SolverResult analyzeResult(const EigenResult& result, const Matrix& H, const Matrix& S) {
    SolverResult analysis;
    analysis.eigen = result;

    if (!result.converged || result.eigenvalues.size() == 0) {
        return analysis;
    }

    analysis.ground_state_energy = result.eigenvalues(0);
    analysis.ground_state_coeffs = result.eigenvectors.col(0);

    // Check normalization: <ψ|S|ψ> should be 1
    Vector c = analysis.ground_state_coeffs;
    analysis.normalization = c.transpose() * S * c;

    // Check energy expectation: <ψ|H|ψ>/<ψ|S|ψ> should equal eigenvalue
    Real H_expect = c.transpose() * H * c;
    analysis.energy_expectation = H_expect / analysis.normalization;

    // Condition number of S
    Eigen::JacobiSVD<Matrix> svd(S);
    Vector sigma = svd.singularValues();
    if (sigma(sigma.size() - 1) > 1e-15) {
        analysis.overlap_condition = sigma(0) / sigma(sigma.size() - 1);
    } else {
        analysis.overlap_condition = std::numeric_limits<Real>::infinity();
    }

    return analysis;
}

} // namespace gem
