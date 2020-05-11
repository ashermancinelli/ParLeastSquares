#pragma once

/*
 * Driver function for calculating the least squares with the
 * Levenberg-Marquardt method
 */
[[nodiscard]]
extern Eigen::VectorXd least_squares(
    Eigen::MatrixXd& S_mat,
    Eigen::MatrixXd& R_back_mat,
    Eigen::MatrixXd& P_mat, 
    Eigen::VectorXd& Keq_constant,
    Eigen::VectorXd& E_Regulation,
    Eigen::VectorXd& log_fcounts,
    Eigen::VectorXd& log_vcounts,
    DiffType dt=DiffType::Analytical,
    const int maxfevs=2e4,
    const double xtol=1e-15);
