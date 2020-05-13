#pragma once
#include <string>
#include <iomanip>

[[nodiscard]]
Eigen::MatrixXd read_mm(const std::string& path);

bool write_mm(const Eigen::MatrixXd& mtx, const std::string& path);

[[nodiscard]]
Eigen::VectorXd read_vector(const std::string& path);

bool write_vector(const Eigen::VectorXd& vec, const std::string& path);
