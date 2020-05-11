#pragma once

[[nodiscard]]
Eigen::MatrixXd read_mm(const std::string& path);

[[nodiscard]]
Eigen::VectorXd read_vector(const std::string& path);
