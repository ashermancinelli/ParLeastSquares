#include <string>
#include <vector>

#include <Eigen/Core>
#include <unsupported/Eigen/SparseExtra>

#include <ParLeastSquares>
#include <test_utils.hpp>
extern std::vector<std::string> search_dirs;

using Eigen::VectorXd;
using Eigen::MatrixXd;

int main(int argc, char** argv)
{
  std::string datadir;
  if (argc > 1)
  {
    datadir = argv[1];
  }
  else
  {
    datadir = "../tests/test_cases/matrix_market/";
  }

  MatrixXd S_mat          = read_mm(datadir + "/S_mat.mtx");
  MatrixXd R_mat          = read_mm(datadir + "/R_mat.mtx");
  MatrixXd P_mat          = read_mm(datadir + "/P_mat.mtx");

  VectorXd Keq_constant   = read_vector(datadir + "/Keq_constant.txt");
  VectorXd state          = read_vector(datadir + "/state.txt");
  VectorXd f_log_counts   = read_vector(datadir + "/f_log_counts.txt");
  VectorXd v_log_counts   = read_vector(datadir + "/v_log_counts.txt");

  VectorXd results        = read_vector(datadir + "/results.txt");

  auto lm = lmder_functor(
      S_mat, R_mat, P_mat,
      Keq_constant, state, f_log_counts);

  VectorXd deriv(S_mat.rows());

  lm(v_log_counts, deriv);

  std::cout << "Got deriv values of:\n"
    << deriv
    << "\n\n";

  return 0;
}
