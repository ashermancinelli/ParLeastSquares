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
    for (int i=1; i<argc; i++)
    {
      search_dirs.push_back(argv[i]);
    }
  }

  auto dir = found_acceptable_dir(search_dirs);

  if (dir)
  {
    datadir = *tmp;
  }
  else
  {
    std::cout << "Could not find valid dir with all needed matrices.\n";
    return 1;
  }

  MatrixXd S_mat = read_mm(datadir + "/S_mat.mtx");
  MatrixXd R_mat = read_mm(datadir + "/R_mat.mtx");
  MatrixXd P_mat = read_mm(datadir + "/P_mat.mtx");
  VectorXd Keq_constant = read_vector(datadir + "/Keq_constant.txt");
  VectorXd state = read_vector(datadir + "/state.txt");
  VectorXd f_log_counts = read_vector(datadir + "/f_log_counts.txt");
  VectorXd v_log_counts = read_vector(datadir + "/v_log_counts.txt");

  VectorXd results = read_vector(datadir + "/results.txt");

  Eigen::VectorXd result = least_squares(
      S_mat,
      R_mat,
      P_mat, 
      Keq_constant,
      state,
      f_log_counts,
      v_log_counts);

  std::cout << "--- Successful least squares call.\n"
    << "--- Result: ";
  // << result;
}
