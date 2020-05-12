#include <TestDriver.hpp>

TestDriver::TestDriver()
{
  const std::optional<std::string> dir =
    found_acceptable_dir(default_search_dirs);
  assert (dir && "Could not find all needed files in given dirs.");
  init(*dir);
}

TestDriver::TestDriver(std::vector<std::string> dirs)
{
  std::optional<std::string> dir = found_acceptable_dir(dirs);
  if (!dir)
  {
    dir = found_acceptable_dir(default_search_dirs);
  }
  assert (dir && "Could not find all needed files in given dirs.");
  init(*dir);
}

TestDriver::TestDriver(const std::string& dir)
{
  assert (all_files_exist(dir) && "Could not find all needed files in given dir");
  init(dir);
}

void TestDriver::init(const std::string& datadir)
{
  S_mat = read_mm(datadir + "/S_mat.mtx");
  R_mat = read_mm(datadir + "/R_mat.mtx");
  P_mat = read_mm(datadir + "/P_mat.mtx");
  results = read_mm(datadir + "/results.mtx");
  Keq_constant = read_vector(datadir + "/Keq_constant.txt");
  state = read_vector(datadir + "/state.txt");
  f_log_counts = read_vector(datadir + "/f_log_counts.txt");
  v_log_counts = read_vector(datadir + "/v_log_counts.txt");
}

void TestDriver::test_minimize(
    DiffType dt,
    const int maxfevs,
    const double xtol)
{
  VectorXd result = least_squares(
      S_mat,
      R_mat,
      P_mat, 
      Keq_constant,
      state,
      f_log_counts,
      v_log_counts,
      dt,
      maxfevs,
      xtol);
  os << "Called least_squares with " << dt
    << " with result:\n" << result << "\n";
}

void TestDriver::test_all()
{
  test_minimize(DiffType::Analytical, 100, 1e-10);
  test_minimize(DiffType::Numerical, 100, 1e-10);
}

void TestDriver::print_summary()
{
  if (fail)
  {
    os << "Tests failed.\n";
  }
  else
  {
    os << "Tests passed.\n";
  }
}
