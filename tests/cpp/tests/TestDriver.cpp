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
  Keq_constant = read_vector(datadir + "/Keq_constant.txt");
  state = read_vector(datadir + "/state.txt");
  f_log_counts = read_vector(datadir + "/f_log_counts.txt");
  v_log_counts = read_vector(datadir + "/v_log_counts.txt");
  results = read_vector(datadir + "/results.txt");
}

void TestDriver::test_interface()
{
}

void TestDriver::test_call()
{
}

void TestDriver::test_minimize_numerical()
{
  VectorXd result = least_squares(
      S_mat,
      R_mat,
      P_mat, 
      Keq_constant,
      state,
      f_log_counts,
      v_log_counts,
      DiffType::Numerical);
  os << "Called least_squares. Got result:\n" << result << "\n";
}

void TestDriver::test_minimize_analytical()
{
  VectorXd result = least_squares(
      S_mat,
      R_mat,
      P_mat, 
      Keq_constant,
      state,
      f_log_counts,
      v_log_counts,
      DiffType::Analytical);
  os << "Called least_squares. Got result:\n" << result << "\n";
}

void TestDriver::test_minimize(DiffType dt)
{
  switch (dt)
  {
    case DiffType::Analytical:
      test_minimize_analytical();
      return;
    case DiffType::Numerical:
      test_minimize_numerical();
      return;
  }
}

void TestDriver::test_all()
{
  test_interface();
  test_call();
  test_minimize(DiffType::Analytical);
  test_minimize(DiffType::Numerical);
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

/*
 * - - - - - - - - - - - - - - - - - - - - - -
 * Utility functions
 * - - - - - - - - - - - - - - - - - - - - - -
 */
inline bool file_exists(const std::string& path)
{
  std::ifstream f(path.c_str());
  return f.good();
}

bool all_files_exist(
    const std::string& dir,
    const std::vector<std::string>& filenames)
{
  for (const auto& filename : filenames)
  {
    if (!file_exists(dir + filename))
      return false;
  }

  return true;
}

std::optional<std::string> found_acceptable_dir(
    const std::vector<std::string>& dirs,
    const std::vector<std::string>& filenames)
{
  for (const auto& dir : dirs)
  {
    std::cout << "Checking dir " << dir << "\n";
    if (all_files_exist(dir, filenames))
    {
      std::cout << "Found all files in dir " << dir << "\n";
      return dir;
    }
  }
  return std::nullopt;
}

/*
 * - - - - - - - - - - - - - - - - - - - - - -
 * \Utility functions
 * - - - - - - - - - - - - - - - - - - - - - -
 */
