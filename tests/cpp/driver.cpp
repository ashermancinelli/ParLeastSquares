/*
 * Asher Mancinelli
 *
 * Pass arbitrary number of directories to be searched for matrix files.
 *
 */

#include <string>
#include <vector>

#include <Eigen/Core>
#include <unsupported/Eigen/SparseExtra>

#include <ParLeastSquares>
#include <TestDriver.hpp>

int main(int argc, char** argv)
{
  std::vector<std::string> datadirs;
  if (argc > 1)
  {
    for (int i=1; i<argc; i++)
    {
      datadirs.push_back(argv[i]);
    }
  }

  TestDriver driver(datadirs);

  driver.test_all();
  driver.print_summary();
}
