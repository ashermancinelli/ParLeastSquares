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

void usage()
{
  std::cout << "Usage:\n"
    "./binary <maxfevs> <xtol> <numerical|analytical> <arbitrary number of paths to data directories>\n";
}

int main(int argc, char** argv)
{
  std::vector<std::string> datadirs;
  int maxfevs;
  double xtol;
  DiffType dt;

  switch (argc)
  {
    case 1:
      maxfevs = 1000;
      xtol = 1e-10;
      dt = DiffType::Numerical;
      break;
    case 2:
    case 3:
      usage();
      return 1;
    default:
      maxfevs = atoi(argv[1]);
      xtol = atof(argv[2]);

      if (strcmp(argv[3], "numerical") == 0)
      {
        dt = DiffType::Numerical;
      }
      else if (strcmp(argv[3], "analytical") == 0)
      {
        dt = DiffType::Analytical;
      }
      else
      {
        usage();
        return 1;
      }

      if (argc > 4)
      {
        for (int i=4; i<argc; i++)
        {
          datadirs.push_back(argv[i]);
        }
      }
  }


  TestDriver driver(datadirs);
  driver.test_minimize(dt, maxfevs, xtol);
}
