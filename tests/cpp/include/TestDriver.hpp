#pragma once
#include <optional>
#include <vector>
#include <string>
#include <string>
#include <iostream>
#include <fstream>
#include <cassert>

#include <ParLeastSquares.hpp>
#include <Eigen/Core>

using Eigen::MatrixXd;
using Eigen::VectorXd;

static std::vector<std::string> default_search_dirs =
{
  "../tests/data/tca_glycolosis",
  "../tests/data/tca_glycolosis_gogat",
  "/qfs/projects/boltzmann/data/matrix_market/tca_glycolosis/",
  "/qfs/projects/boltzmann/data/matrix_market/tca_glycolosis_gogat//",
};

static const std::vector<std::string> default_filenames =
{
  "/S_mat.mtx",
  "/R_mat.mtx",
  "/P_mat.mtx",
  "/Keq_constant.txt",
  "/state.txt",
  "/f_log_counts.txt",
  "/v_log_counts.txt",
};

/*
 * - - - - - - - - - - - - - - - - - - - - - -
 * Utility functions
 * - - - - - - - - - - - - - - - - - - - - - -
 */
inline bool file_exists(const std::string& path);

bool all_files_exist(
    const std::string& dir,
    const std::vector<std::string>& filenames=default_filenames);

std::optional<std::string> found_acceptable_dir(
    const std::vector<std::string>& dirs,
    const std::vector<std::string>& filenames=default_filenames);
/*
 * - - - - - - - - - - - - - - - - - - - - - -
 * \Utility functions
 * - - - - - - - - - - - - - - - - - - - - - -
 */

/*
 * - - - - - - - - - - - - - - - - - - - - - -
 * Class to handle setting up and running tests
 *
 * Finds all files needed for tests and keeps
 * track of failures.
 * - - - - - - - - - - - - - - - - - - - - - -
 */
struct TestDriver
{
  private:
    MatrixXd S_mat;
    MatrixXd R_mat;
    MatrixXd P_mat;
    VectorXd Keq_constant;
    VectorXd state;
    VectorXd f_log_counts;
    VectorXd v_log_counts;
    VectorXd results;
    std::ostream& os = std::cout;
    int fail = 0;

    void init(const std::string& datadir);
    void test_minimize_numerical();
    void test_minimize_analytical();

  public:
    /*
     * It is the responsibility of the constructors to
     * also initialize the matrices and vectors needed
     * for each test.
     */
    TestDriver();
    TestDriver(const std::string&);
    TestDriver(std::vector<std::string>);
    void print_summary();
    void test_interface();
    void test_call();
    void test_minimize(DiffType);
    void test_all();
};
