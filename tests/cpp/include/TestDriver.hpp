#pragma once
#include <optional>
#include <vector>
#include <string>
#include <string>
#include <iostream>
#include <fstream>
#include <cassert>

#include <Eigen/Core>
#include <ParLeastSquares>
#include <test_utilities.hpp>

using Eigen::MatrixXd;
using Eigen::VectorXd;

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
