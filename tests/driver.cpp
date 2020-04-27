#include <string>
#include <vector>

#include <Eigen/Core>
#include <unsupported/Eigen/SparseExtra>

#include "helper_functions.hpp"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Map;

int main(int argc, char** argv)
{
    std::string datadir("test_cases/matrix-market/");
    for (const auto& fn : { "Keq_constant.txt"
                            "P_mat.mtx"
                            "R_mat.mtx"
                            "S_mat.mtx"
                            "f_log_counts.txt"
                            "results.mtx"
                            "f_log_counts.txt"
                            "v_log_counts.txt" })
    {
        std::cout << "---- Searching for " << fn << " in " << datadir << "\n";
        std::ifstream f((datadir + fn).c_str());
        if (!f.good())
        {
            std::cout << "Didn't find " << fn << "\n";
            return 1;
        }
    }
    
    std::string fn = datadir + "P_mat.mtx";
    struct matrix_market* P_mat_raw   = read_mm(fn);
    fn = datadir + "R_mat.mtx";
    struct matrix_market* R_mat_raw   = read_mm(fn);
    fn = datadir + "S_mat.mtx";
    struct matrix_market* S_mat_raw   = read_mm(fn);
    fn = datadir + "results.mtx";
    struct matrix_market* res_mat_raw = read_mm(fn);

    Map<MatrixXd> P_mat;
    Map<MatrixXd> R_mat;
    Map<MatrixXd> S_mat;
    Map<MatrixXd> results;
     
    new (&P_mat)   Map<MatrixXd>(P_mat_raw->data,   P_mat_raw->rows   * P_mat_raw->cols);
    new (&R_mat)   Map<MatrixXd>(R_mat_raw->data,   R_mat_raw->rows   * R_mat_raw->cols);
    new (&S_mat)   Map<MatrixXd>(S_mat_raw->data,   S_mat_raw->rows   * S_mat_raw->cols);
    new (&results) Map<MatrixXd>(res_mat_raw->data, res_mat_raw->rows * res_mat_raw->cols);

    fn = datadir + "Keq_constant.txt";
    std::vector<double> keq       = read_vector(fn);
    fn = datadir + "f_log_counts.txt";
    std::vector<double> f_log     = read_vector(fn);
    fn = datadir + "v_log_counts.txt";
    std::vector<double> v_log     = read_vector(fn);
    fn = datadir + "state.txt";
    std::vector<double> state_raw = read_vector(fn);

    VectorXd Keq_constant(keq.data());
    VectorXd f_log_counts(f_log.data());
    VectorXd v_log_counts(v_log.data());
    VectorXd state(state.data());

    std::cout << "--- Calling least squares\n";

    Eigen::VectorXd result = least_squares(
            S_mat,
            R_mat,
            P_mat, 
            Keq_constant,
            state,
            f_log_counts,
            v_log_counts);

    std::cout << "--- Successful least squares call.\n"
        << "--- Result: "
        << result;
}
