/*
 *
 * All we need to let someone else use this as a
 * package:
 *
 * Basic input S to calculate potential step
 * which metabolites are fixed and which are variable
 *
 */

#include <vector>
#include <iostream>
#include <typeinfo>
#include <thread>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

#include <Eigen/Core>

#include <ParLeastSquares>

namespace py = pybind11;

/*
 *  %%
 *  input must be able to determine optimization routine. We need the follwing variables:
 *  0: state, type:np.array(float), purpose: current enzyme activities
 *  1: v_log_counts, type:np.array(float), purpose:initial guess for optimization
 *  2: f_log_counts, type:np.array(float), purpose:fixed metabolites
 *  3: mu0, type: , purpose: non as of now
 *  4: S_mat, type: np.array(float), purpose: stoichiometric matrix (rxn by matabolites)
 *  5: R_back_mat, type: np.array(float), purpose: reverse stoichiometric matrix (rxn by matabolites)
 *         #note could be calculated from S_mat: R_back_mat = np.where(S_mat<0, S_mat, 0)
 *  6: P_mat, type: np.array(float), purpose: forward stoichiometric matrix (rxn by matabolites), 
 *         # note could be calculated from S_mat: P_mat = np.where(S_mat>0,S_mat,0)
 *  7: Keq_constant, type: np.array(float), purpose: equilibrium constants
 *
 */
void potential_step(
    const int index,
    Eigen::MatrixXd& S_mat,
    Eigen::MatrixXd& R_back_mat,
    Eigen::MatrixXd& P_mat,
    Eigen::VectorXd& Keq_constant,
    Eigen::VectorXd& E_Regulation,
    Eigen::VectorXd& log_fcounts,
    Eigen::VectorXd& log_vcounts,
    Eigen::MatrixXd& returns,
    int tid,
    DiffType dt,
    const int& maxfev,
    const double& xtol)
{
  std::cout << "--- tid<" << tid << "> potential step being calculated...\n";

  Eigen::VectorXd result = least_squares(
      S_mat,
      R_back_mat,
      P_mat, 
      Keq_constant,
      E_Regulation,
      log_fcounts,
      log_vcounts,
      dt,
      maxfev,
      xtol);

  returns.row(tid) = result;
  std::cout << "Returning from tid<" << tid << ">\n";
}


[[nodiscard]] auto dispatch(
    const std::vector<int>& indices,
    Eigen::MatrixXd& S_mat,
    Eigen::MatrixXd& R_back_mat,
    Eigen::MatrixXd& P_mat,
    Eigen::VectorXd& Keq_constant,
    Eigen::VectorXd& E_Regulation,
    Eigen::VectorXd& log_fcounts,
    Eigen::VectorXd& log_vcounts
    ) -> Eigen::MatrixXd
{

  const int n_threads = indices.size();

  /*
   * Returns mxn where:
   * m = num reactions
   * n = num variable metabolites
   */
  Eigen::MatrixXd returns (S_mat.rows(), S_mat.cols()-log_fcounts.size());
  const int maxfev = 2000;
  const double xtol = 1e-15;

  for (int tid=0; tid<n_threads; tid++)
  {
    potential_step(
        indices[tid],
        S_mat,
        R_back_mat,
        P_mat,
        Keq_constant,
        E_Regulation,
        log_fcounts,
        log_vcounts,
        returns,
        tid,
        DiffType::Numerical,
        maxfev,
        xtol);
  }

  /*
   * Not used for now... Sequential solver first!
   *
   std::vector<std::thread> handles(n_threads);
   for (int tid=0; tid<n_threads; tid++)
   {
   handles[tid] = std::thread(
   potential_step,
   indices[tid],
   std::ref(state),
   std::ref(v_log_counts),
   std::ref(f_log_counts),
   std::ref(mu0),
   std::ref(S_mat),
   std::ref(R_back_mat),
   std::ref(P_mat),
   std::ref(Keq_constant),
   returns,
   tid
   );
   }

   for (auto& th : handles) th.join();
   */

  return returns;
}

PYBIND11_MODULE(pstep, m) {
  m.doc() = "Dispatches jobs to calculate potential steps.";
  m.def("dispatch", &dispatch, "Dispatches jobs");
}
