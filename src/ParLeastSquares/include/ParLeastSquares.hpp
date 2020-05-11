#pragma once

#include <cmath>
#include <cassert>
#include <fstream>
#include <iostream>
#include <Eigen/Eigen>
#include <unsupported/Eigen/NonLinearOptimization>

#include <par_ls_defs.hpp>

using namespace Eigen;

static constexpr int CPP98      = 199711;
static constexpr int GCC98      = 199711;
static constexpr int CPP11      = 201103;
static constexpr int GNU11      = 201103;
static constexpr int CPP14      = 201402;
static constexpr int GNU14      = 201402;
static constexpr int CPP1z      = 201500;
static constexpr int CPP17      = 201500;
static constexpr int MY_CPP_STD = __cplusplus;

enum DiffType
{
  Analytical,
  Numerical,
};

struct LMFunctor
{
  typedef double Scalar;
  enum {
    InputsAtCompileTime = Eigen::Dynamic,
    ValuesAtCompileTime = Eigen::Dynamic
  };
  typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
  typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
  typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

  // Number of data points, i.e. values.
  int m;
  // The number of parameters, i.e. inputs.
  int n;

  MatrixXd S;
  MatrixXd R;
  MatrixXd P;
  VectorXd Keq_constant;
  VectorXd E_Regulation;
  VectorXd log_fcounts;

  LMFunctor(
      MatrixXd& _S,
      MatrixXd& _R, 
      MatrixXd& _P,
      VectorXd& _Keq_constant,
      VectorXd& _E_Regulation,
      VectorXd& _log_fcounts):
    m(static_cast<int>(_S.rows())),
    n(static_cast<int>(_S.cols()-_log_fcounts.size())),
    S(_S),
    R(_R),
    P(_P),
    Keq_constant(_Keq_constant),
    E_Regulation(_E_Regulation),
    log_fcounts(_log_fcounts) {}

  // Returns 'm', the number of values.
  inline int values() const { return m; }
  // Returns 'n', the number of inputs.
  inline int inputs() const { return n; }

  int operator()(const VectorXd& log_vcounts, VectorXd& deriv) const;
};

struct LMFunctorAnalytical : public LMFunctor
{
  LMFunctorAnalytical(
      MatrixXd& _S,
      MatrixXd& _R, 
      MatrixXd& _P,
      VectorXd& _Keq_constant,
      VectorXd& _E_Regulation,
      VectorXd& _log_fcounts):
    LMFunctor(_S, _R, _P, _Keq_constant, _E_Regulation, _log_fcounts) {}
  int df(const VectorXd &log_vcounts, MatrixXd &fjac);
};

/*
 * Driver function for calculating the least squares with the
 * Levenberg-Marquardt method
 */
[[nodiscard]]
extern Eigen::VectorXd least_squares(
    Eigen::MatrixXd& S_mat,
    Eigen::MatrixXd& R_back_mat,
    Eigen::MatrixXd& P_mat, 
    Eigen::VectorXd& Keq_constant,
    Eigen::VectorXd& E_Regulation,
    Eigen::VectorXd& log_fcounts,
    Eigen::VectorXd& log_vcounts,
    DiffType dt=DiffType::Analytical,
    const int maxfevs=2e4,
    const double xtol=1e-15);

[[nodiscard]]
Eigen::MatrixXd read_mm(const std::string& path);

[[nodiscard]]
Eigen::VectorXd read_vector(const std::string& path);
