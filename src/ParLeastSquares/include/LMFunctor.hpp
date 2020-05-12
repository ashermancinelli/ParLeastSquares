#pragma once

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

  /* 
   * The number of parameters, i.e. inputs.
   * AKA the number of variable concentrations
   */
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
