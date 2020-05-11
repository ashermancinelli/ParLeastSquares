#pragma once
#include <LMFunctor.hpp>

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
