#include <ParLeastSquares>

/*
 * AKA calc_jacobian
 */
int LMFunctorAnalytical::df(const VectorXd &log_vcounts, MatrixXd &fjac)
{
  //this should be a numerical jacobian
  //Jac is the Jacobian matrix, 
  //an N metabolite time-differential equations by (rows) by 
  //N metabolites derivatives (columns)
  //J_ij = d/dx_i(df_j/dt)

  int nrxns = static_cast<int>(S.rows());
  int nvar = static_cast<int>(log_vcounts.size());//make sure this is length and not 1

  //WARNING Only use to calcualte KQ
  VectorXd log_metabolites(log_vcounts.size() + log_fcounts.size());
  log_metabolites << log_vcounts, log_fcounts;

  VectorXd metabolites = log_metabolites.array().exp();
  VectorXd metabolites_recip = metabolites.array().pow(-1.0);

  //nrxn x metabolite_count <= component product from:  (metabolite_count x 1) * (nrxn x metabolite_count)
  auto MR = metabolites_recip.rowwise().replicate(S.rows()).transpose();
  auto _S = (-1. * S);

  MatrixXd S_recip_metab = MR.array() * _S.array();

  VectorXd log_Q_inv = -1.0 * ( (R * log_metabolites) + (P * log_metabolites));
  VectorXd log_Q = 1.0 * ( (P * log_metabolites) + (R * log_metabolites));

  VectorXd x(nrxns);

  for (int rxn=0; rxn < nrxns; rxn++){
    double Q_inv_rxn = exp(log_Q_inv(rxn));
    double ekq_f = E_Regulation(rxn) * Keq_constant(rxn) * Q_inv_rxn;

    double Q_rxn = exp(log_Q(rxn));
    double ekq_r = E_Regulation(rxn) * pow(Keq_constant(rxn),-1.0) * Q_rxn; 

    x(rxn) = ekq_f + ekq_r;
  }

  //nrxn x metabolite_count <= component (nrxn x 1 ) * (nrxn x metabolite_count), but not sure how cwiseProduce is working. 
  //maybe need to do (x.transpose()).cwiseProduce(S_recip_metab.transpose()).transpose()
  MatrixXd y = x.rowwise().replicate(S.cols()).array() * S_recip_metab.array();

  fjac = ((S.transpose()) * y).block(0, 0, nvar, nvar);

  return 0;
}
