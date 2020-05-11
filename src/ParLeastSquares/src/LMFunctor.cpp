#include <ParLeastSquares>

int LMFunctor::operator()(const VectorXd& log_vcounts, VectorXd& deriv) const
{
  int nrxns = S.rows();
  int nvar = log_vcounts.rows();//make sure this is length and not 1
  int metabolite_count = S.cols();

  VectorXd log_metabolites(log_vcounts.size() + log_fcounts.size());
  log_metabolites << log_vcounts, log_fcounts;

  VectorXd log_Q_inv = -1.0 * ( (R * log_metabolites) + (P * log_metabolites));
  VectorXd log_Q = 1.0 * ( (P * log_metabolites) + (R * log_metabolites));

  VectorXd EKQ_f(nrxns);  //allocate. can break down to one vector but leave as two for clarity right now. 
  VectorXd EKQ_r(nrxns);    

  for (int rxn=0; rxn < nrxns; rxn++){
    double Q_inv_rxn = exp(log_Q_inv(rxn));
    double ekq_f = E_Regulation(rxn) * Keq_constant(rxn) * Q_inv_rxn;

    EKQ_f(rxn) = ekq_f;

    double Q_rxn = exp(log_Q(rxn));
    double ekq_r = E_Regulation(rxn) * pow(Keq_constant(rxn), -1.0) * Q_rxn; 
    EKQ_r(rxn) = ekq_r;
  }

  //printShape("--- S", S);
  // auto _S = S.block(nrxns,nvar); //take all rows (reactions) and only variable columns.

  //(nvar x 1) <=(nvar x nrxns) * (nrxns x 1)
  deriv = (S.topLeftCorner(nrxns, nvar).transpose()) * (EKQ_f - EKQ_r);
  return 0;
}
