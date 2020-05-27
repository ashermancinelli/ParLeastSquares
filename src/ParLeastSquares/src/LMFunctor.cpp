#include <ParLeastSquares>

int LMFunctor::operator()(const VectorXd& log_vcounts, VectorXd& deriv) const
{
  int nrxns = static_cast<int>(S.rows());
  int nvar = static_cast<int>(log_vcounts.size());//make sure this is length and not 1

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

  deriv = (S.topLeftCorner(nrxns, nvar).transpose()) * (EKQ_f - EKQ_r);
  //std::cout << "L2 norm: " << deriv.squaredNorm() << "\n";
  //std::cout<< "log_vcounts " << log_vcounts << "\n";
  //std::cout<< "deriv " << deriv << "\n";
  return 0;
}
