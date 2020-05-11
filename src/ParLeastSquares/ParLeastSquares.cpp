

#include <ParLeastSquares.hpp>

inline void printShape(const std::string& s, const Eigen::MatrixXd& m)
{
  std::cout << "Shape of " << s << ": (" << m.rows() << ", " << m.cols() << ")\n";
}

/*
 * AKA calc_jacobian
 */
int LMFunctorAnalytical::df(const VectorXd &log_vcounts, MatrixXd &fjac)
{
  std::cout << __func__ << " start\n";
  //this should be a numerical jacobian
  //Jac is the Jacobian matrix, 
  //an N metabolite time-differential equations by (rows) by 
  //N metabolites derivatives (columns)
  //J_ij = d/dx_i(df_j/dt)

  int nrxns = S.rows();
  int nvar = log_vcounts.size();//make sure this is length and not 1
  int metabolite_count = S.cols();

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

  std::cout << "df end\n";
  return 0;
}

int LMFunctor::operator()(const VectorXd& log_vcounts, VectorXd& deriv) const
{
  std::cout << __func__ << " start\n";
  //this function should be derivs

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
  std::cout << deriv << "\n\n";
  std::cout << __func__ << " end\n";
  return 0;
}

  [[nodiscard]]
Eigen::VectorXd least_squares_analytical(
    Eigen::MatrixXd& S_mat,
    Eigen::MatrixXd& R_back_mat,
    Eigen::MatrixXd& P_mat, 
    Eigen::VectorXd& Keq_constant,
    Eigen::VectorXd& E_Regulation,
    Eigen::VectorXd& log_fcounts,
    Eigen::VectorXd& log_vcounts,
    const int maxfev,
    const double xtol)
{
  LMFunctorAnalytical functor(
      S_mat,
      R_back_mat,
      P_mat,
      Keq_constant,
      E_Regulation,
      log_fcounts);
  Eigen::LevenbergMarquardt<LMFunctorAnalytical> lm(functor);

  lm.parameters.maxfev = maxfev;
  lm.parameters.xtol = xtol;
  int return_value = lm.minimize(log_vcounts);

  return log_vcounts;
}

  [[nodiscard]]
Eigen::VectorXd least_squares_numerical(
    Eigen::MatrixXd& S_mat,
    Eigen::MatrixXd& R_back_mat,
    Eigen::MatrixXd& P_mat, 
    Eigen::VectorXd& Keq_constant,
    Eigen::VectorXd& E_Regulation,
    Eigen::VectorXd& log_fcounts,
    Eigen::VectorXd& log_vcounts,
    const int maxfev,
    const double xtol)
{
  LMFunctor functor(
      S_mat,
      R_back_mat,
      P_mat,
      Keq_constant,
      E_Regulation,
      log_fcounts);
  Eigen::NumericalDiff<LMFunctor> numDiff(functor);
  Eigen::LevenbergMarquardt<Eigen::NumericalDiff<LMFunctor>, double> lm(numDiff);

  lm.parameters.maxfev = 20000;
  lm.parameters.xtol = 1e-15;

  int return_value = lm.minimize(log_vcounts);

  return log_vcounts;
}

/*
 * Driver function for calculating the least squares with the
 * Levenberg-Marquardt method
 */
  [[nodiscard]]
Eigen::VectorXd least_squares(
    Eigen::MatrixXd& S_mat,
    Eigen::MatrixXd& R_back_mat,
    Eigen::MatrixXd& P_mat, 
    Eigen::VectorXd& Keq_constant,
    Eigen::VectorXd& E_Regulation,
    Eigen::VectorXd& log_fcounts,
    Eigen::VectorXd& log_vcounts,
    DiffType dt,
    const int maxfev,
    const double xtol)
{
  switch (dt)
  {
    case Numerical:
      return least_squares_numerical(
          S_mat,
          R_back_mat,
          P_mat, 
          Keq_constant,
          E_Regulation,
          log_fcounts,
          log_vcounts,
          maxfev,
          xtol);
    case Analytical:
      return least_squares_analytical(
          S_mat,
          R_back_mat,
          P_mat, 
          Keq_constant,
          E_Regulation,
          log_fcounts,
          log_vcounts,
          maxfev,
          xtol);
  };
  assert (false && "Must pass valid differentiation method.");
}

  [[nodiscard]]
Eigen::MatrixXd read_mm(const std::string& path)
{
  std::ifstream file(path);
  // Ignore comments headers
  while (file.peek() == '%') file.ignore(2048, '\n');

  int rows, cols;
  file >> rows >> cols;
  std::cout << "Found dimensions ("
    << rows << ", " << cols << ") for file "
    << path << "\n";

  double d;
  Eigen::MatrixXd mat(rows, cols);
  for (int j=0; j < cols; j++)
    for (int i=0; i < rows; i++)
    {
      file >> d;
      mat(i, j) = d;
    }
  file.close();
  std::cout << mat << "\n";

  return mat;
}

  [[nodiscard]]
Eigen::VectorXd read_vector(const std::string& path)
{
  std::vector<double> vec;
  std::ifstream f(path);
  std::string line;

  while (getline(f, line))
  {
    double d;
    d = std::atof(line.c_str());
    vec.push_back(d);
  }

  std::cout << "Reading vector with size " << vec.size()
    << " from " << path << "\n";
  VectorXd eigen_vec(vec.size());
  for (int i=0; i<vec.size(); i++) eigen_vec(i) = vec[i];
  return eigen_vec;
}
