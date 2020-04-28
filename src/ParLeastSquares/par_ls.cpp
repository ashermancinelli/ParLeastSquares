

#include <ParLeastSquares>

inline void printShape(const std::string& s, const Eigen::MatrixXd& m)
{
    std::cout << "Shape of " << s << ": (" << m.rows() << ", " << m.cols() << ")\n";
}

int lmder_functor::df(const VectorXd &log_vcounts, MatrixXd &fjac)
{
    std::cout << __func__ << " start\n";
    //this should be a numerical jacobian
    //Jac is the Jacobian matrix, 
    //an N metabolite time-differential equations by (rows) by 
    //N metabolites derivatives (columns)
    //J_ij = d/dx_i(df_j/dt)

    int nrxns = S.rows();
    int nvar = log_vcounts.rows();//make sure this is length and not 1
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

    fjac = (S.transpose()) * y;

    std::cout << "df end\n";
    return 0;
}

int lmder_functor::operator()(const VectorXd& log_vcounts, VectorXd& deriv)
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
    std::cout << "log_Q_inv size: " << log_Q_inv.size()
        << "\tlog_Q size: " << log_Q.size() << "\n";

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
    //printShape("deriv", deriv);
    std::cout << __func__ << " end\n";
    return 0;
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
        Eigen::VectorXd& log_vcounts)
{
    lmder_functor functor(
            S_mat,
            R_back_mat,
            P_mat,
            Keq_constant,
            E_Regulation,
            log_fcounts);

    std::cout << "Passed functor initialization\n";

    Eigen::LevenbergMarquardt<lmder_functor> lm(functor);

    std::cout << "Passed eigen thing initialization\n";

    lm.minimize(log_vcounts);

    std::cout << "Passed minimization\n";

    return log_vcounts;
}

void read_mm_scpy(struct matrix_market* mat, std::istream& file)
{
    int num_row, num_col;
    file >> num_row >> num_col;
    mat->rows = num_row;
    mat->cols = num_col;
    const int num_lines = num_col * num_row;
    mat->data = new double[num_lines];
    std::fill(mat->data, mat->data + num_lines, 0.);
    for (int i = 0; i < num_lines; i++)
    {
        double d;
        file >> d;
        mat->data[i] = d;
    }
}

void read_mm_std(struct matrix_market* mat, std::istream& file)
{
    int num_row, num_col, num_lines;
    // Read number of rows and columns
    file >> num_row >> num_col >> num_lines;
    mat->rows = num_row;
    mat->cols = num_col;

    // Create 2D array and fill with zeros
    mat->data = new double[num_row * num_col];      
    double* data = mat->data;

    std::fill(data, data + num_row * num_col, 0.);

    // fill the mat with data
    for (int l = 0; l < num_lines; l++)
    {
        double d;
        int row, col;
        file >> row >> col >> d;
        data[(row - 1) + (col - 1) * num_row] = d;
    }
}

[[nodiscard]]
struct matrix_market* read_mm(std::string& fn, const bool scipy_fmt=true)
{
    std::ifstream file(fn);
    // Ignore comments headers
    while (file.peek() == '%') file.ignore(2048, '\n');

    auto mat = new matrix_market;
    scipy_fmt ? read_mm_scpy(mat, file) : read_mm_std(mat, file);
    file.close();

    return mat;
}

[[nodiscard]]
std::vector<double> read_vector(const std::string& path)
{
    std::vector<double> vec;
    std::ifstream f(path);
    while (!f.eof())
    {
        double d;
        f >> d;
        vec.push_back(d);
    }
    return vec;
}
