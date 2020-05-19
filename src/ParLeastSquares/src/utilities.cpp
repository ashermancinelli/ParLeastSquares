#include <ParLeastSquares>

inline void printShape(const std::string& s, const Eigen::MatrixXd& m)
{
  std::cout << "Shape of " << s << ": (" << m.rows() << ", " << m.cols() << ")\n";
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

bool write_mm(
    const Eigen::MatrixXd& mtx,
    const std::string& path)
{
  std::ofstream file(path);
  if (!file.good())
  {
    return false;
  }

  const int rows = static_cast<int>(mtx.rows()),
    cols = static_cast<int>(mtx.cols());

  file << rows << cols << "\n" << std::setprecision(10);
  for (int j=0; j < cols; j++)
  {
    for (int i=0; i < rows; i++)
    {
      file << mtx(i, j) << "\n";
    }
  }
  return true;
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
  for (unsigned int i=0; i<vec.size(); i++) eigen_vec(i) = vec[i];
  return eigen_vec;
}

bool write_vector(
    const Eigen::VectorXd& vec,
    const std::string& path)
{
  std::ofstream file(path);
  if (!file.good())
  {
    return false;
  }

  file << std::setprecision(10);
  for (int i=0; i<vec.size(); i++)
  {
    file << vec(i) << "\n";
  }
  return true;
}
