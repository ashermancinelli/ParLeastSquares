#include <test_utilities.hpp>

std::ostream& operator<<(std::ostream& out, const DiffType dt)
{

#define GEN_CASE(diff) case diff: out << #diff; break;

  switch (dt)
  {
    GEN_CASE(DiffType::Analytical);
    GEN_CASE(DiffType::Numerical);
  }

#undef GEN_CASE

  return out;
}

inline bool file_exists(const std::string& path)
{
  std::ifstream f(path.c_str());
  return f.good();
}

bool all_files_exist(
    const std::string& dir,
    const std::vector<std::string>& filenames)
{
  for (const auto& filename : filenames)
  {
    if (!file_exists(dir + filename))
      return false;
  }

  return true;
}

std::optional<std::string> found_acceptable_dir(
    const std::vector<std::string>& dirs,
    const std::vector<std::string>& filenames)
{
  for (const auto& dir : dirs)
  {
    std::cout << "Checking dir " << dir << "\n";
    if (all_files_exist(dir, filenames))
    {
      std::cout << "Found all files in dir " << dir << "\n";
      return dir;
    }
  }
  return std::nullopt;
}
