#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <optional>
#include <ParLeastSquares>

static std::vector<std::string> default_search_dirs =
{
  "../tests/data/tca_glycolosis",
  "../tests/data/tca_glycolosis_gogat",
  "/qfs/projects/boltzmann/data/matrix_market/tca_glycolosis/",
  "/qfs/projects/boltzmann/data/matrix_market/tca_glycolosis_gogat//",
};

static const std::vector<std::string> default_filenames =
{
  "/S_mat.mtx",
  "/R_mat.mtx",
  "/P_mat.mtx",
  "/Keq_constant.txt",
  "/state.txt",
  "/f_log_counts.txt",
  "/v_log_counts.txt",
};

inline bool file_exists(const std::string& path);

bool all_files_exist(
    const std::string& dir,
    const std::vector<std::string>& filenames=default_filenames);

std::optional<std::string> found_acceptable_dir(
    const std::vector<std::string>& dirs,
    const std::vector<std::string>& filenames=default_filenames);

std::ostream& operator<<(std::ostream& out, const DiffType dt);
