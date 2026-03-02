namespace io_data {
  Fem
    input_mesh (const std::string);

  std::tuple<std::map<int, Eigen::Vector3d>, std::vector<int>, std::vector<double>>
    input_pmls(const std::string pml_file);

  std::tuple<std::vector<size_t>, std::vector<size_t>>
    input_outputs (const std::string);
}
