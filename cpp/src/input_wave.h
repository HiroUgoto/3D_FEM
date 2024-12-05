namespace input_wave {
  std::tuple<Eigen::VectorXd, double>
    linspace(const double start, const double end, const size_t num);
  Eigen::VectorXd
    ricker(const Eigen::VectorXd tim, const double fp, const double tp, const double amp);
  Eigen::VectorXd
    simple_sin(const Eigen::VectorXd tim, const double fp, const double amp);
  std::tuple<Eigen::VectorXd, Eigen::VectorXd, double>
    input_acc_file(const std::string input_file);
}
