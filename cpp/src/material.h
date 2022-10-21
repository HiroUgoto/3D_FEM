using EV = Eigen::VectorXd ;
using EM = Eigen::MatrixXd ;
using EM3 = Eigen::Matrix3d ;

class Material {
  public:
    size_t id;
    std::string style;
    std::vector<double> param;
    double rmu, rlambda, rho;
    double kh, kv;

    Material();
    Material(size_t id, std::string style, std::vector<double> param);

    void set_init(size_t id, std::string style, std::vector<double> param);
    void print();

  private:
    void set_param();

  public:
    EM mk_d(const size_t dof);
    EM mk_d_spring();
    EM mk_d_slider();

    EM mk_visco(const size_t dof);
    EM3 mk_imp(const size_t dof);
};
