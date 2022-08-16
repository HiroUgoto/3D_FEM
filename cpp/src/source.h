using EV = Eigen::VectorXd ;
using EM = Eigen::MatrixXd ;


class Source {
  public:
    size_t id;
    size_t element_id;
    double dip, width;

    EM dn;
    EV strain_tensor;

    Source(size_t id, double dip, double width, size_t element_id, EM dn);

  private:
    void set_strain_tensor();
};

std::vector<Source> set_source(std::vector<Element> elements, double dip, double width, double sx, double sy, double sz, size_t n);
