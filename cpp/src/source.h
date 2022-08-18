using EV = Eigen::VectorXd ;
using EM = Eigen::MatrixXd ;


class Source {
  public:
    size_t id;
    size_t element_id;
    double strike, dip, rake, length, width;

    EM dn;
    EV strain_tensor;

    Source(size_t id, double strike, double dip, double rake,
            double length, double width, size_t element_id, EM dn);

  private:
    void set_strain_tensor();
};

std::vector<Source> set_source(std::vector<Element> elements,
    double strike, double dip, double rake, double length, double width,
    double sx, double sy, double sz, size_t nl, size_t nw);
