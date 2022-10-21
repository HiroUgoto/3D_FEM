using EV = Eigen::VectorXd ;
using EV3 = Eigen::Vector3d ;
using EM = Eigen::MatrixXd ;
using EM2 = Eigen::Matrix2d ;
using EM3 = Eigen::Matrix3d ;

class Fault {
  public:
    size_t id;
    size_t pelem_id;
    size_t melem_id;
    std::vector<size_t> spring_id;
    std::vector<double> param;

    double p0, tp, tr, dc;
    Element *pelement_p, *melement_p;
    std::vector<size_t> neighbour_elements_id;
    std::vector<EV3> neighbour_elements_xi;
    EM3 R;

    double slip, sliprate, traction, area;
    bool rupture;

    Fault (size_t id, size_t pelem_id, size_t melem_id, std::vector<size_t> spring_id, std::vector<double> param);

  private:
    void set_param();

  public:
    void set_initial_condition0(std::vector<Element>& elements);
    void set_initial_condition1(std::vector<Element>& elements);

  private:
    void find_neighbour_element(std::vector<Element> elements);
    void set_R();

  public:
    void update_friction(double dt);

  private:
    void calc_average_slip(double dt);

  public:
    void calc_traction(std::vector<Element> elements);
    void update_rupture(std::vector<Element>& elements);

  private:
    EV stress_to_traction(EV stress, EV n);

  public:
    void set_spring_kv(std::vector<Element>& elements);
    void set_spring_kvkh(std::vector<Element>& elements);

};
