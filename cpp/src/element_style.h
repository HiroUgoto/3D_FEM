using EV = Eigen::VectorXd ;
using EV3 = Eigen::Vector3d ;
using EM = Eigen::MatrixXd ;

class ElementStyle {
  public:
    size_t dim;
    size_t ng;
    EV xi, w;

    size_t ng_all;
    std::vector<EV> n_list;
    std::vector<EM> dn_list;
    std::vector<double> w_list;

    EV3 center;
    EM dn_center;

    ElementStyle ();
    virtual ~ElementStyle ();

    virtual bool is_inside (EV3& xi);
    virtual EV shape_function_n (double xi, double eta, double zeta);
    virtual EM shape_function_dn (double xi, double eta, double zeta);
};

// ----------------------------------------------------- //
ElementStyle* set_element_style(const std::string style);
void set_gauss_points (const size_t n, EV& xi, EV& w);
void set_tetra_gauss_points (const size_t n, EM& xi, EV& w);

// ----------------------------------------------------- //
class Solid_3d_10Node: public ElementStyle {
  public:
    Solid_3d_10Node ();
    ~Solid_3d_10Node ();
    bool is_inside (EV3& xi);
    EV shape_function_n (double xi, double eta, double zeta);
    EM shape_function_dn (double xi, double eta, double zeta);
};

// ----------------------------------------------------- //
class Solid_3d_8Node: public ElementStyle {
  public:
    Solid_3d_8Node ();
    ~Solid_3d_8Node ();
    bool is_inside (EV3& xi);
    EV shape_function_n (double xi, double eta, double zeta);
    EM shape_function_dn (double xi, double eta, double zeta);
};

// ----------------------------------------------------- //
class Solid_2d_4Node: public ElementStyle {
  public:
    Solid_2d_4Node ();
    ~Solid_2d_4Node ();
    bool is_inside (EV3& xi);
    EV shape_function_n (double xi, double eta, double zeta=0.0);
    EM shape_function_dn (double xi, double eta, double zeta=0.0);
};

// ----------------------------------------------------- //
class Solid_2d_8Node: public ElementStyle {
  public:
    Solid_2d_8Node ();
    ~Solid_2d_8Node ();
    bool is_inside (EV3& xi);
    EV shape_function_n (double xi, double eta, double zeta=0.0);
    EM shape_function_dn (double xi, double eta, double zeta=0.0);
};

// ----------------------------------------------------- //
class Solid_2d_9Node: public ElementStyle {
  public:
    Solid_2d_9Node ();
    ~Solid_2d_9Node ();
    bool is_inside (EV3& xi);
    EV shape_function_n (double xi, double eta, double zeta=0.0);
    EM shape_function_dn (double xi, double eta, double zeta=0.0);
};

// ----------------------------------------------------- //
class Line_1d_2Node: public ElementStyle {
  public:
    Line_1d_2Node ();
    ~Line_1d_2Node ();
    bool is_inside (EV3& xi);
    EV shape_function_n (double xi, double eta=0.0, double zeta=0.0);
    EM shape_function_dn (double xi, double eta=0.0, double zeta=0.0);
};

// ----------------------------------------------------- //
class Line_1d_3Node: public ElementStyle {
  public:
    Line_1d_3Node ();
    ~Line_1d_3Node ();
    bool is_inside (EV3& xi);
    EV shape_function_n (double xi, double eta=0.0, double zeta=0.0);
    EM shape_function_dn (double xi, double eta=0.0, double zeta=0.0);
};

// ----------------------------------------------------- //
class Input_2d_4Node: public Solid_2d_4Node {};
class Input_1d_2Node: public Line_1d_2Node {};
class Input_1d_3Node: public Line_1d_3Node {};

// ----------------------------------------------------- //
class Visco_2d_4Node: public Solid_2d_4Node {};
class Visco_1d_2Node: public Line_1d_2Node {};
class Visco_1d_3Node: public Line_1d_3Node {};

// ----------------------------------------------------- //
class Connect: public Line_1d_2Node {
  public:
    Connect ();
    ~Connect ();
};
