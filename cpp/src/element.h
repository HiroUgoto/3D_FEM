using EV = Eigen::VectorXd ;
using EM = Eigen::MatrixXd ;
using EM2 = Eigen::Matrix2d ;
using EM3 = Eigen::Matrix3d ;

class Element {
  public:
    size_t id;
    int material_id;
    std::string style;
    std::vector<size_t> inode;
    double gravity;

    std::vector<Node*> nodes_p;

    Material material;
    double rho, mass;

    size_t nnode, dim, ng_all, dof, ndof;
    EM xnT;
    EM dn_center;

    EV M_diag, C_diag;
    EM K, C_off_diag;
    EM De, Dv;
    EM3 imp;
    EV force;
    EV strain, stress;

    EM3 R;
    double traction;

    Element (size_t id, std::string style, int material_id, std::vector<size_t> inode);
    void print();

  private:
    void set_style();

  public:
    void set_nodes(std::vector<Node*> nodes_p);
    void set_material(Material* material_p);
    void set_pointer_list();
    void set_xn();

    void mk_local_matrix_init(const size_t dof);
    void mk_local_matrix();
    void mk_local_vector();

    void mk_ku();
    void mk_cv();
    void mk_ku_cv();

    EV mk_u_hstack();

  private:
    EV mk_v_hstack();
    EM mk_u_vstack();

  public:
    void update_inputwave(const EV vel0);
    void mk_source(const EM dn, const EV strain_tensor, const double slip0);
    void mk_T(const EV T);
    void calc_stress();
    EV calc_stress_xi(const EV xi);
    std::tuple<bool, EV3> check_inside(const EV3 x, double margin=0.0);
};

EM mk_m(const EM N);
EM mk_n(const size_t dof, const size_t nnode, const EV n);

EM mk_nqn(const EM N, const EM3 q, const EM3 imp);
std::tuple<double, EM3> mk_q(const size_t dof, const EM xnT, const EM dn);

EM mk_k(const EM B, const EM D);
EM mk_b(const size_t dof, const size_t nnode, const EM dnj);
EM mk_b_T(const size_t dof, const size_t nnode, const EM dnj);

std::tuple<double, EM>  mk_dnj(const EM xnT, const EM dn);
std::tuple<double, EM3> mk_inv_jacobi(const EM xnT, const EM dn);
std::tuple<double, EM3> mk_jacobi(const EM xnT, const EM dn);
