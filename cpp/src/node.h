using EV3 = Eigen::Vector3d ;

class Node {
  public:
    size_t id, dof;
    std::vector<double> xyz;
    std::vector<size_t> freedom;

    EV3 u, um, v;
    EV3 mass, c;
    EV3 force;

    Node ();
    Node (size_t id, std::vector<double> xyz, std::vector<size_t> freedom);

    void print();
    void set_initial_condition();
};
