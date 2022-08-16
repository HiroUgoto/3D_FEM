#include "all.h"
#include <Eigen/Dense>
#include "node.h"
#include "material.h"
#include "element_style.h"
#include "element.h"

using EV = Eigen::VectorXd ;
using EV3 = Eigen::Vector3d ;
using EM = Eigen::MatrixXd ;
using EM2 = Eigen::Matrix2d ;
using EM3 = Eigen::Matrix3d ;

// ------------------------------------------------------------------- //
Element::Element (size_t id, std::string style, int material_id, std::vector<size_t> inode)
  {
    this->id = id;
    this->style = style;
    this->material_id = material_id;
    this->inode = inode;

    this->gravity = 9.8;
    this->nnode = inode.size();

    this->set_style();
  }

// ------------------------------------------------------------------- //
void Element::print() {
    std::cout << this->id << ": ";
    std::cout << this->style << ", ";
    std::cout << this->material_id << ", ";

    for (size_t i = 0 ; i < this->inode.size() ; i++) {
      std::cout << this->inode.at(i) << " ";
    }
    std::cout << "\n";
  }

void Element::set_style() {
    ElementStyle* estyle_p = set_element_style(this->style);

    this->dim = estyle_p->dim;
    this->ng = estyle_p->ng;
    this->xi = estyle_p->xi;
    this->w  = estyle_p->w;

    this->ng_all = estyle_p->ng_all;
    this->n_list  = estyle_p->n_list;
    this->dn_list = estyle_p->dn_list;
    this->w_list = estyle_p->w_list;

    this->dn_center = estyle_p->dn_center;
  }

// ------------------------------------------------------------------- //
void Element::set_nodes(std::vector<Node*> nodes_p) {
    this->nodes_p = nodes_p;
  }

void Element::set_material(Material* material_p) {
    if (material_p == nullptr) {
      this->rho = 0.0;
    } else {
      size_t id = material_p->id;
      std::string style = material_p->style;
      std::vector<double> param = material_p->param;

      this->material.set_init(id,style,param);
      this->rho = this->material.rho;
    }
  }

void Element::set_xn(){
    this->xnT = EM::Zero(3,this->nnode);

    for (size_t inode = 0 ; inode < this->nnode ; inode++ ) {
      Node* node = this->nodes_p[inode];
      this->xnT(0,inode) = node->xyz[0] + node->u[0];
      this->xnT(1,inode) = node->xyz[1] + node->u[1];
      this->xnT(2,inode) = node->xyz[2] + node->u[2];
    }
  }

// ------------------------------------------------------------------- //
void Element::mk_local_matrix_init(const size_t dof){
    this->dof = dof;
    this->ndof = dof*this->nnode;

    this->M_diag = EV::Zero(this->ndof);

    this->K = EM::Zero(this->ndof,this->ndof);
    this->K_diag = EV::Zero(this->ndof);
    this->K_off_diag = EM::Zero(this->ndof,this->ndof);

    this->C = EM::Zero(this->ndof,this->ndof);
    this->C_diag = EV::Zero(this->ndof);
    this->C_off_diag = EM::Zero(this->ndof,this->ndof);

    this->force = EV::Zero(this->ndof);

    if (this->dim == 3) {
      double V = 0.0;
      for (size_t i = 0 ; i < this->ng_all ; i++){
        auto [det, jacobi] = mk_jacobi(this->xnT, this->dn_list[i]);
        V += det * this->w_list[i];
      }
      this->mass = this->rho * V;

    } else if (this->dim == 2) {
      this->imp = this->material.mk_imp(this->dof);
    }
  }

// ------------------------------------------------------------------- //
void Element::mk_local_matrix() {
    if (this->dim == 3) {
      EM M = EM::Zero(this->ndof,this->ndof);

      this->C = EM::Zero(this->ndof,this->ndof);
      this->K = EM::Zero(this->ndof,this->ndof);

      this->De = this->material.mk_d(this->dof);
      this->Dv = this->material.mk_visco(this->dof);

      for (size_t i = 0 ; i < this->ng_all ; i++){
        double detJ;
        EM N, Me, B, K, C;

        auto [det, dnj] = mk_dnj(this->xnT, this->dn_list[i]);

        N = mk_n(this->dof, this->nnode, this->n_list[i]);
        Me = mk_m(N);

        B = mk_b(this->dof, this->nnode, dnj);
        K = mk_k(B, this->De);
        C = mk_k(B, this->Dv);

        detJ = det * this->w_list[i];

        M += Me * detJ;
        this->C += C * detJ;
        this->K += K * detJ;
      }

      double tr_M = M.trace() / this->dof;
      this->M_diag = M.diagonal() * this->mass/tr_M;

      this->K_diag = this->K.diagonal();
      this->K_off_diag = this->K_diag.asDiagonal();
      this->K_off_diag = (this->K) - (this->K_off_diag);

      this->C_diag = this->C.diagonal();
      this->C_off_diag = this->C_diag.asDiagonal();
      this->C_off_diag = (this->C) - (this->C_off_diag);

    } else if (this->dim == 2) {
      if (this->style.find("input") != std::string::npos ||
          this->style.find("visco") != std::string::npos) {
        this->C = EM::Zero(this->ndof,this->ndof);

        for (size_t i = 0 ; i < this->ng_all ; i++){
          double detJ;
          EM N, NqN;

          auto [det, q] = mk_q(this->dof, this->xnT, this->dn_list[i]);

          N = mk_n(this->dof, this->nnode, this->n_list[i]);
          NqN = mk_nqn(N, q, this->imp);

          detJ = det * this->w_list[i];
          this->C += NqN * detJ;
        }

        this->C_diag = this->C.diagonal();
        this->C_off_diag = this->C_diag.asDiagonal();
        this->C_off_diag = (this->C) - (this->C_off_diag);
      }
    }
  }

// ------------------------------------------------------------------- //
void Element::mk_local_vector() {
    if (this->dim == 3) {
      this->force = EV::Zero(this->ndof);

      double V = 0.0;
      for (size_t i = 0 ; i < this->ng_all ; i++){
        double detJ;
        EM N;

        auto [det, jacobi] = mk_jacobi(this->xnT, this->dn_list[i]);
        N = mk_n(this->dof, this->nnode, this->n_list[i]);

        detJ = det * this->w_list[i];

        V += detJ;
        this->force += N.row(2)*detJ * this->gravity;
      }

      this->force *= this->mass / V;
    }
  }

// ------------------------------------------------------------------- //
void Element::mk_ku() {
    EV u, ku;

    u = this->mk_u_hstack();
    ku = this->K * u;

    for (size_t inode = 0 ; inode < this->nnode ; inode++){
      size_t i0 = inode*this->dof;
      for (size_t i = 0 ; i < this->dof ; i++) {
        this->nodes_p[inode]->force(i) += ku(i0+i);
      }
    }
  }

void Element::mk_ku_up() {
    EV u, ku;

    u = this->mk_u_hstack_up();
    ku = this->K * u;

    for (size_t inode = 0 ; inode < this->nnode ; inode++){
      size_t i0 = inode*this->dof;
      for (size_t i = 0 ; i < this->dof ; i++) {
        this->nodes_p[inode]->force(i) += ku(i0+i);
      }
    }

  }

void Element::mk_cv() {
    EV v, cv;

    v = this->mk_v_hstack();
    cv = this->C_off_diag * v;

    for (size_t inode = 0 ; inode < this->nnode ; inode++){
      size_t i0 = inode*this->dof;
      for (size_t i = 0 ; i < this->dof ; i++) {
        this->nodes_p[inode]->force(i) += cv(i0+i);
      }
    }
  }

void Element::mk_ku_cv() {
    EV u, ku;
    EV v, cv;

    u = this->mk_u_hstack();
    v = this->mk_v_hstack();
    ku = this->K * u;
    cv = this->C_off_diag * v;

    for (size_t inode = 0 ; inode < this->nnode ; inode++){
      size_t i0 = inode*this->dof;
      for (size_t i = 0 ; i < this->dof ; i++) {
        this->nodes_p[inode]->force(i) += ku(i0+i) + cv(i0+i);
      }
    }
  }


// ------------------------------------------------------------------- //
EV Element::mk_u_hstack() {
    EV u(this->ndof);

    for (size_t inode = 0 ; inode < this->nnode ; inode++){
      size_t i0 = inode*this->dof;
      Node* node = this->nodes_p[inode];
      for (size_t i = 0 ; i < this->dof ; i++) {
        u(i0+i) = node->u[i];
      }
    }
    return u;
  }

EV Element::mk_v_hstack() {
    EV v(this->ndof);

    for (size_t inode = 0 ; inode < this->nnode ; inode++){
      size_t i0 = inode*this->dof;
      Node* node = this->nodes_p[inode];
      for (size_t i = 0 ; i < this->dof ; i++) {
        v(i0+i) = node->v[i];
      }
    }
    return v;
  }

EV Element::mk_u_hstack_up() {
    EV up(this->ndof);

    for (size_t inode = 0 ; inode < this->nnode ; inode++){
      size_t i0 = inode*this->dof;
      Node* node = this->nodes_p[inode];
      for (size_t i = 0 ; i < this->dof ; i++) {
        up(i0+i) = node->_up[i];
      }
    }
    return up;
  }

EM Element::mk_u_vstack() {
    EM u(this->nnode,this->dof);

    for (size_t inode = 0 ; inode < this->nnode ; inode++){
      Node* node = this->nodes_p[inode];
      for (size_t i = 0 ; i < this->dof ; i++) {
        u(inode,i) = node->u[i];
      }
    }
    return u;
  }


// ------------------------------------------------------------------- //
void Element::update_inputwave(const EV vel0) {
    EV v(this->ndof), cv(this->ndof);

    for (size_t inode = 0 ; inode < this->nnode ; inode++){
      size_t i0 = inode*this->dof;
      for (size_t i = 0 ; i < this->dof ; i++) {
        v(i0+i) = vel0(i);
      }
    }
    cv = this->C * v;

    for (size_t inode = 0 ; inode < this->nnode ; inode++){
      size_t i0 = inode*this->dof;
      for (size_t i = 0 ; i < this->dof ; i++) {
        this->nodes_p[inode]->force(i) -= 2.0*cv(i0+i);
      }
    }
  }

// ------------------------------------------------------------------- //
void Element::update_bodyforce(const EV acc0) {
    this->mk_bodyforce(acc0);

    for (size_t inode = 0 ; inode < this->nnode ; inode++){
      size_t i0 = inode*this->dof;
      for (size_t i = 0 ; i < this->dof ; i++) {
        this->nodes_p[inode]->force(i) -= this->force(i0+i);
      }
    }
  }

void Element::mk_bodyforce(const EV acc0) {
    if (this->dim == 2) {
      this->force = EV::Zero(this->ndof);

      double V = 0.0;
      for (size_t i = 0 ; i < this->ng_all ; i++){
        double detJ;
        EM N;

        auto [det, jacobi] = mk_jacobi(this->xnT, this->dn_list[i]);
        N = mk_n(this->dof, this->nnode, this->n_list[i]);

        detJ = det * this->w_list[i];

        V += detJ;
        this->force += (N.row(0)*acc0[0] + N.row(1)*acc0[1]) *detJ;
      }

      this->force *= this->mass / V;
    }

  }

// ------------------------------------------------------------------- //
void Element::mk_source(const EM dn, const EV strain_tensor, const double slip0) {
  if (this->dim == 3) {
    EM BT;
    EV moment;

    auto [det, dnj] = mk_dnj(this->xnT, dn);
    BT = mk_b_T(this->dof, this->nnode, dnj);
    moment = this->material.rmu * strain_tensor * slip0;
    this->force = BT * moment;
  }
}

// ------------------------------------------------------------------- //
void Element::calc_stress() {
    EV u;
    EM B;

    auto [det, dnj] = mk_dnj(this->xnT, this->dn_center);
    B = mk_b(this->dof, this->nnode, dnj);
    u = this->mk_u_hstack();

    this->strain = B * u;
    this->stress = this->De * this->strain;
  }

// ------------------------------------------------------------------- //
std::tuple<bool, EV>
  Element::check_inside(const EV x) {
    EV xi = EV::Zero(3);
    bool is_inside = false;
    ElementStyle* estyle_p = set_element_style(this->style);

    for (size_t itr=0; itr<20; itr++) {
      EV n = estyle_p->shape_function_n(xi[0],xi[1],xi[2]);
      EM dn = estyle_p->shape_function_dn(xi[0],xi[1],xi[2]);

      EV J_func = this->xnT*n - x;
      auto [det,dJ_func] = mk_jacobi(this->xnT, dn);

      EV r = dJ_func.partialPivLu().solve(J_func);
      if (r.norm() < 1e-8) break;

      xi -= r;
    }

    if ( (-1.0 <= xi[0]) && (xi[0] < 1.0) &&
         (-1.0 <= xi[1]) && (xi[1] < 1.0) &&
         (-1.0 <= xi[2]) && (xi[2] < 1.0) ) {
      is_inside = true;
    }

    return {is_inside, xi};
  }

// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //
EM mk_m(const EM N) {
    EM M;

    M = N.transpose() * N;
    return M;
  }

EM mk_n(const size_t dof, const size_t nnode, const EV n) {
    EM N = EM::Zero(dof,dof*nnode);

    if (dof == 1) {
      for (size_t i = 0; i < nnode; i++){
        N(0,i) = n(i);
      }
    } else if (dof == 2) {
      for (size_t i = 0; i < nnode; i++){
        size_t i0 = 2*i;
        size_t i1 = 2*i+1;

        N(0,i0) = n(i);
        N(1,i1) = n(i);
      }
    } else if (dof == 3) {
      for (size_t i = 0; i < nnode; i++){
        size_t i0 = 3*i;
        size_t i1 = 3*i+1;
        size_t i2 = 3*i+2;

        N(0,i0) = n(i);
        N(1,i1) = n(i);
        N(2,i2) = n(i);
      }
    }

    return N;
  }

// ------------------------------------------------------------------- //
EM mk_nqn(const EM N, const EM q, const EM imp) {
    EM nqn;

    nqn = N.transpose() * q.transpose() * imp * q * N;
    return nqn;
  }

std::tuple<double, EM>
  mk_q(const size_t dof, const EM xnT, const EM dn) {
    EM q, t;
    EV3 t0, t1, n;
    double det;

    t = xnT * dn;

    t0 = t.col(0); t1 = t.col(1);
    n = t0.cross(t1);
    det = n.norm();

    t0 /= t0.norm();
    t1 /= t1.norm();
    n /= det;

    q = EM::Zero(3,3);
    q(0,0) =  n(0); q(0,1) =  n(1); q(0,2) =  n(2);
    q(1,0) = t0(0); q(1,1) = t0(1); q(1,2) = t0(2);
    q(2,0) = t1(0); q(2,1) = t1(1); q(2,2) = t1(2);

    return {det, q};
  }

// ------------------------------------------------------------------- //
EM mk_k(const EM B, const EM D) {
    EM K;

    K = B.transpose() * D * B;
    return K;
  }

EM mk_b(const size_t dof, const size_t nnode, const EM dnj) {
    EM B;

    B = EM::Zero(6,3*nnode);
    for (size_t i = 0; i < nnode; i++){
      size_t i0 = 3*i;
      size_t i1 = 3*i+1;
      size_t i2 = 3*i+2;

      B(0,i0) = dnj(i,0);
      B(1,i1) = dnj(i,1);
      B(2,i2) = dnj(i,2);

      B(3,i0) = dnj(i,1);
      B(3,i1) = dnj(i,0);

      B(4,i1) = dnj(i,2);
      B(4,i2) = dnj(i,1);

      B(5,i0) = dnj(i,2);
      B(5,i2) = dnj(i,0);
    }

    return B;
  }

EM mk_b_T(const size_t dof, const size_t nnode, const EM dnj) {
    EM B;

    B = EM::Zero(3*nnode,6);
    for (size_t i = 0; i < nnode; i++){
      size_t i0 = 3*i;
      size_t i1 = 3*i+1;
      size_t i2 = 3*i+2;

      B(i0,0) = dnj(i,0);
      B(i1,1) = dnj(i,1);
      B(i2,2) = dnj(i,2);

      B(i0,3) = dnj(i,1);
      B(i1,3) = dnj(i,0);

      B(i1,4) = dnj(i,2);
      B(i2,4) = dnj(i,1);

      B(i0,5) = dnj(i,2);
      B(i2,5) = dnj(i,0);
    }

    return B;
  }

// ------------------------------------------------------------------- //
std::tuple<double, EM>
  mk_dnj(const EM xnT, const EM dn) {
    EM dnj;

    auto [det, jacobi_inv] = mk_inv_jacobi(xnT, dn);
    dnj = dn * jacobi_inv;
    return {det, dnj};
  }

std::tuple<double, EM3>
  mk_inv_jacobi(const EM xnT, const EM dn) {
    EM3 jacobi_inv;

    auto [det, jacobi] = mk_jacobi(xnT, dn);
    // jacobi_inv <<  jacobi(1,1), -jacobi(0,1),
    //               -jacobi(1,0),  jacobi(0,0);
    // jacobi_inv /= det;
    jacobi_inv = jacobi.inverse();
    return {det, jacobi_inv};
  }


std::tuple<double, EM3>
  mk_jacobi(const EM xnT, const EM dn) {
    EM3 jacobi = xnT * dn;

    // double det = jacobi(0,0)*jacobi(1,1) - jacobi(0,1)*jacobi(1,0);
    double det = jacobi.determinant();
    return {det, jacobi};
  }
