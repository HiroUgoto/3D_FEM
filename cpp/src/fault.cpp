#include "all.h"
#include <Eigen/Dense>
#include "node.h"
#include "material.h"
#include "element_style.h"
#include "element.h"
#include "fault.h"

using EV = Eigen::VectorXd ;
using EV3 = Eigen::Vector3d ;
using EM = Eigen::MatrixXd ;
using EM2 = Eigen::Matrix2d ;
using EM3 = Eigen::Matrix3d ;


// ------------------------------------------------------------------- //
Fault::Fault (size_t id, size_t pelem_id, size_t melem_id, std::vector<size_t> spring_id, std::vector<double> param)
  {
    this->id = id;
    this->pelem_id = pelem_id;
    this->melem_id = melem_id;
    this->spring_id = spring_id;
    this->param = param;

    this->set_param();
  }

// ------------------------------------------------------------------- //
void Fault::set_param() {
    this->p0 = this->param.at(0);
    this->tp = this->param.at(1);
    this->tr = this->param.at(2);
    this->dc = this->param.at(3);
  }

// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //
void Fault::set_initial_condition0(std::vector<Element>& elements) {
  ElementStyle* estyle_p;
  EV n;

  this->pelement_p = &elements[this->pelem_id];
  this->melement_p = &elements[this->melem_id];

  estyle_p = set_element_style(this->pelement_p->style);
  n = estyle_p->shape_function_n(0.0,0.0,0.0);
  this->Np = mk_n(this->pelement_p->dof,this->pelement_p->nnode,n);
  this->NTp = this->pelement_p->mk_T_init();

  estyle_p = set_element_style(this->melement_p->style);
  n = estyle_p->shape_function_n(0.0,0.0,0.0);
  this->Nm = mk_n(this->melement_p->dof,this->melement_p->nnode,n);
  this->NTm = this->melement_p->mk_T_init();

  this->find_neighbour_element(elements);
  this->set_R();
  this->slip = 0.0;

  this->traction_force = 0.0;
  // this->rupture = true;
  // this->set_spring_kv(elements);
  this->rupture = false;
  this->set_spring_kvkh(elements);

  this->set_spring_c(elements);

  delete estyle_p;
}

void Fault::set_initial_condition1(std::vector<Element>& elements) {
  if (this->p0 > this->tp) {
    this->traction_force = this->tp - this->p0;
    this->rupture = true;
    this->set_spring_kv(elements);
  } else {
    // this->rupture = false;
    // this->set_spring_kvkh(elements);
  }
}

// ------------------------------------------------------------------- //
void Fault::find_neighbour_element(std::vector<Element>& elements) {
  ElementStyle* estyle_p = set_element_style(this->pelement_p->style);
  EV n = estyle_p->shape_function_n(0.0,0.0,0.0);
  EV3 xc = this->pelement_p->xnT*n;
  EM dn, DB;

  for (auto& element : elements) {
    if (element.style.find("3d") != std::string::npos) {
      auto [is_inside, xi] = element.check_inside(xc,0.01);
      if (is_inside) {
        this->neighbour_elements_id.push_back(element.id);

        ElementStyle* estyle_p1 = set_element_style(element.style);
        dn = estyle_p1->shape_function_dn(xi(0),xi(1),xi(2));
        DB = element.calc_stress_xi_init(dn);
        this->neighbour_elements_DB.push_back(DB);

        delete estyle_p1;
      }
    }
  }

  delete estyle_p;
}

// ------------------------------------------------------------------- //
void Fault::set_R() {
  EM dn, t(3,2);
  EV3 t0, t1, n;
  double det;

  t = this->pelement_p->xnT * this->pelement_p->dn_center;
  t0 = t.col(0); t1 = t.col(1);
  n = t0.cross(t1);
  det = n.norm();

  t0 /= t0.norm();
  t1 /= t1.norm();
  n /= det;

  this->R = EM::Zero(3,3);
  this->R(0,0) =  n(0); this->R(0,1) =  n(1); this->R(0,2) =  n(2);
  this->R(1,0) = t0(0); this->R(1,1) = t0(1); this->R(1,2) = t0(2);
  this->R(2,0) = t1(0); this->R(2,1) = t1(1); this->R(2,2) = t1(2);
}

// ------------------------------------------------------------------- //
void Fault::update_time_fault() {
  EV3 Tinput = EV::Zero(3);
  Tinput(1) = this->traction_force;
  EV3 T = this->R.transpose() * Tinput;
  this->pelement_p->mk_T(this->NTp,-T);
  this->melement_p->mk_T(this->NTm, T);
}

// ------------------------------------------------------------------- //
void Fault::update_friction(const double dt) {
  this->calc_average_slip(dt);
  double f = 0.0;
  if (this->rupture) {
    if (this->slip > 0.0) {
      if (this->slip < this->dc) {
        f = this->tp - this->slip*(this->tp - this->tr)/this->dc - this->p0;
      } else {
        f = this->tr - this->p0;
      }
    }
  }

  this->traction_force = f;
}

// ------------------------------------------------------------------- //
void Fault::calc_average_slip(const double dt) {
  EV3 up, um;
  double slip;

  up = this->Np * this->pelement_p->mk_u_hstack();
  um = this->Nm * this->melement_p->mk_u_hstack();

  slip = (this->R * (up-um))(1);
  this->sliprate = (slip-this->slip)/dt;
  this->slip = slip;
}

// ------------------------------------------------------------------- //
void Fault::calc_traction() {
  EV stress, traction;

  this->traction = 0.0;
  size_t num = this->neighbour_elements_id.size();
  for (size_t i=0 ; i<num ; i++) {
      stress = this->neighbour_elements_p[i]->calc_stress_xi(this->neighbour_elements_DB[i]);
      traction = this->stress_to_traction(stress,this->R.row(0));
      this->traction += traction(1);
  }
  this->traction /= num;

}

// ------------------------------------------------------------------- //
void Fault::update_rupture(std::vector<Element>& elements) {
  if (!this->rupture) {
    double t = this->traction + this->p0;
    if (t > this->tp) {
      this->rupture = true;
      this->set_spring_kv(elements);
    }
  }
}

// ------------------------------------------------------------------- //
EV Fault::stress_to_traction(const EV& stress, const EV& n) {
  EM stress_mat = EM::Zero(3,3);

  stress_mat(0,0) = stress(0); stress_mat(0,1) = stress(3); stress_mat(0,2) = stress(5);
  stress_mat(1,0) = stress(3); stress_mat(1,1) = stress(1); stress_mat(1,2) = stress(4);
  stress_mat(2,0) = stress(5); stress_mat(2,1) = stress(4); stress_mat(2,2) = stress(2);

  EV traction = stress_mat * n;
  return traction;
}

// ------------------------------------------------------------------- //
void Fault::set_spring_kv(std::vector<Element>& elements) {
  for (auto& id : this->spring_id) {
    EM D = elements[id].material.mk_d_slider();
    EM R_spring = EM::Zero(6,6);

    R_spring.block(0,0,3,3) = this->R;
    R_spring.block(3,3,3,3) = this->R;
    elements[id].K = R_spring.transpose() * D * R_spring;
  }
}

void Fault::set_spring_kvkh(std::vector<Element>& elements) {
  for (auto& id : this->spring_id) {
    EM D = elements[id].material.mk_d_spring();
    EM R_spring = EM::Zero(6,6);

    R_spring.block(0,0,3,3) = this->R;
    R_spring.block(3,3,3,3) = this->R;
    elements[id].K = R_spring.transpose() * D * R_spring;
  }
}

void Fault::set_spring_c(std::vector<Element>& elements) {
  for (auto& id : this->spring_id) {
    EM C = elements[id].K * 1.e-5;        // reduce oscillation
    elements[id].C_diag = C.diagonal();
    elements[id].C_off_diag = elements[id].C_diag.asDiagonal();
    elements[id].C_off_diag = C - elements[id].C_off_diag;
  }
}
