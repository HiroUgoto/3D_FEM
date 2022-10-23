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
  this->pelement_p = &elements[this->pelem_id];
  this->melement_p = &elements[this->melem_id];

  this->find_neighbour_element(elements);
  this->set_R();
  this->slip = 0.0;

  ElementStyle* estyle_p = set_element_style(this->pelement_p->style);
  std::vector<EM> dn_list = estyle_p->dn_list;
  std::vector<double> w_list = estyle_p->w_list;

  this->area = 0.0;
  for (size_t i = 0 ; i < this->pelement_p->ng_all ; i++){
    auto [det, q] = mk_q(this->pelement_p->dof, this->pelement_p->xnT, dn_list[i]);
    double detJ = det * w_list[i];
    this->area += detJ;
  }

  this->pelement_p->traction = 0.0;
  this->melement_p->traction = 0.0;
  this->rupture = false;
  this->set_spring_kvkh(elements);

  delete estyle_p;
}

void Fault::set_initial_condition1(std::vector<Element>& elements) {
  if (this->p0 > this->tp) {
    this->pelement_p->traction = this->tp - this->p0;
    this->melement_p->traction = this->tp - this->p0;
    this->rupture = true;
    this->set_spring_kv(elements);
  }
}

// ------------------------------------------------------------------- //
void Fault::find_neighbour_element(std::vector<Element> elements) {
  ElementStyle* estyle_p = set_element_style(this->pelement_p->style);
  EV n = estyle_p->shape_function_n(0.0,0.0,0.0);
  EV3 xc = this->pelement_p->xnT*n;

  for (auto& element : elements) {
    if (element.style.find("3d") != std::string::npos) {
      auto [is_inside, xi] = element.check_inside(xc,0.01);
      if (is_inside) {
        this->neighbour_elements_id.push_back(element.id);
        this->neighbour_elements_xi.push_back(xi);
      }
    }
  }

  // if (this->id == 5) {
  //   for (auto& id : this->neighbour_elements_id) {
  //     std::cout << id << " ";
  //   }
  //   std::cout << std::endl;
  //   exit(1);
  // }
  //
  delete estyle_p;
}

// ------------------------------------------------------------------- //
void Fault::set_R() {
  EM dn, t(3,2);
  EV3 t0, t1, n;
  double det;

  ElementStyle* estyle_p = set_element_style(this->pelement_p->style);
  dn = estyle_p->shape_function_dn(0.0,0.0,0.0);

  t = this->pelement_p->xnT*dn;
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

  this->pelement_p->R = this->R;
  this->melement_p->R = this->R;

  delete estyle_p;
}

// ------------------------------------------------------------------- //
void Fault::update_friction(double dt) {
  this->calc_average_slip(dt);
  double f = 0.0;
  if (this->slip < this->dc) {
    f = this->tp - this->slip*(this->tp - this->tr)/this->dc - this->p0;
  } else {
    f = this->tr - this->p0;
  }

  this->pelement_p->traction = f;
  this->melement_p->traction = f;
}

// ------------------------------------------------------------------- //
void Fault::calc_average_slip(double dt) {
  ElementStyle* estyle_p;
  EV n;
  EM Np, Nm;
  EV3 up, um;
  double slip;

  estyle_p = set_element_style(this->pelement_p->style);
  n = estyle_p->shape_function_n(0.0,0.0,0.0);
  Np = mk_n(this->pelement_p->dof,this->pelement_p->nnode,n);
  up = Np * this->pelement_p->mk_u_hstack();

  estyle_p = set_element_style(this->melement_p->style);
  n = estyle_p->shape_function_n(0.0,0.0,0.0);
  Nm = mk_n(this->melement_p->dof,this->melement_p->nnode,n);
  um = Nm * this->melement_p->mk_u_hstack();

  slip = (this->R * (up-um))(1);
  this->sliprate = (slip-this->slip)/dt;
  this->slip = slip;

  delete estyle_p;
}

// ------------------------------------------------------------------- //
void Fault::calc_traction(std::vector<Element> elements) {
  this->traction = 0.0;
  size_t num = this->neighbour_elements_id.size();

  for (size_t i=0 ; i<num ; i++) {
      size_t id = this->neighbour_elements_id[i];
      EV xi = this->neighbour_elements_xi[i];
      EV stress = elements[id].calc_stress_xi(xi);

      EV traction = this->stress_to_traction(stress,this->R.row(0));
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
EV Fault::stress_to_traction(EV stress, EV n) {
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