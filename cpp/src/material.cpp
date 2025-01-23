#include "all.h"
#include <Eigen/Core>
#include "material.h"

using EV = Eigen::VectorXd ;
using EM = Eigen::MatrixXd ;
using EM3 = Eigen::Matrix3d ;

Material::Material () {}
Material::Material (size_t id, std::string style, std::vector<double> param) {
  this->id = id;
  this->style = style;
  this->param = param;

  this->set_param();
}

void Material::set_init(size_t id, std::string style, std::vector<double> param) {
  this->id = id;
  this->style = style;
  this->param = param;

  this->set_param();
}

void Material::print() {
  std::cout << this->id << ": ";
  std::cout << this->style << ", ";
  for (size_t i = 0 ; i < this->param.size() ; i++) {
    std::cout << this->param.at(i) << " ";
  }
  std::cout << "\n";
}

void Material::set_param() {
  if (this->style == "vs_vp_rho") {
    double vs = this->param.at(0);
    double vp = this->param.at(1);
    double rho = this->param.at(2);

    this->rmu = rho*vs*vs;
    this->rlambda = rho*vp*vp - 2.0*this->rmu;
    this->rho = rho;

  } else if (this->style == "nu_vp_rho") {
    double nu = this->param.at(0);
    double vp = this->param.at(1);
    double rho = this->param.at(2);

    this->rmu = rho/2.0*vp*vp*(1.0-2.0*nu)/(1.0-nu);
    this->rlambda = rho*nu*vp*vp/(1.0-nu);
    this->rho = rho;

  } else if (this->style == "nu_vs_rho") {
    double nu = this->param.at(0);
    double vs = this->param.at(1);
    double rho = this->param.at(2);

    this->rmu = rho*vs*vs;
    this->rlambda = 2.0*nu/(1.0-2.0*nu) * this->rmu;
    this->rho = rho;

  } else if (this->style == "vs_vp_rho_Q") {
    double vs = this->param.at(0);
    double vp = this->param.at(1);
    double rho = this->param.at(2);
    double Q = this->param.at(3);
    double fmax = this->param.at(4);

    this->rmu = rho*vs*vs;
    this->rlambda = rho*vp*vp - 2.0*this->rmu;
    this->rho = rho; 

    double wmax = 2.0*M_PI*fmax;
    this->beta = 3.0 / (2.0*Q*wmax);
  }
}

// ------------------------------------------------------------------- //
EM Material::mk_d(const size_t dof) {
    EM D(6,6);

    D = EM::Zero(6,6);
    D(0,0) = this->rlambda + 2.0*this->rmu;
    D(0,1) = this->rlambda;
    D(0,2) = this->rlambda;

    D(1,0) = this->rlambda;
    D(1,1) = this->rlambda + 2.0*this->rmu;
    D(1,2) = this->rlambda;

    D(2,0) = this->rlambda;
    D(2,1) = this->rlambda;
    D(2,2) = this->rlambda + 2.0*this->rmu;

    D(3,3) = this->rmu;
    D(4,4) = this->rmu;
    D(5,5) = this->rmu;

    return D;
  }

// ------------------------------------------------------------------- //
EM Material::mk_visco(const size_t dof) {
    EM D(6,6);
    D = EM::Zero(6,6);

    if (this->style.find("Q") != std::string::npos) {
      double d0 =  4.0/3.0*this->rmu;
      double d1 = -2.0/3.0*this->rmu;

      D(0,0) = d0;
      D(0,1) = d1;
      D(0,2) = d1;

      D(1,0) = d1;
      D(1,1) = d0;
      D(1,2) = d1;

      D(2,0) = d1;
      D(2,1) = d1;
      D(2,2) = d0;

      D(3,3) = this->rmu;
      D(4,4) = this->rmu;
      D(5,5) = this->rmu;

      D *= this->beta;
    }

    return D;
  }

// ------------------------------------------------------------------- //
EM3 Material::mk_imp(const size_t dof) {
    EM3 imp;
    double vs = sqrt(this->rmu/this->rho);
    double vp = sqrt((this->rlambda + 2.0*this->rmu)/this->rho);

    imp = EM::Zero(3,3);
    imp(0,0) = this->rho * vp;
    imp(1,1) = this->rho * vs;
    imp(2,2) = this->rho * vs;

    return imp;
  }
