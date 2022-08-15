#include "all.h"
#include <Eigen/Core>
#include "node.h"
#include "material.h"
#include "element_style.h"
#include "element.h"
#include "source.h"

using EV = Eigen::VectorXd ;
using EM = Eigen::MatrixXd ;

// ------------------------------------------------------------------- //
std::vector<Source> set_source(const std::vector<Element> elements, const double dip, const double width, const double sx, const double sy, const double sz) {
  std::vector<Source> source_list;

  size_t element_id = 0;
  double xi, eta, zeta;
  xi = 0.0; eta = 0.0; zeta = 0.0;

  ElementStyle* estyle_p = set_element_style(elements[element_id].style);
  EM dn = estyle_p->shape_function_dn(xi,eta,zeta);

  // Source source(0,dip,width,element_id,dn);
  Source source(0,dip,width *width,element_id,dn);
  source_list.push_back(source);

  return source_list;
}

// ------------------------------------------------------------------- //
Source::Source(size_t id, double dip, double width, size_t element_id, EM dn)
  {
    this->id = id;
    this->element_id = element_id;
    this->dn = dn;

    this->dip = dip * M_PI/180.0;
    this->width = width;

    this->set_strain_tensor();
  }

  void Source::set_strain_tensor(){
    this->strain_tensor = EV::Zero(6);

    this->strain_tensor(0) = -std::sin(2.0*this->dip) * this->width;
    this->strain_tensor(2) =  std::sin(2.0*this->dip) * this->width;
    this->strain_tensor(5) =  std::cos(2.0*this->dip) * this->width;
  }
