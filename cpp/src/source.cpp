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
std::vector<Source> set_source(std::vector<Element> elements, double dip, double width, double sx, double sy, double sz, size_t n) {
  std::vector<Source> source_list;
  std::vector<double> w(n);

  for (size_t i=0; i<n; i++){
    double w0, w1;
    w0 = width/n* i    - width/2.0;
    w1 = width/n*(i+1) - width/2.0;
    w[i] = (w0+w1)/2.0;
  }

  double dip_rad = dip * M_PI/180.0;
  double dw = width/n  * width;

  EV x = EV::Zero(3);
  for (size_t i=0; i<n; i++) {
    x[0] = sx + w[i]*std::cos(dip_rad);
    x[1] = sy;
    x[2] = sz + w[i]*std::sin(dip_rad);

    for (auto& element : elements) {
      if (element.dim == 3) {
        auto [is_inside,xi] = element.check_inside(x);
        if (is_inside) {
          ElementStyle* estyle_p = set_element_style(elements[element.id].style);
          EM dn = estyle_p->shape_function_dn(xi[0],xi[1],xi[2]);

          Source source(i,dip,dw,element.id,dn);
          source_list.push_back(source);
          break;
        }
      }
    }
  }

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
