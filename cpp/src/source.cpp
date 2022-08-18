#include "all.h"
#include <Eigen/Core>
#include "node.h"
#include "material.h"
#include "element_style.h"
#include "element.h"
#include "source.h"

using EV = Eigen::VectorXd ;
using EM = Eigen::MatrixXd ;
using EV3 = Eigen::Vector3d ;

// ------------------------------------------------------------------- //
std::vector<Source> set_source(std::vector<Element> elements,
    double strike, double dip, double rake, double length, double width,
    double sx, double sy, double sz, size_t nl, size_t nw) {
  std::vector<Source> source_list;
  std::vector<double> l(nl), w(nw);

  for (size_t i=0; i<nl; i++){
    double l0, l1;
    l0 = length/nl* i    - length/2.0;
    l1 = length/nl*(i+1) - length/2.0;
    l[i] = (l0+l1)/2.0;
  }

  for (size_t i=0; i<nw; i++){
    double w0, w1;
    w0 = width/nw* i    - width/2.0;
    w1 = width/nw*(i+1) - width/2.0;
    w[i] = (w0+w1)/2.0;
  }

  double strike_rad = strike * M_PI/180.0;
  double dip_rad = dip * M_PI/180.0;
  double dl = length/nl;
  double dw = width/nw;

  EV3 x = EV::Zero(3);
  size_t id = 0;
  for (size_t i=0; i<nl; i++) {
    for (size_t j=0; j<nw; j++) {
      x[0] = sx + l[i]*std::cos(strike_rad) - w[j]*std::cos(dip_rad)*std::sin(strike_rad);
      x[1] = sy + l[i]*std::sin(strike_rad) + w[j]*std::cos(dip_rad)*std::cos(strike_rad);
      x[2] = sz + w[j]*std::sin(dip_rad);

      for (auto& element : elements) {
        if (element.dim == 3) {
          auto [is_inside,xi] = element.check_inside(x);
          if (is_inside) {
            ElementStyle* estyle_p = set_element_style(elements[element.id].style);
            EM dn = estyle_p->shape_function_dn(xi[0],xi[1],xi[2]);

            Source source(id,strike,dip,rake,dl,dw,element.id,dn);
            source_list.push_back(source);
            id += 1;
            break;
          }
        }
      }
    }
  }

  return source_list;
}

// ------------------------------------------------------------------- //
Source::Source(size_t id, double strike, double dip, double rake,
                double length, double width, size_t element_id, EM dn)
  {
    this->id = id;
    this->element_id = element_id;
    this->dn = dn;

    this->strike = strike * M_PI/180.0;
    this->dip = dip * M_PI/180.0;
    this->rake = rake * M_PI/180.0;

    this->length = length;
    this->width = width;

    this->set_strain_tensor();
  }

  void Source::set_strain_tensor(){
    this->strain_tensor = EV::Zero(6);
    double m0 = this->width * this->length;

    // Mxx
    this->strain_tensor(0) = -m0*( std::sin(this->dip)*std::cos(this->rake)*std::sin(2.0*this->strike)
                                 + std::sin(2.0*this->dip)*std::sin(this->rake)*std::pow(std::sin(this->strike),2.0));
    // Myy
    this->strain_tensor(1) =  m0*( std::sin(this->dip)*std::cos(this->rake)*std::sin(2.0*this->strike)
                                 - std::sin(2.0*this->dip)*std::sin(this->rake)*std::pow(std::cos(this->strike),2.0));
    // Mzz
    this->strain_tensor(2) =  m0*std::sin(2.0*this->dip)*std::sin(this->rake);
    // Mxy
    this->strain_tensor(3) =  m0*( std::sin(this->dip)*std::cos(this->rake)*std::cos(2.0*this->strike)
                                 + std::sin(2.0*this->dip)*std::sin(this->rake)*std::sin(2.0*this->strike)*0.5);
    // Myz
    this->strain_tensor(4) = -m0*( std::cos(this->dip)*std::cos(this->rake)*std::sin(this->strike) \
                                 - std::cos(2.0*this->dip)*std::sin(this->rake)*std::cos(this->strike));
    // Mxz
    this->strain_tensor(5) = -m0*( std::cos(this->dip)*std::cos(this->rake)*std::cos(this->strike) \
                                 + std::cos(2.0*this->dip)*std::sin(this->rake)*std::sin(this->strike));
  }
