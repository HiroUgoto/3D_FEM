#include "all.h"
#include <Eigen/Core>
#include "node.h"

using EV3 = Eigen::Vector3d ;

Node::Node () {}
Node::Node (size_t id, std::vector<double> xyz, std::vector<size_t> freedom)
  {
    this->id = id;
    this->xyz = xyz;
    this->freedom = freedom;
    this->dof = freedom.size();
  }

void Node::print() {
  std::cout << this->id << ": ";
  for (size_t i = 0 ; i < this->xyz.size() ; i++) {
    std::cout << this->xyz.at(i) << " ";
  }
  std::cout << "-- ";
  for (size_t i = 0 ; i < this->freedom.size() ; i++) {
    std::cout << this->freedom.at(i) << " ";
  }
  std::cout << "\n";
}

void Node::set_initial_condition() {
  this->u = EV3::Zero(this->dof);
  this->um = EV3::Zero(this->dof);
  this->v  = EV3::Zero(this->dof);

  this->force = EV3::Zero(this->dof);

  this->mass = EV3::Zero(this->dof);
  this->c    = EV3::Zero(this->dof);
}
