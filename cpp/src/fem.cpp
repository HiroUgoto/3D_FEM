#include "all.h"
#include <Eigen/Core>
#include "node.h"
#include "material.h"
#include "element_style.h"
#include "element.h"
#include "source.h"
#include "fem.h"

using EV = Eigen::VectorXd ;

// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //
Fem::Fem (size_t dof, std::vector<Node> nodes,
                std::vector<Element> elements,
                std::vector<Material> materials)
  {
    this->nnode = nodes.size();
    this->nelem = elements.size();
    this->dof = dof;

    this->nodes = nodes;
    this->elements = elements;
    this->materials = materials;
  }

// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //
void Fem::set_init() {
    this->_set_mesh();
    this->_set_initial_condition();
    this->_set_initial_matrix();
  }

// ----------------------------- //
void Fem::_set_mesh() {
    for (auto& element : this->elements) {

      std::vector<Node*> nodes_p(element.nnode);
      for (size_t i = 0 ; i < element.nnode ; i++ ){
        size_t inode = element.inode[i];

        if (this->nodes[inode].id == inode) {
          nodes_p[i] = &this->nodes[inode];
        } else {
          for (auto& node : this->nodes) {
            if (node.id == inode) {
              nodes_p[i] = &node;
              break;
            }
          }
        }
      }

      element.set_nodes(nodes_p);

      Material* material_p = nullptr;
      if (this->materials[element.material_id].id == element.material_id) {
        material_p = &this->materials[element.material_id];
      } else {
        for (auto& material : this->materials) {
          if (material.id == element.material_id) {
            material_p = &material;
            break;
          }
        }
      }

      element.set_material(material_p);

      if (element.style.find("solid") != std::string::npos) {
        this->solid_elements_p.push_back(&element);
      }
      if (element.style.find("visco") != std::string::npos) {
        this->visco_elements_p.push_back(&element);
      }
      if (element.style.find("input") != std::string::npos) {
        this->input_elements_p.push_back(&element);
      }
      if (element.style.find("connect") != std::string::npos) {
        this->connected_elements_p.push_back(&element);
      }
    }
  }

// ----------------------------- //
void Fem::_set_initial_condition() {
    for (auto& node : this->nodes) {
      node.set_initial_condition();

      size_t total_dof = std::accumulate(node.freedom.begin(),node.freedom.end(),0);
      if (total_dof == dof) {
        this->free_nodes_p.push_back(&node);
      } else {
        this->fixed_nodes_p.push_back(&node);
      }
    }
  }

// ----------------------------- //
void Fem::_set_initial_matrix(){
    for (auto& element : this->elements) {
      element.set_xn();
      element.mk_local_matrix_init(this->dof);
      element.mk_local_matrix();
      element.mk_local_vector();

      size_t id = 0;
      for (size_t inode = 0 ; inode < element.nnode ; inode++) {
        for (size_t i = 0 ; i < this->dof ; i++) {
          element.nodes_p[inode]->mass[i] += element.M_diag[id];
          element.nodes_p[inode]->c[i] += element.C_diag[id];
          id++;
        }
      }

      // if (element.style.find("input") != std::string::npos ||
      //     element.style.find("visco") != std::string::npos) {
      //   size_t id = 0;
      //   for (size_t inode = 0 ; inode < element.nnode ; inode++) {
      //     for (size_t i = 0 ; i < this->dof ; i++) {
      //       element.nodes_p[inode]->c[i] += element.C_diag[id];
      //       id++;
      //     }
      //   }
      // }
    }
  }

// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //
void Fem::set_output(std::tuple<std::vector<size_t>, std::vector<size_t>> outputs) {
    auto [output_node_list, output_element_list] = outputs;

    this->output_nnode = output_node_list.size();
    for (size_t inode = 0 ; inode < this->output_nnode ; inode++) {
      size_t id = output_node_list[inode];
      this->output_nodes_p.push_back(&this->nodes[id]);
    }

    this->output_nelem = output_element_list.size();
    for (size_t ielem = 0 ; ielem < this->output_nelem ; ielem++) {
      size_t id = output_element_list[ielem];
      this->output_elements_p.push_back(&this->elements[id]);
    }
  }

// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //
void Fem::update_init(const double dt) {
    this->dt = dt;
    this->inv_dt2 = 1.0/(2.0*dt);
    this->inv_dtdt = 1.0/(dt*dt);
  }

// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //
void Fem::update_time(const EV acc0, const EV vel0, const bool input_wave) {
    if(input_wave) {
      this->update_time_input(vel0);
    } else {
      exit(1);
    }
  }

// ------------------------------------------------------------------- //
void Fem::update_time_input(const EV vel0) {
    for (auto& node : this->nodes) {
      node.force = EV::Zero(this->dof);
    }

    for (auto& element_p : this->input_elements_p) {
      element_p->update_inputwave(vel0);
    }

    for (auto& element : this->elements) {
      element.mk_ku_cv();
    }

    this->_update_time_set_free_nodes();
    this->_update_time_set_fixed_nodes();
    this->_update_time_set_connected_elements();

    for (auto& element_p : this->output_elements_p) {
      element_p->calc_stress();
    }
  }

// ------------------------------------------------------------------- //
void Fem::update_time_source(const std::vector<Source> sources, const double slip0) {
    for (auto& node : this->nodes) {
      node.force = EV::Zero(this->dof);
    }

    this->_update_time_source(sources,slip0);

    // for (auto& element : this->elements) {
    //   element.mk_ku_cv();
    // }

    for (auto& element_p : this->solid_elements_p) {
      element_p->mk_ku();
    }
    for (auto& element_p : this->visco_elements_p) {
      element_p->mk_ku_cv();
    }

    this->_update_time_set_free_nodes();
    this->_update_time_set_fixed_nodes();
    this->_update_time_set_connected_elements();

    for (auto& element_p : this->output_elements_p) {
      element_p->calc_stress();
    }
  }

void Fem::_update_time_source(const std::vector<Source> sources, const double slip0) {
  for (auto& source : sources){
    size_t id = source.element_id;

    this->elements[id].mk_source(source.dn,source.strain_tensor,slip0);

    for (size_t inode = 0 ; inode < this->elements[id].nnode ; inode++){
      size_t i0 = inode*this->dof;
      for (size_t i = 0 ; i < this->dof ; i++) {
        this->elements[id].nodes_p[inode]->force(i) -= this->elements[id].force[i0+i];
      }
    }
  }
}

// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //
void Fem::_update_time_set_free_nodes() {
    for (auto& node_p : this->free_nodes_p) {
      EV u = node_p->u;
      for (size_t i = 0 ; i < node_p->dof ; i++) {
        double inv_mc = 1.0 / (node_p->mass[i] + 0.5*this->dt*node_p->c[i]);
        double mass_inv_mc = node_p->mass[i] * inv_mc;
        double c_inv_mc = node_p->c[i] * inv_mc * 0.5*this->dt;
        double dtdt_inv_mc = this->dt*this->dt*inv_mc;

        node_p->u[i] = mass_inv_mc*(2.0*u[i] - node_p->um[i])
                + c_inv_mc*node_p->um[i] - dtdt_inv_mc*node_p->force[i];

      }
      node_p->v = (node_p->u - node_p->um) * this->inv_dt2;
      node_p->um = u;
    }
  }

void Fem::_update_time_set_fixed_nodes() {
    for (auto& node_p : this->fixed_nodes_p) {
      EV u = node_p->u;
      for (size_t i = 0 ; i < node_p->dof ; i++) {
        if (node_p->freedom[i] == 0) {
          node_p->u[i] = 0.0;
        } else {
          double inv_mc = 1.0 / (node_p->mass[i] + 0.5*this->dt*node_p->c[i]);
          double mass_inv_mc = node_p->mass[i] * inv_mc;
          double c_inv_mc = node_p->c[i] * inv_mc * 0.5*this->dt;
          double dtdt_inv_mc = this->dt*this->dt*inv_mc;

          node_p->u[i] = mass_inv_mc*(2.0*u[i] - node_p->um[i])
                  + c_inv_mc*node_p->um[i] - dtdt_inv_mc*node_p->force[i];
        }
      }
      node_p->v = (node_p->u - node_p->um) * this->inv_dt2;
      node_p->um = u;
    }
  }

void Fem::_update_time_set_connected_elements() {
    for (auto& element_p : this->connected_elements_p) {
      EV u = EV::Zero(element_p->dof);
      for (size_t inode = 0 ; inode < element_p->nnode ; inode++) {
        u += element_p->nodes_p[inode]->u;
      }
      for (size_t inode = 0 ; inode < element_p->nnode ; inode++) {
        element_p->nodes_p[inode]->u = u/element_p->nnode;
      }
    }
  }
