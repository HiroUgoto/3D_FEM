#include "all.h"
#include <Eigen/Core>
#include "node.h"
#include "material.h"
#include "element_style.h"
#include "element.h"
#include "fault.h"
#include "source.h"
#include "fem.h"
#include "io_data.h"

// ------------------------------------------------------------------- //
Fem io_data::input_mesh (const std::string mesh_file) {
    size_t nnode, nelem, nfault, nmaterial, dof;
    std::string line;

    // Open file //
    std::ifstream f(mesh_file);

    // Read header //
    std::getline(f, line);
    std::istringstream iss(line);
    iss >> nnode >> nelem >> nfault >> nmaterial >> dof;

    // Read nodes //
    std::vector<Node> nodes;
    for (size_t inode = 0 ; inode < nnode ; inode++) {
      size_t id;
      std::vector<double> xyz(3);
      std::vector<size_t> freedom(dof);

      std::getline(f, line);
      if (line[line.size()-1] == '\n') line.erase(line.size()-1);
      if (line[line.size()-1] == '\r') line.erase(line.size()-1);
      // std::cout << line + "\n";

      size_t s = line.find_first_not_of(' ');
      size_t e = line.find_last_not_of(' ');
      std::istringstream iss(line.substr(s,e-s+1));
      iss >> id;
      for (size_t i = 0 ; i < 3 ; i++) {
        iss >> xyz.at(i);
      }
      for (size_t i = 0 ; i < dof ; i++) {
        iss >> freedom.at(i);
      }

      Node node(id,xyz,freedom);
      nodes.push_back(node);
    }

    // Read elements //
    std::vector<Element> elements;
    for (size_t ielem = 0 ; ielem < nelem ; ielem++) {
      size_t id;
      int material_id;
      std::string style;
      std::vector<size_t> inode;

      std::getline(f, line);
      if (line[line.size()-1] == '\n') line.erase(line.size()-1);
      if (line[line.size()-1] == '\r') line.erase(line.size()-1);
      // std::cout << line + "\n";

      size_t s = line.find_first_not_of(' ');
      size_t e = line.find_last_not_of(' ');
      std::istringstream iss(line.substr(s,e-s+1));
      iss >> id >> style >> material_id ;
      while(!iss.eof()) {
        size_t in;
        iss >> in;
        inode.push_back(in);
      }
      Element element(id,style,material_id,inode);
      elements.push_back(element);
    }

    // Read faults //
    std::vector<Fault> faults;
    for (size_t ifault = 0 ; ifault < nfault ; ifault++) {
      size_t id;
      size_t pelem_id, melem_id;
      std::vector<size_t> spring_id;
      std::vector<double> param;

      std::getline(f, line);
      if (line[line.size()-1] == '\n') line.erase(line.size()-1);
      if (line[line.size()-1] == '\r') line.erase(line.size()-1);
      // std::cout << line + "\n";

      size_t s = line.find_first_not_of(' ');
      size_t e = line.find_last_not_of(' ');
      std::istringstream iss(line.substr(s,e-s+1));
      iss >> id >> pelem_id >> melem_id ;
      for (size_t i=0 ; i<4 ; i++) {
        size_t sid;
        iss >> sid;
        spring_id.push_back(sid);
      }
      while(!iss.eof()) {
        double ip;
        iss >> ip;
        param.push_back(ip);
      }

      Fault fault(id,pelem_id,melem_id,spring_id,param);
      faults.push_back(fault);
    }

    // Read materials //
    std::vector<Material> materials;
    for (size_t imaterial = 0 ; imaterial < nmaterial ; imaterial++) {
      size_t id;
      std::string style;
      std::vector<double> param;

      std::getline(f, line);
      if (line[line.size()-1] == '\n') line.erase(line.size()-1);
      if (line[line.size()-1] == '\r') line.erase(line.size()-1);
      // std::cout << line + "\n";

      size_t s = line.find_first_not_of(' ');
      size_t e = line.find_last_not_of(' ');
      std::istringstream iss(line.substr(s,e-s+1));
      iss >> id >> style ;
      while(!iss.eof()) {
        double ip;
        iss >> ip;
        param.push_back(ip);
      }

      Material material(id,style,param);
      materials.push_back(material);
    }

    Fem fem(dof,nodes,elements,faults,materials);
    return fem;
  }

// ------------------------------------------------------------------- //
std::tuple<std::vector<size_t>, std::vector<size_t>, std::vector<size_t>>
  io_data::input_outputs (const std::string output_file) {
    size_t nnode, nelem, nfault;
    std::string line;

    // Open file //
    std::ifstream f(output_file);

    // Read header //
    std::getline(f, line);
    std::istringstream iss(line);
    iss >> nnode >> nelem >> nfault;

    // Read nodes //
    std::vector<size_t> nodes;
    for (size_t inode = 0 ; inode < nnode ; inode++) {
      size_t id;

      std::getline(f, line);
      // std::cout << line + "\n";
      std::istringstream iss(line);
      iss >> id;

      nodes.push_back(id);
    }

    std::vector<size_t> elements;
    for (size_t ielem = 0 ; ielem < nelem ; ielem++) {
      size_t id;

      std::getline(f, line);
      // std::cout << line + "\n";
      std::istringstream iss(line);
      iss >> id;

      elements.push_back(id);
    }

    std::vector<size_t> faults;
    for (size_t ifault = 0 ; ifault < nfault ; ifault++) {
      size_t id;

      std::getline(f, line);
      // std::cout << line + "\n";
      std::istringstream iss(line);
      iss >> id;

      faults.push_back(id);
    }


    return {nodes, elements, faults};
  }
