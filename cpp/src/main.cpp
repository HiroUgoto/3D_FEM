#include "all.h"
#include <Eigen/Core>
#include "node.h"
#include "material.h"
#include "element_style.h"
#include "element.h"
#include "source.h"
#include "fem.h"
#include "io_data.h"
#include "input_wave.h"

using EV = Eigen::VectorXd ;
using EM = Eigen::MatrixXd ;

int main() {

  clock_t start = clock();

  // ----- Input FEM Mesh ----- //
  Fem fem = io_data::input_mesh("input/mesh.in");
  auto outputs = io_data::input_outputs("input/output.in");
  std::string output_dir = "result/";

  // ----- FEM Set up ----- //
  fem.set_init();
  fem.set_output(outputs);
  // exit(1);

  // ----- Define input wave ----- //
  size_t fsamp = 100;

  double fp = 0.5;
  double duration = 6.0;

  EV slip_rate;
  auto [tim, dt] = input_wave::linspace(0,duration,(int)(fsamp*duration));
  slip_rate = input_wave::ricker(tim,fp,1.0/fp,1.0);
  size_t ntim = tim.size();


  double strike = 270.0;
  double dip = 30.0;
  double rake = 90.0;

  double length = 1000.0;
  double width = 1000.0;
  double sx = 2500.0;
  double sy = 2500.0;
  double sz = 2500.0;

  auto sources = set_source(fem.elements,strike,dip,rake,length,width,sx,sy,sz,5,5);

  // std::ofstream f0(output_dir + "input.acc");
  // for (size_t it = 0 ; it < ntim ; it++) {
  //   f0 << tim(it) ;
  //   f0 << " " << wave_acc(it) ;
  //   f0 << "\n";
  // }
  // f0.close();
  // exit(1);

  // ----- Prepare time solver ----- //
  fem.update_init(dt);

  EM output_dispx = EM::Zero(ntim,fem.output_nnode);
  EM output_dispz = EM::Zero(ntim,fem.output_nnode);
  EM output_velx = EM::Zero(ntim,fem.output_nnode);
  EM output_velz = EM::Zero(ntim,fem.output_nnode);
  EM output_accx = EM::Zero(ntim,fem.output_nnode);
  EM output_accz = EM::Zero(ntim,fem.output_nnode);

  // ----- time iteration ----- //
  double slip0 = 0.0;

  for (size_t it = 0 ; it < ntim ; it++) {
    slip0 += slip_rate[it]*dt;

    fem.update_time_source(sources,slip0);

    for (size_t i = 0 ; i < fem.output_nnode ; i++) {
      Node* node_p = fem.output_nodes_p[i];
      output_dispx(it,i) = node_p->u(0);
      output_dispz(it,i) = node_p->u(1);
      output_velx(it,i) = node_p->v(0);
      output_velz(it,i) = node_p->v(1);
      output_accx(it,i) = node_p->a(0);
      output_accz(it,i) = node_p->a(1);
    }

    if (it%20 == 0) {
      std::cout << it << " t= " << it*dt << " ";
      std::cout << output_dispx(it,0) << "\n";
    }
  }

  clock_t end = clock();
  std::cout << "elapsed_time: " << (double)(end - start) / CLOCKS_PER_SEC << "[sec]\n";

  // --- Write output file --- //
  std::ofstream fd(output_dir + "output_x.disp");
  for (size_t it = 0 ; it < ntim ; it++) {
    fd << tim(it) ;
    for (size_t i = 0 ; i < fem.output_nnode ; i++) {
      fd << " " << output_dispx(it,i);
    }
    fd << "\n";
  }
  fd.close();

  std::ofstream fv(output_dir + "output_x.vel");
  for (size_t it = 0 ; it < ntim ; it++) {
    fv << tim(it) ;
    for (size_t i = 0 ; i < fem.output_nnode ; i++) {
      fv << " " << output_velx(it,i);
    }
    fv << "\n";
  }
  fv.close();

  std::ofstream fa(output_dir + "output_x.acc");
  for (size_t it = 0 ; it < ntim ; it++) {
    fa << tim(it) ;
    for (size_t i = 0 ; i < fem.output_nnode ; i++) {
      fa << " " << output_accx(it,i);
    }
    fa << "\n";
  }
  fa.close();

}
