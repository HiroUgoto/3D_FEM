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

  // ----- Define source ----- //
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

  // std::ofstream f0(output_dir + "slip_rate.dat");
  // for (size_t it = 0 ; it < ntim ; it++) {
  //   f0 << tim(it) ;
  //   f0 << " " << slip_rate(it) ;
  //   f0 << "\n";
  // }
  // f0.close();
  // exit(1);

  // ----- Prepare time solver ----- //
  fem.update_init(dt);

  EM output_dispx = EM::Zero(ntim,fem.output_nnode);
  EM output_dispy = EM::Zero(ntim,fem.output_nnode);
  EM output_dispz = EM::Zero(ntim,fem.output_nnode);
  EM output_velx = EM::Zero(ntim,fem.output_nnode);
  EM output_vely = EM::Zero(ntim,fem.output_nnode);
  EM output_velz = EM::Zero(ntim,fem.output_nnode);
  EM output_accx = EM::Zero(ntim,fem.output_nnode);
  EM output_accy = EM::Zero(ntim,fem.output_nnode);
  EM output_accz = EM::Zero(ntim,fem.output_nnode);

  // ----- time iteration ----- //
  double slip0 = 0.0;

  for (size_t it = 0 ; it < ntim ; it++) {
    slip0 += slip_rate[it]*dt;

    fem.update_time_source(sources,slip0);

    for (size_t i = 0 ; i < fem.output_nnode ; i++) {
      Node* node_p = fem.output_nodes_p[i];
      output_dispx(it,i) = node_p->u(0);
      output_dispy(it,i) = node_p->u(1);
      output_dispz(it,i) = node_p->u(2);
      output_velx(it,i) = node_p->v(0);
      output_vely(it,i) = node_p->v(1);
      output_velz(it,i) = node_p->v(2);
      output_accx(it,i) = node_p->a(0);
      output_accy(it,i) = node_p->a(1);
      output_accz(it,i) = node_p->a(2);
    }

    if (it%40 == 0) {
      std::cout << it << " t= " << it*dt << " ";
      std::cout << output_dispx(it,0) << "\n";
    }
  }

  clock_t end = clock();
  std::cout << "elapsed_time: " << (double)(end - start) / CLOCKS_PER_SEC << "[sec]\n";

  // --- Write output file --- //
  std::ofstream fdx(output_dir + "output_x.disp");
  std::ofstream fdy(output_dir + "output_y.disp");
  for (size_t it = 0 ; it < ntim ; it++) {
    fdx << tim(it) ;
    fdy << tim(it) ;
    for (size_t i = 0 ; i < fem.output_nnode ; i++) {
      fdx << " " << output_dispx(it,i);
      fdy << " " << output_dispy(it,i);
    }
    fdx << "\n";
    fdy << "\n";
  }
  fdx.close();
  fdy.close();

  std::ofstream fvx(output_dir + "output_x.vel");
  std::ofstream fvy(output_dir + "output_y.vel");
  for (size_t it = 0 ; it < ntim ; it++) {
    fvx << tim(it) ;
    fvy << tim(it) ;
    for (size_t i = 0 ; i < fem.output_nnode ; i++) {
      fvx << " " << output_velx(it,i);
      fvy << " " << output_vely(it,i);
    }
    fvx << "\n";
    fvy << "\n";
  }
  fvx.close();
  fvy.close();

  std::ofstream fax(output_dir + "output_x.acc");
  std::ofstream fay(output_dir + "output_y.acc");
  for (size_t it = 0 ; it < ntim ; it++) {
    fax << tim(it) ;
    fay << tim(it) ;
    for (size_t i = 0 ; i < fem.output_nnode ; i++) {
      fax << " " << output_accx(it,i);
      fay << " " << output_accy(it,i);
    }
    fax << "\n";
    fay << "\n";
  }
  fax.close();
  fay.close();

}
