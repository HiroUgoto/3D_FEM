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

  auto [tim, dt] = input_wave::linspace(0,duration,(int)(fsamp*duration));
  size_t ntim = tim.size();
  EV slip_rate(ntim);
  slip_rate = input_wave::ricker(tim,fp,1.0/fp,1.0);


  double strike = 270.0;
  double dip = 30.0;
  double rake = 90.0;

  double length = 1000.0;
  double width = 1000.0;
  double sx = 2500.0;
  double sy = 2500.0;
  double sz = 2500.0;

  auto sources = set_source(fem.elements,strike,dip,rake,length,width,sx,sy,sz,2,2);

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

  EM output_velx(ntim,fem.output_nnode);
  EM output_vely(ntim,fem.output_nnode);
  EM output_velz(ntim,fem.output_nnode);

  // ----- time iteration ----- //
  double slip0 = 0.0;

  for (size_t it = 0 ; it < ntim ; it++) {
    slip0 += slip_rate[it]*dt;

    fem.update_time_source(sources,slip0);

    for (size_t i = 0 ; i < fem.output_nnode ; i++) {
      Node* node_p = fem.output_nodes_p[i];
      output_velx(it,i) = node_p->v(0);
      output_vely(it,i) = node_p->v(1);
      output_velz(it,i) = node_p->v(2);
    }

    if (it%40 == 0) {
      std::cout << it << " t= " << it*dt << " ";
      std::cout << output_velx(it,0) << "\n";
    }
  }

  clock_t end = clock();
  std::cout << "elapsed_time: " << (double)(end - start) / CLOCKS_PER_SEC << "[sec]\n";

  // --- Write output file --- //
  std::ofstream fvx(output_dir + "output_x.vel");
  std::ofstream fvy(output_dir + "output_y.vel");
  std::ofstream fvz(output_dir + "output_z.vel");
  for (size_t it = 0 ; it < ntim ; it++) {
    fvx << tim(it) ;
    fvy << tim(it) ;
    fvz << tim(it) ;
    for (size_t i = 0 ; i < fem.output_nnode ; i++) {
      fvx << " " << output_velx(it,i);
      fvy << " " << output_vely(it,i);
      fvz << " " << output_velz(it,i);
    }
    fvx << "\n";
    fvy << "\n";
    fvz << "\n";
  }
  fvx.close();
  fvy.close();
  fvz.close();
}
