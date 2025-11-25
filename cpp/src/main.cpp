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

  // ----- Define source ----- //
  // size_t fsamp = 100;
  // double fp = 0.5;
  // double duration = 6.0;

  // auto [tim, dt] = input_wave::linspace(0,duration,(int)(fsamp*duration));
  // size_t ntim = tim.size();
  // EV wave_acc(ntim);
  // wave_acc = input_wave::ricker(tim,fp,1.0/fp,1.0);

  // ---- Read input wave //
  // auto [tim,wave_acc,dt] = input_wave::input_acc_file("input/scaled_input_acc.txt");
  // size_t ntim = tim.size();

  // double polarity = 45;  // [deg] N[XX]E
  // EV wave_accx(ntim);  
  // EV wave_accy(ntim); 

  // double polarity_rad = polarity * M_PI/180.0;
  // for (size_t it = 0 ; it < ntim ; it++) {
  //   wave_accx[it] = wave_acc[it] * std::cos(polarity_rad); 
  //   wave_accy[it] = wave_acc[it] * std::sin(polarity_rad); 
  // }

  // std::ofstream fa(output_dir + "input.acc");
  // std::ofstream fv(output_dir + "input.vel");
  // double velx = 0.0;
  // double vely = 0.0;
  // for (size_t it = 0 ; it < ntim ; it++) {
  //   fa << tim(it) ;
  //   fa << " " << wave_accx[it] ;
  //   fa << " " << wave_accy[it] ;
  //   fa << "\n";

  //   velx += wave_accx[it]*dt;
  //   vely += wave_accy[it]*dt;
  //   fv << tim(it) ;
  //   fv << " " << velx ;
  //   fv << " " << vely ;
  //   fv << "\n";
  // }
  // fa.close();
  // fv.close();
  // // exit(1);

  // ----- Define EQ source ----- //
  size_t fsamp = 1000;
  double fp = 2.0;
  double duration = 2.5;

  auto [tim, dt] = input_wave::linspace(0,duration,(int)(fsamp*duration));
  size_t ntim = tim.size();
  EV input_stf(ntim);
  input_stf = input_wave::ricker(tim,fp,1.5/fp,1.0);

  // auto [tim,input_stf,dt] = input_wave::input_acc_file("input/scaled_input_STF.txt");
  // size_t ntim = tim.size();

  double strike = 0.0;
  double dip = 45.0;
  double rake = 90.0;

  // double length = 1000.0;
  // double width = 1000.0;
  double sx = 250.0;
  double sy = 250.0;
  double sz = 250.0;
  double mw = 3.0;

  double rmu = 1000*1000*2100;
  double m0 = pow(10, 1.5*mw+9.1);

  double length = sqrt(m0/rmu);
  double width = length;

  auto sources = set_source(fem.elements,strike,dip,rake,length,width,sx,sy,sz,1,1);

  std::ofstream f0(output_dir + "input_stf.dat");
  for (size_t it = 0 ; it < ntim ; it++) {
    f0 << tim(it) ;
    f0 << " " << input_stf(it) ;
    f0 << "\n";
  }
  f0.close();
  // exit(1);

  // ----- Prepare time solver ----- //
  fem.update_init(dt);

  EM output_velx(ntim,fem.output_nnode);
  EM output_vely(ntim,fem.output_nnode);
  EM output_velz(ntim,fem.output_nnode);

  EM output_dispx(ntim,fem.output_nnode);
  EM output_dispy(ntim,fem.output_nnode);
  EM output_dispz(ntim,fem.output_nnode);

  // ----- time iteration ----- //
  // EV vel0(3);
  // vel0[0] = 0.0; vel0[1] = 0.0; vel0[2] = 0.0; 

  for (size_t it = 0 ; it < ntim ; it++) {
    // vel0[0] += wave_accx[it]*dt;
    // vel0[1] += wave_accy[it]*dt;
    // fem.update_time_input(vel0);

    fem.update_time_source(sources,input_stf[it]);

    for (size_t i = 0 ; i < fem.output_nnode ; i++) {
      Node* node_p = fem.output_nodes_p[i];
      output_velx(it,i) = node_p->v(0);
      output_vely(it,i) = node_p->v(1);
      output_velz(it,i) = node_p->v(2);
      output_dispx(it,i) = node_p->u(0);
      output_dispy(it,i) = node_p->u(1);
      output_dispz(it,i) = node_p->u(2);
    }

    if (it%40 == 0) {
      std::cout << it << " t= " << it*dt << " ";
      std::cout << output_vely(it,0) << "\n";
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

  std::ofstream fdx(output_dir + "output_x.disp");
  std::ofstream fdy(output_dir + "output_y.disp");
  std::ofstream fdz(output_dir + "output_z.disp");
  for (size_t it = 0 ; it < ntim ; it++) {
    fdx << tim(it) ;
    fdy << tim(it) ;
    fdz << tim(it) ;
    for (size_t i = 0 ; i < fem.output_nnode ; i++) {
      fdx << " " << output_dispx(it,i);
      fdy << " " << output_dispy(it,i);
      fdz << " " << output_dispz(it,i);
    }
    fdx << "\n";
    fdy << "\n";
    fdz << "\n";
  }
  fdx.close();
  fdy.close();
  fdz.close();
}
