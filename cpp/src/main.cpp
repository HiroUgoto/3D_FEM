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
  size_t fsamp = 500;
  double duration = 4.0;

  auto [tim, dt] = input_wave::linspace(0,duration,(int)(fsamp*duration));
  size_t ntim = tim.size();

  // ----- Fault setup ----- //
  fem.set_initial_fault();

  // ----- Prepare time solver ----- //
  fem.update_init(dt);

  EM output_velx(ntim,fem.output_nnode);
  EM output_vely(ntim,fem.output_nnode);
  EM output_velz(ntim,fem.output_nnode);

  EM output_slip(ntim,fem.output_nfault);
  EM output_sliprate(ntim,fem.output_nfault);
  EM output_traction(ntim,fem.output_nfault);

  // ----- time iteration ----- //
  for (size_t it = 0 ; it < ntim ; it++) {
    fem.update_time_dynamic_fault(tim(it));

    for (size_t i = 0 ; i < fem.output_nnode ; i++) {
      Node* node_p = fem.output_nodes_p[i];
      output_velx(it,i) = node_p->v(0);
      output_vely(it,i) = node_p->v(1);
      output_velz(it,i) = node_p->v(2);
    }

    for (size_t i = 0 ; i < fem.output_nfault ; i++) {
      Fault* fault_p = fem.output_faults_p[i];
      output_slip(it,i) = fault_p->slip;
      output_sliprate(it,i) = fault_p->sliprate;
      output_traction(it,i) = fault_p->traction + fault_p->p0;
    }

    if (it%100 == 0) {
      std::cout << it << " t= " << tim(it) << std::endl;
      std::cout << "  sliprate: " ;
      std::cout << output_sliprate(it,0) << " ";
      std::cout << output_sliprate(it,1) << " ";
      std::cout << output_sliprate(it,2) << " ";
      std::cout << output_sliprate(it,3) << std::endl;

      std::cout << "  traction: " ;
      std::cout << output_traction(it,0) << " ";
      std::cout << output_traction(it,1) << " ";
      std::cout << output_traction(it,2) << " ";
      std::cout << output_traction(it,3) << std::endl;
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

  std::ofstream fs(output_dir + "output_slip.dat");
  std::ofstream fsr(output_dir + "output_sliprate.dat");
  std::ofstream ft(output_dir + "output_traction.dat");
  for (size_t it = 0 ; it < ntim ; it++) {
    fs  << tim(it) ;
    fsr << tim(it) ;
    ft  << tim(it) ;
    for (size_t i = 0 ; i < fem.output_nfault ; i++) {
      fs  << " " << output_slip(it,i);
      fsr << " " << output_sliprate(it,i);
      ft  << " " << output_traction(it,i);
    }
    fs  << "\n";
    fsr << "\n";
    ft  << "\n";
  }
  fs.close();
  fsr.close();
  ft.close();
}
