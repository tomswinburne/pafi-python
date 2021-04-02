#ifndef CSIM_H
#define CSIM_H

#include "LAMMPSSimulator.hpp"

class LAMMPS_TI_Simulator : public LAMMPSSimulator {

public:
  LAMMPS_TI_Simulator(MPI_Comm &instance_comm, Parser &p, Holder &h, int t)
  : LAMMPSSimulator(instance_comm, p, h, t) {
    /*
      anything extra to initialize ?
    */
  };

  void constrained_average() {
    double *lmp_ptr;
    auto v = parser->split_line(parser->configuration["SampleFixes"]);
    int nfixes = v.size()/2;
    LAMMPSSimulator::run_script("PreSample");
    std::string cmd = "run "+parser->configuration["SampleSteps"];
    LAMMPSSimulator::run_commands(cmd);
    // need to have one for each worker....
    for(int j=0;j<nfixes;j++) {
      int f_s = std::stoi(v[2*j+1]);
      for(int i=0;i<f_s;i++) {
        lmp_ptr = (double *) lammps_extract_fix(lmp,&*v[2*j].begin(),0,1,i,0);
        results["f_"+v[2*j]+"_"+std::to_string(i)] = *lmp_ptr;
        lammps_free(lmp_ptr);
      }
    }
    LAMMPSSimulator::run_script("PostSample");
  };

  void reset() {
    LAMMPSSimulator::reset();
  };

  void sample(Holder params, double *dev) {
    reset();

    // change potential....
    LAMMPSSimulator::run_script("Input");
    // set up natoms etc
    LAMMPSSimulator::fill_lammps_vectors();

    LAMMPSSimulator::run_script("Input");

    LAMMPSSimulator::sample(params,dev);
  };
};


#endif
