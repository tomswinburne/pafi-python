#ifndef CSIM_H
#define CSIM_H
#include "LAMMPSSimulator.hpp"

/*
  Simulator performs all the .... simulation
  See comments below for details
*/
class CustomSimulator : public LAMMPSSimulator {

public:
  CustomSimulator(MPI_Comm &instance_comm, Parser &p, Holder &h, int t)
  : LAMMPSSimulator(instance_comm, p, h, t){
    /*
      anything extra to initialize ?
    */
  };

  void reset() {
    LAMMPSSimulator::reset();
  };

  void constrained_average() {
    /*
      The main LAMMPSSimulator::sample run fills
      std::map<std::string,double> results with:
      0: MaxDev : total |X_0-X_t|
      1: MaxJump : max per-atom entry of |X_0-X_t|
      2: MinEnergy : potential energy following in-plane minimization
      3: Valid : considered valid so far: MaxJump< MaxJump in config file
      4: aveF : time averaged PAFI free energy gradient
      6: avePsi : time averaged PAFI MEP,MFEP overlap
      7: dXTangent : projection of [X_t-X_0] on plane normal. should be 0 ~ 1e-15
      8: postT : equipartition temperature following sample run
      9: preT : equipartition temperature before sample run
      10: stdF : Time variance of PAFI free energy gradient
    */

    double *lmp_ptr;
    auto v = parser->split_line(parser->configuration["SampleFixes"]);

    int nfixes = v.size()/2;
    LAMMPSSimulator::run_script("PreSample");
    std::string cmd = "run "+parser->configuration["SampleSteps"];
    LAMMPSSimulator::run_commands(cmd);

    // need to have one for each worker....

    for(int j=0;j<nfixes;j++) {

      int f_s = std::stoi(v[2*j+1]); // number of fields from SampleFixes

      if(f_s==0) continue;

      if(f_s==1) { // returns a scalar
        lmp_ptr = (double *) lammps_extract_fix(lmp,&*v[2*j].begin(),0,0,0,0);
        results["f_"+v[2*j]+"_0"] = *lmp_ptr;
        lammps_free(lmp_ptr);
      } else { // returns a vector
        for(int i=0;i<f_s;i++) {
          lmp_ptr = (double *) lammps_extract_fix(lmp,&*v[2*j].begin(),0,1,i,0);
          results["f_"+v[2*j]+"_"+std::to_string(i)] = *lmp_ptr;
          lammps_free(lmp_ptr);
        }
      }
      
    }
    LAMMPSSimulator::run_script("PostSample");
  };
};


#endif
