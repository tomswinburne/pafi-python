#ifndef CGA_H
#define CGA_H

#include "GenericGatherer.hpp"

/*
  Gatherer collates simulation results and writes output files
  See comments below for details
*/
class CustomGatherer : public GenericGatherer {
public:
  CustomGatherer(Parser &p, int _nW, int di, int _rank) :
    GenericGatherer(p, _nW, di, _rank){
    field_width=15; // narrower print out (more fields)
    has_dV_fix = false;
    has_lambda_sweep = false;
  };

  void screen_output_header(bool end=true) {
    /*
      Header line for screen output
      GenericGatherer::screen_output_header only outputs for rank==0,
      but this function will output for any rank
      boolean argument : linebreak at end of output.
      If false, other arguments can be appended.

      Standard PAFI equivalent:
      GenericGatherer::screen_output_header(true);

      Here, we additionally print:
       - "Lambda" sweep parameter
       - average result of "SampleFixes", named dV

      We check that these fields exist, using find()
      PAFI will still run if they don't exist,
      but the .csv output may not be readable!
    */
    GenericGatherer::screen_output_header(false);
    if(rank>0) return; // only print for rank==0

    // check if Lambda is in sweep variables
    if(params.find("Lambda")!=params.end())
      has_lambda_sweep = true;

    // Check if "dV" in SampleFixes
    if(parser->configuration.find("SampleFixes")!=parser->configuration.end()) {
      auto v = parser->split_line(parser->configuration["SampleFixes"]);
      for(int nf=0;nf<v.size()/2;nf++) if(v[2*nf]=="dV") has_dV_fix=true;
    }

    if(has_lambda_sweep) std::cout<<std::setw(field_width)<<"Lambda";
    if(has_dV_fix) std::cout<<std::setw(field_width)<<"av(dV)";
    std::cout<<std::endl;
  };

  void screen_output_line(bool end=true) {
    /*
      Header line for screen output of a single sampling run
      GenericGatherer::screen_output_line only outputs for rank==0,
      but this function will output for any rank
      boolean argument : linebreak at end of output.
      If false, other arguments can be appended.

      Standard PAFI equivalent:
      GenericGatherer::screen_output_line(true);

      See above for description of this usage case

     */
    GenericGatherer::screen_output_line(false);
    if(rank>0) return; // only print for rank==0

    if(has_lambda_sweep) std::cout<<std::setw(field_width)<<params["Lambda"];

    if(has_dV_fix) std::cout<<std::setw(field_width)<<ens_results["f_dV"].first;

    std::cout<<std::endl;
  };
private:
  bool has_dV_fix, has_lambda_sweep;
};

#endif
