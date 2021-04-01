#ifndef UTILS_H
#define UTILS_H


#include <math.h>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <array>

//#include <Eigen/Dense> // For Supercell and Hessian

/*
typedef std::map<std::string,double> SimParam;

class DataGatherer {
public:
  DataGatherer(Parser &p,std::vector<double> pathway_r) {
    params = &p;
    dsize=0;
    ens_data==NULL;

  }
  Parser *params;
  int dsize=0;
  double *ens_data;
  int dsize,nWorkers;

  std::map<std::string,std::pair<std::vector<double>> all_results;

};


*/


class DataGatherer {
public:
  DataGatherer() {
    dsize=0;
    ens_data==NULL;
  }

  virtual int initialize(Parser &p, std::string dump_file,int _nWorkers) {
    if(raw.is_open()) raw.close();
    nWorkers = _nWorkers;
    params = &p;
    raw.open(dump_file.c_str(),std::ofstream::out);
    if(raw.is_open()) return 1;
    std::cout<<"Could not open "<<dump_file<<"! EXIT"<<std::endl;
    return 0;
  };

  virtual void prepare(std::map<std::string,double> &results) {

    if(!raw.is_open()) return;
    int i=1,j=0;

    raw<<"# 0: r ";
    for(auto res: results) {
      raw<<i++<<": "<<res.first<<"  ";
      all_results[res.first] = *(new std::vector<double>);
    }
    raw<<std::endl;

    dsize = results.size();
    if(ens_data==NULL) ens_data = new double[dsize*2+1];


    // ensemble average
    for(j=0;j<2*dsize+1;j++) ens_data[j] = 0.0;
  };

  int ensemble_average(double r, int *valid, double *all_data) {
    if(!raw.is_open()) return 0;

    int total_valid=0,i,j;

    // raw output
    raw<<r<<" ";
    for(i=0;i<nWorkers*dsize;i++) raw<<all_data[i]<<" ";
    raw<<std::endl;


    // all results

    for(i=0;i<nWorkers;i++) {
      sample_r.push_back(r);
      j=0; for(auto &res: all_results)
        res.second.push_back(all_data[i*dsize+j++]);
    }

    // unnormalize for sum
    for(j=0;j<dsize;j++) {
      ens_data[j+dsize] *= ens_data[2*dsize];
      ens_data[j+dsize] += ens_data[j] * ens_data[j];
    }
    for(j=0;j<2*dsize;j++) ens_data[j] *= ens_data[2*dsize];

    // add more data
    for(i=0;i<nWorkers;i++) if(valid[i]==1) {
      total_valid++;
      ens_data[2*dsize] += 1.0;
      for(j=0;j<dsize;j++) {
        ens_data[j]+=all_data[i*dsize+j];
        ens_data[j+dsize] += all_data[i*dsize+j] * all_data[i*dsize+j];
      }
    }

    // renormalize for access
    if(ens_data[2*dsize]>0.5) {
      for(j=0;j<2*dsize;j++) ens_data[j] /= ens_data[2*dsize];
      for(j=0;j<dsize;j++) {
        ens_data[j+dsize] -= ens_data[j] * ens_data[j];
        ens_data[j+dsize] /= ens_data[2*dsize]; // N^2 for aves-of-aves
      }
    }
    #ifdef VERBOSE
    std::cout<<"END DataGatherer:: ensemble()"<<std::endl;
    #endif
    return total_valid;

  };

  void next() {
    // wipe ens_data
    for(int j=0;j<2*dsize+1;j++) ens_data[j] = 0.0;
  };

  void close() {
    if(raw.is_open()) raw.close();
  };


Parser *params;
std::ofstream raw;
double *ens_data;
int dsize,nWorkers;
std::map<std::string,std::vector<double>> all_results;
std::vector<double> sample_r;

};
/*

class SampleParameters {
public:
  SampleParameters(Parser &p,std::vector<double> pathway_r) {




  sample_r.clear();
  double dr = 0.1;
  if (params->nPlanes>1)
    dr = (params->stopr-params->startr)/(double)(params->nPlanes-1);

  if(params->spline_path and not params->match_planes) {
    for (double r = params->startr; r <= params->stopr+0.5*dr; r += dr )
      sample_r.push_back(r);
  } else {
    for(auto r: pathway_r) if(r>=0.0 && r<=1.0) sample_r.push_back(r);
  }


  Parser *params;
  std::string dump_suffix, dev_file, dump_file, raw_file;
  std::map<std::string,double> sim_param; // simulation parameters

  std::vector<double> sample_r; // obligatory
  std::vector<double> temperature_r; // obligatory
  void determine_r(std::vector )

};
*/

#endif
