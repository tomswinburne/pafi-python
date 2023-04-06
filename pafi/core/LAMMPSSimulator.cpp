#include "LAMMPSSimulator.hpp"

LAMMPSSimulator::LAMMPSSimulator (MPI_Comm &instance_comm, Parser &p, Holder &h, int t)
  : GenericSimulator::GenericSimulator(instance_comm, p, h, t) {
  // set up LAMMPS
  char str1[32];
  char **lmparg = new char*[5];
  lmparg[0] = NULL; // required placeholder for program name
  lmparg[1] = (char *) "-screen";
  lmparg[2] = (char *) "none";
  lmparg[3] = (char *) "-log";
  if(parser->loglammps) {
    sprintf(str1,"log.lammps.%d",tag);
    lmparg[4] = str1;
  } else lmparg[4] = (char *) "none";

  lmp = new LAMMPS(5,lmparg,*comm);

  // check version number and existence of fix-pafi and quit if otherwise
  if(!check_lammps_compatibility()) return;

  // reset the system
  reset();

  // run input script
  run_script("Input");

  natoms=0;
  species=NULL;
  q=NULL;
  image=NULL;
  id=NULL;
  lt=NULL; // for scattering

  // these won't change size even if we rereun "input"
  fill_lammps_vectors();

  simulator_name = "LAMMPSSimulator";

};

bool LAMMPSSimulator::check_lammps_compatibility() {
  // check version number and existence of fix-pafi
  int lammps_release_int = lammps_version(lmp); // YYYYMMDD

  pafi_package = "USER-MISC";

  if(lammps_release_int<20201101) {
    if(local_rank==0)
      std::cout<<"LAMMPSSimulator(): Require version > 28July2020!"<<std::endl;
    return false;
  }
  if(lammps_release_int>=20210728) pafi_package = "EXTRA-FIX";

  #ifdef VERBOSE
  if(local_rank==0)
    std::cout<<"LAMMPSSimulator(): Searching for "<<pafi_package<<std::endl;
  #endif

  has_pafi = (bool)lammps_config_has_package((char *) pafi_package.c_str());
  #ifdef VERBOSE
  if(local_rank==0)
    std::cout<<"LAMMPSSimulator(): has_pafi: "<<has_pafi<<std::endl;
  #endif

  return has_pafi;

};

void LAMMPSSimulator::reset() {
  run_commands("clear");
  made_fix=false;
  made_compute=false;
};

void LAMMPSSimulator::fill_lammps_vectors() {

  #ifdef VERBOSE
  if(local_rank==0) std::cout<<"LAMMPSSimulator(): Ran input script"<<std::endl;
  #endif
  int new_natoms = *((int *) lammps_extract_global(lmp,(char *) "natoms"));
  if(natoms>0 and new_natoms!=natoms and local_rank==0) {
    std::cout<<"LAMMPSSimulator(): Atom count changed on reload!!"<<std::endl;
  }
  natoms = new_natoms;
  #ifdef VERBOSE
  if(local_rank==0) std::cout<<"LAMMPSSimulator(): natoms: "<<natoms<<std::endl;
  #endif

  /* possible starting point for parallel IO- not really required here.
  nlocal = (3*natoms) / local_size;
  if(local_rank==local_size-1) nlocal += (3*natoms)%local_size;
  offset =  nlocal*local_rank;
  */
  nlocal = 3*natoms;
  offset = 0;

  if(id==NULL) id = new int[natoms];
  gather("id",1,id);
  #ifdef VERBOSE
  if(local_rank==0) std::cout<<"LAMMPSSimulator(): gathered id"<<std::endl;
  #endif

  // get cell info TODO - allow to vary? can always redo....
  pbc.load(getCellData());

  // Get type / image info
  if(species==NULL) species = new int[natoms];
  gather("type",1,species);
  #ifdef VERBOSE
  if(local_rank==0) std::cout<<"LAMMPSSimulator(): gathered type"<<std::endl;
  #endif

  s_flag=true;

  if(image==NULL) image = new int[natoms];
  gather("image",1,image);
  #ifdef VERBOSE
  if(local_rank==0) std::cout<<"LAMMPSSimulator(): gathered image"<<std::endl;
  #endif

  if(x==NULL) x = new double[3*natoms]; // mainly for declaration
  gather("x",3,x);
  #ifdef VERBOSE
  if(local_rank==0) std::cout<<"LAMMPSSimulator(): gathered x"<<std::endl;
  #endif

  if(lt==NULL) lt = new double[natoms];

  if(!has_pafi and local_rank==0) {
    std::cout<<"PAFI Error: missing "<<pafi_package<<" package in LAMMPS"<<std::endl;
    if(error_count>0) std::cout<<last_error()<<std::endl;
  }
};

/*
  Load xyz configuration from data file and put in vector
*/
void LAMMPSSimulator::load_config(std::string file_string, double *x) {
  made_fix=false;
  made_compute=false;
  std::string cmd = "delete_atoms group all\nread_data ";
  cmd += file_string;
  cmd += " add merge";
  run_commands(cmd);
  gather("x",3,x);
  if(error_count>0 && local_rank==0) std::cout<<last_error()<<std::endl;
};

/*
  Parse and run script from configuration file
*/
void LAMMPSSimulator::run_script(std::string sn){
  std::vector<std::string> strv = parser->Script(sn);
  if(strv.size()>0) {
    strv.push_back("run 0");
    run_commands(strv);
  }
};

/*
  Parse and run script from string with linebreaks
*/
void LAMMPSSimulator::run_commands(std::vector<std::string> strv) {
  #ifdef VERBOSE
  if(local_rank==0) std::cout<<"LAMMPSSimulator.run_commands(): "<<std::endl;
  #endif
  std::vector<std::string> strv_parsed;
  for(auto s:strv) {
    parser->insert_params(s,*params);
    strv_parsed.push_back(s);
    #ifdef VERBOSE
    if(local_rank==0) std::cout<<s<<std::endl;
    #endif
    lammps_command((void *)lmp,(char *)s.c_str());
    log_error(strv_parsed);
  }
  //MPI_Barrier(*comm);

};

void LAMMPSSimulator::run_commands(std::string strv) {
  #ifdef VERBOSE
  if(local_rank==0) std::cout<<"LAMMPSSimulator.run_commands():\n"<<strv<<std::endl;
  #endif
  run_commands(parser->split_lines(strv));
};


void LAMMPSSimulator::log_error(std::string lc) {
  if(bool(lammps_has_error(lmp))) {
    error_count++;
    char error_message[2048];
    int error_type = lammps_get_last_error_message(lmp,error_message,2048);
    if ((error_count==1) or (error_message!=last_error_message) )
      if(local_rank==0) {
        std::cout<<"LAMMPSSimulator.log_error():\n\t"<<error_message<<std::endl;
        std::cout<<"\tFrom command: "<<lc<<std::endl;
      }
    last_error_message = error_message;
    last_command = lc+"\n";
  }
};

void LAMMPSSimulator::log_error(std::vector<std::string> lc) {
  if(bool(lammps_has_error(lmp))) {
    error_count++;
    char error_message[2048];
    int error_type = lammps_get_last_error_message(lmp,error_message,2048);
    last_command = "";
    for(auto s:lc) last_command += s+"\n";
    if ((error_count==1) or (error_message!=last_error_message) )
      if(local_rank==0) {
        std::cout<<"LAMMPSSimulator.log_error():\n\t"<<error_message<<std::endl;
        std::cout<<"\tFrom commands: \n";
        for(auto s:lc) std::cout<<"\t"<<s<<std::endl;
      }
    last_error_message = error_message;
  }
};

// over load for type
void LAMMPSSimulator::gather(std::string name, int c, int *v){
  lammps_gather(lmp,(char *)name.c_str(),0,c,v);
  std::string lc = "lammps_gather("+name+",int,"+std::to_string(c)+")";
  log_error(lc);
};

// over load for type
void LAMMPSSimulator::gather(std::string name, int c, double *v){
  lammps_gather(lmp,(char *)name.c_str(),1,c,v);
  std::string lc = "lammps_gather("+name+",double,"+std::to_string(c)+")";
  log_error(lc);
};

// over load for type
void LAMMPSSimulator::scatter(std::string name, int c, int *v){
  lammps_scatter(lmp,(char *)name.c_str(),0,c,v);
  std::string lc="lammps_scatter("+name+",int,"+std::to_string(c)+")";
  log_error(lc);
};

// over load for type
void LAMMPSSimulator::scatter(std::string name, int c, double * v){
  lammps_scatter(lmp,(char *)name.c_str(),1,c,v);
  std::string lc="lammps_scatter("+name+",double,"+std::to_string(c)+")";
  log_error(lc);
};

std::string LAMMPSSimulator::last_error(){
  std::string res = "\nworker "+std::to_string(tag)+" had ";
  res += std::to_string(error_count)+" errors\n\tlast message:\n";
  res += last_error_message+"\n\tfrom commands:\n"+last_command+"\n";
  return res;
};

/*
  Fill configuration, path, tangent and tangent gradient. Return tangent norm
*/
void LAMMPSSimulator::populate(double r, double &norm_mag, double T) {

  rescale_cell(T); // updates scale vector


  std::string cmd;
  double ncom[]={0.,0.,0.};
  norm_mag=0.;
  // check for __pafipath fix and compute
  if (!made_fix) {
    #ifdef VERBOSE
    if(local_rank==0)
      std::cout<<"LAMMPSSimulator.populate(): making __pafipath fix"<<std::endl;
    #endif
    cmd = "fix __pafipath all property/atom ";
    cmd += "d_ux d_uy d_uz d_nx d_ny d_nz d_dnx d_dny d_dnz\nrun 0";
    run_commands(cmd);
    made_fix=true;
  }


  std::string xyz[3]; xyz[0]="x"; xyz[1]="y"; xyz[2]="z";
  for(int i=0;i<natoms;i++) for(int j=0;j<3;j++) x[3*i+j] = pathway[3*i+j].deriv(0,r) * scale[j];
  scatter("x",3,x);

  for(int j=0;j<3;j++) {
    #ifdef VERBOSE
    if(local_rank==0)
      std::cout<<"LAMMPSSimulator.populate(): Scattering d_u"+xyz[j]<<std::endl;
    #endif
    for(int i=0;i<natoms;i++)  lt[i] = pathway[3*i+j].deriv(0,r) * scale[j];
    scatter("d_u"+xyz[j],1,lt);
  }

  for(int i=0;i<natoms;i++) for(int j=0;j<3;j++) x[3*i+j] = pathway[3*i+j].deriv(1,r) * scale[j];

  // Center of mass projection and tangent length
  for(int i=0;i<3*natoms;i++) ncom[i%3] += x[i]/(double)natoms;
  for(int i=0;i<3*natoms;i++) x[i] -= ncom[i%3];
  for(int i=0;i<3*natoms;i++) norm_mag += x[i]*x[i];
  norm_mag = sqrt(norm_mag);
  for(int i=0;i<3*natoms;i++) x[i] /= norm_mag;
  for(int j=0;j<3;j++) {
    for(int i=0;i<natoms;i++)  lt[i] = x[3*i+j];
    scatter("d_n"+xyz[j],1,lt);
    for(int i=0;i<natoms;i++)  lt[i] = pathway[3*i+j].deriv(2,r) * scale[j] / norm_mag / norm_mag;
    scatter("d_dn"+xyz[j],1,lt);
  }

  run_commands("run 0");

  if(!made_compute) {
    #ifdef VERBOSE
    if(local_rank==0)
      std::cout<<
      "LAMMPSSimulator.populate(): making __pafipath compute"<<std::endl;
    #endif
    cmd = "compute __pafipath all property/atom ";
    cmd += "d_ux d_uy d_uz d_nx d_ny d_nz d_dnx d_dny d_dnz\nrun 0";
    run_commands(cmd);
    made_compute=true;
  }
};

/* Rescale simulation cell */
void LAMMPSSimulator::rescale_cell(double T) {
  #ifdef VERBOSE
  if(local_rank==0)
    std::cout<<
    "LAMMPSSimulator.rescale_cell(): T = "<<T<<std::endl;
  #endif
  double newscale[3];
  std::string cmd, ssx,ssy,ssz;
  expansion(T,newscale);
  ssx = std::to_string(newscale[0]/scale[0]);
  ssy = std::to_string(newscale[1]/scale[1]);
  ssz = std::to_string(newscale[2]/scale[2]);
  cmd ="change_box all x scale "+ssx+" y scale "+ssy+" z scale "+ssz+"\n";
  cmd += "run 0";
  run_commands(cmd);
  scale[0]=newscale[0]; scale[1]=newscale[1]; scale[2]=newscale[2];
  #ifdef VERBOSE
  if(local_rank==0)
    std::cout<<
    "END LAMMPSSimulator.rescale_cell(): T = "<<T<<std::endl;
  #endif

};

/*
  Main sample run. Results vector should have thermalization temperature,
  sample temperature <f>, <f^2>, <psi> and <x-u>.n, and max_jump
*/


void LAMMPSSimulator::sample(Holder params, double *dev) {

  error_count = 0;
  last_error_message="";
  results.clear();
  double a_disp=0.0,max_disp = 0.0;
  std::string cmd;
  double norm_mag, sampleT, dm;
  double *lmp_ptr;
  if(params.find("ReactionCoordinate")==params.end()) {
    if(local_rank==0)
      std::cout<<"LAMMPSSimulator: No ReactionCoordinate! Exiting!"<<std::endl;
    return;
  }
  double r = params["ReactionCoordinate"];
  if(params.find("Temperature")==params.end()) {
    if(local_rank==0)
      std::cout<<"LAMMPSSimulator: No Temperature! Exiting!"<<std::endl;
    return;
  }
  double T = params["Temperature"];

  int overdamped = 0;
  if(parser->configuration.find("OverDamped")==parser->configuration.end()) {
    if(local_rank==0)
      std::cout<<"LAMMPSSimulator: No OverDamped! Defaulting to 0"<<std::endl;
    parser->configuration["OverDamped"] = std::to_string(0);
  }
  overdamped = std::stoi(parser->configuration["OverDamped"]);


  populate(r,norm_mag,0.0);
  run_script("PreRun");
  populate(r,norm_mag,T);

  // pafi fix
  parser->configuration["Temperature"] = std::to_string(T);
  cmd += "fix hp all pafi __pafipath "+parser->configuration["Temperature"]+" ";
  cmd += parser->configuration["Friction"]+" ";
  cmd += parser->seed_str()+" overdamped ";
  cmd += parser->configuration["OverDamped"]+" com 1\nrun 0";
  run_commands(cmd);

  if(parser->preMin) {
    #ifdef VERBOSE
    if(local_rank==0)
      std::cout<<"LAMMPSSimulator.sample(): minimizing"<<std::endl;
    #endif
    cmd = "min_style fire\n minimize 0 0.001 ";
    cmd += parser->configuration["MinSteps"]+" "+parser->configuration["MinSteps"];
    run_commands(cmd);
  }

  run_script("PreTherm");

  // time average
  MinEnergy = getEnergy();
  results["MinEnergy"] = MinEnergy;

  lmp_ptr = (double *) lammps_extract_fix(lmp,(char *)"hp",0,1,4,0);
  refP = *lmp_ptr;
  lammps_free(lmp_ptr);

  cmd = "reset_timestep 0\n";
  cmd += "fix ae all ave/time 1 "+parser->configuration["ThermWindow"]+" ";
  cmd += parser->configuration["ThermSteps"]+" ";
  if(overdamped==1) cmd += "c_pe\n";
  else cmd += "c_thermo_temp\n";
  cmd += "run "+parser->configuration["ThermSteps"];
  run_commands(cmd);

  // pre temperature
  lmp_ptr = (double *) lammps_extract_fix(lmp,(char *)"ae",0,0,0,0);
  if(overdamped==1) sampleT = (*lmp_ptr-MinEnergy)/natoms/1.5/BOLTZ;
  else sampleT = *lmp_ptr;
  lammps_free(lmp_ptr);

  results["preT"] = sampleT;
  run_commands("unfix ae\nrun 0");


  // time averages for sampling TODO: groupname for ave/atom
  std::string SampleSteps = parser->configuration["SampleSteps"];
  cmd = "reset_timestep 0\n";
  cmd += "fix ae all ave/time 1 "+SampleSteps+" "+SampleSteps+" ";
  if(overdamped==1) cmd += "c_pe\n";
  else cmd += "c_thermo_temp\n";
  if(parser->postDump)
    cmd += "fix ap all ave/atom 1 "+SampleSteps+" "+SampleSteps+" x y z\n";

  cmd += "fix af all ave/time 1 "+SampleSteps+" "+SampleSteps;
  cmd += " f_hp[1] f_hp[2] f_hp[3] f_hp[4]\n";
  run_commands(cmd);

  // run SampleSteps and add to results
  constrained_average();

  lmp_ptr = (double *) lammps_extract_fix(lmp,(char *)"ae",0,0,0,0);
  if(overdamped==1) sampleT = (*lmp_ptr-MinEnergy)/natoms/1.5/BOLTZ;
  else sampleT = *lmp_ptr;
  lammps_free(lmp_ptr);
  results["postT"] = sampleT;

  lmp_ptr = (double *) lammps_extract_fix(lmp,(char *)"af",0,1,0,0);
  results["aveF"] = *lmp_ptr * norm_mag * -1.0;
  lammps_free(lmp_ptr);

  lmp_ptr = (double *) lammps_extract_fix(lmp,(char *)"af",0,1,1,0);
  results["stdF"] = *lmp_ptr * norm_mag * norm_mag;
  results["stdF"] -= results["aveF"] * results["aveF"];
  lammps_free(lmp_ptr);

  lmp_ptr = (double *) lammps_extract_fix(lmp,(char *)"af",0,1,2,0);
  results["avePsi"] = *lmp_ptr;
  lammps_free(lmp_ptr);

  lmp_ptr = (double *) lammps_extract_fix(lmp,(char *)"af",0,1,3,0);
  results["dXTangent"] = *lmp_ptr;
  lammps_free(lmp_ptr);

  // post minmization - max jump
  cmd = "min_style fire\n minimize 0 0.001 ";
  cmd += parser->configuration["MinSteps"]+" "+parser->configuration["MinSteps"];
  run_commands(cmd);
  gather("x",3,dev);
  for(int i=0;i<3*natoms;i++) dev[i] -= path(i,r,0,scale[i%3]);
  pbc.wrap(dev,3*natoms);
  for(int i=0;i<natoms;i++) {
    a_disp = 0.0;
    for(int j=0;j<3;j++) a_disp += dev[3*i+j]*dev[3*i+j];
    max_disp = std::max(a_disp,max_disp);
    if(!parser->postDump) for(int j=0;j<3;j++) dev[3*i+j] = 0.0;
  }
  results["MaxJump"] = sqrt(max_disp);

  results["Valid"] = double(bool(results["MaxJump"]<parser->maxjump_thresh));

  // deviation average
  //run_commands("reset_timestep "+SampleSteps); // for fix calculation
  if(parser->postDump) {
    gather("f_ap",3,dev);
    for(int i=0;i<3*natoms;i++) dev[i] = dev[i]/scale[i%3]-path(i,r,0,1.0);
    pbc.wrap(dev,3*natoms);
    dm = 0.0;
    for(int i=0;i<3*natoms;i++) {
      dev[i] *= scale[i%3];
      dm = std::max(dm,sqrt(dev[i] * dev[i]));
    }
    results["MaxDev"] = dm;
    run_commands("unfix ap");
  } else results["MaxDev"] = sqrt(max_disp);

  // reset
  run_commands("unfix ae\nunfix af\nunfix hp");

  // Stress Fixes
  run_script("PostRun");

  // rescale back
  populate(r,norm_mag,0.0);

};

double LAMMPSSimulator::getEnergy() {
  run_commands("run 0");
  double * lmpE = (double *) lammps_extract_compute(lmp,(char *) "pe",0,0);
  if (lmpE == NULL) {
    std::string cmd = "compute pe all pe\nvariable pe equal pe\nrun 0";
    run_commands(cmd);
    lmpE = (double *) lammps_extract_compute(lmp,(char *) "pe",0,0);
  }
  double baseE = *lmpE;
  return baseE;
};


double LAMMPSSimulator::getForceEnergy(double *f) {
  run_commands("run 0");
  double * lmpE = (double *) lammps_extract_compute(lmp,(char *) "pe",0,0);
  if (lmpE == NULL) {
    std::string cmd = "compute pe all pe\nvariable pe equal pe\nrun 0";
    run_commands(cmd);
    lmpE = (double *) lammps_extract_compute(lmp,(char *) "pe",0,0);
  }
  gather("f",3,f);
  double baseE = *lmpE;
  return baseE;
};

double LAMMPSSimulator::get_fix(std::string fixid,int type, int index) {
  double * lmp_val = (double *) lammps_extract_fix(lmp,(char *)fixid.c_str(),0,type,index,0);
  double val = *lmp_val;
  return val;
}

// Fill 9D array with Lx, Ly, Lz, xy, xz, yz, then periodicity in x, y, z
std::array<double,9> LAMMPSSimulator::getCellData() {
  std::array<double,9> cell;
  for(int i=0; i<9; i++) cell[i] = 0.;
  double *lmp_ptr;
  double *boxlo = (double *) lammps_extract_global(lmp,(char *) "boxlo");
  double *boxhi = (double *) lammps_extract_global(lmp,(char *) "boxhi");
  int *pv = (int *) lammps_extract_global(lmp,(char *) "periodicity");

  std::array<std::string,3> od;
  od[0]="xy"; od[1]="xz"; od[2]="yz";
  for(int i=0; i<3; i++) {
    cell[i] = boxhi[i]-boxlo[i];
    lmp_ptr = (double *) lammps_extract_global(lmp,(char *)od[i].c_str());
    cell[3+i] = *lmp_ptr;
    cell[6+i] = pv[i];
  }
  return cell;
};

void LAMMPSSimulator::constrained_average() {
  std::string cmd = "run "+parser->configuration["SampleSteps"];
  run_commands(cmd);
};



void LAMMPSSimulator::close() {
  GenericSimulator::close();
  delete [] id;
  delete [] species;
  delete [] image;
  delete [] lt;
  lammps_close(lmp);
};

std::string LAMMPSSimulator::header(double mass=55.85) {
  std::string res="LAMMPS dump file\n\n";
  res += std::to_string(natoms)+" atoms\n";
  res += "1 atom types\n\n"; // TODO species
  res += "0. "+std::to_string(pbc.cell[0][0]*scale[0])+" xlo xhi\n";
  res += "0. "+std::to_string(pbc.cell[1][1]*scale[1])+" ylo yhi\n";
  res += "0. "+std::to_string(pbc.cell[2][2]*scale[2])+" zlo zhi\n";
  res += std::to_string(pbc.cell[0][1]*scale[1])+" ";
  res += std::to_string(pbc.cell[0][2]*scale[2])+" ";
  res += std::to_string(pbc.cell[1][2]*scale[2])+" xy xz yz\n";
  //res += std::to_string(pbc.cell[2][2])
  res += "\nMasses\n\n 1 "+std::to_string(mass)+"\n Atoms\n\n";
  return res;
};

void LAMMPSSimulator::lammps_dump_path(std::string fn, double r) {

  std::ofstream out;
  if(local_rank==0){
    out.open(fn.c_str(),std::ofstream::out);
    out<<header();

    double ncom[]={0.,0.,0.};
    double c,nm=0.;

    for(int i=0;i<natoms;i++) for(int j=0;j<3;j++) {
      ncom[j] += pathway[3*i+j].deriv(1,r)*scale[j] / natoms;
    }

    for(int i=0;i<natoms;i++) for(int j=0;j<3;j++) {
      c = pathway[3*i+j].deriv(1,r)*scale[j]-ncom[j];
      nm += c * c;
    }
    nm = sqrt(nm);

    for(int i=0;i<natoms;i++) {
      out<<i+1<<" 1 "; // TODO multispecies
      for(int j=0;j<3;j++) out<<pathway[3*i+j](r)*scale[j]<<" "; // x y z
      out<<std::endl;
    }
    out<<"\nPafiPath\n"<<std::endl;
    for(int i=0;i<natoms;i++) {
      out<<i+1<<" "; // TODO multispecies
      for(int j=0;j<3;j++) out<<pathway[3*i+j](r)*scale[j]<<" "; // path
      for(int j=0;j<3;j++) out<<(pathway[3*i+j].deriv(1,r)*scale[j]-ncom[j])/nm<<" ";
      for(int j=0;j<3;j++) out<<pathway[3*i+j].deriv(2,r)*scale[j]/nm/nm<<" ";
      out<<std::endl;
    }
    out.close();
  }
};

void LAMMPSSimulator::lammps_write_dev(std::string fn, double r, double *dev, double *dev_sq) {
  std::ofstream out;
  out.open(fn.c_str(),std::ofstream::out);
  out<<header();
  for(int i=0;i<natoms;i++){
    out<<i+1<<" 1 "; // TODO multispecies
    for(int j=0;j<3;j++) out<<pathway[3*i+j](r)*scale[j]<<" "; // x y z
    for(int j=0;j<3;j++) out<<dev[3*i+j]<<" "; // <dev>
    for(int j=0;j<3;j++) out<<sqrt(dev_sq[3*i+j]-dev[3*i+j]*dev[3*i+j])<<" ";
    out<<std::endl;
  }
  out.close();
};
