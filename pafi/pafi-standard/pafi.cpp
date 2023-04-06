#include "pafi.hpp"

int main(int narg, char **arg) {
MPI_Init(&narg,&arg);
  MPI_Comm world=MPI_COMM_WORLD;
  run<Simulator,Gatherer>(world,"./config.xml");
  MPI_Finalize();
};
