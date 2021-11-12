#include "bem3_aw_b1.h"

int main(int argc,char **argv)
{
  DMDA md;

  read_dmda(argc,argv,&md);
  print_dmda(&md);
  initialize_dmda(&md);
  solve_bieq_dmda(&md);
  dat_write_dmda(argv[3],&md);
  finalize_dmda(&md);
  return 0;
}
