void MR3_get_board_version(int *version, int *freq)
{
  printf("** MR3_get_board_version is called **\n");
  *version=999;
  *freq=999;
  return;
}


void MR3_rscale_and_pack(int n, double (*x)[3], int *atype,
			 longlong (*rj)[2])
{
  printf("** MR3_rscale_and_pack is called **\n");
}


void MR3_rscale_and_pack_old(int n, double (*x)[3], int *atype,
			     longlong (*rj)[2])
{
  printf("** MR3_rscale_and_pack_old is called **\n");
}


void MR3SetTable_emu(char *filename, int tblno, int flag)
{
  printf("** MR3SetTable_emu is called **\n");
}


int MR3_get_total_num_chips(void)
{
  printf("** MR3_get_total_num_chips is called **\n");
  return 0;
}


void MR3print_chip_temperature(int flag)
{
  printf("** MR3print_chip_temperature is called **\n");
}


int MR3_get_retry_count(void)
{
  printf("** MR3_get_retry_count is called **\n");
  return 0;
}


void m3_free_unit(M3_UNIT *unit)
{
  printf("** m3_free_unit is called **\n");
  vg_exit(1);
}


M3_UNIT *m3_allocate_periodic_unit(const char *fname, int mode, double size, int dum)
{
  printf("** m3_allocate_periodic_unit is called **\n");
  vg_exit(1);
  return NULL;
}


void m3_set_positions(M3_UNIT *unit, double (*pos)[3], int n)
{
  printf("** m3_set_positions is called **\n");
  vg_exit(1);
}


void m3_set_charges(M3_UNIT *unit, const double *q, int n)
{
  printf("** m3_set_charges is called **\n");
  vg_exit(1);
}


void m3_set_cells(M3_UNIT *unit, M3_CELL *cell, int ncell)
{
  printf("** m3_set_cells is called **\n");
  vg_exit(1);
}


void m3_calculate_forces(M3_UNIT *unit, double (*xi)[3], int ni, double (*force)[3])
{
  printf("** m3_calculate_forces is called **\n");
  vg_exit(1);
}


void m3_calculate_potentials(M3_UNIT *unit, double (*xi)[3], int ni, double *eng)
{
  printf("** m3_calculate_potentials is called **\n");
  vg_exit(1);
}


void m3_set_types(M3_UNIT *unit, const int *atype, int n)
{
  printf("** m3_set_types is called **\n");
  vg_exit(1);
}


void m3_set_pipeline_types(M3_UNIT *unit, const int *atype, int n)
{
  printf("** m3_set_pipeline_types is called **\n");
  vg_exit(1);
}


void m3_set_charge_matrix(M3_UNIT *unit, const double *gscale, int nati, int natj)
{
  printf("** m3_set_charge_matrix is called **\n");
  vg_exit(1);
}


void m3_set_rscale_matrix(M3_UNIT *unit, const double *rscale, int nati, int natj)
{
  printf("** m3_set_rscale_matrix is called **\n");
  vg_exit(1);
}





