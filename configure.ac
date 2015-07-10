AC_INIT([ExaFMM], [1.0])
AC_CONFIG_SRCDIR([examples/serial.cxx])
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE([foreign subdir-objects -Wall -Werror])
m4_include([m4/ax_prog_cxx_mpi.m4])
AX_PROG_CXX_MPI
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile examples/Makefile])
AC_OUTPUT
