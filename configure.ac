AC_INIT([exaFMM], [1.0])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_FILES([Makefile examples/Makefile])
AC_CONFIG_HEADER(config.h)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([examples/fmm.cxx])
AM_INIT_AUTOMAKE([foreign subdir-objects -Wall -Werror])
AC_PROG_CC
AC_PROG_CXX
AC_PROG_F77
AX_CC_MAXOPT
AX_CXX_MAXOPT
AX_F77_MAXOPT

# OpenMP
AC_ARG_ENABLE(openmp, [AC_HELP_STRING([--enable-openmp],[use OpenMP])], enable_openmp=$e\
nableval, enable_openmp=no)
if test "$enable_openmp" = "yes"; then
   AC_DEFINE(HAVE_OPENMP,1,[Use OpenMP])
   AX_OPENMP([], [AC_MSG_ERROR([don't know how to enable OpenMP])])
fi
AM_CONDITIONAL(OPENMP, test "$enable_openmp" = "yes")

# Intel TBB
AC_ARG_WITH(tbb, [AC_HELP_STRING([--with-tbb],[use Intel TBB])], with_tbb=$withval, with_tbb=no)
if test "$with_tbb" = "yes"; then
        AC_DEFINE(WITH_TBB,1,[Use Intel TBB])
fi
AM_CONDITIONAL(WITH_TBB, test "$with_tbb" = "yes")

# Intel Cilk
AC_ARG_WITH(cilk, [AC_HELP_STRING([--with-cilk],[use Intel Cilk])], with_cilk=$withval, with_cilk=no)
if test "$with_cilk" = "yes"; then
        AC_DEFINE(WITH_CILK,1,[Use Intel Cilk])
fi
AM_CONDITIONAL(WITH_CILK, test "$with_cilk" = "yes")

# MassiveThreads
AC_ARG_WITH(mthread, [AC_HELP_STRING([--with-mthread],[use MassiveThreads])], with_mthread=$withval, with_mthread=no)
if test "$with_mthread" = "yes"; then
        AC_DEFINE(WITH_MTHREAD,1,[Use MassiveThreads])
fi
AM_CONDITIONAL(WITH_MTHREAD, test "$with_mthread" = "yes")

# QThreads
AC_ARG_WITH(qthread, [AC_HELP_STRING([--with-qthread],[use QThreads])], with_qthread=$withval, with_qthread=no)
if test "$with_qthread" = "yes"; then
        AC_DEFINE(WITH_QTHREAD,1,[Use QThreads])
fi
AM_CONDITIONAL(WITH_QTHREAD, test "$with_qthread" = "yes")

# Single precision
AC_ARG_ENABLE(single, [AC_HELP_STRING([--enable-single],[compile exaFMM in single precision])], ok=$enableval, ok=no)
AC_ARG_ENABLE(float, [AC_HELP_STRING([--enable-float],[synonym for --enable-single])], ok=$enableval)
if test "$ok" = "yes"; then
        AC_DEFINE(EXAFMM_SINGLE,1,[Define to compile in single precision.])
        PRECISION=s
else
        PRECISION=d
fi
AM_CONDITIONAL(EXAFMM_SINGLE, test "$ok" = "yes")
AC_SUBST(PRECISION)

# MPI
AC_ARG_ENABLE(mpi, [AC_HELP_STRING([--enable-mpi],[compile exaFMM MPI library])], enable_mpi=$enableval, enable_mpi=no)
if test "$enable_mpi" = "yes"; then
   ACX_MPI([],[AC_MSG_ERROR([could not find mpi library for --enable-mpi])])
   AC_CHECK_PROG(MPIRUN, mpirun, mpirun)
   AC_SUBST(MPIRUN)
   CXX=$MPICXX
fi
AM_CONDITIONAL(MPI, test "$enable_mpi" = "yes")

# SSE3
AC_ARG_ENABLE(sse3, [AC_HELP_STRING([--enable-sse3],[enable SSE/SSE3 optimizations])], have_sse3=$enableval, have_s\
se3=no)
if test "$have_sse3" = "yes"; then
        AC_DEFINE(HAVE_SSE3,1,[Define to enable SSE/SSE3 optimizations.])
        if test "$PRECISION" != "d" -a "$PRECISION" != "s"; then
                AC_MSG_ERROR([SSE3 requires single or double precision])
        fi
	if test "x$SSE3_CXXFLAGS" = x; then
   	   AX_CHECK_COMPILER_FLAGS(-$flag, [SSE3_CXXFLAGS="-msse3"],
   				   	   [AC_MSG_ERROR([Need a version of gcc with -msse3])])
        fi
fi
AM_CONDITIONAL(HAVE_SSE3, test "$have_sse3" = "yes")
AC_SUBST(SSE3_CXXFLAGS)

# AVX
AC_ARG_ENABLE(avx, [AC_HELP_STRING([--enable-avx],[enable AVX optimizations])], have_avx=$enableval, have_avx=no)
if test "$have_avx" = "yes"; then
   AC_DEFINE(HAVE_AVX,1,[Define to enable AVX optimizations.])
        if test "$PRECISION" != "d" -a "$PRECISION" != "s"; then
                AC_MSG_ERROR([AVX requires single or double precision])
        fi
	if test "x$AVX_CXXFLAGS" = x; then
   	  AX_CHECK_COMPILER_FLAGS(-mavx, [AVX_CXXFLAGS="-mavx"],
                		  	 [AC_MSG_ERROR([Need a version of gcc with -mavx])])
	fi
fi
AM_CONDITIONAL(HAVE_AVX, test "$have_avx" = "yes")
AC_SUBST(AVX_CXXFLAGS)

# Neon
AC_ARG_ENABLE(neon, [AC_HELP_STRING([--enable-neon],[enable ARM NEON optimizations])], have_neon=$enableval, have_n\
eon=no)
if test "$have_neon" = "yes"; then
        AC_DEFINE(HAVE_NEON,1,[Define to enable ARM NEON optimizations.])
        if test "$PRECISION" != "s"; then
                AC_MSG_ERROR([NEON requires single precision])
        fi
	if test "x$NEON_CXXFLAGS" = x; then
   	  AX_CHECK_COMPILER_FLAGS(-mfpu=neon, [NEON_CXXFLAGS="-mfpu=neon"],
            				      [AC_MSG_ERROR([Need a version of gcc with -mfpu=neon])])
        fi
fi
AM_CONDITIONAL(HAVE_NEON, test "$have_neon" = "yes")
AC_SUBST(NEON_CXXFLAGS)

# Disable SIMD
AC_ARG_ENABLE(simd, [AC_HELP_STRING([--disable-simd],[disable SIMD optimizations])], use_simd=$enableval, use_simd=yes)
if test "$use_simd" = "no"; then
        AC_DEFINE(USE_SIMD,0,[Disable SIMD optimizations.])
fi
AM_CONDITIONAL(USE_SIMD, test "$use_simd" = "yes")

# Weights for partitioning
AC_ARG_ENABLE(weight, [AC_HELP_STRING([--enable-weight],[enable weight for partitioning])], use_weight=$enableval, use_weight=no)
if test "$use_weight" = "yes"; then
        AC_DEFINE(USE_WEIGHT,1,[Enable weight for partitioning.])
fi
AM_CONDITIONAL(USE_WEIGHT, test "$use_weight" = "yes")

# Kahan summation
AC_ARG_ENABLE(kahan, [AC_HELP_STRING([--enable-kahan],[enable Kahan summation])], use_kahan=$enableval, use_kahan=no)
if test "$use_kahan" = "yes"; then
        AC_DEFINE(USE_KAHAN,1,[Enable Kahan summation.])
fi
AM_CONDITIONAL(USE_KAHAN, test "$use_kahan" = "yes")

# PAPI
AC_ARG_ENABLE(papi, [AC_HELP_STRING([--enable-papi],[enable PAPI performance counter])], use_papi=$enableval, use_papi=no)
if test "$use_papi" = "yes"; then
        AC_DEFINE(USE_PAPI,1,[Enable PAPI performance counter.])
fi
AM_CONDITIONAL(USE_PAPI, test "$use_papi" = "yes")

AC_OUTPUT
