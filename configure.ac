AC_INIT([exaFMM], [1.0])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_FILES([Makefile examples/Makefile])
AC_CONFIG_HEADER(config.h)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([examples/serial.cxx])
AM_INIT_AUTOMAKE([foreign subdir-objects -Wall -Werror])
AC_PROG_CC
AC_PROG_CXX

# Single precision
AC_ARG_ENABLE(single, [AC_HELP_STRING([--enable-single],[compile exaFMM in single precision])], ok=$enableval, ok=no)
AC_ARG_ENABLE(float, [AC_HELP_STRING([--enable-float],[synonym for --enable-single])], ok=$enableval)
if test "$ok" = "yes"; then
        AC_DEFINE(EXAFMM_SINGLE,1,[Define to compile in single precision.])
        PRECISION=s
else
        PRECISION=d
fi
AM_CONDITIONAL(SINGLE, test "$ok" = "yes")
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

# SSE
AC_ARG_ENABLE(sse, [AC_HELP_STRING([--enable-sse],[enable SSE optimizations])], have_sse=$enableval, have_sse=no)
if test "$have_sse" = "yes"; then
        if test "$PRECISION" != "s"; then
                AC_MSG_ERROR([SSE requires single precision])
        fi
fi

# SSE3
AC_ARG_ENABLE(sse3, [AC_HELP_STRING([--enable-sse3],[enable SSE/SSE3 optimizations])], have_sse3=$enableval, have_s\
se2=no)
if test "$have_sse" = "yes"; then have_sse3=yes; fi
if test "$have_sse3" = "yes"; then
        AC_DEFINE(HAVE_SSE3,1,[Define to enable SSE/SSE3 optimizations.])
        if test "$PRECISION" != "d" -a "$PRECISION" != "s"; then
                AC_MSG_ERROR([SSE3 requires single or double precision])
        fi
fi
AM_CONDITIONAL(HAVE_SSE3, test "$have_sse3" = "yes")

# AVX
AC_ARG_ENABLE(avx, [AC_HELP_STRING([--enable-avx],[enable AVX optimizations])], have_avx=$enableval, have_avx=no)
if test "$have_avx" = "yes"; then
   AC_DEFINE(HAVE_AVX,1,[Define to enable AVX optimizations.])
        if test "$PRECISION" != "d" -a "$PRECISION" != "s"; then
                AC_MSG_ERROR([AVX requires single or double precision])
        fi
fi
AM_CONDITIONAL(HAVE_AVX, test "$have_avx" = "yes")

# Altivec
AC_ARG_ENABLE(altivec, [AC_HELP_STRING([--enable-altivec],[enable Altivec optimizations])], have_altivec=$enableval\
, have_altivec=no)
if test "$have_altivec" = "yes"; then
        AC_DEFINE(HAVE_ALTIVEC,1,[Define to enable Altivec optimizations.])
        if test "$PRECISION" != "s"; then
	   AC_MSG_ERROR([Altivec requires single precision])
        fi
fi
AM_CONDITIONAL(HAVE_ALTIVEC, test "$have_altivec" = "yes")

# Neon
AC_ARG_ENABLE(neon, [AC_HELP_STRING([--enable-neon],[enable ARM NEON optimizations])], have_neon=$enableval, have_n\
eon=no)
if test "$have_neon" = "yes"; then
        AC_DEFINE(HAVE_NEON,1,[Define to enable ARM NEON optimizations.])
        if test "$PRECISION" != "s"; then
                AC_MSG_ERROR([NEON requires single precision])
        fi
fi
AM_CONDITIONAL(HAVE_NEON, test "$have_neon" = "yes")

# Disable SIMD
AC_ARG_ENABLE(simd, [AC_HELP_STRING([--disable-simd],[disable SIMD optimizations])], use_simd=$enableval, use_simd=yes)
if test "$use_simd" = "no"; then
        AC_DEFINE(USE_SIMD,0,[Disable SIMD optimizations.])
fi
AM_CONDITIONAL(USE_SIMD, test "$use_simd" = "yes")

AX_CXX_MAXOPT

case "${ax_cv_c_compiler_vendor}" in
   intel) # Stop icc from defining __GNUC__, except on MacOS where this fails                                       
        case "${host_os}" in
            *darwin*) ;; # icc -no-gcc fails to compile some system headers                                         
            *)
               AX_CHECK_COMPILER_FLAGS([-no-gcc], [CXX="$CXX -no-gcc"])
               ;;
        esac
        ;;

   hp) # must (sometimes) manually increase cpp limits to handle fftw3.h                                            
        AX_CHECK_COMPILER_FLAGS([-Wp,-H128000],
                                [CXX="$CXX -Wp,-H128000"])
        ;;

   portland) # -Masmkeyword required for asm("") cycle counters                                                     
        AX_CHECK_COMPILER_FLAGS([-Masmkeyword],
                                [CXX="$CXX -Masmkeyword"])
        ;;
esac

# SIMD CXXFLAGS
if test "$have_sse3" = "yes" -a "x$SSE3_CXXFLAGS" = x; then
   AX_CHECK_COMPILER_FLAGS(-$flag, [SSE3_CXXFLAGS="-msse3"],
   				   [AC_MSG_ERROR([Need a version of gcc with -$flag])])
fi

if test "$have_avx" = "yes" -a "x$AVX_CXXFLAGS" = x; then
   AX_CHECK_COMPILER_FLAGS(-mavx, [AVX_CXXFLAGS="-mavx"],
                		  [AC_MSG_ERROR([Need a version of gcc with -mavx])])
fi

if test "$have_altivec" = "yes" -a "x$ALTIVEC_CXXFLAGS" = x; then
   # -DFAKE__VEC__ is a workaround because gcc-3.3 does not                                                
   # #define __VEC__ with -maltivec.                                                                       
   AX_CHECK_COMPILER_FLAGS(-faltivec, [ALTIVEC_CXXFLAGS="-faltivec"],
                		      [AX_CHECK_COMPILER_FLAGS(-maltivec -mabi=altivec,
                    		      [ALTIVEC_CXXFLAGS="-maltivec -mabi=altivec -DFAKE__VEC__"],
                    		      [AX_CHECK_COMPILER_FLAGS(-fvec, [ALTIVEC_CXXFLAGS="-fvec"],
                        	      [AC_MSG_ERROR([Need a version of gcc with -maltivec])])])])
fi

if test "$have_neon" = "yes" -a "x$NEON_CXXFLAGS" = x; then
   AX_CHECK_COMPILER_FLAGS(-mfpu=neon, [NEON_CXXFLAGS="-mfpu=neon"],
                		       [AC_MSG_ERROR([Need a version of gcc with -mfpu=neon])])
fi

AC_SUBST(SSE3_CXXFLAGS)
AC_SUBST(AVX_CXXFLAGS)
AC_SUBST(ALTIVEC_CXXFLAGS)
AC_SUBST(NEON_CXXFLAGS)

AC_OUTPUT
