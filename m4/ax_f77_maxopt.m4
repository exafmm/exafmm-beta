dnl @synopsis AX_CXX_MAXOPT
dnl @summary turn on optimization flags for the C compiler
dnl @category C
dnl
dnl Try to turn on "good" C optimization flags for various compilers
dnl and architectures, for some definition of "good".  (In our case,
dnl good for FFTW and hopefully for other scientific codes.  Modify
dnl as needed.)
dnl
dnl The user can override the flags by setting the CXXFLAGS environment
dnl variable.  The user can also specify --enable-portable-binary in
dnl order to disable any optimization flags that might result in
dnl a binary that only runs on the host architecture.
dnl
dnl Note also that the flags assume that ANSI C aliasing rules are
dnl followed by the code (e.g. for gcxx's -fstrict-aliasing), and that
dnl floating-point computations can be re-ordered as needed.
dnl
dnl Requires macros: AX_CHECK_COMPILER_FLAGS, AX_COMPILER_VENDOR,
dnl                  AX_GCC_ARCHFLAG, AX_GCC_X86_CPUID
dnl
dnl @version 2005-05-30
dnl @license GPLWithACException
dnl @author Steven G. Johnson <stevenj@alum.mit.edu> and Matteo Frigo.
AC_DEFUN([AX_F77_MAXOPT],
[
AC_REQUIRE([AC_PROG_F77])
AC_REQUIRE([AC_CANONICAL_HOST])
AC_LANG_PUSH([Fortran 77])

AX_COMPILER_VENDOR

AC_ARG_ENABLE(portable-binary, [AC_HELP_STRING([--enable-portable-binary], [disable compiler optimizations that would produce unportable binaries])],
			       acx_maxopt_portable=$withval, acx_maxopt_portable=yes)

# Try to determine "good" native compiler flags if none specified via FFLAGS
if test "$ac_test_FFLAGS" != "set"; then
  FFLAGS=""
  case $ax_cv_c_compiler_vendor in
    dec) FFLAGS="-newc -w0 -O5 -ansi_alias -ansi_args -fp_reorder -tune host"
     if test "x$acx_maxopt_portable" = xno; then
           FFLAGS="$FFLAGS -arch host"
         fi;;

    sun) FFLAGS="-native -fast -xO5 -dalign"
     if test "x$acx_maxopt_portable" = xyes; then
        FFLAGS="$FFLAGS -xarch=generic"
         fi;;

    hp)  FFLAGS="+Oall +Optrs_ansi +DSnative"
     if test "x$acx_maxopt_portable" = xyes; then
        FFLAGS="$FFLAGS +DAportable"
	 fi;;

    ibm) if test "x$acx_maxopt_portable" = xno; then
           xlc_opt="-qarch=auto -qtune=auto"
	    else
           xlc_opt="-qtune=auto"
	    fi
         AX_CHECK_COMPILER_FLAGS($xlc_opt,
		FFLAGS="-O3 -qansialias -w $xlc_opt",
               [FFLAGS="-O3 -qansialias -w"
                echo "******************************************************"
                echo "*  You seem to have the IBM  C compiler.  It is      *"
                echo "*  recommended for best performance that you use:    *"
                echo "*                                                    *"
                echo "*    FFLAGS=-O3 -qarch=xxx -qtune=xxx -qansialias -w *"
                echo "*                      ^^^        ^^^                *"
                echo "*  where xxx is pwr2, pwr3, 604, or whatever kind of *"
                echo "*  CPU you have.  (Set the FFLAGS environment var.   *"
                echo "*  and re-run configure.)  For more info, man cxx.    *"
                echo "******************************************************"])
         ;;

    intel) FFLAGS="-O3 -ansi_alias"
    if test "x$acx_maxopt_portable" = xno; then
      icc_archflag=unknown
        icc_flags=""
	  case $host_cpu in
	      i686*|x86_64*)
              # icc accepts gcc assembly syntax, so these should work:
	            AX_GCC_X86_CPUID(0)
              AX_GCC_X86_CPUID(1)
	            case $ax_cv_gcc_x86_cpuid_0 in # see AX_GCC_ARCHFLAG
                *:756e6547:*:*) # Intel
                  case $ax_cv_gcc_x86_cpuid_1 in
                    106a?:*:*:*) icc_flags="-xSSE4.2";;
                    1067?:*:*:*) icc_flags="-xS";;
                    *6a?:*[[234]]:*:*|*6[[789b]]?:*:*:*) icc_flags="-xK";;
                    *6f?:*:*:*) icc_flags="-xT";;
                    *f3[[347]]:*:*:*|*f4[[1347a]]:*:*:*) icc_flags="-xP -xN -xW -xK";;
                    *f??:*:*:*) icc_flags="-xN -xW -xK";;
                  esac ;;
              esac ;;
          esac
          if test "x$icc_flags" != x; then
            for flag in $icc_flags; do
              AX_CHECK_COMPILER_FLAGS($flag, [icc_archflag=$flag; break])
            done
          fi
          AC_MSG_CHECKING([for icc architecture flag])
	    AC_MSG_RESULT($icc_archflag)
          if test "x$icc_archflag" != xunknown; then
            FFLAGS="$FFLAGS $icc_archflag"
          fi
        fi
	;;

    gnu)
     # default optimization flags for gcc on all systems
     FFLAGS="-O3 -fomit-frame-pointer"

     # -malign-double for x86 systems
     AX_CHECK_COMPILER_FLAGS(-malign-double, FFLAGS="$FFLAGS -malign-double")

     #  -fstrict-aliasing for gcc-2.95+
     AX_CHECK_COMPILER_FLAGS(-fstrict-aliasing,
     FFLAGS="$FFLAGS -fstrict-aliasing")

     # note that we enable "unsafe" fp optimization with other compilers, too
     AX_CHECK_COMPILER_FLAGS(-ffast-math, FFLAGS="$FFLAGS -ffast-math")

     AX_GCC_ARCHFLAG($acx_maxopt_portable)
     ;;
  esac

  if test -z "$FFLAGS"; then
  echo ""
  echo "**********************************************************"
        echo "* WARNING: Don't know the best FFLAGS for this system  *"
        echo "* Use ./configure FFLAGS=... to specify your own flags *"
	echo "* (otherwise, a default of FFLAGS=-O3 will be used)    *"
	echo "**********************************************************"
	echo ""
        FFLAGS="-O3"
  fi

  AX_CHECK_COMPILER_FLAGS($FFLAGS, [], [
  echo ""
        echo "**********************************************************"
        echo "* WARNING: The guessed FFLAGS don't seem to work with  *"
        echo "* your compiler.                                         *"
        echo "* Use ./configure FFLAGS=... to specify your own flags *"
        echo "**********************************************************"
        echo ""
        FFLAGS=""
  ])

fi
AC_LANG_POP([Fortran 77])
])
