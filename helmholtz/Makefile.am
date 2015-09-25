include ../Makefile.am.include

bin_PROGRAMS = helmholtz_spherical
helmholtz_spherical_SOURCES = fmm.cxx kernel.cxx
helmholtz_spherical_CPPFLAGS = $(AM_CPPFLAGS) -DHelmholtz -DSpherical -DEXPANSION=40 -DFP64
