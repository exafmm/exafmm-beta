# ========================================================================================
#  Originally http://www.gnu.org/software/autoconf-archive/ax_compiler_flags_cxxflags.html
# ========================================================================================
#
# SYNOPSIS
#
#   AX_COMPILER_FLAGS()
#
# DESCRIPTION
#
#   Add warning flags for given compiler to VARIABLE, which defaults to
#   ax_compiler_c/cxx/fcflags.  VARIABLE is AC_SUBST-ed by
#   this macro, but must be manually added to the CFLAGS/CXXFLAGS/FCFLAGS
#   variable for each target in the code base.
#
#   This macro depends on the environment set up by AX_COMPILER_FLAGS.
#
# LICENSE
#
#   Copyright (c) 2015 David King <amigadave@amigadave.com>
#   Copyright (c) 2016 Rio Yokota <rioyokota@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.  This file is offered as-is, without any
#   warranty.

#serial 8

AC_DEFUN([AX_COMPILER_FLAGS],[
    AC_REQUIRE([AC_PROG_SED])

    # We need to turn warnings to errors for AX_APPEND_COMPILE_FLAGS to be able to work.
    # Different compilers require different flags for this:
    # GNU: -Werror
    # Intel: -diag-error warn
    # Clang: -Werror=unknown-warning-option
    AX_CHECK_COMPILE_FLAG([-Werror=unknown-warning-option],[
        ax_compiler_flags_test="-Werror=unknown-warning-option"
    ],[
        ax_compiler_flags_test="-Werror"
    ])
    AX_CHECK_COMPILE_FLAG([-diag-error warn],[
        ax_compiler_flags_test="-diag-error warn"
    ])

    # Base flags
    AX_APPEND_COMPILE_FLAGS([dnl
        -fno-strict-aliasing dnl
    ],ax_compiler_[]_AC_LANG_ABBREV[]flags,[$ax_compiler_flags_test])

    AS_IF([test "$ax_enable_compile_warnings" != "no"],[
        # "yes" flags
        AX_APPEND_COMPILE_FLAGS([dnl
            -Wall dnl
            -Wextra dnl
            -Wundef dnl
            -Wwrite-strings dnl
            -Wpointer-arith dnl
            -Wmissing-declarations dnl
            -Wredundant-decls dnl
            -Wno-unused-parameter dnl
            -Wno-missing-field-initializers dnl
            -Wformat=2 dnl
            -Wcast-align dnl
            -Wformat-nonliteral dnl
            -Wformat-security dnl
            -Wsign-compare dnl
            -Wstrict-aliasing dnl
            -Wshadow dnl
            -Winline dnl
            -Wpacked dnl
            -Wmissing-format-attribute dnl
            -Wmissing-noreturn dnl
            -Winit-self dnl
            -Wredundant-decls dnl
            -Wmissing-include-dirs dnl
            -Wunused-but-set-variable dnl
            -Warray-bounds dnl
            -Wreturn-type dnl
            -Wno-overloaded-virtual dnl
            -Wswitch-enum dnl
            -Wswitch-default dnl
        ],ax_compiler_[]_AC_LANG_ABBREV[]flags,[$ax_compiler_flags_test])
    ])

    # In the flags below, when disabling specific flags, always add *both*
    # -Wno-foo and -Wno-error=foo. This fixes the situation where (for example)
    # we enable -Werror, disable a flag, and a build bot passes C/CXX/FCFLAGS=-Wall,
    # which effectively turns that flag back on again as an error.
    for flag in $ax_compiler_[]_AC_LANG_ABBREV[]flags; do
        AS_CASE([$flag],
                [-Wno-*=*],[],
                [-Wno-*],[
                    AX_APPEND_COMPILE_FLAGS([-Wno-error=$(AS_ECHO([$flag]) | $SED 's/^-Wno-//')],
                                            ax_compiler_[]_AC_LANG_ABBREV[]flags,
                                            [$ax_compiler_flags_test])
                ])
    done

    # Substitute the variables
    AC_SUBST(ax_compiler_[]_AC_LANG_ABBREV[]flags)
])dnl AX_COMPILER_FLAGS
