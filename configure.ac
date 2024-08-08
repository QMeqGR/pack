AC_INIT([PACK], m4_esyscmd([./scripts/get_version | tr -d '\n']))
AC_CONFIG_SRCDIR(src/main.c)
AM_INIT_AUTOMAKE([-Wall no-define foreign])


#AC_SUBST([CONFIG_STATUS_DEPENDENCIES], ['$(top_srcdir)/src/main.c'])

AC_DEFINE([HAVE_MPI], [], [Enable MPI])

m4_include(acx_mpi.m4)

AC_ARG_WITH(mpi, [AC_HELP_STRING([--with-mpi], 
                 [enable MPI support [default=yes]])], 
                 [case "${withval}" in 
                       yes|no) with_mpi=$withval;; 
                       *)
                       AC_MSG_ERROR(bad value ${withval} for --with-mpi);; 
                  esac], 
                 [with_mpi=yes])

AC_PROG_CC

if test "x$with_mpi" = "xyes"; then 
        ACX_MPI([], [AC_MSG_ERROR([could not find mpi library for --with-mpi])]) 
        #AX_MPI
        AC_DEFINE(HAVE_MPI) 
        MPI_CLIBS=`mpicc --showme:link` 
        MPI_CFLAGS=`mpicc --showme:compile` 
        AC_SUBST(MPI_CLIBS) 
        AC_SUBST(MPI_CFLAGS) 
        CC="$MPICC"
else 
  MPICC="$CC" 
  #MPICXX="$CXX" 
  AC_SUBST(MPICC) 
  #AC_SUBST(MPICXX) 
fi

LIBS="$MPILIBS $LIBS"

# Checks for libraries.
# Link math lib
AC_CHECK_LIB(m, pow)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([memory.h stdlib.h string.h wchar.h wctype.h math.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_STRUCT_TM

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memmove memset strcasecmp strncasecmp strstr])

AC_CONFIG_HEADERS([src/config.h])

AM_PROG_CC_C_O

AC_CONFIG_FILES([
  Makefile
  src/Makefile
])
AC_OUTPUT
