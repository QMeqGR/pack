## PACK - a Monte Carlo structure searching program

## Place generated object files (.o) into the same directory as their source
## files, in order to avoid collisions when non-recursive make is used.
AUTOMAKE_OPTIONS = subdir-objects

## Additional flags to pass to aclocal when it is invoked automatically at
## make time. The ${ACLOCAL_FLAGS} variable is picked up from the environment
## to provide a way for the user to supply additional arguments.
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

## Define executables which will be installed into the
## directory named by the predefined variable $(bindir).
bin_PROGRAMS = pack ewald

## Define the list of source files for the targets. The file extension
## .cc is recognized by Automake, and causes it to produce rules which invoke
## the C++ compiler to produce an object file (.o) from each source file. The
## header files (.h) do not result in object files by themselves, but will be
## included in distribution archives of the project.
pack_SOURCES    = packlib.c packlib.h main.c amoeba_ehm.c amotry_ehm.c printing.c metropolis.c \
                  en_rep.c dsm.c simplex.c autoadjust_pCnt.c init_parms.c init_simplex.c \
                  simp_trans.c init_estat.c init_energ.c compare.c tot_energy.c \
                  hard_dist.c init_cell.c matrix.c make_obj.c np_force_mv.c \
                  rescale_object.c lat_pathological.c translate.c rotate_object.c dist_fcns.c \
                  cell_basis.c ewald_sums.c chk_bounds.c rezone.c pack_frac.c \
                  pressure.c autoadjust_rejpct.c force_npbc.c cent_of_mass.c conversions.c \
                  wang_landau.c center_np.c lnfmod.c wall_check.c global_defs.h global_vars.h \
                  min-map.c chk_constr.c

ewald_SOURCES   = ewald.c matrix.c packlib.c ewald_sums.c \
                  cell_basis.c rezone.c

#AM_CFLAGS += $(MPI_CFLAGS)
#LINK += $(MPI_CLIBS)

AM_CFLAGS = -I$(top_srcdir)/src/
