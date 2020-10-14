dnl $Id$

dnl Process this file with autoconf to produce a configure script.
dnl
dnl usage : aclocal -I config/ && autoconf && ./configure && make 

dnl then  : VAR=VAL ./configure

dnl In some configurations (AIX) the next line is needed:
dnl MPIFC=mpxlf95 ./configure

dnl then  : ./configure VAR=VAL 
dnl then  : ./configure --help=short
dnl then  : ./configure --help

dnl the PSBLAS modules get this task difficult to accomplish!
dnl SEE  : --module-path --include-path

dnl NOTE : There is no cross compilation support.

dnl NOTE : missing ifort and kl* library handling..
dnl NOTE : odd configurations like ifc + gcc still await in the mist of the unknown


###############################################################################
###############################################################################
# 
#  This script is used by the PSBLAS to determine the compilers, linkers, and
# libraries to build its libraries executable code.
#  Its behaviour is driven on the compiler it finds or it is dictated to work
# with.
#
###############################################################################
###############################################################################

# NOTE: the literal for version (the second argument to AC_INIT should be a literal!)
AC_INIT([AMG4PSBLAS],1.0.0, [https://github.com/sfilippone/amg4psblas/issues])

# VERSION is the file containing the PSBLAS version code
# FIXME
amg4psblas_cv_version="1.0.0"

# A sample source file
AC_CONFIG_SRCDIR([mlprec/amg_prec_type.f90])
# Our custom M4 macros are in the 'config' directory
AC_CONFIG_MACRO_DIR([config])
AC_MSG_NOTICE([
--------------------------------------------------------------------------------
	Welcome to the $PACKAGE_NAME $amg4psblas_cv_version configure Script.

	This creates Make.inc, but if you read carefully the
	documentation, you can make your own by hand for your needs.

        ./configure  --with-psblas=/path/to/psblas
	See ./configure --help=short fore more info.
--------------------------------------------------------------------------------
		])

###############################################################################
# 			FLAGS and LIBS user customization
###############################################################################

dnl NOTE : no spaces before the comma, and no brackets before the second argument!
PAC_ARG_WITH_PSBLAS
PSBLAS_DIR="$pac_cv_psblas_dir";
PSBLAS_INCDIR="$pac_cv_psblas_incdir";
PSBLAS_MODDIR="$pac_cv_psblas_moddir";
PSBLAS_LIBDIR="$pac_cv_psblas_libdir";
AC_MSG_CHECKING([for PSBLAS install dir])
if test "X$PSBLAS_DIR" != "X" ; then 
   case $PSBLAS_DIR in 
   	  /*) ;; 
	   *)  AC_MSG_ERROR([The PSBLAS installation dir must be an absolute pathname
 specified with --with-psblas=/path/to/psblas])
 esac
  if test ! -d "$PSBLAS_DIR" ; then 
    AC_MSG_ERROR([Could not find PSBLAS build dir $PSBLAS_DIR!])	
  fi
  AC_MSG_RESULT([$PSBLAS_DIR])
  pac_cv_status_file="$PSBLAS_INCDIR/Make.inc.psblas"
  if test ! -f "$pac_cv_status_file" ; then 
    AC_MSG_NOTICE([Could not find an installation  in $PSBLAS_DIR.])
    pac_cv_status_file="NONE";
  fi
else
  pac_cv_status_file="NONE";
fi	

dnl. $pac_cv_status_file
AC_MSG_NOTICE([Loaded $pac_cv_status_file $FC $MPIFC $BLACS_LIBS])
AM_INIT_AUTOMAKE
dnl Specify required version of autoconf.
AC_PREREQ(2.59)
#
# Installation. 
#
#
AC_PROG_INSTALL

AC_MSG_CHECKING([where to install])
case $prefix in  
   \/* )   eval "INSTALL_DIR=$prefix";;
   * ) eval "INSTALL_DIR=/usr/local/amg4psblas";;
esac
case $libdir in 
   \/* )   eval "INSTALL_LIBDIR=$libdir";;
   * ) eval "INSTALL_LIBDIR=$INSTALL_DIR/lib";;
esac
case $includedir in 
   \/* )   eval "INSTALL_INCLUDEDIR=$includedir";;
   * ) eval "INSTALL_INCLUDEDIR=$INSTALL_DIR/include";;
esac
INSTALL_MODULESDIR=$INSTALL_DIR/modules
case $docsdir in 
   \/* )   eval "INSTALL_DOCSDIR=$docsdir";;
   * ) eval "INSTALL_DOCSDIR=$INSTALL_DIR/docs";;
esac
case $samplesdir in 
   \/* )   eval "INSTALL_SAMPLESDIR=$samplesdir";;
   * ) eval "INSTALL_SAMPLESDIR=$INSTALL_DIR/samples";;
esac
AC_MSG_RESULT([$INSTALL_DIR $INSTALL_INCLUDEDIR $INSTALL_MODULESDIR $INSTALL_LIBDIR $INSTALL_DOCSDIR $INSTALL_SAMPLESDIR])

dnl
dnl We set our own FC flags, ignore those from AC_PROG_FC but not those from the
dnl environment variable. Same for C
dnl 
save_FCFLAGS="$FCFLAGS";
AC_PROG_FC([ftn xlf2003_r xlf2003 xlf95_r xlf95 xlf90 xlf pgf95 pgf90 ifort ifc  nagfor gfortran])
FCFLAGS="$save_FCFLAGS";
save_CFLAGS="$CFLAGS";
AC_PROG_CC([xlc pgcc icc gcc cc ])
CFLAGS="$save_CFLAGS";
dnl AC_PROG_CXX

dnl AC_PROG_F90 doesn't exist, at the time of writing this !
dnl AC_PROG_F90

# Sanity checks, although redundant (useful when debugging this configure.ac)!
if test "X$FC" == "X" ; then
	AC_MSG_ERROR([Problem : No Fortran compiler specified nor found!])
fi

if test "X$CC" == "X" ; then
	AC_MSG_ERROR([Problem : No C compiler specified nor found!])
fi
if eval "$FC -qversion 2>&1 | grep XL 2>/dev/null" ; then
	# Some configurations of the XLF want "-WF," prepended to -D.. flags.
	# TODO : discover the exact conditions when the usage of -WF is needed.
	amg_cv_define_prepend="-WF,"
        if eval "$MPIFC -qversion 2>&1 | grep -e\"Version: 10\.\" 2>/dev/null"; then
        	FDEFINES="$amg_cv_define_prepend-DXLF_10 $FDEFINES"		  
        fi

 # Note : there could be problems with old xlf compiler versions ( <10.1 )
 # since (as far as it is known to us) -WF, is not used in earlier versions.
 # More problems could be undocumented yet.
fi 
###############################################################################
#			Suitable MPI compilers detection
###############################################################################
# Note: Someday we will contemplate a fake MPI - configured version of PSBLAS
###############################################################################
# First check whether the user required our serial (fake) mpi.
PAC_ARG_SERIAL_MPI

#Note : we miss the name of the Intel C compiler
if test x"$pac_cv_serial_mpi" == x"yes" ; then
   FAKEMPI="fakempi.o";
   MPIFC="$FC";
   MPICC="$CC";
else 
AC_LANG([C])
if test "X$MPICC" = "X" ; then
    # This is our MPICC compiler preference: it will override ACX_MPI's first try.
    AC_CHECK_PROGS([MPICC],[mpxlc mpiicc mpcc mpicc cc])
fi
ACX_MPI([], [AC_MSG_ERROR([[Cannot find any suitable MPI implementation for C]])])


AC_LANG([Fortran])

if test "X$MPIFC" = "X" ; then
    # This is our MPIFC compiler preference: it will override ACX_MPI's first try.
    AC_CHECK_PROGS([MPIFC],[mpxlf2003_r mpxlf2003 mpxlf95_r mpxlf90  mpiifort mpf95 mpf90 mpifort mpif95 mpif90 ftn ])
fi

ACX_MPI([], [AC_MSG_ERROR([[Cannot find any suitable MPI implementation for Fortran]])])

FC="$MPIFC" ;
CC="$MPICC";
fi

AC_LANG([C])

dnl Now on, MPIFC should be set, and MPICC

###############################################################################
# Sanity checks, although redundant (useful when debugging this configure.ac)!
###############################################################################

if test "X$MPIFC" == "X" ; then
	AC_MSG_ERROR([Problem : No MPI Fortran compiler specified nor found!])
fi

if test "X$MPICC" == "X" ; then
	AC_MSG_ERROR([Problem : No MPI C compiler specified nor found!])
fi

###############################################################################
# 			FLAGS and LIBS user customization
###############################################################################

dnl NOTE : no spaces before the comma, and no brackets before the second argument!
PAC_ARG_WITH_FLAGS(ccopt,CCOPT)
PAC_ARG_WITH_FLAGS(fcopt,FCOPT)
PAC_ARG_WITH_LIBS
PAC_ARG_WITH_FLAGS(clibs,CLIBS)
PAC_ARG_WITH_FLAGS(flibs,FLIBS)

dnl candidates for removal:
PAC_ARG_WITH_FLAGS(library-path,LIBRARYPATH)
PAC_ARG_WITH_FLAGS(include-path,INCLUDEPATH)
PAC_ARG_WITH_FLAGS(module-path,MODULE_PATH)

# we just gave the user the chance to append values to these variables
PAC_ARG_WITH_EXTRA_LIBS

###############################################################################
# Sanity checks, although redundant (useful when debugging this configure.ac)!
###############################################################################

###############################################################################
#	Compiler identification (sadly, it is necessary)
###############################################################################
psblas_cv_fc=""

dnl Do we use gfortran & co ? Compiler identification.
dnl NOTE : in /autoconf/autoconf/fortran.m4 there are plenty of better tests!

PAC_CHECK_HAVE_GFORTRAN(
	[psblas_cv_fc="gcc"],
)

PAC_CHECK_HAVE_CRAYFTN(
	[psblas_cv_fc="cray"],
)

if test x"$psblas_cv_fc" == "x" ; then
 if eval "$MPIFC -qversion 2>&1 | grep XL 2>/dev/null" ; then
	psblas_cv_fc="xlf"
	# Some configurations of the XLF want "-WF," prepended to -D.. flags.
	# TODO : discover the exact conditions when the usage of -WF is needed.
	psblas_cv_define_prepend="-WF,"
        if eval "$MPIFC -qversion 2>&1 | grep -e\"Version: 10\.\" 2>/dev/null"; then
        	FDEFINES="$psblas_cv_define_prepend-DXLF_10 $FDEFINES"		  
        fi

 # Note : there could be problems with old xlf compiler versions ( <10.1 )
 # since (as far as it is known to us) -WF, is not used in earlier versions.
 # More problems could be undocumented yet.
 elif eval "$MPIFC -V  2>&1 | grep Sun 2>/dev/null" ; then
   # Sun compiler detection
  
   psblas_cv_fc="sun"
 elif eval "$MPIFC -V  2>&1 | grep Portland 2>/dev/null" ; then
   # Portland group compiler detection
   
   psblas_cv_fc="pg"
 elif eval "$MPIFC -V 2>&1 | grep Intel.*Fortran.*Compiler 2>/dev/null" ; then
   # Intel compiler identification
    
   psblas_cv_fc="ifc"
 elif eval "$MPIFC -v  2>&1 | grep NAG 2>/dev/null" ; then
   psblas_cv_fc="nag"
   FC="$MPIFC"
 else
   psblas_cv_fc=""
  # unsupported MPI Fortran compiler
     AC_MSG_NOTICE([[Unknown Fortran compiler, proceeding with fingers crossed !]])
 fi
fi
if test "X$psblas_cv_fc" == "Xgcc" ; then
PAC_HAVE_MODERN_GFORTRAN(
	[],
        [AC_MSG_ERROR([Bailing out.])]
)
fi


###############################################################################
#	Linking, symbol mangling, and misc tests
###############################################################################

# Note : This is functional to Make.inc rules and structure (see below).
AC_LANG([C])
AC_CHECK_SIZEOF(void *)
# Define for platforms with 64 bit (void *     ) pointers
if test X"$ac_cv_sizeof_void_p" == X"8" ; then
  CDEFINES="-DPtr64Bits $CDEFINES"
fi
AC_LANG([Fortran])
__AC_FC_NAME_MANGLING
if test "X$psblas_cv_fc" == X"pg" ; then
   FC=$save_FC
fi
AC_LANG([C])
dnl AC_MSG_NOTICE([Fortran name mangling: $ac_cv_fc_mangling])
[pac_fc_case=${ac_cv_fc_mangling%%,*}]
[pac_fc_under=${ac_cv_fc_mangling#*,}]
[pac_fc_sec_under=${pac_fc_under#*,}]
[pac_fc_sec_under=${pac_fc_sec_under# }]
[pac_fc_under=${pac_fc_under%%,*}]
[pac_fc_under=${pac_fc_under# }]
AC_MSG_CHECKING([defines for C/Fortran name interfaces])
if test "x$pac_fc_case" == "xlower case"; then
   if test "x$pac_fc_under" == "xunderscore"; then
      if test "x$pac_fc_sec_under" == "xno extra underscore"; then 
          pac_f_c_names="-DLowerUnderscore"
      elif test "x$pac_fc_sec_under" == "xextra underscore"; then
          pac_f_c_names="-DLowerDoubleUnderscore"
      else
          pac_f_c_names="-DUNKNOWN"
dnl  AC_MSG_NOTICE([Fortran name mangling extra underscore unknown case])	
      fi
   elif test "x$pac_fc_under" == "xno underscore"; then
      pac_f_c_names="-DLowerCase"
   else		       
      pac_f_c_names="-DUNKNOWN"
dnl  AC_MSG_NOTICE([Fortran name mangling underscore unknown case])	
   fi
elif test "x$pac_fc_case" == "xupper case"; then
   if test "x$pac_fc_under" == "xunderscore"; then
      if test "x$pac_fc_sec_under" == "xno extra underscore"; then 
          pac_f_c_names="-DUpperUnderscore"
      elif test "x$pac_fc_sec_under" == "xextra underscore"; then
          pac_f_c_names="-DUpperDoubleUnderscore"
      else
          pac_f_c_names="-DUNKNOWN"
dnl  AC_MSG_NOTICE([Fortran name mangling extra underscore unknown case])	
      fi
   elif test "x$pac_fc_under" == "xno underscore"; then
      pac_f_c_names="-DUpperCase"
   else		       
      pac_f_c_names="-DUNKNOWN"
dnl  AC_MSG_NOTICE([Fortran name mangling underscore unknown case])	
   fi
dnl  AC_MSG_NOTICE([Fortran name mangling UPPERCASE not handled])	
else
   pac_f_c_names="-DUNKNOWN"   
dnl   AC_MSG_NOTICE([Fortran name mangling unknown case])	
fi
CDEFINES="$pac_f_c_names $CDEFINES"

AC_MSG_RESULT([ $pac_f_c_names ])

###############################################################################
# 			Make.inc generation logic
###############################################################################
# Honor CFLAGS if they were specified explicitly, but --with-ccopt take precedence
if test "X$CCOPT" == "X" ; then
  CCOPT="$CFLAGS";
fi	
if test "X$CCOPT" == "X" ; then
 if test "X$psblas_cv_fc" == "Xgcc" ; then
 # note that no space should be placed around the equality symbol in assignements
 # Note : 'native' is valid _only_ on GCC/x86  (32/64 bits)
   CCOPT="-O3  $CCOPT"

 elif test "X$psblas_cv_fc" == X"xlf" ; then
   # XL compiler : consider using -qarch=auto
     CCOPT="-O3 -qarch=auto $CCOPT"
 elif test "X$psblas_cv_fc" == X"ifc" ; then
   # other compilers ..
     CCOPT="-O3 $CCOPT"
 elif test "X$psblas_cv_fc" == X"pg" ; then
   # other compilers ..
     CCOPT="-fast $CCOPT"
  # NOTE : PG & Sun use -fast instead -O3
 elif test "X$psblas_cv_fc" == X"sun" ; then
   # other compilers ..
     CCOPT="-fast $CCOPT"
 elif test "X$psblas_cv_fc" == X"cray" ; then
     CCOPT="-O3 $CCOPT"
     MPICC="cc"
 elif test "X$psblas_cv_fc" == X"nag" ; then
   # using GCC in conjunction with NAG.
     CCOPT="-O2"
 else 
     CCOPT="-O2 $CCOPT"
 fi
fi
#CFLAGS="${CCOPT}"

# Honor FCFLAGS if they were specified explicitly, but --with-fcopt take precedence
if test "X$FCOPT" == "X" ; then
  FCOPT="$FCFLAGS";
fi	
if test "X$FCOPT" == "X" ; then
 if test "X$psblas_cv_fc" == "Xgcc" ; then
 # note that no space should be placed around the equality symbol in assignations
 # Note : 'native' is valid _only_ on GCC/x86  (32/64 bits)
   FCOPT="-O3  $FCOPT"
 elif test "X$psblas_cv_fc" == X"xlf" ; then
   # XL compiler : consider using -qarch=auto
     FCOPT="-O3 -qarch=auto  -qlanglvl=extended -qxlf2003=polymorphic:autorealloc $FCOPT"
     FCFLAGS="-qhalt=e -qlanglvl=extended -qxlf2003=polymorphic:autorealloc $FCFLAGS"
 elif test "X$psblas_cv_fc" == X"ifc" ; then
   # other compilers ..
     FCOPT="-O3 $FCOPT"
 elif test "X$psblas_cv_fc" == X"pg" ; then
   # other compilers ..
     FCOPT="-fast $FCOPT"
  # NOTE : PG & Sun use -fast instead -O3
 elif test "X$psblas_cv_fc" == X"sun" ; then
   # other compilers ..
     FCOPT="-fast $FCOPT"
 elif test "X$psblas_cv_fc" == X"cray" ; then
     FCOPT="-O3 -em $FCOPT"
 elif test "X$psblas_cv_fc" == X"nag" ; then
   # NAG  compiler ..
     FCOPT="-O2 "
  # NOTE : PG & Sun use -fast instead -O3
 else 
     FCOPT="-O2 $FCOPT"
 fi
fi
if test "X$psblas_cv_fc" == X"nag" ; then
   # Add needed options 
     FCOPT="$FCOPT -dcfuns -f2003 -wmismatch=mpi_scatterv,mpi_alltoallv,mpi_gatherv,mpi_allgatherv"
     EXTRA_OPT="-mismatch_all"
fi


# COPT,FCOPT are aliases for CFLAGS,FCFLAGS .

##############################################################################
# Compilers variables selection 
##############################################################################
FC=${FC}
CC=${CC}
MPCC=${MPICC}


##############################################################################
# Choice of our compilers, needed by Make.inc
##############################################################################
if test "X$psblas_cv_fc" == X"cray" 
then
   MODEXT=".mod"
   FMFLAG="-I"
   FIFLAG="-I"
   BASEMODNAME=PSB_BASE_MOD
   PRECMODNAME=PSB_PREC_MOD
   METHDMODNAME=PSB_KRYLOV_MOD
   UTILMODNAME=PSB_UTIL_MOD

else 
   AX_F90_MODULE_EXTENSION
   AX_F90_MODULE_FLAG
   MODEXT=".$ax_cv_f90_modext"
   FMFLAG="${ax_cv_f90_modflag%%[ ]*}"
   FIFLAG=-I
   BASEMODNAME=psb_base_mod
   PRECMODNAME=psb_prec_mod
   METHDMODNAME=psb_krylov_mod
   UTILMODNAME=psb_util_mod
fi


##############################################################################
# Choice of our compilers, needed by Make.inc
##############################################################################
if test "X$FLINK" == "X" ; then
   FLINK=${MPF90}
fi

PAC_FORTRAN_HAVE_PSBLAS([AC_MSG_RESULT([yes.])],
	[AC_MSG_ERROR([no.  Could not find working version of PSBLAS.])])

PAC_FORTRAN_PSBLAS_VERSION()

if test "x$pac_cv_psblas_major" == "xunknown"; then
  AC_MSG_ERROR([PSBLAS version major "$pac_cv_psblas_major".])
fi
if test "x$pac_cv_psblas_minor" == "xunknown"; then
   AC_MSG_ERROR([PSBLAS version minor "$pac_cv_psblas_minor".])
fi
if test "x$pac_cv_psblas_patchlevel" == "xunknown"; then
   AC_MSG_ERROR([PSBLAS patchlevel "$pac_cv_psblas_patchlevel".])
fi
if (( $pac_cv_psblas_major < 3 )) ||
   (  (( $pac_cv_psblas_major == 3 )) && (( $pac_cv_psblas_minor < 7 ))) ; then
  AC_MSG_ERROR([I need at least PSBLAS version 3.7.]) 
else
  AC_MSG_NOTICE([Am configuring with PSBLAS version $pac_cv_psblas_major.$pac_cv_psblas_minor.$pac_cv_psblas_patchlevel.])
fi
PAC_FORTRAN_PSBLAS_INTEGER_SIZES()
AC_MSG_NOTICE([PSBLAS size of LPK "$pac_cv_psblas_lpk".])
###############################################################################
#	Parachute rules for ar and ranlib ... (could cause problems)
###############################################################################

if test "X$AR" == "X" ; then
  AR="ar"
fi

if test "X$RANLIB" == "X" ; then
  RANLIB="ranlib"
fi

# This should be portable
AR="${AR} -cur"


###############################################################################
# NOTE :
# Missing stuff : 
# 		In the case the detected fortran compiler is ifort, icc or gcc
#		should be valid options.
#		The same for pg (Portland Group compilers).
###############################################################################


###############################################################################
# Custom test : do we have a module or include for MPI Fortran interface?
if test x"$pac_cv_serial_mpi" == x"yes" ; then
   FDEFINES="$psblas_cv_define_prepend-DSERIAL_MPI $psblas_cv_define_prepend-DMPI_MOD $FDEFINES";
else 
    PAC_FORTRAN_CHECK_HAVE_MPI_MOD(
 	  [FDEFINES="$psblas_cv_define_prepend-DMPI_MOD $FDEFINES"],
	  [FDEFINES="$psblas_cv_define_prepend-DMPI_H $FDEFINES"])	  
fi

#
# Tests for support of various Fortran features; some of them are critical,
# some optional
#

PAC_FORTRAN_TEST_EXTENDS( 
	[],
        [AC_MSG_ERROR([Sorry, cannot build PSBLAS without support for EXTENDS. 
 Please get a Fortran compiler that supports it, e.g. GNU Fortran 4.8.])]
)

PAC_FORTRAN_TEST_CLASS_TBP( 
	[],
        [AC_MSG_ERROR([Sorry, cannot build PSBLAS without support for CLASS and type bound procedures. 
 Please get a Fortran compiler that supports them, e.g. GNU Fortran 4.8.])]
)

PAC_FORTRAN_TEST_SOURCE( 
	[],
        [AC_MSG_ERROR([Sorry, cannot build PSBLAS without support for SOURCE= allocation. 
 Please get a Fortran compiler that supports it, e.g. GNU Fortran 4.8.])]
)

PAC_FORTRAN_HAVE_MOVE_ALLOC(
	[],
        [AC_MSG_ERROR([Sorry, cannot build PSBLAS without support for MOVE_ALLOC. 
 Please get a Fortran compiler that supports it, e.g. GNU Fortran 4.8.])]
)

PAC_FORTRAN_TEST_ISO_C_BIND( 
	[],
        [AC_MSG_ERROR([Sorry, cannot build PSBLAS without support for ISO_C_BINDING. 
 Please get a Fortran compiler that supports it, e.g. GNU Fortran 4.8.])]
)

PAC_FORTRAN_TEST_ISO_FORTRAN_ENV(
	[],
        [AC_MSG_ERROR([Sorry, cannot build PSBLAS without support for ISO_FORTRAN_ENV. 
 Please get a Fortran compiler that supports it, e.g. GNU Fortran 4.8.])]
)

PAC_FORTRAN_TEST_VOLATILE(
	[],
        [AC_MSG_ERROR([Sorry, cannot build PSBLAS without support for	VOLATILE])]	
)

PAC_FORTRAN_TEST_FINAL( 
	[],
        [AC_MSG_ERROR([Sorry, cannot build PSBLAS without support for FINAL])]	
)


#
# Optional features
#


PAC_FORTRAN_TEST_GENERICS( 
	[],
	[FDEFINES="$psblas_cv_define_prepend-DHAVE_BUGGY_GENERICS $FDEFINES"]
)

PAC_FORTRAN_TEST_FLUSH( 
	[FDEFINES="$psblas_cv_define_prepend-DHAVE_FLUSH_STMT $FDEFINES"],
)


PAC_FORTRAN_TEST_SAME_TYPE( 
	[FDEFINES="$psblas_cv_define_prepend-DHAVE_SAME_TYPE_AS $FDEFINES"],
)

PAC_FORTRAN_TEST_EXTENDS_TYPE( 
	[FDEFINES="$psblas_cv_define_prepend-DHAVE_EXTENDS_TYPE_OF $FDEFINES"],
)

PAC_FORTRAN_TEST_MOLD( 
	[FDEFINES="$psblas_cv_define_prepend-DHAVE_MOLD $FDEFINES"],
)

###############################################################################
#		BLAS  library presence checks
###############################################################################

# Note : The libmkl.a (Intel Math Kernel Library) library could be used, too.
# It is sufficient to specify it as -lmkl in the CLIBS or FLIBS or LIBS 
# and specify its path adjusting -L/path in CFLAGS.

# Right now it is a matter of user's taste when linking custom applications.
# But PSBLAS examples could take advantage of these libraries, too.
AC_LANG([Fortran])
PAC_BLAS([], [AC_MSG_ERROR([[Cannot find BLAS library, specify a path using --with-blas=DIR/LIB (for example --with-blas=/usr/path/lib/libcxml.a)]])])
PAC_LAPACK(
[FDEFINES="$psblas_cv_define_prepend-DHAVE_LAPACK $FDEFINES"],
)
AC_LANG([C])
###############################################################################
#		Auxiliary packages
###############################################################################


PAC_CHECK_MUMPS
#
# For the time being:
#    1. Disable MUMPS when IPK=8 (would need a check on whether MUMPS
#    		      	   	  was compiled with 8 bytes)
#    2. Enable even with LPK=8, internally it will check if
#    	the problem size fits into 4 bytes, very likely since we
#    	are mostly using MUMPS at coarse level.
#
if test "x$amg4psblas_cv_have_mumps" == "xyes" ; then
   if test "x$pac_cv_psblas_ipk" == "x8" ; then
      AC_MSG_NOTICE([PSBLAS defines PSB_IPK_ as $pac_cv_psblas_ipk. MUMPS interfacing disabled. ])
      MUMPS_FLAGS="";
      MUMPS_LIBS="";
      amg4psblas_cv_have_mumps=no;
   fi
fi   
if test "x$amg4psblas_cv_have_mumps" == "xyes" ; then
   if test "x$pac_cv_psblas_lpk" == "x8" ; then
      AC_MSG_NOTICE([PSBLAS defines PSB_LPK_ as $pac_cv_psblas_lpk. MUMPS interfacing will fail when called in global mode on very large matrices. ])
   fi	
   MUMPS_FLAGS="-DHave_MUMPS_ $MUMPS_INCLUDES"
   if test "x$pac_mumps_fmods_ok" == "xyes" ; then 
      FDEFINES="$amg_cv_define_prepend-DHAVE_MUMPS_ $amg_cv_define_prepend-DHAVE_MUMPS_MODULES_ $MUMPS_MODULES $FDEFINES"
   elif test "x$pac_mumps_fincs_ok" == "xyes" ; then
      FDEFINES="$amg_cv_define_prepend-DHAVE_MUMPS_ $amg_cv_define_prepend-DHAVE_MUMPS_INCLUDES_ $MUMPS_FINCLUDES $FDEFINES"
   else
      # This should not happen
      MUMPS_FLAGS=""
      MUMPS_LIBS=""
   fi
else
   MUMPS_FLAGS=""
   MUMPS_LIBS=""
fi

PAC_CHECK_UMFPACK
if test "x$amg4psblas_cv_have_umfpack" == "xyes" ; then 
   UMF_FLAGS="-DHave_UMF_ $UMF_INCLUDES"
   FDEFINES="$amg_cv_define_prepend-DHAVE_UMF_ $FDEFINES"
else
   UMF_FLAGS=""
fi

PAC_CHECK_SUPERLU
if test "x$amg4psblas_cv_have_superlu" == "xyes" ; then 
   SLU_FLAGS="-DHave_SLU_ -DSLU_VERSION_$pac_slu_version $SLU_INCLUDES"
   FDEFINES="$amg_cv_define_prepend-DHAVE_SLU_ $FDEFINES"
else
   SLU_FLAGS=""
fi

PAC_CHECK_SUPERLUDIST()

if test "x$amg4psblas_cv_have_superludist" == "xyes" ; then
   pac_sludist_version="$amg4psblas_cv_superludist_major";      
   if (($amg4psblas_cv_superludist_major==6)); then
      if (($amg4psblas_cv_superludist_minor>=3)); then
        pac_sludist_version="63";
      fi
   fi	
   SLUDIST_FLAGS=""	       
   SLUDIST_FLAGS="-DHave_SLUDist_ -DSLUD_VERSION_$pac_sludist_version $SLUDIST_INCLUDES"
   FDEFINES="$amg_cv_define_prepend-DHAVE_SLUDIST_ $FDEFINES"
else
   SLUDIST_FLAGS=""
fi

##############################################
FINCLUDES="$PSBLAS_INCLUDES"

if test "x$pac_cv_status_file" != "xNONE"; then
   PSBLAS_INSTALL_MAKEINC='include $(PSBLAS_INCDIR)/Make.inc.psblas';
else
   PSBLAS_INSTALL_MAKEINC="";
fi

if test "x$pac_cv_status_file" != "xNONE"; then
  COMPILERULES='';
else
COMPILERULES='
FLINK=$(MPFC)

# These should be portable rules, arent they?
.c.o:
	$(CC) $(CCOPT) $(CINCLUDES) $(CDEFINES) -c $< -o $@
.f90.o:
	$(FC) $(FCOPT) $(FINCLUDES) -c $< -o $@
.F90.o:
	$(FC) $(FCOPT) $(FINCLUDES) $(FDEFINES) -c $< -o $@'

COMPILERULES="
# The following ones are the variables used by the PSBLAS make scripts.

FC=$FC
CC=$CC
FCOPT=$FCOPT
CCOPT=$CCOPT
FMFLAG=$FMFLAG
FIFLAG=$FIFLAG
EXTRA_OPT=$EXTRA_OPT

BASEMODNAME=$BASEMODNAME
PRECMODNAME=$PRECMODNAME
METHDMODNAME=$METHDMODNAME
UTILMODNAME=$UTILMODNAME

LIBS=$LIBS
BLAS=$BLAS_LIBS

# These three should be always set!
MPFC=$MPFC
MPCC=$MPCC

AR=$AR
RANLIB=$RANLIB
PSBLDLIBS=$LAPACK $BLAS_LIBS $METIS_LIB $LIBS

$COMPILERULES"
fi 


###############################################################################
# Variable substitutions : the Make.inc.in will have these @VARIABLES@
# substituted.

AC_SUBST(PSBLAS_DIR)
AC_SUBST(PSBLAS_INCDIR)
AC_SUBST(PSBLAS_MODDIR)
AC_SUBST(PSBLAS_LIBDIR)
AC_SUBST(PSBLAS_INCLUDES)
AC_SUBST(PSBLAS_INSTALL_MAKEINC)
AC_SUBST(PSBLAS_LIBS)

AC_SUBST(INSTALL)
AC_SUBST(INSTALL_DATA)
AC_SUBST(INSTALL_DIR)
AC_SUBST(INSTALL_LIBDIR)
AC_SUBST(INSTALL_INCLUDEDIR)
AC_SUBST(INSTALL_MODULESDIR)
AC_SUBST(INSTALL_DOCSDIR)
AC_SUBST(INSTALL_SAMPLESDIR)

AC_SUBST(EXTRA_LIBS)
AC_SUBST(MUMPS_FLAGS)
AC_SUBST(MUMPS_LIBS)
AC_SUBST(SLU_FLAGS)
AC_SUBST(SLU_LIBS)
AC_SUBST(UMF_FLAGS)
AC_SUBST(UMF_LIBS)
AC_SUBST(SLUDIST_FLAGS)
AC_SUBST(SLUDIST_LIBS)
AC_SUBST(FDEFINES)
AC_SUBST(MODEXT)
AC_SUBST(COMPILERULES)
###############################################################################
# the following files will be created by Automake

AC_CONFIG_FILES([Make.inc])
AC_OUTPUT()

###############################################################################

dnl Please note that brackets around variable identifiers are absolutely needed for compatibility..
AC_MSG_NOTICE([
	${PACKAGE_NAME} ${amg4psblas_cv_version} has been configured as follows:

	PSBLAS library        : ${PSBLAS_DIR}
	MUMPS detected        : ${amg4psblas_cv_have_mumps}
	SuperLU detected      : ${amg4psblas_cv_have_superlu}
	SuperLU_Dist detected : ${amg4psblas_cv_have_superludist}
	UMFPack detected      : ${amg4psblas_cv_have_umfpack}

	If you are satisfied, run 'make' to build ${PACKAGE_NAME} and its documentation; otherwise
	type  ./configure --help=short for a complete list of configure options specific to ${PACKAGE_NAME}.
dnl	To install the program and its documentation, run 'make install' if you are root,
dnl	or run 'su -c "make install"' if you are not root.
])

###############################################################################

