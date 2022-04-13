dnl
dnl $Id$
dnl
dnl 20080206
dnl M4 macros for the PSBLAS library and useful for packages using PSBLAS.
dnl

dnl @synopsis PAC_CHECK_LIBS
dnl
dnl Tries to detect the presence of a specific function among various libraries, using AC_CHECK_LIB
dnl repeatedly on the specified libraries.
dnl 
dnl Example use:
dnl
dnl PAC_CHECK_LIBS([atlas blas],
dnl		[dgemm],
dnl		[have_dgemm=yes],
dnl		[have_dgemm=no])
dnl 
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl
dnl 20080211 modified slighty from original.
AC_DEFUN([PAC_CHECK_LIBS],
[
 pac_check_libs_ok=no
 [for pac_check_libs_f in $2 
 do ]
 [for pac_check_libs_l in $1 
 do ]
    if test x"$pac_check_libs_ok" == xno ; then
     AC_CHECK_LIB([$pac_check_libs_l],[$pac_check_libs_f], [pac_check_libs_ok=yes; pac_check_libs_LIBS="-l$pac_check_libs_l"],[],[$5])
    fi
  done
  done
 # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
 [ if test x"$pac_check_libs_ok" = xyes ; then
	$3
 else
        pac_check_libs_ok=no
        $4
 fi
 ]
])dnl 


dnl @synopsis PAC_FORTRAN_HAVE_MOVE_ALLOC( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program with move_alloc (a Fortran 2003 function).
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl
AC_DEFUN([PAC_FORTRAN_HAVE_MOVE_ALLOC],
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([support for Fortran MOVE_ALLOC intrinsic])
 AC_LANG_PUSH([Fortran])
 ac_ext='f90';
 AC_COMPILE_IFELSE([ program test_move_alloc
		       integer, allocatable :: a(:), b(:)
		       allocate(a(3))
		       call move_alloc(a, b)
		       print *, allocated(a), allocated(b)
		       print *, b
		     end program test_move_alloc],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AC_FD_CC
		     cat conftest.$ac_ext >&AC_FD_CC
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])



dnl @synopsis PAC_CHECK_HAVE_CRAYFTN( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will check if MPIFC is $FC.
dnl The check will proceed by compiling a small Fortran program
dnl containing the _CRAYFTN macro, which should be defined in the
dnl gfortran compiled programs.
dnl
dnl On pass, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_CHECK_HAVE_CRAYFTN,
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([for Cray Fortran])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='F90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
           program main
#ifdef _CRAYFTN 
              print *, "Cray FTN!"
#else
        this program will fail
#endif
           end],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AC_FD_CC
		     cat conftest.$ac_ext >&AC_FD_CC
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])



dnl @synopsis PAC_CHECK_HAVE_GFORTRAN( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will check if MPIFC is $FC.
dnl The check will proceed by compiling a small Fortran program
dnl containing the __GNUC__ macro, which should be defined in the
dnl gfortran compiled programs.
dnl
dnl On pass, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl
AC_DEFUN(PAC_CHECK_HAVE_GFORTRAN,
[AC_MSG_CHECKING([for GNU Fortran])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='F90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
           program main
#ifdef __GNUC__ 
              print *, "GCC!"
#else
        this program will fail
#endif
           end],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AC_FD_CC
		     cat conftest.$ac_ext >&AC_FD_CC
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])


dnl @synopsis PAC_HAVE_MODERN_GFORTRAN( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will check if the GNU fortran version is suitable for PSBLAS.
dnl If yes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl
dnl Note : Will use MPIFC; if unset, will use '$FC'.
dnl 
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl
AC_DEFUN(PAC_HAVE_MODERN_GFORTRAN,
 [AC_MSG_CHECKING([for recent GNU Fortran])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='F90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
           program main
#if ( __GNUC__ >= 4 && __GNUC_MINOR__ >= 8 ) || ( __GNUC__ > 4 )
              print *, "ok"
#else
        this program will fail
#endif
           end],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])
		     AC_MSG_NOTICE([Sorry, we require GNU Fortran version 4.8.4 or later.])
		     echo "configure: failed program was:" >&AC_FD_CC
		     cat conftest.$ac_ext >&AC_FD_CC
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])


dnl @synopsis PAC_ARG_WITH_FLAGS(lcase_name, UCASE_NAME)
dnl
dnl Test for --with-lcase_name="compiler/loader flags".  if defined, prepend 
dnl flags to standard UCASE_NAME definition.
dnl
dnl Use this macro to facilitate additional special flags that should be
dnl passed on to the preprocessor/compilers/loader.
dnl
dnl NOTE : Renamed after TAC_ARG_WITH_FLAGS as in the Trilinos-8.0.4 package.
dnl 
dnl NOTE : This macro works in a way the user should invoke
dnl         --with-flags=...
dnl	   only once, otherwise the first one will take effect.
dnl
dnl Example use:
dnl 
dnl PAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
dnl 
dnl tests for --with-cxxflags and pre-pends to CXXFLAGS
dnl 
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl @notes  Michele Martone <michele.martone@uniroma2.it>
dnl
AC_DEFUN([PAC_ARG_WITH_FLAGS],
[
AC_MSG_CHECKING([whether additional [$2] flags should be added (should be invoked only once)])
dnl AC_MSG_CHECKING([whether additional [$2] flags should be added])
AC_ARG_WITH($1,
AC_HELP_STRING([--with-$1], 
[additional [$2] flags to be added: will prepend to [$2]]),
[
$2="${withval} ${$2}"
AC_MSG_RESULT([$2 = ${$2}])
],
AC_MSG_RESULT(no)
)
])


dnl @synopsis PAC_ARG_WITH_LIBS
dnl
dnl Test for --with-libs="name(s)".
dnl 
dnl Prepends the specified name(s) to the list of libraries to link 
dnl with.  
dnl
dnl note: Renamed after PAC_ARG_WITH_LIBS as in the Trilinos package.
dnl
dnl Example use:
dnl
dnl PAC_ARG_WITH_LIBS
dnl 
dnl tests for --with-libs and pre-pends to LIBS
dnl
dnl @author Jim Willenbring <jmwille@sandia.gov>
dnl
AC_DEFUN([PAC_ARG_WITH_LIBS],
[
AC_MSG_CHECKING([whether additional libraries are needed])
AC_ARG_WITH(libs,
AC_HELP_STRING([--with-libs], 
[List additional link flags  here.  For example, --with-libs=-lspecial_system_lib
or --with-libs=-L/path/to/libs]),
[
LIBS="${withval} ${LIBS}"
AC_MSG_RESULT([LIBS = ${LIBS}])
],
AC_MSG_RESULT(no)
)
]
)


dnl @synopsis PAC_ARG_WITH_EXTRA_LIBS
dnl
dnl Test for --with-extra-libs="name(s)".
dnl 
dnl Appends the specified name(s) to the list of libraries to link 
dnl with.  
dnl
dnl note: Renamed after PAC_ARG_WITH_EXTRA_LIBS as in the Trilinos package.
dnl
dnl Example use:
dnl
dnl PAC_ARG_WITH_EXTRA_LIBS
dnl 
dnl tests for --with-extra-libs and pre-pends to LIBS
dnl
dnl @author Jim Willenbring <jmwille@sandia.gov>
dnl
AC_DEFUN([PAC_ARG_WITH_EXTRA_LIBS],
[
AC_MSG_CHECKING([whether additional libraries are needed])
AC_ARG_WITH(extra-libs,
AC_HELP_STRING([--with-extra-libs], 
[List additional link flags  here.  For example, --with-extra-libs=-lspecial_system_lib
or --with-extra-libs=-L/path/to/libs]),
[
EXTRA_LIBS="${withval}"
AC_MSG_RESULT([EXTRA_LIBS = ${EXTRA_LIBS}])
],
AC_MSG_RESULT(no)
)
]
)

dnl @synopsis PAC_ARG_WITH_PSBLAS
dnl
dnl Test for --with-psblas="pathname".
dnl 
dnl Defines the path to PSBLAS build dir.
dnl
dnl note: Renamed after PAC_ARG_WITH_LIBS as in the Trilinos package.
dnl
dnl Example use:
dnl
dnl PAC_ARG_WIT_PSBLAS
dnl 
dnl tests for --with-psblas and pre-pends to PSBLAS_PATH
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN([PAC_ARG_WITH_PSBLAS],
[
AC_ARG_WITH(psblas,
AC_HELP_STRING([--with-psblas=DIR], [The install directory for PSBLAS, for example,
 --with-psblas=/opt/packages/psblas-3.5]),
[pac_cv_psblas_dir=$withval],
[pac_cv_psblas_dir=''])
AC_ARG_WITH(psblas-incdir, AC_HELP_STRING([--with-psblas-incdir=DIR], [Specify the directory for PSBLAS C includes.]),
        [pac_cv_psblas_incdir=$withval],
        [pac_cv_psblas_incdir=''])
AC_ARG_WITH(psblas-moddir, AC_HELP_STRING([--with-psblas-moddir=DIR], [Specify the directory for PSBLAS Fortran modules.]),
        [pac_cv_psblas_moddir=$withval],
        [pac_cv_psblas_moddir=''])
AC_ARG_WITH(psblas-libdir, AC_HELP_STRING([--with-psblas-libdir=DIR], [Specify the directory for PSBLAS library.]),
        [pac_cv_psblas_libdir=$withval],
        [pac_cv_psblas_libdir=''])
if test x"$pac_cv_psblas_incdir" == "x" ; then
   pac_cv_psblas_incdir="$pac_cv_psblas_dir/include";
fi
if test x"$pac_cv_psblas_moddir" == "x" ; then
   pac_cv_psblas_moddir="$pac_cv_psblas_dir/modules";
fi
if test x"$pac_cv_psblas_libdir" == "x" ; then
   pac_cv_psblas_libdir="$pac_cv_psblas_dir/lib";
fi

])

dnl @synopsis PAC_FORTRAN_HAVE_PSBLAS( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program using the PSBLAS library
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl
AC_DEFUN(PAC_FORTRAN_HAVE_PSBLAS,
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([for working installation of PSBLAS])
 AC_LANG_PUSH([Fortran])
 ac_objext='o'
 ac_ext='f90'
 ac_fc="${MPIFC-$FC}";
 save_FCFLAGS="$FCFLAGS";
dnl FCFLAGS=" $FMFLAG$PSBLAS_DIR/include $save_FCFLAGS"
dnl save_FCFLAGS="$FCFLAGS";
 save_LDFLAGS="$LDFLAGS";
if test "x$pac_cv_psblas_incdir" != "x"; then 
dnl  AC_MSG_NOTICE([psblas include dir $pac_cv_psblas_incdir])
 PSBLAS_INCLUDES="$FMFLAG$pac_cv_psblas_incdir"
elif test "x$pac_cv_psblas_dir" != "x"; then 
dnl AC_MSG_NOTICE([psblas dir $pac_cv_psblas_dir])
 PSBLAS_INCLUDES="$FMFLAG$pac_cv_psblas_dir/include"
fi
if test "x$pac_cv_psblas_moddir" != "x"; then 
dnl  AC_MSG_NOTICE([psblas modules dir $pac_cv_psblas_moddir])
 PSBLAS_INCLUDES="$FMFLAG$pac_cv_psblas_moddir $PSBLAS_INCLUDES"
elif test "x$pac_cv_psblas_dir" != "x"; then 
dnl AC_MSG_NOTICE([psblas dir $pac_cv_psblas_dir])
 PSBLAS_INCLUDES="$FMFLAG$pac_cv_psblas_dir/modules $PSBLAS_INCLUDES"
fi
 FCFLAGS=" $PSBLAS_INCLUDES $save_FCFLAGS"
AC_COMPILE_IFELSE([
		    program test
		    use psb_base_mod
		    end program test],
		   [ifelse([$1], , :, [ $1   ])],
		   [  ifelse([$2], , , [ $2 ])])
AC_LANG_POP([Fortran])
rm -f conftest*])


dnl @synopsis PAC_FORTRAN_PSBLAS_VERSION( )
dnl
dnl Will try to compile, link and run  a program using the PSBLAS library. \
dnl  Checks for version major,  minor and patchlevel
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_FORTRAN_PSBLAS_VERSION,
[AC_MSG_CHECKING([for version of PSBLAS])
AC_LANG_PUSH([Fortran])	  
ac_exeext=''
ac_objext='o'
ac_ext='f90'
save_FCFLAGS=$FCFLAGS;
FCFLAGS=" $PSBLAS_INCLUDES $save_FCFLAGS"
save_LDFLAGS=$LDFLAGS;
## if test "x$pac_cv_psblas_libdir" != "x"; then 
## dnl AC_MSG_NOTICE([psblas lib dir $pac_cv_psblas_libdir])
##  PSBLAS_LIBS="-L$pac_cv_psblas_libdir"
## elif test "x$pac_cv_psblas_dir" != "x"; then 
## dnl AC_MSG_NOTICE([psblas dir $pac_cv_psblas_dir])
##  PSBLAS_LIBS="-L$pac_cv_psblas_dir/lib"
## fi
PSBLAS_LIBS="-lpsb_krylov -lpsb_prec -lpsb_util -lpsb_base -L$PSBLAS_LIBDIR"
LDFLAGS=" $PSBLAS_LIBS $save_LDFLAGS"

dnl ac_compile='${MPIFC-$FC} -c -o conftest${ac_objext} $FMFLAG$PSBLAS_DIR/include $FMFLAG$PSBLAS_DIR/lib conftest.$ac_ext  1>&5'
dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FCFLAGS $LDFLAGS conftest.$ac_ext $FMFLAG$PSBLAS_DIR/include -L$PSBLAS_DIR/lib -lpsb_base $LIBS 1>&5'
ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FCFLAGS conftest.$ac_ext $LDFLAGS $LIBS  1>&5'
dnl Warning : square brackets are EVIL!

AC_LINK_IFELSE(
    [AC_LANG_SOURCE([[
		program test
		use psb_base_mod, only : psb_version_major_
		print *,psb_version_major_
		end program test]])],
	       [pac_cv_psblas_major=`./conftest${ac_exeext} | sed 's/^ *//'`],
	       [pac_cv_psblas_major="unknown"])
  
AC_LINK_IFELSE(
    [AC_LANG_SOURCE([[
		program test
		use psb_base_mod, only : psb_version_minor_
		print *,psb_version_minor_
		end program test]])],
	       [pac_cv_psblas_minor=`./conftest${ac_exeext} | sed 's/^ *//'`],
	       [pac_cv_psblas_minor="unknown"])
  
AC_LINK_IFELSE(
    [AC_LANG_SOURCE([[
		program test
		use psb_base_mod, only : psb_patchlevel_
		print *,psb_patchlevel_
		end program test]])],
	       [pac_cv_psblas_patchlevel=`./conftest${ac_exeext} | sed 's/^ *//'`],
	       [pac_cv_psblas_patchlevel="unknown"])
LDFLAGS="$save_LDFLAGS";
FCFLAGS="$save_FCFLAGS";

AC_MSG_RESULT([Done])
AC_LANG_POP([Fortran])])


dnl @synopsis PAC_FORTRAN_PSBLAS_INTEGER_SIZES( )
dnl
dnl Will try to compile, link and run  a program using the PSBLAS library. \
dnl  Checks for size of integers
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_FORTRAN_PSBLAS_INTEGER_SIZES,
[AC_MSG_CHECKING([for PSBLAS integer sizes])
AC_LANG_PUSH([Fortran])	  
ac_exeext=''
ac_objext='o'
ac_ext='f90'
save_FCFLAGS=$FCFLAGS;
FCFLAGS=" $PSBLAS_INCLUDES $save_FCFLAGS"
save_LDFLAGS=$LDFLAGS;
if test "x$pac_cv_psblas_libdir" != "x"; then 
dnl AC_MSG_NOTICE([psblas lib dir $pac_cv_psblas_libdir])
 PSBLAS_LIBS="-L$pac_cv_psblas_libdir"
elif test "x$pac_cv_psblas_dir" != "x"; then 
dnl AC_MSG_NOTICE([psblas dir $pac_cv_psblas_dir])
 PSBLAS_LIBS="-L$pac_cv_psblas_dir/lib"
fi
if test "x$pac_cv_psblas_moddir" != "x"; then 
dnl  AC_MSG_NOTICE([psblas modules dir $pac_cv_psblas_moddir])
 PSBLAS_INCLUDES="$FMFLAG$pac_cv_psblas_moddir $PSBLAS_INCLUDES"
elif test "x$pac_cv_psblas_dir" != "x"; then 
dnl AC_MSG_NOTICE([psblas dir $pac_cv_psblas_dir])
 PSBLAS_INCLUDES="$FMFLAG$pac_cv_psblas_dir/modules $PSBLAS_INCLUDES"
fi
FCFLAGS=" $PSBLAS_INCLUDES $save_FCFLAGS"
PSBLAS_LIBS="-lpsb_krylov -lpsb_prec -lpsb_util -lpsb_base $PSBLAS_LIBS"
LDFLAGS=" $PSBLAS_LIBS $save_LDFLAGS"

dnl ac_compile='${MPIFC-$FC} -c -o conftest${ac_objext} $FMFLAG$PSBLAS_DIR/include $FMFLAG$PSBLAS_DIR/lib conftest.$ac_ext  1>&5'
dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FCFLAGS $LDFLAGS conftest.$ac_ext $FMFLAG$PSBLAS_DIR/include -L$PSBLAS_DIR/lib -lpsb_base $LIBS 1>&5'
dnl Warning : square brackets are EVIL!

AC_LINK_IFELSE([
		program test
		use iso_fortran_env
		use psb_base_mod
		if (psb_ipk_ == int32) then
		  print *, '4'
		elseif (psb_ipk_ == int64) then
		  print *, '8'
		else
		  print *, '-1'
		end if
		end program test],
	       [pac_cv_psblas_ipk=`./conftest${ac_exeext} | sed 's/^ *//'`],
	       [pac_cv_psblas_ipk="unknown"])
  
AC_LINK_IFELSE([
		program test
		use iso_fortran_env
		use psb_base_mod
		if (psb_lpk_ == int32) then
		  print *, '4'
		elseif (psb_lpk_ == int64) then
		  print *, '8'
		else
		  print *, '-1'
		end if
		end program test],
	       [pac_cv_psblas_lpk=`./conftest${ac_exeext} | sed 's/^ *//'`],
	       [pac_cv_psblas_lpk="unknown"])
LDFLAGS="$save_LDFLAGS";
FCFLAGS="$save_FCFLAGS";

AC_MSG_RESULT([Done])
AC_LANG_POP([Fortran])])



dnl @synopsis PAC_CHECK_UMFPACK
dnl
dnl Will try to find the UMFPACK library and headers.
dnl
dnl Will use $CC
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_CHECK_UMFPACK,
[AC_ARG_WITH(umfpack, AC_HELP_STRING([--with-umfpack=LIBNAME], [Specify the library name for UMFPACK and its support libraries. 
Default: "-lumfpack -lamd"]),
        [amg4psblas_cv_umfpack=$withval],
        [amg4psblas_cv_umfpack='-lumfpack -lamd'])
AC_ARG_WITH(umfpackdir, AC_HELP_STRING([--with-umfpackdir=DIR], [Specify the directory for UMFPACK library and includes.]),
        [amg4psblas_cv_umfpackdir=$withval],
        [amg4psblas_cv_umfpackdir=''])
AC_ARG_WITH(umfpackincdir, AC_HELP_STRING([--with-umfpackincdir=DIR], [Specify the directory for UMFPACK includes.]),
        [amg4psblas_cv_umfpackincdir=$withval],
        [amg4psblas_cv_umfpackincdir=''])
AC_ARG_WITH(umfpacklibdir, AC_HELP_STRING([--with-umfpacklibdir=DIR], [Specify the directory for UMFPACK library.]),
        [amg4psblas_cv_umfpacklibdir=$withval],
        [amg4psblas_cv_umfpacklibdir=''])

AC_LANG_PUSH([C])
save_LIBS="$LIBS"
save_CPPFLAGS="$CPPFLAGS"
if test "x$amg4psblas_cv_umfpackincdir" != "x"; then 
 AC_MSG_NOTICE([umfp include dir $amg4psblas_cv_umfpackincdir])
 UMF_INCLUDES="-I$amg4psblas_cv_umfpackincdir"
elif test "x$amg4psblas_cv_umfpackdir" != "x"; then 
 AC_MSG_NOTICE([umfp dir $amg4psblas_cv_umfpackdir])
 UMF_INCLUDES="-I$amg4psblas_cv_umfpackdir" 
fi
CPPFLAGS="$UMF_INCLUDES $CPPFLAGS"
AC_CHECK_HEADER([umfpack.h],
 [pac_umf_header_ok=yes],
 [pac_umf_header_ok=no; UMF_INCLUDES=""])
if test "x$amg4psblas_cv_umfpacklibdir" != "x"; then 
   LIBS="-L$amg4psblas_cv_umfpacklibdir $LIBS $EXTRA_LIBS"
   UMF_LIBDIR="-L$amg4psblas_cv_umfpacklibdir"
elif test "x$amg4psblas_cv_umfpackdir" != "x"; then 
   LIBS="-L$amg4psblas_cv_umfpackdir $LIBS $EXTRA_LIBS"
   UMF_LIBDIR="-L$amg4psblas_cv_umfpackdir"
fi
if test "x$pac_umf_header_ok" == "xno" ; then
  unset ac_cv_header_umfpack_h
  UMF_INCLUDES="-I$amg4psblas_cv_umfpackdir" 
  CPPFLAGS="$UMF_INCLUDES $SAVE_CPPFLAGS"
  AC_CHECK_HEADER([umfpack.h],
  [pac_umf_header_ok=yes],
  [pac_umf_header_ok=no; UMF_INCLUDES=""])
fi
if test "x$pac_umf_header_ok" == "xno" ; then
dnl Maybe Include or include subdirs? 
  unset ac_cv_header_umfpack_h
  UMF_INCLUDES="-I$amg4psblas_cv_umfpackdir/include -I$amg4psblas_cv_umfpackdir/Include "
  CPPFLAGS="$UMF_INCLUDES $SAVE_CPPFLAGS"

  AC_MSG_CHECKING([for umfpack_di_symbolic in $UMF_INCLUDES])
  AC_CHECK_HEADER([umfpack.h],
    [pac_umf_header_ok=yes],
    [pac_umf_header_ok=no; UMF_INCLUDES=""])
fi
if test "x$pac_umf_header_ok" == "xno" ; then
dnl Maybe new structure with UMFPACK UFconfig AMD? 
   unset ac_cv_header_umfpack_h
   UMF_INCLUDES="-I$amg4psblas_cv_umfpackdir/UFconfig -I$amg4psblas_cv_umfpackdir/UMFPACK/Include -I$amg4psblas_cv_umfpackdir/AMD/Include"
   CPPFLAGS="$UMF_INCLUDES $SAVE_CPPFLAGS"
   AC_CHECK_HEADER([umfpack.h],
     [pac_umf_header_ok=yes],
     [pac_umf_header_ok=no; UMF_INCLUDES=""])
fi


if test "x$pac_umf_header_ok" == "xyes" ; then 
      UMF_LIBS="$amg4psblas_cv_umfpack $UMF_LIBDIR"
      LIBS="$UMF_LIBS -lm $LIBS $EXTRA_LIBS";
      AC_MSG_CHECKING([for umfpack_di_symbolic in $UMF_LIBS])
      AC_TRY_LINK_FUNC(umfpack_di_symbolic, 
       [amg4psblas_cv_have_umfpack=yes;pac_umf_lib_ok=yes; ],
       [amg4psblas_cv_have_umfpack=no;pac_umf_lib_ok=no; UMF_LIBS=""])
      AC_MSG_RESULT($pac_umf_lib_ok)
     if test "x$pac_umf_lib_ok" == "xno" ; then 
        dnl Maybe Lib or lib? 
        UMF_LIBDIR="-L$amg4psblas_cv_umfpackdir/Lib -L$amg4psblas_cv_umfpackdir/lib"
        UMF_LIBS="$amg4psblas_cv_umfpack $UMF_LIBDIR -lm $SAVE_LIBS  $EXTRA_LIBS"
        LIBS="$UMF_LIBS"
        
      AC_MSG_CHECKING([for umfpack_di_symbolic in $UMF_LIBS])
      AC_TRY_LINK_FUNC(umfpack_di_symbolic, 
       [amg4psblas_cv_have_umfpack=yes;pac_umf_lib_ok=yes; ],
       [amg4psblas_cv_have_umfpack=no;pac_umf_lib_ok=no; UMF_LIBS=""])
      AC_MSG_RESULT($pac_umf_lib_ok)
     fi
     if test "x$pac_umf_lib_ok" == "xno" ; then 
        dnl Maybe UMFPACK/Lib? 
        UMF_LIBDIR="-L$amg4psblas_cv_umfpackdir/AMD/Lib -L$amg4psblas_cv_umfpackdir/UMFPACK/Lib"
        UMF_LIBS="$amg4psblas_cv_umfpack $UMF_LIBDIR -lm $SAVE_LIBS $EXTRA_LIBS"
             LIBS="$UMF_LIBS"
      AC_MSG_CHECKING([for umfpack_di_symbolic in $UMF_LIBS])
      AC_TRY_LINK_FUNC(umfpack_di_symbolic, 
       [amg4psblas_cv_have_umfpack=yes;pac_umf_lib_ok=yes; ],
       [amg4psblas_cv_have_umfpack=no;pac_umf_lib_ok=no; UMF_LIBS=""])
      AC_MSG_RESULT($pac_umf_lib_ok)
     fi
fi
LIBS="$SAVE_LIBS";
CPPFLAGS="$SAVE_CPPFLAGS";
AC_LANG_POP([C])
])dnl 

dnl @synopsis PAC_CHECK_SUPERLU
dnl
dnl Will try to find the SUPERLU library and headers.
dnl
dnl Will use $CC
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_CHECK_SUPERLU,
[AC_ARG_WITH(superlu, AC_HELP_STRING([--with-superlu=LIBNAME], [Specify the library name for SUPERLU library.
Default: "-lsuperlu"]),
        [amg4psblas_cv_superlu=$withval],
        [amg4psblas_cv_superlu='-lsuperlu'])
AC_ARG_WITH(superludir, AC_HELP_STRING([--with-superludir=DIR], [Specify the directory for SUPERLU library and includes.]),
        [amg4psblas_cv_superludir=$withval],
        [amg4psblas_cv_superludir=''])
AC_ARG_WITH(superluincdir, AC_HELP_STRING([--with-superluincdir=DIR], [Specify the directory for SUPERLU includes.]),
        [amg4psblas_cv_superluincdir=$withval],
        [amg4psblas_cv_superluincdir=''])
AC_ARG_WITH(superlulibdir, AC_HELP_STRING([--with-superlulibdir=DIR], [Specify the directory for SUPERLU library.]),
        [amg4psblas_cv_superlulibdir=$withval],
        [amg4psblas_cv_superlulibdir=''])
AC_LANG_PUSH([C])
save_LIBS="$LIBS"
save_CPPFLAGS="$CPPFLAGS"
if test "x$amg4psblas_cv_superluincdir" != "x"; then 
 AC_MSG_NOTICE([slu include dir $amg4psblas_cv_superluincdir])
 SLU_INCLUDES="-I$amg4psblas_cv_superluincdir"
elif test "x$amg4psblas_cv_superludir" != "x"; then 
 AC_MSG_NOTICE([slu dir $amg4psblas_cv_superludir])
 SLU_INCLUDES="-I$amg4psblas_cv_superludir"
fi
if test "x$amg4psblas_cv_superlulibdir" != "x"; then 
   SLU_LIBS="-L$amg4psblas_cv_superlulibdir"
elif test "x$amg4psblas_cv_superludir" != "x"; then 
   SLU_LIBS="-L$amg4psblas_cv_superludir"
fi

LIBS="$SLU_LIBS $LIBS"
CPPFLAGS="$SLU_INCLUDES $save_CPPFLAGS"
AC_CHECK_HEADERS([slu_ddefs.h],
		[pac_slu_header_ok=yes],
		[pac_slu_header_ok=no; SLU_INCLUDES=""])
if test "x$pac_slu_header_ok" == "xno" ; then 
dnl Maybe include subdirs? 
  unset ac_cv_header_slu_ddefs_h
  SLU_INCLUDES="-I$amg4psblas_cv_superludir/include "
  CPPFLAGS="$SLU_INCLUDES $save_CPPFLAGS"

  AC_CHECK_HEADERS([slu_ddefs.h],
		 [pac_slu_header_ok=yes],
		 [pac_slu_header_ok=no; SLU_INCLUDES=""])
fi
if test "x$pac_slu_header_ok" == "xno" ; then 
dnl Maybe Include subdirs? 
  unset ac_cv_header_slu_ddefs_h
  SLU_INCLUDES="-I$amg4psblas_cv_superludir/Include "
  CPPFLAGS="$SLU_INCLUDES $save_CPPFLAGS"

  AC_CHECK_HEADERS([slu_ddefs.h],
		 [pac_slu_header_ok=yes],
		 [pac_slu_header_ok=no; SLU_INCLUDES=""])
fi

if test "x$pac_slu_header_ok" == "xyes" ; then 
 SLU_LIBS="$amg4psblas_cv_superlu $SLU_LIBS"
 LIBS="$SLU_LIBS -lm $save_LIBS";
 AC_MSG_CHECKING([for superlu_malloc in $SLU_LIBS])
 AC_TRY_LINK_FUNC(superlu_malloc, 
		  [amg4psblas_cv_have_superlu=yes;pac_slu_lib_ok=yes;],
		  [amg4psblas_cv_have_superlu=no;pac_slu_lib_ok=no; SLU_LIBS=""; ])
 if test "x$pac_slu_lib_ok" == "xno" ; then 
    dnl Maybe lib?
    SLU_LIBS="$amg4psblas_cv_superlu -L$amg4psblas_cv_superludir/lib";
    LIBS="$SLU_LIBS -lm $save_LIBS";
    AC_TRY_LINK_FUNC(superlu_malloc, 
		     [amg4psblas_cv_have_superlu=yes;pac_slu_lib_ok=yes;],
		     [amg4psblas_cv_have_superlu=no;pac_slu_lib_ok=no; SLU_LIBS=""; ])
 fi
 if test "x$pac_slu_lib_ok" == "xno" ; then 
    dnl Maybe lib64?
    SLU_LIBS="$amg4psblas_cv_superlu -L$amg4psblas_cv_superludir/lib64";
    LIBS="$SLU_LIBS -lm $save_LIBS";
    AC_TRY_LINK_FUNC(superlu_malloc, 
		     [amg4psblas_cv_have_superlu=yes;pac_slu_lib_ok=yes;],
		     [amg4psblas_cv_have_superlu=no;pac_slu_lib_ok=no; SLU_LIBS=""; ])
 fi
 AC_MSG_RESULT($pac_slu_lib_ok)
fi
if test "x$pac_slu_header_ok" == "xyes" ; then 
   AC_MSG_CHECKING([for superlu version 5])
   AC_LANG_PUSH([C])
   AC_COMPILE_IFELSE(
       [AC_LANG_SOURCE([[#include "slu_ddefs.h"
			 int testdslu()
			 { SuperMatrix AC, *L, *U;
			   int *perm_r, *perm_c,  *etree,  panel_size, permc_spec, relax, info;
			   superlu_options_t options;   SuperLUStat_t stat;
			   GlobalLU_t Glu;   
			   dgstrf(&options, &AC, relax, panel_size, etree,
				  NULL, 0, perm_c, perm_r, L, U, &Glu, &stat, &info);               
			   
			 }]])],
       [ AC_MSG_RESULT([yes]);      pac_slu_version="5";],
       [ AC_MSG_RESULT([no]);      pac_slu_version="4";])
   AC_LANG_POP([C])
fi   

LIBS="$save_LIBS";
CPPFLAGS="$save_CPPFLAGS";
AC_LANG_POP([C])
])dnl 

dnl @synopsis PAC_CHECK_SUPERLU_Dist
dnl
dnl Will try to find the SUPERLU_Dist library and headers.
dnl
dnl Will use $CC
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_CHECK_SUPERLUDIST,
[AC_ARG_WITH(superludist, AC_HELP_STRING([--with-superludist=LIBNAME], [Specify the libname for SUPERLUDIST library. Requires you also specify SuperLU. Default: "-lsuperlu_dist"]),
        [amg4psblas_cv_superludist=$withval],
        [amg4psblas_cv_superludist='-lsuperlu_dist'])
AC_ARG_WITH(superludistdir, AC_HELP_STRING([--with-superludistdir=DIR], [Specify the directory for SUPERLUDIST library and includes.]),
        [amg4psblas_cv_superludistdir=$withval],
        [amg4psblas_cv_superludistdir=''])

AC_ARG_WITH(superludistincdir, AC_HELP_STRING([--with-superludistincdir=DIR], [Specify the directory for SUPERLUDIST includes.]),
        [amg4psblas_cv_superludistincdir=$withval],
        [amg4psblas_cv_superludistincdir=''])

AC_ARG_WITH(superludistlibdir, AC_HELP_STRING([--with-superludistlibdir=DIR], [Specify the directory for SUPERLUDIST library.]),
        [amg4psblas_cv_superludistlibdir=$withval],
        [amg4psblas_cv_superludistlibdir=''])

AC_LANG_PUSH([C])
save_LIBS="$LIBS"
save_CPPFLAGS="$CPPFLAGS"
save_CC="$CC"
CC=${MPICC}
CPP="${CC} -E"
if test "x$amg4psblas_cv_superludistincdir" != "x"; then 
 AC_MSG_NOTICE([sludist dir $amg4psblas_cv_superludistincdir]) 
 SLUDIST_INCLUDES="-I$amg4psblas_cv_superludistincdir"
elif test "x$amg4psblas_cv_superludistdir" != "x"; then 
 AC_MSG_NOTICE([sludist dir $amg4psblas_cv_superludistdir]) 
 SLUDIST_INCLUDES="-I$amg4psblas_cv_superludistdir"
fi
if test "x$amg4psblas_cv_superludistlibdir" != "x"; then 
   SLUDIST_LIBS="-L$amg4psblas_cv_superludistlibdir"
elif test "x$amg4psblas_cv_superludistdir" != "x"; then 
   SLUDIST_LIBS="-L$amg4psblas_cv_superludistdir"
fi

LIBS="$SLUDIST_LIBS $save_LIBS"
CPPFLAGS="$SLUDIST_INCLUDES $save_CPPFLAGS"

AC_CHECK_HEADERS([superlu_ddefs.h],
 [pac_sludist_header_ok=yes],
 [pac_sludist_header_ok=no; SLUDIST_INCLUDES=""])
if test "x$pac_sludist_header_ok" == "xno" ; then 
dnl Maybe include subdirs? 
  unset ac_cv_header_superlu_ddefs_h
  SLUDIST_INCLUDES="-I$amg4psblas_cv_superludistdir/include"
  CPPFLAGS="$SLUDIST_INCLUDES $save_CPPFLAGS"

 AC_CHECK_HEADERS([superlu_ddefs.h],
		 [pac_sludist_header_ok=yes],
		 [pac_sludist_header_ok=no; SLUDIST_INCLUDES=""; SLUDIST_LIBS=""; ])
fi
if test "x$pac_sludist_header_ok" == "xno" ; then 
dnl Maybe Include subdirs? 
  unset ac_cv_header_superlu_ddefs_h
  SLUDIST_INCLUDES="-I$amg4psblas_cv_superludistdir/Include"
  CPPFLAGS="$SLUDIST_INCLUDES $save_CPPFLAGS"

 AC_CHECK_HEADERS([superlu_ddefs.h],
		 [pac_sludist_header_ok=yes],
		 [pac_sludist_header_ok=no; SLUDIST_INCLUDES=""; SLUDIST_LIBS=""; ])
fi

if test "x$pac_sludist_header_ok" == "xyes" ; then 
      SLUDIST_LIBS="$amg4psblas_cv_superludist $SLUDIST_LIBS"
      LIBS="$SLUDIST_LIBS -lm $save_LIBS";
      AC_MSG_CHECKING([for superlu_malloc_dist in $SLUDIST_LIBS])
      AC_TRY_LINK_FUNC(superlu_malloc_dist, 
       [amg4psblas_cv_have_superludist=yes;pac_sludist_lib_ok=yes;],
       [amg4psblas_cv_have_superludist=no;pac_sludist_lib_ok=no; 
          SLUDIST_LIBS=""; ])
  if test "x$pac_sludist_lib_ok" == "xno" ; then 
     dnl Maybe lib?
     SLUDIST_LIBS="$amg4psblas_cv_superludist  -L$amg4psblas_cv_superludistdir/lib";
     LIBS="$SLUDIST_LIBS -lm $save_LIBS";
     AC_MSG_CHECKING([for superlu_malloc_dist in $SLUDIST_LIBS])
     AC_TRY_LINK_FUNC(superlu_malloc_dist, 
		     [amg4psblas_cv_have_superludist=yes;pac_sludist_lib_ok=yes;],
		     [amg4psblas_cv_have_superludist=no;pac_sludist_lib_ok=no; 
		      SLUDIST_LIBS="";])
  fi
  if test "x$pac_sludist_lib_ok" == "xno" ; then 
     dnl Maybe lib64?
     SLUDIST_LIBS="$amg4psblas_cv_superludist  -L$amg4psblas_cv_superludistdir/lib64";
     LIBS="$SLUDIST_LIBS -lm $save_LIBS";
     AC_MSG_CHECKING([for superlu_malloc_dist in $SLUDIST_LIBS])
     AC_TRY_LINK_FUNC(superlu_malloc_dist, 
		     [amg4psblas_cv_have_superludist=yes;pac_sludist_lib_ok=yes;],
		     [amg4psblas_cv_have_superludist=no;pac_sludist_lib_ok=no; 
		      SLUDIST_LIBS="";SLUDIST_INCLUDES=""])
  fi
  AC_MSG_RESULT([$pac_sludist_lib_ok $SLUDIST_LIBS])
fi

if test "x$pac_sludist_lib_ok" == "xyes" ; then
    
  AC_LANG_PUSH([C])
  ac_cc=${MPICC-$CC}
  ac_exeext="";
  CPPFLAGS="$SLUDIST_INCLUDES $save_CPPFLAGS"
  LIBS="$SLUDIST_LIBS -lm $save_LIBS";
  AC_LINK_IFELSE(
      [AC_LANG_SOURCE([[#include <superlu_defs.h>
		  #include <stdio.h>
		  
		  void main()
		  { int i=SUPERLU_DIST_MAJOR_VERSION;      
		    printf("%d\n",i);
		    } ]])],
	       [amg4psblas_cv_superludist_major=`./conftest${ac_exeext} | sed 's/^ *//'`],
	       [amg4psblas_cv_superludist_major="unknown"])
  AC_LINK_IFELSE(
      [AC_LANG_SOURCE([[#include <superlu_defs.h>
		  #include <stdio.h>
		  
		  void main()
		  { int i=SUPERLU_DIST_MINOR_VERSION;      
		    printf("%d\n",i);
		    }]])],
	       [amg4psblas_cv_superludist_minor=`./conftest${ac_exeext} | sed 's/^ *//'`],
	       [amg4psblas_cv_superludist_minor="unknown"])
  AC_LANG_POP([C])
  if test "x$amg4psblas_cv_superludist_major" == "xunknown" ; then
     AC_MSG_CHECKING([for superlu_dist version 4])
     AC_LANG_PUSH([C])
     ac_cc=${MPICC-$CC}
     AC_COMPILE_IFELSE(
         [AC_LANG_SOURCE([[   #include "superlu_ddefs.h" 
			    int testdslud()
			    {  LUstruct_t *LUstruct;
			       int n; 
			       LUstructInit(n, LUstruct);     
			    }]])],
	 [ AC_MSG_RESULT([yes]);   amg4psblas_cv_superludist_major="4"; amg4psblas_cv_superludist_minor="";],
	 [ AC_MSG_RESULT([no]);    amg4psblas_cv_superludist_major="3"; amg4psblas_cv_superludist_minor="";])
     AC_LANG_POP([C])
     if test "x$amg4psblas_cv_superludist_major" == "x4" ; then
	AC_MSG_CHECKING([for superlu_dist version 5])
	AC_LANG_PUSH([C])
	ac_cc=${MPICC-$CC}
	AC_COMPILE_IFELSE(
            [AC_LANG_SOURCE([[   #include "superlu_ddefs.h"
				 int testdslud()
				 {     superlu_dist_options_t options;
				       int n;
				       set_default_options_dist(&options);
				 }]])],
	    [ AC_MSG_RESULT([yes]);  amg4psblas_cv_superludist_major="5";  amg4psblas_cv_superludist_minor="";],
	    [ AC_MSG_RESULT([no]);   amg4psblas_cv_superludist_major="4";  amg4psblas_cv_superludist_minor="";])
	AC_LANG_POP([C])
     fi
   fi
   AC_MSG_NOTICE([SuperLU_dist version $amg4psblas_cv_superludist_major.$amg4psblas_cv_superludist_minor.])

 else
    SLUDIST_LIBS="";
    SLUDIST_INCLUDES="";
 fi
 LIBS="$save_LIBS";
 CPPFLAGS="$save_CPPFLAGS";
 CC="$save_CC";
 AC_LANG_POP([C])
])dnl 

dnl @synopsis PAC_CHECK_MUMPS
dnl
dnl Will try to find the MUMPS library and headers.
dnl
dnl Will use $CC
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_CHECK_MUMPS,
[AC_ARG_WITH(mumps, AC_HELP_STRING([--with-mumps=LIBNAME], [Specify the libname for MUMPS. Default: autodetect with minimum "-lmumps_common -lpord"]),
        [amg4psblas_cv_mumps=$withval],
        [amg4psblas_cv_mumps='-lsmumps -ldmumps -lcmumps -lzmumps -lmumps_common -lpord'])
 AC_ARG_WITH(mumpsdir, AC_HELP_STRING([--with-mumpsdir=DIR], [Specify the directory for MUMPS library and includes. Note: you will need to add auxiliary libraries with --extra-libs; this depends on how MUMPS was configured and installed, at a minimum you will need SCALAPACK and BLAS]),
        [amg4psblas_cv_mumpsdir=$withval],
        [amg4psblas_cv_mumpsdir=''])

AC_ARG_WITH(mumpsincdir, AC_HELP_STRING([--with-mumpsincdir=DIR], [Specify the directory for MUMPS includes.]),
        [amg4psblas_cv_mumpsincdir=$withval],
        [amg4psblas_cv_mumpsincdir=''])

AC_ARG_WITH(mumpsmoddir, AC_HELP_STRING([--with-mumpsmoddir=DIR], [Specify the directory for MUMPS Fortran modules.]),
        [amg4psblas_cv_mumpsmoddir=$withval],
        [amg4psblas_cv_mumpsmoddir=''])

AC_ARG_WITH(mumpslibdir, AC_HELP_STRING([--with-mumpslibdir=DIR], [Specify the directory for MUMPS library.]),
        [amg4psblas_cv_mumpslibdir=$withval],
        [amg4psblas_cv_mumpslibdir=''])

AC_LANG_PUSH([Fortran])
save_LIBS="$LIBS"
save_FC="$FC"
FC=${MPIFC}
if test "x$amg4psblas_cv_mumpsincdir" != "x"; then 
 AC_MSG_NOTICE([mumps dir $amg4psblas_cv_mumpsincdir]) 
 MUMPS_INCLUDES="-I$amg4psblas_cv_mumpsincdir"
elif test "x$amg4psblas_cv_mumpsdir" != "x"; then 
 AC_MSG_NOTICE([mumps dir $amg4psblas_cv_mumpsdir]) 
 MUMPS_INCLUDES="-I$amg4psblas_cv_mumpsdir"
fi
if test "x$amg4psblas_cv_mumpsmoddir" != "x"; then 
 AC_MSG_NOTICE([mumps dir $amg4psblas_cv_mumpsmoddir]) 
 MUMPS_MODULES="$FMFLAG$amg4psblas_cv_mumpsmoddir"
elif test "x$amg4psblas_cv_mumpsdir" != "x"; then 
 AC_MSG_NOTICE([mumps dir $amg4psblas_cv_mumpsdir]) 
 MUMPS_MODULES="$FMFLAG$amg4psblas_cv_mumpsdir"
fi
if test "x$amg4psblas_cv_mumpslibdir" != "x"; then 
   MUMPS_LIBS="-L$amg4psblas_cv_mumpslibdir"
elif test "x$amg4psblas_cv_mumpsdir" != "x"; then 
   MUMPS_LIBS="-L$amg4psblas_cv_mumpsdir"
fi

LIBS="$MUMPS_LIBS $save_LIBS $EXTRA_LIBS"
CPPFLAGS="$MUMPS_INCLUDES $save_CPPFLAGS"

AC_LANG_PUSH([C])
AC_CHECK_HEADER([dmumps_c.h],
 [pac_mumps_header_ok=yes],
 [pac_mumps_header_ok=no; MUMPS_INCLUDES=""])
if test "x$pac_mumps_header_ok" == "xno" ; then 
   dnl Maybe Include or include subdirs?
   unset ac_cv_header_dmumps_c_h
   MUMPS_INCLUDES="-I$amg4psblas_cv_mumpsdir/include"
   CPPFLAGS="$MUMPS_INCLUDES $save_CPPFLAGS"
   AC_CHECK_HEADER([dmumps_c.h],
		   [pac_mumps_header_ok=yes],
		   [pac_mumps_header_ok=no; MUMPS_INCLUDES=""])
fi
if test "x$pac_mumps_header_ok" == "xno" ; then 
   dnl Maybe Include or include subdirs?
   unset ac_cv_header_dmumps_c_h
   MUMPS_INCLUDES="-I$amg4psblas_cv_mumpsdir/Include"
   CPPFLAGS="$MUMPS_INCLUDES $save_CPPFLAGS"
   AC_CHECK_HEADER([dmumps_c.h],
		   [pac_mumps_header_ok=yes],
		   [pac_mumps_header_ok=no; MUMPS_INCLUDES=""])
fi
   
AC_LANG_POP([C])
ac_objext='o'
ac_ext='f90'
ac_fc="${MPIFC-$FC}";
save_FCFLAGS="$FCFLAGS";
FCFLAGS="$MUMPS_MODULES $save_FCFLAGS"
AC_COMPILE_IFELSE([
		    program test
		    use dmumps_struc_def
		    end program test],
		   [pac_mumps_fmods_ok=yes; amg4psblas_cv_mumpmoddir="$MUMPS_MODULES";],
		   [pac_mumps_fmods_ok=no; MUMPS_MODULES=""])
if test "x$pac_mumps_fmods_ok" == "xno" ; then 
   dnl Maybe Include or include subdirs? 
   MUMPS_MODULES="$FMFLAG$amg4psblas_cv_mumpsdir/include"
   FCFLAGS="$MUMPS_MODULES $save_FCFLAGS"
   
   AC_COMPILE_IFELSE([
		      program test
		      use dmumps_struc_def
		      end program test],
		     [pac_mumps_fmods_ok=yes amg4psblas_cv_mumpsmoddir="$MUMPS_MODULES";],
		     [pac_mumps_fmods_ok=no; MUMPS_MODULES=""])
fi
if test "x$pac_mumps_fmods_ok" == "xno" ; then 
   dnl Maybe Include or include subdirs? 
   MUMPS_MODULES="$FMFLAG$amg4psblas_cv_mumpsdir/Include"
   FCFLAGS="$MUMPS_MODULES $save_FCFLAGS"
   
   AC_COMPILE_IFELSE([
		      program test
		      use dmumps_struc_def
		      end program test],
		     [pac_mumps_fmods_ok=yes amg4psblas_cv_mumpsmoddir="$MUMPS_MODULES";],
		     [pac_mumps_fmods_ok=no; MUMPS_MODULES=""])
fi
if test "x$pac_mumps_fmods_ok" == "xno" ; then 
   dnl Maybe Modules or modules subdirs? 
   MUMPS_MODULES="$FMFLAG$amg4psblas_cv_mumpsdir/modules"
   FCFLAGS="$MUMPS_MODULES $save_FCFLAGS"
   
   AC_COMPILE_IFELSE([
		      program test
		      use dmumps_struc_def
		      end program test],
		     [pac_mumps_fmods_ok=yes amg4psblas_cv_mumpsmoddir="$MUMPS_MODULES";],
		     [pac_mumps_fmods_ok=no; MUMPS_MODULES=""])
fi
if test "x$pac_mumps_fmods_ok" == "xno" ; then 
   dnl Maybe Modules or modules subdirs? 
   MUMPS_MODULES="$FMFLAG$amg4psblas_cv_mumpsdir/Modules"
   FCFLAGS="$MUMPS_MODULES $save_FCFLAGS"
   
   AC_COMPILE_IFELSE([
		      program test
		      use dmumps_struc_def
		      end program test],
		     [pac_mumps_fmods_ok=yes amg4psblas_cv_mumpsmoddir="$MUMPS_MODULES";],
		     [pac_mumps_fmods_ok=no; MUMPS_MODULES=""])
fi

if test "x$pac_mumps_fmods_ok" == "xno" ; then
    dnl Maybe only Xmumps.h?
    if test "x$amg4psblas_cv_mumpsincdir" != "x"; then 
	MUMPS_FINCLUDES="$FIFLAG$amg4psblas_cv_mumpsincdir"
    elif test "x$amg4psblas_cv_mumpsdir" != "x"; then
	MUMPS_FINCLUDES="$FIFLAG$amg4psblas_cv_mumpsincdir"
    fi
   
    FCFLAGS="$MUMPS_FINCLUDES $save_FCFLAGS"
    AC_COMPILE_IFELSE([
			 program test
			 include 'dmumps_struc.h'
			 end program test],
		      [pac_mumps_fincs_ok=yes; amg4psblas_cv_mumpmoddir="$MUMPS_FINCLUDES";],
		      [pac_mumps_fincs_ok=no; MUMPS_FINCLUDES=""])

    if test "x$pac_mumps_fincs_ok" == "xno" ; then 
	dnl Maybe Include or include subdirs? 
	MUMPS_FINCLUDES="$FIFLAG$amg4psblas_cv_mumpsdir/include"
	FCFLAGS="$MUMPS_FINCLUDES $save_FCFLAGS"
	AC_COMPILE_IFELSE([
			     program test
			     include 'dmumps_struc.h'
			     end program test],
			  [pac_mumps_fincs_ok=yes; amg4psblas_cv_mumpmoddir="$MUMPS_FINCLUDES";],
			  [pac_mumps_fincs_ok=no; MUMPS_FINCLUDES=""])
    fi
    if test "x$pac_mumps_fincs_ok" == "xno" ; then 
	dnl Maybe Include or include subdirs? 
	MUMPS_FINCLUDES="$FIFLAG$amg4psblas_cv_mumpsdir/Include"
	FCFLAGS="$MUMPS_FINCLUDES $save_FCFLAGS"
	AC_COMPILE_IFELSE([
			     program test
			     include 'dmumps_struc.h'
			     end program test],
			  [pac_mumps_fincs_ok=yes; amg4psblas_cv_mumpmoddir="$MUMPS_FINCLUDES";],
			  [pac_mumps_fincs_ok=no; MUMPS_FINCLUDES=""])
    fi
    if test "x$pac_mumps_fincs_ok" == "xno" ; then 
	dnl Maybe Modules or modules subdirs? 
	MUMPS_FINCLUDES="$FIFLAG$amg4psblas_cv_mumpsdir/modules"
	FCFLAGS="$MUMPS_FINCLUDES $save_FCFLAGS"
	AC_COMPILE_IFELSE([
			     program test
			     include 'dmumps_struc.h'
			     end program test],
			  [pac_mumps_fincs_ok=yes; amg4psblas_cv_mumpmoddir="$MUMPS_FINCLUDES";],
			  [pac_mumps_fincs_ok=no; MUMPS_FINCLUDES=""])
    fi
    if test "x$pac_mumps_fincs_ok" == "xno" ; then 
	dnl Maybe Modules or modules subdirs? 
	MUMPS_FINCLUDES="$FIFLAG$amg4psblas_cv_mumpsdir/Modules"
	FCFLAGS="$MUMPS_FINCLUDES $save_FCFLAGS"
	AC_COMPILE_IFELSE([
			     program test
			     include 'dmumps_struc.h'
			     end program test],
			  [pac_mumps_fincs_ok=yes; amg4psblas_cv_mumpmoddir="$MUMPS_FINCLUDES";],
			  [pac_mumps_fincs_ok=no; MUMPS_FINCLUDES=""])
    fi
    
fi   

if test "x$pac_mumps_fmods_ok" == "xyes" || test "x$pac_mumps_fincs_ok" == "xyes" ; then  
    MUMPS_LIBS="$amg4psblas_cv_mumps $MUMPS_LIBS"
    LIBS="$MUMPS_LIBS  $save_LIBS  $EXTRA_LIBS";
    AC_MSG_CHECKING([for dmumps in $MUMPS_LIBS])
			AC_TRY_LINK_FUNC(dmumps, 
					 [amg4psblas_cv_have_mumps=yes;pac_mumps_lib_ok=yes;],
					 [amg4psblas_cv_have_mumps=no;pac_mumps_lib_ok=no; 
					  MUMPS_LIBS=""; ])
  if test "x$pac_mumps_lib_ok" == "xno" ; then 
     dnl Maybe lib?
     MUMPS_LIBS="$amg4psblas_cv_mumps  -L$amg4psblas_cv_mumpsdir/lib";
     LIBS="$MUMPS_LIBS  $save_LIBS  $EXTRA_LIBS";
     AC_TRY_LINK_FUNC(dmumps, 
		     [amg4psblas_cv_have_mumps=yes;pac_mumps_lib_ok=yes;],
		     [amg4psblas_cv_have_mumps=no;pac_mumps_lib_ok=no; 
		      MUMPS_LIBS="";MUMPS_MODULES=""])
  fi
  AC_MSG_RESULT($pac_mumps_lib_ok)
fi
		      
 LIBS="$save_LIBS";
 CPPFLAGS="$save_CPPFLAGS";
 FCFLAGS="$save_FCFLAGS";
 FC="$save_FC";
 AC_LANG_POP([Fortran])
])dnl 



dnl @synopsis PAC_ARG_SERIAL_MPI
dnl
dnl Test for --enable-serial
dnl 
dnl 
dnl
dnl Example use:
dnl
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN([PAC_ARG_SERIAL_MPI],
[AC_MSG_CHECKING([whether we want serial  mpi stubs])
AC_ARG_ENABLE(serial,
AC_HELP_STRING([--enable-serial], 
[Specify whether to enable a fake mpi library to run in serial mode. ]),
[
pac_cv_serial_mpi="yes";
]
dnl ,
dnl [pac_cv_serial_mpi="no";]
)
if test x"$pac_cv_serial_mpi" == x"yes" ; then
   AC_MSG_RESULT([yes.])
else
 pac_cv_serial_mpi="no";
 AC_MSG_RESULT([no.])
fi
]
)

dnl @synopsis PAC_FORTRAN_CHECK_HAVE_MPI_MOD( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will determine if the fortran compiler MPIFC needs to include mpi.h or needs
dnl to use the mpi module.
dnl
dnl If yes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl 
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl Modified Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_FORTRAN_CHECK_HAVE_MPI_MOD,
 [AC_MSG_CHECKING([for Fortran MPI mod])
  AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='F90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
           program test
             use mpi
           end program test],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AC_FD_CC
		     cat conftest.$ac_ext >&AC_FD_CC
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])

dnl @synopsis PAC_FORTRAN_TEST_EXTENDS( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program checking the EXTENDS Fortran support.
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
AC_DEFUN(PAC_FORTRAN_TEST_EXTENDS,
ac_exeext=''
ac_ext='f90'
ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FCFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([support for Fortran EXTENDS])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='F90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
program conftest
  type foo
    integer :: i
  end type foo
  type, extends(foo) :: bar
    integer j
  end type bar 
  type(bar) :: barvar
end program conftest],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AC_FD_CC
		     cat conftest.$ac_ext >&AC_FD_CC
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])


dnl @synopsis PAC_ARG_LONG_INTEGERS
dnl
dnl Test for --enable-long-integers
dnl 
dnl 
dnl
dnl Example use:
dnl
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN([PAC_ARG_LONG_INTEGERS],
[
AC_MSG_CHECKING([whether we want long (8 bytes) integers])
AC_ARG_ENABLE(long-integers,
AC_HELP_STRING([--enable-long-integers], 
[Specify usage of 64 bits integers. ]),
[
pac_cv_long_integers="yes";
]
dnl ,
dnl [pac_cv_long_integers="no";]
)
if test x"$pac_cv_long_integers" == x"yes" ; then
   AC_MSG_RESULT([yes.])
else
 pac_cv_long_integers="no";
 AC_MSG_RESULT([no.])
fi
]
)

dnl @synopsis PAC_FORTRAN_TEST_CLASS_TBP( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program checking the TBP Fortran support.
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
AC_DEFUN(PAC_FORTRAN_TEST_CLASS_TBP,
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([support for Fortran CLASS TBP])
AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='F90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
module conftest_mod
  type foo
    integer :: i 
  contains
    procedure, pass(a) :: doit
    procedure, pass(a) :: getit
  end type foo

  private doit,getit
contains
  subroutine  doit(a) 
    class(foo) :: a
    
    a%i = 1
    write(*,*) 'FOO%DOIT base version'
  end subroutine doit
  function getit(a) result(res)
    class(foo) :: a
    integer :: res

    res = a%i
  end function getit

end module conftest_mod
program conftest
  use conftest_mod
  type(foo) :: foovar
end program conftest],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AC_FD_CC
		     cat conftest.$ac_ext >&AC_FD_CC
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])


dnl @synopsis PAC_FORTRAN_TEST_SOURCE( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program checking the SOURCE=  Fortran support.
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
AC_DEFUN(PAC_FORTRAN_TEST_SOURCE,
[AC_MSG_CHECKING([support for Fortran SOURCE= allocation])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='f90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
program xtt
  type foo
    integer :: i
  end type foo
  type, extends(foo) :: new_foo
    integer :: j
  end type new_foo
  class(foo), allocatable  :: fooab
  type(new_foo) :: nfv 
  integer :: info

  allocate(fooab, source=nfv, stat=info)

end program xtt],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AC_FD_CC
		     cat conftest.$ac_ext >&AC_FD_CC
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])

dnl @synopsis PAC_FORTRAN_HAVE_MOVE_ALLOC( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program with move_alloc (a Fortran 2003 function).
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl
AC_DEFUN([PAC_FORTRAN_HAVE_MOVE_ALLOC],
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([support for Fortran MOVE_ALLOC intrinsic])
 AC_LANG_PUSH([Fortran])
 ac_ext='f90';
 AC_COMPILE_IFELSE([ program test_move_alloc
		       integer, allocatable :: a(:), b(:)
		       allocate(a(3))
		       call move_alloc(a, b)
		       print *, allocated(a), allocated(b)
		       print *, b
		     end program test_move_alloc],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AC_FD_CC
		     cat conftest.$ac_ext >&AC_FD_CC
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])



dnl @synopsis PAC_FORTRAN_TEST_ISO_C_BIND( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program checking the ISO C Binding  Fortran support.
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
AC_DEFUN(PAC_FORTRAN_TEST_ISO_C_BIND,
[AC_MSG_CHECKING([support for Fortran ISO_C_BINDING module])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='f90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
program conftest
  use iso_c_binding
end program conftest],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AC_FD_CC
		     cat conftest.$ac_ext >&AC_FD_CC
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])



dnl @synopsis PAC_FORTRAN_TEST_VOLATILE( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program checking the VOLATILE Fortran support.
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
AC_DEFUN(PAC_FORTRAN_TEST_VOLATILE,
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([support for Fortran VOLATILE])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='F90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
program conftest
  integer, volatile :: i, j
end program conftest],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AC_FD_CC
		     cat conftest.$ac_ext >&AC_FD_CC
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])

dnl @synopsis PAC_MAKE_IS_GNUMAKE
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
define(PAC_MAKE_IS_GNUMAKE,[
AC_MSG_CHECKING(for gnumake)
MAKE=${MAKE:-make}

if $MAKE --version 2>&1 | grep -e"GNU Make" >/dev/null; then 
    AC_MSG_RESULT(yes)
    psblas_make_gnumake='yes'
else
    AC_MSG_RESULT(no)
    psblas_make_gnumake='no'
fi
])dnl


dnl @synopsis PAC_FORTRAN_TEST_GENERICS( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile a program checking the GENERIC Fortran support.
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
AC_DEFUN(PAC_FORTRAN_TEST_GENERICS,
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([test GENERIC interfaces])
AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='F90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
module conftest

  interface foo 
    subroutine i_sub_foo(v)
      integer, intent(inout) :: v(:)
    end subroutine i_sub_foo
  end interface foo

  interface bar
    procedure i_sub_foo
  end interface bar

end module conftest],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AC_FD_CC
		     cat conftest.$ac_ext >&AC_FD_CC
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])

dnl @synopsis PAC_FORTRAN_TEST_FLUSH( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program checking the FLUSH Fortran support.
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
AC_DEFUN(PAC_FORTRAN_TEST_FLUSH,
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([support for Fortran FLUSH statement])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='f90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
program conftest
   integer :: iunit=10
   open(10)
   write(10,*) 'Test '
   flush(10)
   close(10)
end program conftest],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AC_FD_CC
		     cat conftest.$ac_ext >&AC_FD_CC
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])



dnl @synopsis PAC_FORTRAN_TEST_ISO_FORTRAN_ENV( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will determine if the fortran compiler MPIFC supports ISO_FORTRAN_ENV
dnl
dnl If yes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl 
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_FORTRAN_TEST_ISO_FORTRAN_ENV,
[AC_MSG_CHECKING([support for ISO_FORTRAN_ENV])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='f90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
           program test
             use iso_fortran_env
           end program test],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AC_FD_CC
		     cat conftest.$ac_ext >&AC_FD_CC
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])


dnl @synopsis PAC_FORTRAN_TEST_FINAL( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program checking the FINAL Fortran support.
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
AC_DEFUN(PAC_FORTRAN_TEST_FINAL,
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([support for Fortran FINAL])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='f90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
module conftest_mod
  type foo
    integer :: i 
  contains
    final  :: destroy_foo
  end type foo

  private destroy_foo
contains
  subroutine destroy_foo(a)
    type(foo) :: a
     ! Just a test
  end subroutine destroy_foo
end module conftest_mod
program conftest
  use conftest_mod
  type(foo) :: foovar
end program conftest],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AC_FD_CC
		     cat conftest.$ac_ext >&AC_FD_CC
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])

dnl @synopsis PAC_FORTRAN_TEST_SAME_TYPE( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program checking the SAME_TYPE_AS Fortran support.
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
AC_DEFUN(PAC_FORTRAN_TEST_SAME_TYPE,
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([support for Fortran SAME_TYPE_AS])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='f90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
program stt
  type foo
    integer :: i
  end type foo
  type, extends(foo) :: new_foo
    integer :: j
  end type new_foo
  type(foo) :: foov
  type(new_foo) :: nfv1, nfv2

    
  write(*,*) 'foov == nfv1? ', same_type_as(foov,nfv1)
  write(*,*) 'nfv2 == nfv1? ', same_type_as(nfv2,nfv1)
end program stt],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AC_FD_CC
		     cat conftest.$ac_ext >&AC_FD_CC
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])

dnl @synopsis PAC_FORTRAN_TEST_EXTENDS_TYPE( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program checking the EXTENDS_TYPE_OF Fortran support.
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
AC_DEFUN(PAC_FORTRAN_TEST_EXTENDS_TYPE,
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([support for Fortran EXTENDS_TYPE_OF])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='f90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
program xtt
  type foo
    integer :: i
  end type foo
  type, extends(foo) :: new_foo
    integer :: j
  end type new_foo
  type(foo) :: foov
  type(new_foo) :: nfv1, nfv2

  write(*,*) 'nfv1 extends foov? ', extends_type_of(nfv1,foov)
end program xtt],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AC_FD_CC
		     cat conftest.$ac_ext >&AC_FD_CC
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])


dnl @synopsis PAC_FORTRAN_TEST_MOLD( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program checking the MOLD=  Fortran support.
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
AC_DEFUN(PAC_FORTRAN_TEST_MOLD,
[AC_MSG_CHECKING([support for Fortran MOLD= allocation])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='f90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
program xtt
  type foo
    integer :: i
  end type foo
  type, extends(foo) :: new_foo
    integer :: j
  end type new_foo
  class(foo), allocatable  :: fooab
  type(new_foo) :: nfv 
  integer :: info

  allocate(fooab, mold=nfv, stat=info)

end program xtt],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AC_FD_CC
		     cat conftest.$ac_ext >&AC_FD_CC
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])


dnl @synopsis PAC_BLAS([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl modified from ACX_BLAS([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the BLAS
dnl linear-algebra interface (see http://www.netlib.org/blas/). On
dnl success, it sets the BLAS_LIBS output variable to hold the
dnl requisite library linkages.
dnl
dnl To link with BLAS, you should link with:
dnl
dnl 	$BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order. FLIBS is the output variable of the
dnl AC_F77_LIBRARY_LDFLAGS macro (called if necessary by ACX_BLAS), and
dnl is sometimes necessary in order to link with F77 libraries. Users
dnl will also need to use AC_F77_DUMMY_MAIN (see the autoconf manual),
dnl for the same reason.
dnl
dnl Many libraries are searched for, from ATLAS to CXML to ESSL. The
dnl user may also use --with-blas=<lib> in order to use some specific
dnl BLAS library <lib>. In order to link successfully, however, be
dnl aware that you will probably need to use the same Fortran compiler
dnl (which can be set via the F77 env. var.) as was used to compile the
dnl BLAS library.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a BLAS
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands to
dnl run it if it is not found. If ACTION-IF-FOUND is not specified, the
dnl default action will define HAVE_BLAS.
dnl
dnl This macro requires autoconf 2.50 or later.
dnl
dnl @category InstalledPackages
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl @version 2001-12-13
dnl @license GPLWithACException
dnl
dnl modified by salvatore.filippone@uniroma2.it
dnl
dnl shifted check for ESSL as it was generating erroneous results on
dnl AIX SP5. 
dnl Modified with new name to handle Fortran compilers (such as NAG) 
dnl for which the linking MUST be done with the compiler (i.e.: 
dnl trying to link the Fortran version of the BLAS with the C compiler 
dnl would fail even when linking in the compiler's library)

AC_DEFUN([PAC_BLAS], [
AC_PREREQ(2.50)
dnl AC_REQUIRE([AC_FC_LIBRARY_LDFLAGS])
pac_blas_ok=no 

AC_ARG_WITH(blas,
	[AC_HELP_STRING([--with-blas=<lib>], [use BLAS library <lib>])])
case $with_blas in
	yes | "") ;;
	no) pac_blas_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) BLAS_LIBS="$with_blas" ;;
	*) BLAS_LIBS="-l$with_blas" ;;
esac
AC_ARG_WITH(blasdir,
	[AC_HELP_STRING([--with-blasdir=<dir>], [search for BLAS library in <dir>])])
case $with_blasdir in
  "") ;;
      *) if test -d $with_blasdir; then 
	    BLAS_LIBDIR="-L$with_blasdir";
	    fi ;;
esac
# Get fortran linker names of BLAS functions to check for.
#AC_FC_FUNC(sgemm)
#AC_FC_FUNC(dgemm)

pac_blas_save_LIBS="$LIBS"
#LIBS="$LIBS $FLIBS"
AC_LANG([Fortran])

# First, check BLAS_LIBS environment variable
if test $pac_blas_ok = no; then
if test "x$BLAS_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $BLAS_LIBDIR $LIBS"
	AC_MSG_CHECKING([for sgemm in $BLAS_LIBS])
	AC_TRY_LINK_FUNC(sgemm, [pac_blas_ok=yes], [BLAS_LIBS=""])
	AC_MSG_RESULT($pac_blas_ok)
	LIBS="$save_LIBS"
fi
fi

LIBS="$BLAS_LIBDIR $save_LIBS "
# BLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
if test $pac_blas_ok = no; then
	AC_LANG([C])
	AC_CHECK_LIB(atlas, ATL_xerbla,
		[AC_LANG([Fortran])
		 AC_CHECK_LIB(f77blas, sgemm,
		[AC_LANG([C])
		 AC_CHECK_LIB(cblas, cblas_dgemm,
			[pac_blas_ok=yes
			 BLAS_LIBS="-lcblas -lf77blas -latlas $BLAS_LIBDIR"],
			[], [-lf77blas -latlas])],
			[], [-latlas])])

fi
if test $pac_blas_ok = no; then
	AC_LANG([C])
	AC_CHECK_LIB(satlas, ATL_xerbla,
		[AC_LANG([Fortran])
		 AC_CHECK_LIB(satlas, sgemm,
		[AC_LANG([C])
		 AC_CHECK_LIB(satlas, cblas_dgemm,
			[pac_blas_ok=yes
			 BLAS_LIBS="-lsatlas $BLAS_LIBDIR"],
			[], [-lsatlas])],
			[], [-lsatlas])])

fi

# BLAS in PhiPACK libraries? (requires generic BLAS lib, too)
if test $pac_blas_ok = no; then
        AC_LANG([Fortran])
	AC_CHECK_LIB(blas, sgemm,
		[AC_CHECK_LIB(dgemm, dgemm,
		[AC_CHECK_LIB(sgemm, sgemm,
			[pac_blas_ok=yes; BLAS_LIBS="-lsgemm -ldgemm -lblas $BLAS_LIBDIR"],
			[], [-lblas])],
			[], [-lblas])])
fi


# BLAS in OpenBLAS? 
if test $pac_blas_ok = no; then
  AC_LANG([Fortran])
  AC_CHECK_LIB(openblas, sgemm, [pac_blas_ok=yes;BLAS_LIBS="-lopenblas $BLAS_LIBDIR"])
fi
# BLAS in Alpha CXML library? 
if test $pac_blas_ok = no; then
	AC_CHECK_LIB(cxml, sgemm, [pac_blas_ok=yes;BLAS_LIBS="-lcxml $BLAS_LIBDIR"])
fi

# BLAS in Alpha DXML library? (now called CXML, see above)
if test $pac_blas_ok = no; then
	AC_CHECK_LIB(dxml, sgemm, [pac_blas_ok=yes;BLAS_LIBS="-ldxml $BLAS_LIBDIR"])

fi

# BLAS in Sun Performance library?
if test $pac_blas_ok = no; then
	if test "x$GCC" != xyes; then # only works with Sun CC
		AC_CHECK_LIB(sunmath, acosp,
			[AC_CHECK_LIB(sunperf, sgemm,
        			[BLAS_LIBS="-xlic_lib=sunperf -lsunmath $BLAS_LIBDIR"
                                 pac_blas_ok=yes],[],[-lsunmath])])

	fi
fi

# BLAS in SCSL library?  (SGI/Cray Scientific Library)
if test $pac_blas_ok = no; then
	AC_CHECK_LIB(scs, sgemm, [pac_blas_ok=yes; BLAS_LIBS="-lscs $BLAS_LIBDIR"])
fi

# BLAS in SGIMATH library?
if test $pac_blas_ok = no; then
	AC_CHECK_LIB(complib.sgimath, $sgemm,
		     [pac_blas_ok=yes; BLAS_LIBS="-lcomplib.sgimath $BLAS_LIBDIR"])
fi

# BLAS in IBM ESSL library? (requires generic BLAS lib, too)
if test $pac_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm,
		[AC_CHECK_LIB(essl, sgemm,
			[pac_blas_ok=yes; BLAS_LIBS="-lessl -lblas $BLAS_LIBDIR"],
			[], [-lblas $FLIBS])])
	fi
# BLAS in generic BLAS library? 
if test $pac_blas_ok = no; then
  AC_LANG([Fortran])
  AC_CHECK_LIB(blas, sgemm, , [pac_blas_ok=yes;BLAS_LIBS="-lblas $BLAS_LIBDIR"])
fi
	
# BLAS linked to by default?  (happens on some supercomputers)
if test $pac_blas_ok = no; then
	AC_TRY_LINK_FUNC(sgemm, [pac_blas_ok=yes], [BLAS_LIBS=""])
dnl	AC_CHECK_FUNC(sgemm, [pac_blas_ok=yes])
fi

# Generic BLAS library?
if test $pac_blas_ok = no; then
  AC_LANG([Fortran])
  AC_CHECK_LIB(blas, sgemm, [pac_blas_ok=yes; BLAS_LIBS="-lblas $BLAS_LIBDIR"])
fi

dnl AC_SUBST(BLAS_LIBS)

LIBS="$pac_blas_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$pac_blas_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_BLAS,1,[Define if you have a BLAS library.]),[$1])
        :
else
        pac_blas_ok=no
        $2
fi
])dnl PAC_BLAS



dnl @synopsis PAC_LAPACK([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl @synopsis ACX_LAPACK([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the LAPACK
dnl linear-algebra interface (see http://www.netlib.org/lapack/). On
dnl success, it sets the LAPACK_LIBS output variable to hold the
dnl requisite library linkages.
dnl
dnl To link with LAPACK, you should link with:
dnl
dnl     $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order. BLAS_LIBS is the output variable of the ACX_BLAS
dnl macro, called automatically. FLIBS is the output variable of the
dnl AC_F77_LIBRARY_LDFLAGS macro (called if necessary by ACX_BLAS), and
dnl is sometimes necessary in order to link with F77 libraries. Users
dnl will also need to use AC_F77_DUMMY_MAIN (see the autoconf manual),
dnl for the same reason.
dnl
dnl The user may also use --with-lapack=<lib> in order to use some
dnl specific LAPACK library <lib>. In order to link successfully,
dnl however, be aware that you will probably need to use the same
dnl Fortran compiler (which can be set via the F77 env. var.) as was
dnl used to compile the LAPACK and BLAS libraries.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a LAPACK
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands to
dnl run it if it is not found. If ACTION-IF-FOUND is not specified, the
dnl default action will define HAVE_LAPACK.
dnl
dnl @category InstalledPackages
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl @version 2002-03-12
dnl @license GPLWithACException
dnl modified by salvatore.filippone@uniroma2.it
dnl shifted check for ESSL as it was generating erroneous results on
dnl AIX SP5. 
dnl Modified with new name to handle Fortran compilers (such as NAG) 
dnl for which the linking MUST be done with the compiler (i.e.: 
dnl trying to link the Fortran version of the BLAS with the C compiler 
dnl would fail even when linking in the compiler's library)

AC_DEFUN([PAC_LAPACK], [
AC_REQUIRE([PAC_BLAS])
pac_lapack_ok=no

AC_ARG_WITH(lapack,
        [AC_HELP_STRING([--with-lapack=<lib>], [use LAPACK library <lib>])])
case $with_lapack in
        yes | "") ;;
        no) pac_lapack_ok=disable ;;
        -* | */* | *.a | *.so | *.so.* | *.o) LAPACK_LIBS="$with_lapack" ;;
        *) LAPACK_LIBS="-l$with_lapack" ;;
esac

# Get fortran linker name of LAPACK function to check for.
#AC_FC_FUNC(cheev)

# We cannot use LAPACK if BLAS is not found
if test "x$pac_blas_ok" != xyes; then
        pac_lapack_ok=noblas
fi

# First, check LAPACK_LIBS environment variable
if test "x$LAPACK_LIBS" != x; then
        save_LIBS="$LIBS"; LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
        AC_MSG_CHECKING([for cheev in $LAPACK_LIBS])
	AC_LANG([Fortran])
	dnl Warning : square brackets are EVIL!
	cat > conftest.$ac_ext <<EOF
        program test_cheev 
          call cheev
        end 
EOF
	if AC_TRY_EVAL(ac_link) && test -s conftest${ac_exeext}; then
	  pac_lapack_ok=yes
	  AC_MSG_RESULT([yes])	
	else
	  AC_MSG_RESULT([no])	
	  echo "configure: failed program was:" >&AC_FD_CC
	  cat conftest.$ac_ext >&AC_FD_CC
	fi 
	rm -f conftest*
        LIBS="$save_LIBS"
        if test pac_lapack_ok = no; then
                LAPACK_LIBS=""
        fi
        AC_LANG([C])
fi

# LAPACK linked to by default?  (is sometimes included in BLAS lib)
if test $pac_lapack_ok = no; then
        save_LIBS="$LIBS"; LIBS="$LIBS $BLAS_LIBS $FLIBS"
        AC_MSG_CHECKING([for cheev in default libs])
	AC_LANG([Fortran])
	dnl Warning : square brackets are EVIL!
	cat > conftest.$ac_ext <<EOF
        program test_cheev 
          call cheev
        end 
EOF
	if AC_TRY_EVAL(ac_link) && test -s conftest${ac_exeext}; then
	  pac_lapack_ok=yes
	  AC_MSG_RESULT([yes])	
	else
	  AC_MSG_RESULT([no])	
	  echo "configure: failed program was:" >&AC_FD_CC
	  cat conftest.$ac_ext >&AC_FD_CC
	fi 
	rm -f conftest*
        LIBS="$save_LIBS"
        AC_LANG([C])
fi

# Generic LAPACK library?
for lapack in lapack lapack_rs6k; do
        if test $pac_lapack_ok = no; then
                save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
		AC_LANG([Fortran])
		AC_CHECK_LIB($lapack, cheev,
                    [pac_lapack_ok=yes; LAPACK_LIBS="-l$lapack"], [], [$FLIBS])
		AC_LANG([C])
                LIBS="$save_LIBS"
        fi
done

AC_SUBST(LAPACK_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$pac_lapack_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_LAPACK,1,[Define if you have LAPACK library.]),[$1])
        :
else
        pac_lapack_ok=no
        $2
fi
])dnl PAC_LAPACK


dnl @synopsis PAC_ARG_WITH_IPK
dnl
dnl Test for --with-ipk 
dnl 
dnl 
dnl
dnl Example use: --with-ipk=4
dnl
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN([PAC_ARG_WITH_IPK],
[
AC_MSG_CHECKING([what size in bytes we want for local indices and data])
AC_ARG_WITH(ipk,
	    AC_HELP_STRING([--with-ipk=<bytes>], 
			   [Specify the size in bytes for local indices and data, default 4 bytes. ]),
	    [pac_cv_ipk_size=$withval;],
	    [pac_cv_ipk_size=4;]
	   )
if test x"$pac_cv_ipk_size" == x"4"  || test  x"$pac_cv_ipk_size" == x"8" ; then
   AC_MSG_RESULT([Size: $pac_cv_ipk_size.])
else
  AC_MSG_RESULT([Unsupported value for IPK: $pac_cv_ipk_size, defaulting to 4.])
  pac_cv_ipk_size=4;
fi
]
)

dnl @synopsis PAC_ARG_WITH_LPK
dnl
dnl Test for --with-lpk 
dnl 
dnl 
dnl
dnl Example use: --with-lpk=8
dnl
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN([PAC_ARG_WITH_LPK],
[
 AC_MSG_CHECKING([what size in bytes we want for global indices and data])
 AC_ARG_WITH(lpk,
	     AC_HELP_STRING([--with-lpk=<bytes>], 
			    [Specify the size in bytes for global indices and data, default 8 bytes. ]),
	     [pac_cv_lpk_size=$withval;],
	     [pac_cv_lpk_size=8;]
	    )
if test x"$pac_cv_lpk_size" == x"4" || test x"$pac_cv_lpk_size" == x"8"; then
  AC_MSG_RESULT([Size: $pac_cv_lpk_size.])
else
  AC_MSG_RESULT([Unsupported value for LPK: $pac_cv_lpk_size, defaulting to 8.])
  pac_cv_lpk_size=8;
fi
]
)


dnl @synopsis PAC_CHECK_METIS
dnl
dnl Will try to find the METIS library and headers.
dnl
dnl Will use $CC
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_CHECK_METIS,
[AC_ARG_WITH(metis, AC_HELP_STRING([--with-metis=LIBNAME], [Specify the library name for METIS library. 
Default: "-lmetis"]),
        [psblas_cv_metis=$withval],
        [psblas_cv_metis='-lmetis'])
AC_ARG_WITH(metisincfile, AC_HELP_STRING([--with-metisincfile=DIR], [Specify the name  for METIS include file.]),
        [psblas_cv_metisincfile=$withval],
        [psblas_cv_metisincfile='metis.h'])
AC_ARG_WITH(metisdir, AC_HELP_STRING([--with-metisdir=DIR], [Specify the directory for METIS library and includes.]),
        [psblas_cv_metisdir=$withval],
        [psblas_cv_metisdir=''])
AC_ARG_WITH(metisincdir, AC_HELP_STRING([--with-metisincdir=DIR], [Specify the directory for METIS includes.]),
        [psblas_cv_metisincdir=$withval],
        [psblas_cv_metisincdir=''])
AC_ARG_WITH(metislibdir, AC_HELP_STRING([--with-metislibdir=DIR], [Specify the directory for METIS library.]),
        [psblas_cv_metislibdir=$withval],
        [psblas_cv_metislibdir=''])

AC_LANG([C])
SAVE_LIBS="$LIBS"
SAVE_CPPFLAGS="$CPPFLAGS"
if test "x$psblas_cv_metisdir" != "x"; then 
   METIS_LIBDIR="-L$psblas_cv_metisdir"
   LIBS="-L$psblas_cv_metisdir $LIBS"
   METIS_INCLUDES="-I$psblas_cv_metisdir"
   CPPFLAGS="$METIS_INCLUDES $CPPFLAGS"
fi
if test "x$psblas_cv_metisincdir" != "x"; then 
   METIS_INCLUDES="-I$psblas_cv_metisincdir"
   CPPFLAGS="$METIS_INCLUDES $CPPFLAGS"
fi
if test "x$psblas_cv_metislibdir" != "x"; then 
   LIBS="-L$psblas_cv_metislibdir $LIBS"
   METIS_LIBDIR="-L$psblas_cv_metislibdir"
fi

AC_MSG_NOTICE([metis dir $psblas_cv_metisdir])
AC_CHECK_HEADERS([limits.h $psblas_cv_metisincfile],
 [pac_metis_header_ok=yes],
 [pac_metis_header_ok=no; METIS_INCLUDES=""])
if test "x$pac_metis_header_ok" == "xno" ; then
dnl Maybe Include or include subdirs? 
  unset ac_cv_header_metis_h
  METIS_INCLUDES="-I$psblas_cv_metisdir/include -I$psblas_cv_metisdir/Include "
  CPPFLAGS="$METIS_INCLUDES $SAVE_CPPFLAGS"

 AC_MSG_CHECKING([for metis_h in $METIS_INCLUDES])
 AC_CHECK_HEADERS([limits.h  $psblas_cv_metisincfile],
    [pac_metis_header_ok=yes],
    [pac_metis_header_ok=no; METIS_INCLUDES=""])
fi
if test "x$pac_metis_header_ok" == "xno" ; then
dnl Maybe new structure with METIS UFconfig METIS? 
   unset ac_cv_header_metis_h
   METIS_INCLUDES="-I$psblas_cv_metisdir/UFconfig -I$psblas_cv_metisdir/METIS/Include -I$psblas_cv_metisdir/METIS/Include"
   CPPFLAGS="$METIS_INCLUDES $SAVE_CPPFLAGS"
   AC_CHECK_HEADERS([limits.h  $psblas_cv_metisincfile],
     [pac_metis_header_ok=yes],
     [pac_metis_header_ok=no; METIS_INCLUDES=""])
fi

if test "x$pac_metis_header_ok" == "xyes" ; then
   AC_LANG_PUSH([C])
   AC_MSG_CHECKING([for METIS integer size])
   AC_LINK_IFELSE([AC_LANG_SOURCE(
	#include <stdio.h>
	#include "$psblas_cv_metisincfile"
        void main(){
		    printf("%d\n",IDXTYPEWIDTH);
		    }
	       )],
	       [pac_cv_metis_idx=`./conftest${ac_exeext} | sed 's/^ *//'`],
	       [pac_cv_metis_idx="unknown"])
      AC_MSG_RESULT($pac_cv_metis_idx)

   AC_LANG_POP()
fi

if test "x$pac_metis_header_ok" == "xyes" ; then
   AC_LANG_PUSH([C])
   AC_MSG_CHECKING([for METIS real size])
   AC_LINK_IFELSE([AC_LANG_SOURCE(
	#include <stdio.h>
	#include "$psblas_cv_metisincfile"
        void main(){
		    printf("%d\n",REALTYPEWIDTH);
		    }
	       )],
	       [pac_cv_metis_real=`./conftest${ac_exeext} | sed 's/^ *//'`],
	       [pac_cv_metis_real="unknown"])
      AC_MSG_RESULT($pac_cv_metis_real)

   AC_LANG_POP()
fi

if test "x$pac_metis_header_ok" = "xyes" ; then 
      psblas_cv_metis_includes="$METIS_INCLUDES"
      METIS_LIBS="$psblas_cv_metis $METIS_LIBDIR"
      LIBS="$METIS_LIBS -lm $LIBS";
      AC_MSG_CHECKING([for METIS_PartGraphKway in $METIS_LIBS])
      AC_TRY_LINK_FUNC(METIS_PartGraphKway, 
       [psblas_cv_have_metis=yes;pac_metis_lib_ok=yes; ],
       [psblas_cv_have_metis=no;pac_metis_lib_ok=no; METIS_LIBS=""])
      AC_MSG_RESULT($pac_metis_lib_ok)
     if test "x$pac_metis_lib_ok" = "xno" ; then 
        dnl Maybe Lib or lib? 
        METIS_LIBDIR="-L$psblas_cv_metisdir/Lib -L$psblas_cv_metisdir/lib"
        METIS_LIBS="$psblas_cv_metis $METIS_LIBDIR"
        LIBS="$METIS_LIBS -lm $SAVE_LIBS"
        
      AC_MSG_CHECKING([for METIS_PartGraphKway in $METIS_LIBS])
      AC_TRY_LINK_FUNC(METIS_PartGraphKway, 
       [psblas_cv_have_metis=yes;pac_metis_lib_ok=yes; ],
       [psblas_cv_have_metis=no;pac_metis_lib_ok=no; METIS_LIBS=""])
      AC_MSG_RESULT($pac_metis_lib_ok)
     fi

     if test "x$pac_metis_lib_ok" = "xno" ; then 
        dnl Maybe METIS/Lib? 
        METIS_LIBDIR="-L$psblas_cv_metisdir/METIS/Lib -L$psblas_cv_metisdir/METIS/Lib"
        METIS_LIBS="$psblas_cv_metis $METIS_LIBDIR"
        LIBS="$METIS_LIBS -lm $SAVE_LIBS"
      AC_MSG_CHECKING([for METIS_PartGraphKway in $METIS_LIBS])
      AC_TRY_LINK_FUNC(METIS_PartGraphKway, 
       [psblas_cv_have_metis=yes;pac_metis_lib_ok="yes"; ],
       [psblas_cv_have_metis=no;pac_metis_lib_ok="no"; METIS_LIBS=""])
      AC_MSG_RESULT($pac_metis_lib_ok)
      fi
 fi
dnl AC_MSG_NOTICE([ metis lib ok $pac_metis_lib_ok])

 if test "x$pac_metis_lib_ok" = "xyes" ; then 
      AC_MSG_CHECKING([for METIS_SetDefaultOptions in $LIBS])
      AC_TRY_LINK_FUNC(METIS_SetDefaultOptions, 
        [psblas_cv_have_metis=yes;pac_metis_lib_ok=yes; ],
        [psblas_cv_have_metis=no;pac_metis_lib_ok="no. Unusable METIS version, sorry."; METIS_LIBS=""
      ])
      AC_MSG_RESULT($pac_metis_lib_ok)

fi

LIBS="$SAVE_LIBS";
CPPFLAGS="$SAVE_CPPFLAGS";
])dnl 


