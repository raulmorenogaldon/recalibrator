#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([recalibrator], [0.1])

AC_CONFIG_SRCDIR([src/recalibrate.c])
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_AUX_DIR([config])

AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_AWK
AC_PROG_CC
AC_PROG_MAKE_SET
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_LIBTOOL
AC_SUBST(LIBTOOL_DEPS)

# Checks for libraries.
AC_CHECK_LIB([z], [gzread])
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([xml2],[xmlStrPrintf])
AC_CHECK_LIB([pthread],[pthread_mutex_init])
AC_CHECK_LIB([cprops],[cp_mempool_alloc])
AC_CHECK_LIB([openmp],[omp_num_threads])
AC_CHECK_LIB([argtable2],[arg_lit0])

AC_LIB_SAMTOOLS
AC_LIB_COMPBIO


# Checks for header files.
AC_HEADER_STDC
AC_HEADER_TIME
AC_CHECK_HEADERS([stdio.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_UINT32_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([pow strdup strrchr])


AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_FILES([\
				Makefile	\
				src/Makefile	\
				src/common/Makefile	\
				src/aux_library/Makefile	\
				src/bam_recal/Makefile	\
				])
AC_OUTPUT
