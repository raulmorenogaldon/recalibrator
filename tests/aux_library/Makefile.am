ACLOCAL_AMFLAGS = -I m4

TESTS = check_aux_library
check_PROGRAMS 				= check_aux_library
check_aux_library_SOURCES 	= check_aux_library.c
check_aux_library_CFLAGS 	= $(CHECK_CFLAGS)
check_aux_library_CPPFLAGS 	= -I$(top_builddir)/src/aux_library -I$(top_builddir)/src/common
check_aux_library_LDADD 	= $(top_builddir)/src/librecal.la $(CHECK_LIBS)


CLEANFILES = *~