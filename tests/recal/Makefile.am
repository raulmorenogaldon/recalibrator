ACLOCAL_AMFLAGS = -I m4

TESTS = check_recal
check_PROGRAMS 			= check_recal
check_recal_SOURCES 	= check_recal.c
check_recal_CFLAGS 		= $(CHECK_CFLAGS)
check_recal_CPPFLAGS 	= -I$(top_builddir)/include
check_recal_LDADD 		= $(top_builddir)/src/librecal.la $(top_builddir)/src/timestats/libtimestats.la $(CHECK_LIBS)


CLEANFILES = *~