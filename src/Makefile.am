ACLOCAL_AMFLAGS = -I m4

SUBDIRS = common aux_library bam_recal

INCLUDES = 		-I$(top_builddir)/src/common \
				-I$(top_builddir)/src/aux_library \
				-I$(top_builddir)/src/bam_recal 

lib_LTLIBRARIES = librecal.la
librecal_la_SOURCES = dummy.c
librecal_la_LIBADD = 	$(top_builddir)/src/common/libcommon.la	\
						$(top_builddir)/src/aux_library/librecal_aux.la	\
						$(top_builddir)/src/bam_recal/librecal_bam.la


bin_PROGRAMS = recal
recal_SOURCES = recalibrate.c
recal_LDADD = librecal.la

CLEANFILES = *~