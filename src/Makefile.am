ACLOCAL_AMFLAGS = -I m4

SUBDIRS = aux_library bam_recal

lib_LTLIBRARIES = librecal.la
librecal_la_SOURCES = dummy.c
librecal_la_LIBADD = 	$(src_builddir)/aux_library/librecal_aux.la	\
						$(src_builddir)/bam_recal/librecal_bam.la


bin_PROGRAMS = recal
recal_SOURCES = recalibrate.c
recal_LDADD = librecal.la