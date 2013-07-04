ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src tests
EXTRA_DIST = AUTHORS ChangeLog NEWS README THANKS

DISTCHECK_CONFIGURE_FLAGS	= 	--with-samtools=/home/rmoreno/ext/samtools-0.1.18/ 	\
								--with-compbio=/home/rmoreno/workspace/compbio-hpg/	\
								--prefix=/home/rmoreno/recalibrator

CLEANFILES = *~