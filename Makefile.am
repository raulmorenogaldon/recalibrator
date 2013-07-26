ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src tests
EXTRA_DIST = AUTHORS ChangeLog NEWS README THANKS

#DISTCHECK_CONFIGURE_FLAGS	= 	--with-samtools=/home/rmoreno/ext/samtools-0.1.18/ 	\
#								--with-compbio=/home/rmoreno/workspace/compbio-hpg/	\
#								--prefix=/home/rmoreno/recalibrator

DISTCHECK_CONFIGURE_FLAGS	= 	--with-opencb=/home/rmoreno/workspace/OpenCB/	\
								--prefix=/home/rmoreno/recalibrator \
								CFLAGS="-O3 -std=c99 -Wall -D_REENTRANT -D_XOPEN_SOURCE=600 -fPIC"

CLEANFILES = *~