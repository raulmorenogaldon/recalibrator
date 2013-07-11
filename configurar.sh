#!/bin/sh

CFLAGS="-O3 -std=c99 -Wall -D_REENTRANT -D_XOPEN_SOURCE=500 -fPIC"		\
./configure 	--with-samtools=/home/rmoreno/ext/samtools-0.1.18/ 	\
		--with-compbio=/home/rmoreno/workspace/compbio-hpg/	\
		--prefix=/home/rmoreno/recalibrator
