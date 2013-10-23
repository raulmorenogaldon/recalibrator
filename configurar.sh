#!/bin/sh

./configure 	--prefix=/home/rmoreno/recalibrator \
				--with-opencb=/home/rmoreno/workspace/OpenCB/ \
				CFLAGS="-O0 -g -pg -std=c99 -Wall -Winline -save-temps -D_REENTRANT -D_XOPEN_SOURCE=600 -fPIC -fopenmp"
				
#./configure 	--prefix=/home/rmoreno/recalibrator \
#				--with-opencb=/home/rmoreno/workspace/OpenCB/ \
#				CFLAGS="-O3 -std=c99 -Wall -D_REENTRANT -D_XOPEN_SOURCE=600 -fPIC"	

#./configure 	--prefix=/home/rmoreno/recalibrator \
#				--with-samtools=/home/rmoreno/ext/samtools-0.1.18/ \
#				--with-compbio=/home/rmoreno/workspace/compbio-hpg/ \
#				CFLAGS="-O0 -g -std=c99 -Wall -D_REENTRANT -D_XOPEN_SOURCE=600 -fPIC"		
				