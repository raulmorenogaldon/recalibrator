CC = gcc
compiler = gcc
#compiler = mpicc
#linker = tau_cc.sh
#linker = mpicc
linker = gcc

CFLAGS = -O3 -ansi -Wall -std=c99 -D_REENTRANT -D_XOPEN_SOURCE=500
CFLAGS_NOSSE2 = -O3 -ansi -Wall -std=c99 -D_REENTRANT -D_XOPEN_SOURCE=500 -mfpmath=387 -mno-sse2
CFLAGS_SSE2 = -O3 -ansi -Wall -std=c99 -D_REENTRANT -D_XOPEN_SOURCE=500 -mfpmath=387 -msse2
CFLAGS_DEBUG = -O0 -ansi -Wall -ggdb -std=c99 -D_REENTRANT -D_XOPEN_SOURCE=500

# Main folders
SRC_DIR = $(PWD)/src
INC_DIR = $(PWD)/include
COMPBIO_LIB_DIR = $(PWD)/../compbio-hpg
LIB_DIR = $(PWD)/lib
BIN_DIR = $(PWD)/bin

SAMTOOLS_DIR = $(HOME)/ext/samtools-0.1.18
CPROPS_DIR = $(HOME)/ext/cprops/include
LIBXML_DIR = $(HOME)/ext/libxml2/include/libxml2
LIBXML_LIB_DIR = $(HOME)/ext/libxml2/lib
CPROPS_LIB_DIR = $(HOME)/ext/cprops/lib
BIOINFO_LIBS_DIR = $(COMPBIO_LIB_DIR)/bioinfo-libs
COMMON_LIBS_DIR = $(COMPBIO_LIB_DIR)/common-libs

CONTAINERS_DIR = $(COMMON_LIBS_DIR)/containers
COMMONS_DIR = $(COMMON_LIBS_DIR)/commons
BIOFORMATS_DIR = $(BIOINFO_LIBS_DIR)/bioformats
BWT_DIR = $(BIOINFO_LIBS_DIR)/aligners/bwt

# Include and lib folders
INCLUDES = -I . -I ./include -I $(LIB_DIR) -I $(COMPBIO_LIB_DIR) -I $(BIOINFO_LIBS_DIR) -I $(COMMON_LIBS_DIR) -I $(INC_DIR) -I $(SAMTOOLS_DIR) -I $(CPROPS_DIR) -I $(LIBXML_DIR) -I $(BWT_DIR) -I /usr/include/libxml2 -I /usr/local/include
LIBS = -L$(LIB_DIR) -L$(COMPBIO_LIB_DIR) -L$(SAMTOOLS_DIR) -L$(CPROPS_LIB_DIR) -L$(LIBXML_LIB_DIR) -lm -lbam -lcurl -lxml2 -lpthread -lcprops -fopenmp -largtable2 

# Object file dependencies
MISC_OBJS = $(COMMONS_DIR)/*.o $(CONTAINERS_DIR)/*.o $(BIOFORMATS_DIR)/bam-sam/*.o $(BWT_DIR)/*.o 

# Project source files
#HPG_BAM_FILES = bamaux.c recal.c timestats.c bampool.c
HPG_BAM_FILES = bamaux.c recal.c timestats.c

# Project object files
#HPG_BAM_OBJS = bamaux.o recal.o timestats.o bampool.o
HPG_BAM_OBJS = bamaux.o recal.o timestats.o

ALL_OBJS = $(HPG_BAM_OBJS) $(MISC_OBJS)

# Targets
all: compile-dependencies recal recal-debug 

recal: compile-dependencies
	cd $(SRC_DIR) &&                                                         \
	$(compiler) $(CFLAGS) -c recalibrate.c $(HPG_BAM_FILES) $(INCLUDES) $(LIBS) &&    \
	$(linker) $(CFLAGS) -o $(BIN_DIR)/$@ recalibrate.o $(ALL_OBJS) $(INCLUDES) $(LIBS)
	
recal-debug: compile-dependencies
	cd $(SRC_DIR) &&                                                         \
	$(compiler) $(CFLAGS_DEBUG) -c recalibrate.c $(HPG_BAM_FILES) $(INCLUDES) $(LIBS) -DDEBUG &&    \
	$(linker) $(CFLAGS_DEBUG) -o $(BIN_DIR)/$@ recalibrate.o $(ALL_OBJS) $(INCLUDES) $(LIBS)

compile-dependencies: bam-dependencies
	cd $(COMMONS_DIR) && make compiler=$(compiler) &&           \
	cd $(CONTAINERS_DIR) && make compiler=$(compiler) 

bam-dependencies:
	cd $(BIOFORMATS_DIR)/bam-sam &&  \
        $(compiler) $(CFLAGS) -c -o $(BIOFORMATS_DIR)/bam-sam/alignment.o $(BIOFORMATS_DIR)/bam-sam/alignment.c $(INCLUDES) $(LIBS) && \
		$(compiler) $(CFLAGS) -c -o $(BIOFORMATS_DIR)/bam-sam/bam_file.o $(BIOFORMATS_DIR)/bam-sam/bam_file.c $(INCLUDES) $(LIBS) && \
		$(compiler) $(CFLAGS) -c -o $(BWT_DIR)/genome.o $(BWT_DIR)/genome.c $(INCLUDES) $(LIBS)
		
		
		
		
		
		
recal-sse2: bam-dependencies compile-dependencies
	cd $(SRC_DIR) &&                                                         \
	$(compiler) $(CFLAGS_SSE2) -c recalibrate.c $(HPG_BAM_FILES) $(INCLUDES) $(LIBS) &&    \
	$(linker) $(CFLAGS_SSE2) -o $(BIN_DIR)/$@ recalibrate.o $(ALL_OBJS) $(INCLUDES) $(LIBS)
		
		
		
recal-nosse2: bam-dependencies compile-dependencies
	cd $(SRC_DIR) &&                                                         \
	$(compiler) $(CFLAGS_NOSSE2) -c recalibrate.c $(HPG_BAM_FILES) $(INCLUDES) $(LIBS) &&    \
	$(linker) $(CFLAGS_NOSSE2) -o $(BIN_DIR)/$@ recalibrate.o $(ALL_OBJS) $(INCLUDES) $(LIBS)



clean:
	-rm -f $(SRC_DIR)/*~ $(SRC_DIR)/\#*\# $(SRC_DIR)/*.o 
	-rm -f $(COMMONS_DIR)/*.o
	-rm -f $(CONTAINERS_DIR)/*.o
	-rm -f $(BIOFORMATS_DIR)/bam-sam/*.o
	-rm -f $(BIN_DIR)/*
	-rm -f $(BIN_DIR)/*
