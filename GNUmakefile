#                    Copyright (C) 2015, UChicago Argonne, LLC
#                               All Rights Reserved
# 
#                               Generic IO (ANL-15-066)
#                     Hal Finkel, Argonne National Laboratory
# 
#                              OPEN SOURCE LICENSE
# 
# Under the terms of Contract No. DE-AC02-06CH11357 with UChicago Argonne,
# LLC, the U.S. Government retains certain rights in this software.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
#   1. Redistributions of source code must retain the above copyright notice,
#      this list of conditions and the following disclaimer.
# 
#   2. Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in the
#      documentation and/or other materials provided with the distribution.
# 
#   3. Neither the names of UChicago Argonne, LLC or the Department of Energy
#      nor the names of its contributors may be used to endorse or promote
#      products derived from this software without specific prior written
#      permission.
# 
# *****************************************************************************
# 
#                                  DISCLAIMER
# THE SOFTWARE IS SUPPLIED “AS IS” WITHOUT WARRANTY OF ANY KIND.  NEITHER THE
# UNTED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR
# UCHICAGO ARGONNE, LLC, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY,
# EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE
# ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, DATA, APPARATUS,
# PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
# PRIVATELY OWNED RIGHTS.
# 
# *****************************************************************************

CC = gcc
CXX = g++

MPICC = mpicc
MPICXX = mpicxx

all: fe-progs mpi-progs
sql: fe-sqlite

FEDIR = frontend
FE_CFLAGS := -g -fPIC -O3 -fopenmp
FE_CPPFLAGS := -Ithirdparty/blosc -Ithirdparty/sqlite -DGENERICIO_NO_MPI

MPIDIR = mpi
MPI_CFLAGS := -g -O3
MPI_CPPFLAGS := -Ithirdparty/blosc 

$(FEDIR):
	mkdir -p $(FEDIR)

$(FEDIR)/%.o: thirdparty/blosc/%.c | $(FEDIR)
	$(CC) $(FE_CFLAGS) $(FE_CPPFLAGS) -c -o $@ $<

$(FEDIR)/%.o: %.cxx | $(FEDIR)
	$(CXX) $(FE_CFLAGS) $(FE_CPPFLAGS) -c -o $@ $<

$(FEDIR)/GenericIOPrint: $(FEDIR)/GenericIOPrint.o $(FEDIR)/GenericIO.o  $(FEDIR)/blosc.o $(FEDIR)/blosclz.o $(FEDIR)/shuffle.o
	$(CXX) $(FE_CFLAGS) -o $@ $^ 

$(FEDIR)/GenericIOVerify: $(FEDIR)/GenericIOVerify.o $(FEDIR)/GenericIO.o $(FEDIR)/blosc.o $(FEDIR)/blosclz.o $(FEDIR)/shuffle.o
	$(CXX) $(FE_CFLAGS) -o $@ $^ 

FE_UNAME := $(shell uname -s)
ifeq ($(FE_UNAME),Darwin)
FE_SHARED := -bundle
else
FE_SHARED := -shared
endif
$(FEDIR)/GenericIOSQLite.so: $(FEDIR)/GenericIOSQLite.o $(FEDIR)/GenericIO.o
	$(CXX) $(FE_CFLAGS) $(FE_SHARED) -o $@ $^

SQLITE_CPPFLAGS := -DSQLITE_ENABLE_COLUMN_METADATA=1 -DSQLITE_DISABLE_DIRSYNC=1 -DSQLITE_ENABLE_FTS3=3 -DSQLITE_ENABLE_RTREE=1 -DSQLITE_ENABLE_UNLOCK_NOTIFY=1 -DSQLITE_ENABLE_LOAD_EXTENSION=1 -DHAVE_READLINE=1

$(FEDIR)/sqbuild:
	mkdir -p $(FEDIR)/sqbuild

$(FEDIR)/sqbuild/%.o: thirdparty/sqlite/%.c | $(FEDIR)/sqbuild
	$(CC) $(FE_CFLAGS) $(FE_CPPFLAGS) $(SQLITE_CPPFLAGS) -c -o $@ $<

$(FEDIR)/sqlite3: $(FEDIR)/sqbuild/sqlite3.o $(FEDIR)/sqbuild/shell.o
	$(CC) $(FE_CFLAGS) -o $@ $^ -pthread -lreadline -lrt -ldl

$(MPIDIR):
	mkdir -p $(MPIDIR)

$(MPIDIR)/%.o: thirdparty/blosc/%.c | $(MPIDIR)
	$(MPICC) $(MPI_CFLAGS) $(MPI_CPPFLAGS) -c -o $@ $<

$(MPIDIR)/%.o: %.cxx | $(MPIDIR)
	$(MPICXX) $(MPI_CFLAGS) $(MPI_CPPFLAGS) -c -o $@ $<

$(MPIDIR)/GenericIOPrint: $(MPIDIR)/GenericIOPrint.o $(MPIDIR)/GenericIO.o $(MPIDIR)/blosc.o $(MPIDIR)/blosclz.o $(MPIDIR)/shuffle.o
	$(MPICXX) $(MPI_CFLAGS) -o $@ $^ 

$(MPIDIR)/GenericIOVerify: $(MPIDIR)/GenericIOVerify.o $(MPIDIR)/GenericIO.o $(MPIDIR)/blosc.o $(MPIDIR)/blosclz.o $(MPIDIR)/shuffle.o
	$(MPICXX) $(MPI_CFLAGS) -o $@ $^ 

frontend-progs: $(FEDIR)/GenericIOPrint $(FEDIR)/GenericIOVerify
fe-progs: frontend-progs

mpi-progs: $(MPIDIR)/GenericIOPrint $(MPIDIR)/GenericIOVerify

frontend-sqlite: $(FEDIR)/GenericIOSQLite.so $(FEDIR)/sqlite3
fe-sqlite: frontend-sqlite

clean:
	rm -rf frontend mpi

