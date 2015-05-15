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

