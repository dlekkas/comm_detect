# Simple Makefile

OUTDIR = bin
OBJDIR = obj
SRCDIR = src
TSTDIR = tests
INCLUDEDIR = include

CC = /cluster/apps/gcc/6.3.0/bin/g++
MPICC = mpicxx
CFLAGS = -c -std=c++14 -Wall -O3 -g
BOOSTFLAGS = -lpthread
OPENMPFLAGS = -fopenmp

all: $(OUTDIR)/plp_test $(OUTDIR)/plm_test $(OUTDIR)/plm_test_mpi $(OUTDIR)/modularity_test

$(OUTDIR)/plp_test: $(OBJDIR)/plp_test.o $(OBJDIR)/graph.o $(OBJDIR)/plp.o
	$(CC) $(OBJDIR)/plp_test.o $(OBJDIR)/graph.o $(OBJDIR)/plp.o -o $@ $(BOOSTFLAGS) $(OPENMPFLAGS)

$(OUTDIR)/plm_test: $(OBJDIR)/plm_test.o $(OBJDIR)/graph.o $(OBJDIR)/plm.o
	$(CC) $(OBJDIR)/plm_test.o $(OBJDIR)/graph.o $(OBJDIR)/plm.o -o $@ $(BOOSTFLAGS) $(OPENMPFLAGS)

$(OUTDIR)/plm_test_mpi: $(OBJDIR)/plm_test_mpi.o $(OBJDIR)/plm_mpi.o # $(OBJDIR)/graph.o
	$(MPICC) $(OBJDIR)/plm_test_mpi.o $(OBJDIR)/graph.o $(OBJDIR)/plm_mpi.o -o $@

$(OUTDIR)/modularity_test: $(OBJDIR)/modularity_test.o $(OBJDIR)/graph.o
	$(CC) $(OBJDIR)/modularity_test.o $(OBJDIR)/graph.o -o $@ $(BOOSTFLAGS) $(OPENMPFLAGS)

$(OBJDIR)/plp_test.o: $(TSTDIR)/plp_test.cc $(INCLUDEDIR)/graph.h $(INCLUDEDIR)/plp.h
	$(CC) -o $@ $(CFLAGS) $(OPENMPFLAGS) $(TSTDIR)/plp_test.cc

$(OBJDIR)/plm_test.o: $(TSTDIR)/plm_test.cc $(INCLUDEDIR)/graph.h $(INCLUDEDIR)/plm.h
	$(CC) -o $@ $(CFLAGS) $(OPENMPFLAGS) $(TSTDIR)/plm_test.cc

$(OBJDIR)/plm_test_mpi.o: $(TSTDIR)/plm_test_mpi.cc $(INCLUDEDIR)/graph.h $(INCLUDEDIR)/plm.h
	$(MPICC) -o $@ $(CFLAGS) $(TSTDIR)/plm_test_mpi.cc

$(OBJDIR)/modularity_test.o: $(TSTDIR)/modularity_test.cc $(INCLUDEDIR)/network.h $(INCLUDEDIR)/modularity.h
	$(CC) -o $@ $(CFLAGS) $(OPENMPFLAGS) $(TSTDIR)/modularity_test.cc

$(OBJDIR)/graph.o: $(SRCDIR)/graph.cc $(INCLUDEDIR)/graph.h
	$(CC) -o $@ $(CFLAGS) $(SRCDIR)/graph.cc

$(OBJDIR)/plp.o: $(SRCDIR)/plp.cc $(INCLUDEDIR)/plp.h
	$(CC) -o $@ $(CFLAGS) $(OPENMPFLAGS) $(SRCDIR)/plp.cc

$(OBJDIR)/plm.o: $(SRCDIR)/plm.cc $(INCLUDEDIR)/plm.h $(INCLUDEDIR)/modularity.h
	$(CC) -o $@ $(CFLAGS) $(OPENMPFLAGS) $(SRCDIR)/plm.cc

$(OBJDIR)/plm_mpi.o: $(SRCDIR)/plm_mpi.cc $(INCLUDEDIR)/plm.h $(INCLUDEDIR)/modularity.h
	$(MPICC) -o $@ $(CFLAGS) $(SRCDIR)/plm_mpi.cc

clean:
	rm -f $(OBJDIR)/*.o $(OUTDIR)/*
