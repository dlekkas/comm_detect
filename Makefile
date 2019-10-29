# Simple Makefile

OUTDIR = bin
OBJDIR = obj
SRCDIR = src
TSTDIR = tests
INCLUDEDIR = include

CC = g++
CFLAGS = -c -std=c++14 -Wall -O3 -g -fopenmp
BOOSTFLAGS = -lpthread
OPENMPFLAGS = -fopenmp

all: $(OUTDIR)/plp_test $(OUTDIR)/modularity_test

$(OUTDIR)/plp_test: $(OBJDIR)/plp_test.o $(OBJDIR)/graph.o $(OBJDIR)/plp.o
	$(CC) $(OBJDIR)/plp_test.o $(OBJDIR)/graph.o $(OBJDIR)/plp.o -o $@ $(BOOSTFLAGS) $(OPENMPFLAGS)

$(OUTDIR)/modularity_test: $(OBJDIR)/modularity_test.o $(OBJDIR)/graph.o
	$(CC) $(OBJDIR)/modularity_test.o $(OBJDIR)/graph.o -o $@ $(BOOSTFLAGS)

$(OBJDIR)/plp_test.o: $(TSTDIR)/plp_test.cc $(INCLUDEDIR)/graph.h $(INCLUDEDIR)/plp.h
	$(CC) -o $@ $(CFLAGS) $(TSTDIR)/plp_test.cc

$(OBJDIR)/modularity_test.o: $(TSTDIR)/modularity_test.cc $(INCLUDEDIR)/network.h $(INCLUDEDIR)/modularity.h
	$(CC) -o $@ $(CFLAGS) $(TSTDIR)/modularity_test.cc

$(OBJDIR)/graph.o: $(SRCDIR)/graph.cc $(INCLUDEDIR)/graph.h
	$(CC) -o $@ $(CFLAGS) $(SRCDIR)/graph.cc

$(OBJDIR)/plp.o: $(SRCDIR)/plp.cc $(INCLUDEDIR)/plp.h
	$(CC) -o $@ $(CFLAGS) $(SRCDIR)/plp.cc

$(OBJDIR)/plm.o: $(SRCDIR)/plm.cc $(INCLUDEDIR)/plm.h
	$(CC) -o $@ $(CFLAGS) $(SRCDIR)/plm.cc

clean:
	rm -f $(OBJDIR)/*.o $(OUTDIR)/*
