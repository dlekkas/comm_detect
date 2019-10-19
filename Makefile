# Simple Makefile

OUTDIR = bin
OBJDIR = obj
SRCDIR = src
TSTDIR = tests
INCLUDEDIR = include

CC = g++
CFLAGS = -c -std=c++14 -Wall -O3 -g
BOOSTFLAGS = -lpthread

all: $(OUTDIR)/plp_test

$(OUTDIR)/plp_test: $(OBJDIR)/plp_test.o $(OBJDIR)/graph.o
	$(CC) $(OBJDIR)/plp_test.o $(OBJDIR)/graph.o -o $@ $(BOOSTFLAGS)

$(OBJDIR)/plp_test.o: $(TSTDIR)/plp_test.cc $(INCLUDEDIR)/graph.h
	$(CC) -o $@ $(CFLAGS) $(TSTDIR)/plp_test.cc

$(OBJDIR)/graph.o: $(SRCDIR)/graph.cc $(INCLUDEDIR)/graph.h
	$(CC) -o $@ $(CFLAGS) $(SRCDIR)/graph.cc

clean:
	rm -f $(OBJDIR)/*.o $(OUTDIR)/*
