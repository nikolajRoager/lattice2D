#I don't like literal tabs
.RECIPEPREFIX := ;

CC = g++

IDIR =include
ODIR=obj
SRCDIR=src
OUTDIR=bin

CXXFLAGS=-Wall -I$(IDIR)

LIBS = -std=c++17 -O2 -larmadillo -O2

#First test to see if this works
default:  eigenstates.pdf

example_states.pdf: example_states.tsv
;gnuplot plot.gpi

eigenstates.pdf : eigenstates.tsv
;gnuplot plot1.gpi


log.txt eigenstates.tsv: $(OUTDIR)/lattice_states.exe
;$(OUTDIR)/lattice_states.exe 4 6 8  1> log.txt 2>  eigenstates.tsv

example_log.txt example_states.tsv: $(OUTDIR)/debug_showstates.exe
;$(OUTDIR)/debug_showstates.exe 4 6 8  1> example_log.txt 2>  example_states.tsv


$(OUTDIR)/lattice_states.exe: $(ODIR)/lattice_states.o $(ODIR)/generate_states.o
;$(CC) -o $@ $^ $(CXXFLAGS) $(LIBS)

$(ODIR)/lattice_states.o: $(SRCDIR)/lattice_states.cpp $(IDIR)/generate_states.hpp
;$(CC) -c -o $@ $< $(CXXFLAGS) $(LIBS)

$(ODIR)/generate_states.o: $(SRCDIR)/generate_states.cpp $(IDIR)/generate_states.hpp
;$(CC) -c -o $@ $< $(CXXFLAGS) $(LIBS)

$(ODIR)/debug_showstates.o: $(SRCDIR)/debug_showstates.cpp $(IDIR)/generate_states.hpp
;$(CC) -c -o $@ $< $(CXXFLAGS) $(LIBS)

$(OUTDIR)/debug_showstates.exe: $(ODIR)/debug_showstates.o $(ODIR)/generate_states.o
;$(CC) -o $@ $^ $(CXXFLAGS) $(LIBS)

.PHONY: clean

clean:
;rm -f $(ODIR)/*.o
;rm -f $(OUTDIR)/*.exe
;rm -f log*
;rm -f *.tsv
;rm -f *.pdf
