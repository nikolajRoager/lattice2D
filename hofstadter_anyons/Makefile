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
default:  densities.pdf


densities.pdf : densities.tsv
;gnuplot plot.gpi


log.txt densities.tsv: $(OUTDIR)/two_anyons.exe Makefile
;$(OUTDIR)/two_anyons.exe 6 6 1 1 4 4  1> log.txt 2>  densities.tsv



$(OUTDIR)/two_anyons.exe: $(ODIR)/two_anyons.o $(ODIR)/get_state.o $(ODIR)/generate_states.o $(ODIR)/approx.o
;$(CC) -o $@ $^ $(CXXFLAGS) $(LIBS)

$(ODIR)/two_anyons.o: $(SRCDIR)/two_anyons.cpp $(IDIR)/generate_states.hpp $(IDIR)/approx.hpp
;$(CC) -c -o $@ $< $(CXXFLAGS) $(LIBS)

$(ODIR)/generate_states.o: $(SRCDIR)/generate_states.cpp $(IDIR)/generate_states.hpp  $(IDIR)/approx.hpp
;$(CC) -c -o $@ $< $(CXXFLAGS) $(LIBS)

$(ODIR)/get_state.o: $(SRCDIR)/get_state.cpp $(IDIR)/get_state.hpp  $(IDIR)/approx.hpp
;$(CC) -c -o $@ $< $(CXXFLAGS) $(LIBS)

$(ODIR)/approx.o: $(SRCDIR)/approx.cpp $(IDIR)/approx.hpp
;$(CC) -c -o $@ $< $(CXXFLAGS) $(LIBS)


.PHONY: clean

clean:
;rm -f $(ODIR)/*.o
;rm -f $(OUTDIR)/*.exe
;rm -f log*
;rm -f *.tsv
;rm -f *.pdf
