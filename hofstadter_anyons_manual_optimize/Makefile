#I don't like literal tabs
.RECIPEPREFIX := ;

CC = g++

IDIR =include
ODIR=obj
SRCDIR=src
OUTDIR=bin

CXXFLAGS=-Wall -I$(IDIR)

LIBS = -std=c++17 -O2 -larmadillo -O2 -lncurses

default:  q.pdf log0.txt log1.txt log2.txt log3.txt log4.txt log5.txt
;cat log0.txt
;cat log1.txt
;cat log2.txt
;cat log3.txt


$(ODIR)/minimize.o: $(SRCDIR)/minimize.cpp $(IDIR)/minimize.hpp
;$(CC) -c -o $@ $< $(CXXFLAGS) $(LIBS)


q.pdf : q0.tsv q1.tsv q2.tsv q3.tsv
;gnuplot plot.gpi


log0.txt q0.tsv: $(OUTDIR)/two_anyons.exe Makefile V0.tsv
;$(OUTDIR)/two_anyons.exe 6 6 1 1 5 5 V0.tsv 1> log0.txt 2>  q0.tsv


log1.txt q1.tsv: $(OUTDIR)/two_anyons.exe Makefile V1.tsv
;$(OUTDIR)/two_anyons.exe 6 6 1 1 5 5 V1.tsv 1> log1.txt 2>  q1.tsv

log2.txt q2.tsv: $(OUTDIR)/two_anyons.exe Makefile V2.tsv
;$(OUTDIR)/two_anyons.exe 6 6 1 1 5 5 V2.tsv 1> log2.txt 2>  q2.tsv


log3.txt q3.tsv: $(OUTDIR)/two_anyons.exe Makefile V3.tsv
;$(OUTDIR)/two_anyons.exe 6 6 1 1 5 5 V3.tsv 1> log3.txt 2>  q3.tsv



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
;rm -f q*.tsv
;rm -f  *.pdf
