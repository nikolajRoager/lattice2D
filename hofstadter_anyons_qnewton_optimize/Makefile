#I don't like literal tabs
.RECIPEPREFIX := ;

CC = g++

IDIR =include
ODIR=obj
SRCDIR=src
OUTDIR=bin

CXXFLAGS=-Wall -I$(IDIR)

LIBS = -std=c++17 -O2 -larmadillo  -lpthread

default:  q.pdf log.txt
;cat log.txt


$(OUTDIR)/minimum_failure.exe: $(SRCDIR)/minimum_failure.cpp
;$(CC) -o $@ $< $(CXXFLAGS) $(LIBS)

$(ODIR)/minimize.o: $(SRCDIR)/minimize.cpp $(IDIR)/minimize.hpp
;$(CC) -c -o $@ $< $(CXXFLAGS) $(LIBS)


q.pdf : q.tsv
;gnuplot plot.gpi


log.txt q.tsv: $(OUTDIR)/two_anyons.exe Makefile V0.tsv
;$(OUTDIR)/two_anyons.exe 6 6 anyon_location.tsv V0.tsv 1> log.txt 2>  q.tsv


$(OUTDIR)/two_anyons.exe: $(ODIR)/two_anyons.o $(ODIR)/get_state.o $(ODIR)/generate_states.o $(ODIR)/approx.o $(ODIR)/minimize.o
;$(CC) -o $@ $^ $(CXXFLAGS) $(LIBS)

$(ODIR)/two_anyons.o: $(SRCDIR)/two_anyons.cpp $(IDIR)/generate_states.hpp $(IDIR)/approx.hpp  $(IDIR)/minimize.hpp
;$(CC) -c -o $@ $< $(CXXFLAGS) $(LIBS)

$(ODIR)/generate_states.o: $(SRCDIR)/generate_states.cpp $(IDIR)/generate_states.hpp  $(IDIR)/approx.hpp
;$(CC) -c -o $@ $< $(CXXFLAGS) $(LIBS)

$(ODIR)/get_state.o: $(SRCDIR)/get_state.cpp $(IDIR)/get_state.hpp  $(IDIR)/approx.hpp
;$(CC) -c -o $@ $< $(CXXFLAGS) $(LIBS)

$(ODIR)/approx.o: $(SRCDIR)/approx.cpp $(IDIR)/approx.hpp
;$(CC) -c -o $@ $< $(CXXFLAGS) $(LIBS)

test_log.txt: $(OUTDIR)/test_minimize.exe
;$(OUTDIR)/test_minimize.exe -v > test_log.txt

$(OUTDIR)/test_minimize.exe: $(ODIR)/test_minimize.o $(ODIR)/minimize.o
;$(CC) -o $@ $^ $(CXXFLAGS) $(LIBS)


$(ODIR)/test_minimize.o: $(SRCDIR)/test_minimize.cpp $(IDIR)/minimize.hpp
;$(CC) -c -o $@ $< $(CXXFLAGS) $(LIBS)
.PHONY: clean

clean:
;rm -f $(ODIR)/*.o
;rm -f $(OUTDIR)/*.exe
;rm -f log.txt
;rm -f q*.tsv
;rm -f  *.pdf
