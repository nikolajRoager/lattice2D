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
default:  log.txt


log.txt: $(OUTDIR)/test.exe
;$(OUTDIR)/test.exe > log.txt

$(OUTDIR)/test.exe: $(ODIR)/main.o $(ODIR)/minimize.o
;$(CC) -o $@ $^ $(CXXFLAGS) $(LIBS)

$(ODIR)/main.o: $(SRCDIR)/main.cpp $(IDIR)/minimize.hpp
;$(CC) -c -o $@ $< $(CXXFLAGS) $(LIBS)


$(ODIR)/minimize.o: $(SRCDIR)/minimize.cpp $(IDIR)/minimize.hpp
;$(CC) -c -o $@ $< $(CXXFLAGS) $(LIBS)

.PHONY: clean

clean:
;rm -f $(ODIR)/*.o
;rm -f $(OUTDIR)/*.exe
;rm -f log*
;rm -f *.tsv
;rm -f *.pdf
