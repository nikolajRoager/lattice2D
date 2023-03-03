.RECIPEPREFIX := ;

CSC = mcs
RUN = mono

all : outA.txt
;cat outA.txt
;gnuplot plot.gpi

outA.txt rosenbrock.tsv himmelblau.tsv : partA.exe
;$(RUN) $< > outA.txt


partA.exe:main.cs minimizer.dll matrix.dll vector.dll
;$(CSC) -target:exe -out:$@ -reference:vector.dll -reference:matrix.dll  -reference:minimizer.dll $<

minimizer.dll : minimizer.cs matrix.dll vector.dll
;$(CSC) $< -target:library -reference:matrix.dll -reference:vector.dll  -out:$@

matrix.dll : matrix.cs vector.dll
;$(CSC) $< -target:library -reference:vector.dll  -out:$@

vector.dll : vector.cs
;$(CSC) $< -target:library  -out:$@

.PHONY: clean

clean:
;rm -f *.dll *.exe log* *.png out*
