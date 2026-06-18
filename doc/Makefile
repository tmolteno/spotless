PDF = spotless.pdf
SRC = spotless.tex
BIB = spotless.bib

all: $(PDF)

$(PDF): $(SRC) $(BIB)
	pdflatex $(SRC)
	bibtex spotless
	pdflatex $(SRC)
	pdflatex $(SRC)

clean:
	rm -f *.aux *.log *.out *.bbl *.blg *.pdf *.toc
