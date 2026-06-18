PDF = spotless.pdf

all: $(PDF)

$(PDF): spotless.tex
	pdflatex spotless.tex
	pdflatex spotless.tex

clean:
	rm -f *.aux *.log *.out *.bbl *.blg *.pdf
