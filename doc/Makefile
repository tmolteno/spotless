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

SPOTLESS_ARGS=--ms ../test_data/test.ms --SVG --healpix --fov 170 --res 20arcmin
images:
	uv run spotless ${SPOTLESS_ARGS} --title 'spotless'
	uv run spotless ${SPOTLESS_ARGS} --multimodel --title 'multimodel'
