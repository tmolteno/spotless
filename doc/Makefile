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
	rm -f gridless.png residual.png spotless.png spotless_direct.png
	rm -f multimodel_gridless.png multimodel_residual.png multimodel_spotless.png multimodel_direct.png

SPOTLESS_ARGS=--ms ../test_data/test.ms --SVG --healpix --fov 170 --res 20arcmin
images:
	uv run spotless ${SPOTLESS_ARGS} --title 'spotless' --log spotless_run.txt
	uv run spotless ${SPOTLESS_ARGS} --multimodel --title 'multimodel' --log multimodel_run.txt

pngs:
	for svg in *.svg; do \
		inkscape --export-type=png --export-width=1024 "$$svg" \
			--export-filename="$${svg%.svg}.png"; \
	done
	$(MAKE) links

links:
	@for pair in \
		"spotless_gridless_*:gridless" \
		"spotless_residual_*:residual" \
		"spotless_spotless_2*:spotless" \
		"spotless_spotless_direct_*:spotless_direct" \
		"multimodel_gridless_*:multimodel_gridless" \
		"multimodel_residual_*:multimodel_residual" \
		"multimodel_spotless_2*:multimodel_spotless" \
		"multimodel_spotless_direct_*:multimodel_direct"; do \
		glob=$${pair%%:*}; dst=$${pair##*:}.png; \
		for src in $$glob.png; do \
			if [ -f "$$src" ]; then \
				ln -sf "$$src" "$$dst"; \
				echo "  $$dst -> $$src"; \
			fi; \
		done; \
	done
