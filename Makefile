
sync:
	uv sync

develop:
	uv pip install -e .

test:
	uv run pytest

lint:
	uv run flake8 spotless --count --exit-zero --max-complexity=10 --max-line-length=127 --statistics

clean:
	rm -rf dist/ build/
	rm -rf spotless.egg-info

build:
	uv build

TART_ARGS=--ms test_data/test.ms --healpix --fov 160deg --res 120arcmin --debug
#TART_ARGS=--file test_data/test_data.json --healpix --fov 160deg --res 30arcmin
tart:
	uv run spotless  ${TART_ARGS} --HDF tart.h5 --SVG --title tart

ms:
	uv run spotless ${TART_ARGS} --multimodel --HDF ms.hdf --SVG --title ms

disko:
	uv run disko ${TART_ARGS} --tikhonov --alpha 2 --SVG  --HDF disko.hdf --title disko

draw:
	uv run disko_draw ms.hdf --show-sources --SVG ms.svg

big_ms:
	tart2ms --hdf ~/Downloads/tart_obs/tart_obs/vis_2022-10-21_10_19_46.661434.hdf \
		--rephase obs-midpoint --single-field \
		--ms test_data/test.ms --clobber
