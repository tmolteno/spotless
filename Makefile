VENVDIR=~/.tartvenv

develop: venv
	${VENV}/pip3 install -e .

test: venv
	${VENV}/python3 -m pytest

lint:
	flake8 spotless --count --exit-zero --max-complexity=10 --max-line-length=127 --statistics

clean:
	rm -rf dist/ build/
	rm -rf spotless.egg-info
	
test_upload:
	python3 setup.py sdist
	twine upload --repository testpypi dist/*

upload:
	python3 setup.py sdist
	twine upload --repository pypi dist/*

TART_ARGS=--ms test_data/tart.ms --healpix --fov 180deg --res 60arcmin
#TART_ARGS=--file test_data/test_data.json --healpix --fov 160deg --res 30arcmin
tart:
	${VENV}/spotless  ${TART_ARGS} --HDF tart.h5 --SVG --title tart

ms:
	${VENV}/spotless ${TART_ARGS} --ms test_data/test.ms --multimodel --HDF ms.hdf --SVG --title ms

disko:
	${VENV}/disko ${TART_ARGS} --ms test_data/test.ms --tikhonov --alpha 2 --SVG  --HDF disko.hdf --title disko

draw:
	${VENV}/disko_draw ms.hdf --show-sources --SVG ms.svg
	
include Makefile.venv
