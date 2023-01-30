test:
	pytest3

develop:
	python3 setup.py develop --user

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

TART_ARGS=--ms test_data/test.ms --nvis 8000 --healpix --fov 160deg --res 120arcmin
#TART_ARGS=--file test_data/test_data.json --healpix --fov 160deg --res 30arcmin
tart:
	spotless  ${TART_ARGS} --SVG --title tart

ms:
	spotless ${TART_ARGS} --ms test_data/test.ms --multimodel --HDF ms.hdf --SVG --title ms

disko:
	disko ${TART_ARGS} --ms test_data/test.ms --tikhonov --alpha 2 --SVG  --HDF disko.hdf --title disko

draw:
	disko_draw ms.hdf --show-sources --SVG ms.svg

big_ms:
	tart2ms --hdf ~/Downloads/tart_obs/tart_obs/vis_2022-10-21_10_19_46.661434.hdf \
		--rephase obs-midpoint --single-field \
		--ms test_data/test.ms --clobber
