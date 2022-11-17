test:
	python3 setup.py test

develop:
	pip3 install -e .

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

tart:
	spotless --file test_data/test_data.json --multimodel --healpix --fov 160deg --res 30arcmin --SVG

ms:
	spotless --file test_data/test_data.json --healpix --fov 160deg --res 60arcmin --HDF mytest.hdf

draw:
	disko_draw mytest.hdf --show-sources --SVG mytest.svg
