test:
	python3 setup.py test

develop:
	pip3 install -e .

lint:
	pylint --extension-pkg-whitelist=numpy --ignored-modules=numpy --extension-pkg-whitelist=astropy spotless

	
test_upload:
	python3 setup.py sdist
	twine upload --repository testpypi dist/*

upload:
	python3 setup.py sdist
	twine upload --repository pypi dist/*

tart:
	spotless --file doc/publication_data.json --display --PNG
