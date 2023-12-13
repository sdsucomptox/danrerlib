.PHONY: clean html deploy

clean:
	@echo "Cleaning generated files"
	cd sphinx && make clean
	rm -rf docs/*

html:
	@echo "Building HTML documentation"
	cd sphinx && make html
	cp -a sphinx/_build/html/. docs/

deploy: clean html
	@echo "Deploying documentation"
	cd docs && touch .nojekyll