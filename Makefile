.PHONY: clean clean-build clean-pyc clean-test test help docs deps build version-check tag release
.DEFAULT_GOAL := help

define BROWSER_PYSCRIPT
import os, webbrowser, sys

from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

BROWSER := python -c "$$BROWSER_PYSCRIPT"

VERSION := $(shell poetry version -s)
TAG := v$(VERSION)


clean: clean-build clean-pyc clean-test ## remove all build, test, coverage and Python artifacts

clean-build: ## remove build artifacts
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-pyc: ## remove Python file artifacts
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test: ## remove test and coverage artifacts
# 	rm -fr .tox/
	rm -f .coverage
	rm -fr htmlcov/
	rm -fr .pytest_cache

test: ## run tests quickly with the default Python
	poetry run pytest --doctest-modules src/scanpex/ tests/

help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

docs: ## generate Sphinx HTML documentation, including API docs
	rm -f docs/modules.rst
	poetry run sphinx-apidoc -o docs/ src/scanpex
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
	$(BROWSER) docs/_build/html/index.html

deps: ## export dependencies
	poetry export --with dev -f requirements.txt -o ./docs/requirements.txt

synclib: ## initiate automated sync dependencies
	python synclib.py

build: clean ## build the package
	poetry build

version-check: build ## build check with twine
	poetry run twine check dist/*

tag: ## create tags
	@echo "Current version is $(VERSION)"
	@echo "Creating git tag: $(TAG)"
	git tag $(TAG)

release: clean-test test tag ## release current version
	@echo "Pushing tag $(TAG) to origin..."
	git push origin $(TAG)
