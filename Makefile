.PHONY: \
	clean clean-build clean-pyc clean-test \
	test help docs build version-check \
	tag release publish \
	patch minor major _bump-release


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
	rm -f .coverage
	rm -fr htmlcov/
	rm -fr .pytest_cache


test: ## run tests quickly with the default Python
	poetry run pytest --doctest-modules src/scanpex/ tests/


help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)


docs: ## generate Sphinx HTML documentation, including API docs
	rm -f docs/modules.rst docs/scanpex*.rst
	poetry run sphinx-apidoc \
		--implicit-namespaces \
		--module-first \
		-o docs/ src/scanpex
	poetry run $(MAKE) -C docs clean
	poetry run $(MAKE) -C docs html
	$(BROWSER) docs/_build/html/index.html


synclib: ## initiate automated sync dependencies
	python synclib.py


build: clean ## build the package
	poetry build


version-check: build ## build and validate distribution artifacts
	poetry run twine check dist/*


patch: ## bump patch version, validate, commit, tag, and push
	@$(MAKE) _bump-release LEVEL=patch


minor: ## bump minor version, validate, commit, tag, and push
	@$(MAKE) _bump-release LEVEL=minor


major: ## bump major version, validate, commit, tag, and push
	@$(MAKE) _bump-release LEVEL=major


_bump-release:
	@set -e; \
	\
	if [ -n "$$(git status --porcelain)" ]; then \
		echo "ERROR: working tree is not clean."; \
		echo "Commit or stash changes before releasing."; \
		exit 1; \
	fi; \
	\
	OLD_VER=$$(poetry version -s); \
	echo "Current version: $$OLD_VER"; \
	\
	poetry version $(LEVEL); \
	NEW_VER=$$(poetry version -s); \
	echo "New version: $$NEW_VER"; \
	\
	perl -i -pe \
		's/"__version":\s*"[^"]*"/"__version": "'$$NEW_VER'"/' \
		cookiecutter.json; \
	\
	echo "Synchronizing lock file..."; \
	poetry lock; \
	\
	echo "Checking Poetry configuration..."; \
	poetry check --lock; \
	\
	echo "Running tests..."; \
	$(MAKE) test; \
	\
	echo "Building and validating distribution..."; \
	$(MAKE) version-check; \
	\
	echo "Committing version bump..."; \
	git add pyproject.toml poetry.lock cookiecutter.json; \
	git commit -m ":wrench: $(LEVEL) $$OLD_VER -> $$NEW_VER"; \
	\
	echo "Creating tag v$$NEW_VER..."; \
	git tag -a v$$NEW_VER -m "v$$NEW_VER"; \
	\
	echo "Pushing main and tag..."; \
	git push origin main; \
	git push origin v$$NEW_VER; \
	\
	echo "Release preparation completed: v$$NEW_VER"


publish: ## publish the current built distribution to PyPI
	@echo "Publishing version $$(poetry version -s) to PyPI..."
	poetry publish


release: publish ## alias for publishing the validated distribution
