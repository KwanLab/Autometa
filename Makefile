.PHONY: clean black create_environment install image docs clean unit_test unit_test_data unit_test_data_download unit_test_data_build unit_test_wip unit_test_entrypoints

#################################################################################
# GLOBALS                                                                       #
#################################################################################

PROJECT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
PROJECT_NAME = autometa
PYTHON_INTERPRETER = python3
# This was retrieved from https://drive.google.com/file/d/1bSlPldaq3C6Cf9Y5Rm7iwtUDcjxAaeEk/view?usp=sharing
TEST_DATA_FILEID = 1bSlPldaq3C6Cf9Y5Rm7iwtUDcjxAaeEk

ifeq (,$(shell which conda))
HAS_CONDA=False
else
HAS_CONDA=True
endif

#################################################################################
# COMMANDS                                                                      #
#################################################################################

## Delete all compiled Python files
clean:
	find . -type f -name "*.py[co]" -exec rm -r {} +
	find . -type d -name "__pycache__" -exec rm -r {} +
	find . -type d -name "htmlcov" -exec rm -r {} +
	find . -type d -name "Autometa.egg-info" -exec rm -r {} +
	find . -type d -name "dist" -exec rm -r {} +
	find . -type d -name "build" -exec rm -r {} +

## Apply black formatting
black:
	black --exclude autometa/validation autometa

## Set up python interpreter environment
create_environment: requirements.txt
ifeq (True,$(HAS_CONDA))
		@echo ">>> Detected conda, creating conda environment."
ifeq (3,$(findstring 3,$(PYTHON_INTERPRETER)))
	conda create -c conda-forge -c bioconda --name $(PROJECT_NAME) python=3 --file=requirements.txt
else
	conda create -c conda-forge -c bioconda --name $(PROJECT_NAME) python=2.7 --file=requirements.txt
endif
	@echo ">>> New conda env created. Activate with:\nsource activate $(PROJECT_NAME)"
else
	$(PYTHON_INTERPRETER) -m pip install -q virtualenv virtualenvwrapper
	@echo ">>> Installing virtualenvwrapper if not already installed.\nMake sure the following lines are in shell startup file\n\
	export WORKON_HOME=$$HOME/.virtualenvs\nexport PROJECT_HOME=$$HOME/Devel\nsource /usr/local/bin/virtualenvwrapper.sh\n"
	@bash -c "source `which virtualenvwrapper.sh`;mkvirtualenv $(PROJECT_NAME) --python=$(PYTHON_INTERPRETER)"
	@echo ">>> New virtualenv created. Activate with:\nworkon $(PROJECT_NAME)"
endif

#################################################################################
# PROJECT RULES                                                                 #
#################################################################################

## Install autometa from source
install:
	python setup.py install

## Install dependencies for test environment
test_environment: tests/requirements.txt
	python -m pip install --requirement=tests/requirements.txt

## Build docker image from Dockerfile (auto-taggged as jason-c-kwan/autometa:<current-branch>)
image: Dockerfile
	docker build . -t jason-c-kwan/autometa:`git branch --show-current`

## Build documentation for autometa.readthedocs.io
docs:
	make clean html -C docs
	@echo "docs built. Open docs/build/html/index.html to view"

## Download test_data.json for unit testing
unit_test_data_download:
	gdown --id $(TEST_DATA_FILEID) -O tests/data/test_data.json

## Build test_data.json file for unit testing (requires all files from https://drive.google.com/open?id=189C6do0Xw-X813gspsafR9r8m-YfbhTS be downloaded into tests/data/)
unit_test_data_build: tests/data/records.fna
	python make_test_data.py

## Run all unit tests
unit_test: tests/data/test_data.json test_environment
	python -m pytest --durations=0 --cov=autometa --emoji --cov-report=html tests

## Run unit tests marked with WIP
unit_test_wip: tests/data/test_data.json test_environment
	python -m pytest -m "wip" --durations=0 --cov=autometa --emoji --cov-report=html tests

## Run unit tests marked with entrypoint
unit_test_entrypoints: tests/data/test_data.json test_environment
	python -m pytest -m "entrypoint" --durations=0 --cov=autometa --emoji --cov-report=html tests


#################################################################################
# Self Documenting Commands                                                     #
#################################################################################

.DEFAULT_GOAL := help

# Inspired by <http://marmelab.com/blog/2016/02/29/auto-documented-makefile.html>
# sed script explained:
# /^##/:
# 	* save line in hold space
# 	* purge line
# 	* Loop:
# 		* append newline + line to hold space
# 		* go to next line
# 		* if line starts with doc comment, strip comment character off and loop
# 	* remove target prerequisites
# 	* append hold space (+ newline) to line
# 	* replace newline plus comments by `---`
# 	* print line
# Separate expressions are necessary because labels cannot be delimited by
# semicolon; see <http://stackoverflow.com/a/11799865/1968>
.PHONY: help
help:
	@echo "$$(tput bold)Available rules:$$(tput sgr0)"
	@echo
	@sed -n -e "/^## / { \
		h; \
		s/.*//; \
		:doc" \
		-e "H; \
		n; \
		s/^## //; \
		t doc" \
		-e "s/:.*//; \
		G; \
		s/\\n## /---/; \
		s/\\n/ /g; \
		p; \
	}" ${MAKEFILE_LIST} \
	| LC_ALL='C' sort --ignore-case \
	| awk -F '---' \
		-v ncol=$$(tput cols) \
		-v indent=19 \
		-v col_on="$$(tput setaf 6)" \
		-v col_off="$$(tput sgr0)" \
	'{ \
		printf "%s%*s%s ", col_on, -indent, $$1, col_off; \
		n = split($$2, words, " "); \
		line_length = ncol - indent; \
		for (i = 1; i <= n; i++) { \
			line_length -= length(words[i]) + 1; \
			if (line_length <= 0) { \
				line_length = ncol - indent - length(words[i]) - 1; \
				printf "\n%*s ", -indent, " "; \
			} \
			printf "%s ", words[i]; \
		} \
		printf "\n"; \
	}' \
	| more $(shell test $(shell uname) = Darwin && echo '--no-init --raw-control-chars')
