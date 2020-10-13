
hello:
	@echo "Please inspect Makefile for list of commands"

docs:
	make clean html -C docs && open docs/build/html/index.html

clean:
	@echo "Removing everything under 'htmlcov'..."
	@rm -rf htmlcov && make clean -C docs

test: tests/data/test_data.json
	python -m pytest --durations=0 --cov=autometa --emoji --cov-report html

test-wip: tests/data/test_data.json
	python -m pytest -m "wip" --durations=0 --cov=autometa --emoji --cov-report html

test-entrypoints: tests/data/test_data.json
	python -m pytest -m "entrypoint" --durations=0 --cov=autometa --emoji --cov-report html

.PHONY: hello docs clean test test-wip test-entrypoints
