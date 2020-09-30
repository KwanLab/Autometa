
hello:
	@echo "Please inspect Makefile for list of commands"

docs:
	make clean html -C docs && open docs/build/html/index.html

clean:
	@echo "Removing everything under 'htmlcov'..."
	@rm -rf htmlcov && make clean -C docs


test-wip: tests/data/test_data.json
	python -m pytest -m "wip" --durations=0 --cov=autometa --emoji --cov-report html

test: tests/data/test_data.json
	python -m pytest --durations=0 --cov=autometa --emoji --cov-report html

.PHONY: test test-wip hello clean develop docs
