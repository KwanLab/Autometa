
hello:
	@echo "Please inspect Makefile for list of commands"


test-wip: tests/test_data.json
	pytest -m "wip" --durations=0 --cov=autometa --emoji --cov-report html

test: tests/test_data.json
	pytest --durations=0 --cov=autometa --emoji --cov-report html

.Phony:
	test test-wip hello
