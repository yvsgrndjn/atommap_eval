.PHONY: test lint install clean dist

install:
	pip install -e ".[dev]"

test:
	pytest -v

lint:
	ruff check src tests

format:
	black src tests
	ruff check src tests --fix

clean:
	rm -rf build dist *.egg-info .pytest_cache .mypy_cache __pycache__

dist:
	python -m build

