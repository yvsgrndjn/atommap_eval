.PHONY: test lint install clean dist clean-dev

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


clean-dev:
	@echo "🧹 Cleaning build artifacts, caches, and bytecode..."
	rm -rf __pycache__ */__pycache__ */*/__pycache__
	rm -rf *.egg-info */*.egg-info */*/*.egg-info
	rm -rf .pytest_cache .mypy_cache .ruff_cache .coverage
	rm -rf dist build
	@echo "✅ Clean complete."