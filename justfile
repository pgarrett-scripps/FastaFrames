# Install the package with dev dependencies
install:
    pip install -e ".[dev]"

# Run all checks (lint, typecheck, test)
check: lint typecheck test

# Lint with ruff
lint:
    ruff check src tests

# Auto-fix lint issues
fix:
    ruff check --fix src tests

# Format code
fmt:
    ruff format src tests

# Check formatting without modifying
fmt-check:
    ruff format --check src tests

# Type check with ty
typecheck:
    ty check src

# Run tests
test *ARGS:
    python -m pytest tests {{ARGS}}

# Run tests with verbose output
test-v:
    python -m pytest tests -v

# Build the package
build:
    python -m build

# Clean build artifacts
clean:
    rm -rf build dist src/*.egg-info .ruff_cache .pytest_cache
