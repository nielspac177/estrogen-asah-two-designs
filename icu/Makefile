# Reproducibility entrypoint. All commands run inside the uv-managed env.
# UV_LINK_MODE=copy is safe everywhere and required on exFAT/network volumes.
export UV_LINK_MODE = copy

UV := uv run

.PHONY: setup test lint fmt all clean

setup:              ## create/refresh the environment
	uv sync --extra dev

test:               ## run unit + synthetic-fixture tests (no PhysioNet data)
	$(UV) pytest

lint:               ## static checks
	$(UV) ruff check .

fmt:                ## autoformat
	$(UV) ruff check --fix .
	$(UV) ruff format .

all:                ## run the full pipeline in order (synthetic unless config/paths.yaml is set)
	@for f in $$(ls scripts/run_*.py 2>/dev/null | sort); do \
		echo "=== $$f ==="; $(UV) python "$$f" || exit 1; \
	done

clean:              ## remove derived outputs
	rm -rf outputs/* .pytest_cache .ruff_cache
