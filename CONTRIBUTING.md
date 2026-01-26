# Contributing to pyLocusZoom

Thank you for your interest in contributing to pyLocusZoom!

## Development Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/michael-denyer/pyLocusZoom.git
   cd pyLocusZoom
   ```

2. Install dependencies with uv:
   ```bash
   uv sync --all-extras
   ```

3. Run tests:
   ```bash
   uv run python -m pytest tests/ -v
   ```

4. Run linting:
   ```bash
   uv run ruff check src/
   uv run ruff format --check src/ tests/
   ```

## Code Style

- Follow [PEP 8](https://peps.python.org/pep-0008/) guidelines
- Use [ruff](https://github.com/astral-sh/ruff) for linting and formatting
- Maximum line length: 88 characters
- Use Google-style docstrings

### Docstring Example

```python
def function_name(param1: str, param2: int = 10) -> bool:
    """Short one-line description.

    Longer description if needed.

    Args:
        param1: Description of first parameter.
        param2: Description of second parameter.

    Returns:
        Description of return value.

    Raises:
        ValueError: When param1 is empty.
    """
```

## Pull Request Process

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/your-feature`
3. Make your changes
4. Run tests and linting:
   ```bash
   uv run python -m pytest tests/ -v
   uv run ruff check src/
   uv run ruff format src/ tests/
   ```
5. Commit with a descriptive message
6. Push and create a pull request

## Testing

- Write tests for new functionality
- Use pytest fixtures from `tests/conftest.py`
- Mock external dependencies (PLINK, network calls)
- Aim for test coverage of new code

### Running Tests

```bash
# All tests
uv run python -m pytest tests/ -v

# Specific test file
uv run python -m pytest tests/test_plotter.py -v

# With coverage
uv run python -m pytest tests/ --cov=pylocuszoom --cov-report=html
```

## Architecture

See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for project structure.

### Key Modules

| Module | Purpose |
|--------|---------|
| `plotter.py` | Main LocusZoomPlotter class |
| `backends/` | Rendering backends (matplotlib, plotly, bokeh) |
| `ld.py` | PLINK LD calculation |
| `gene_track.py` | Gene/exon visualization |
| `recombination.py` | Recombination map handling |
| `eqtl.py` | eQTL data support |

## Reporting Issues

- Use GitHub Issues
- Include Python version, OS, and package versions
- Provide a minimal reproducible example
- Include full error traceback

## License

By contributing, you agree that your contributions will be licensed under the GPL-3.0-or-later license.
