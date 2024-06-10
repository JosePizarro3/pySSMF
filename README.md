![](https://coveralls.io/repos/github/JosePizarro3/pySSMF/badge.svg?branch=develop)

# pySSMF

## Installation steps

In order to install `pySSMF`, you will need a machine with Python3.9 and create a virtual environment.

First, open a terminal, clone, and enter in the folder of this project:
```sh
git clone https://github.com/JosePizarro3/pySSMF.git
cd pySSMF
```

Create the virtual environment specifying the version of Python3.9:
```sh
python3.9 -m venv .pyenv
source .pyenv/bin/activate
```

We recommend using `uv` for pip installing all packages:
```sh
pip install uv
```

Make sure to have pip upgraded to its latest version:
```sh
uv pip install --upgrade pip
```

Now, install the dependencies in editable mode (this is done by adding the flag `-e` and it is useful when debugging the code):
```sh
uv pip install -e '.[dev]'
```

The `[dev]` option will allow to install the optional-dependencies specified in the configuration `pyproject.toml` file.

### Run the tests

You can run local tests using the `pytest` package:

```sh
python -m pytest -sv tests
```

where the `-s` and `-v` options toggle the output verbosity.

Our CI/CD pipeline produces a more comprehensive test report using `coverage` and `coveralls` packages. We suggest you to generate your own coverage reports locally by doing:

```sh
uv pip install coverage coveralls
python -m pytest --cov=src tests
```

### Run linting and auto-formatting

We use [Ruff](https://docs.astral.sh/ruff/) for auto-formatting our Python modules. This package is included as part of the `[dev]` dependencies of the project, and can be run in the terminal:

```sh
ruff check .
ruff format . --check
```

Ruff auto-formatting is also a part of the GitHub workflow actions. Make sure that before you make a Pull Request, `ruff format . --check` runs in your local without any errors otherwise the workflow action will fail.

### Using VSCode for debugging

I recommend to use VSCode as your visual editor for debugging. The project comes with a `settings.json` file which already takes care of formatting everytime you save changes in a file (using [Ruff](https://docs.astral.sh/ruff/)).

You can also create your own personal debugger `launch.json`. As an example, this is my current debug file:
```json
{
    "version": "0.2.0",
    "configurations": [
        {
            "name": "pySSMF python",
            "type": "debugpy",
            "request": "launch",
            "program": "${workspaceFolder}/debug.py",
            "python": "${workspaceFolder}/.pyenv/bin/python",
            "console": "integratedTerminal",
            "justMyCode": true,
            "env": {
                "PYTHONDONTWRITEBYTECODE": "1"
            }
        },
        {
            "name": "pySSMF pytest",
            "type": "debugpy",
            "request": "launch",
            "module": "pytest",
            "args": [
                "-sv",
                "${workspaceFolder}/tests"
            ],
            "python": "${workspaceFolder}/.pyenv/bin/python",
            "console": "integratedTerminal",
            "justMyCode": true,
            "env": {
                "PYTHONDONTWRITEBYTECODE": "1"
            }
        },
    ]
}
```

Feel free to modify the relative paths under `"program"` or `"args"` respectively.
