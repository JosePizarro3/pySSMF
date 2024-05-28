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

Make sure to have pip upgraded to its latest version:
```sh
pip install --upgrade pip
```

Now, install the dependencies in editable mode (this is done by adding the flag `-e` and it is useful when debugging the code):
```sh
pip install -e '.[dev]'
```

The `[dev]` option will allow to install the optional-dependencies specified in the configuration `pyproject.toml` file.

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
