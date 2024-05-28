# SSMF

The Slave-Spin Mean-Field (SSMF) approach is a methodology based on solving the Hubbard-Kanamori model based on tight-bindings to obtain the orbital **quasiparticle weights** and **occupations**.

It is based on the Brinkmann-Rice picture, so that it can trace the evolution of these quantities in terms of U, J<sub>H</sub>, T, total nominal filling; but once the Mott phase is reached (i.e., all quasiparticle weights are 0), then there is not more information. Thus it is a good approach to study how strongly correlated are certain metals: weakly correlated if the quasiparticle weights are 0.7 - 1.0, strongly correlated if the quasiparticle weights are lower than 0.4 - 0.5. It can also reproduce the **orbital-selective Mott transitions**, meaning some orbital quasiparticles weights are zero, while others are finite.



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
