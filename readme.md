# Parallel Least Squares

Optimizes multiple inputs in parallel using the Levenberg-Marquadt algorithm.

# Requirements

- CMake >= 3.12
- Eigen3
- python3
- numpy

# Installation

Preferably using a virtual environment:
```console
$ python3 -m venv my_venv
$ source my_venv/bin/activate
$ python -m pip install --user --editable .
```

Upon installation of requirements, you should now be able to run tests under 
the `tests` directory. Note that you may have to add the project root to your
linker path like so:
```console
$ LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/path/to/project" python -c 'import pstep'
```

This is because an internal shared library must be installed on your system's path
to be found by python. Actually installing the module (not with --editable) will take
care of this.
