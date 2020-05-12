# Parallel Least Squares

Optimizes multiple inputs in parallel using the Levenberg-Marquadt algorithm.

# Requirements

- CMake >= 3.12
- Eigen3
- python3
- numpy

# Building

You may have to specify directories of dependencies, like so:

```console
$ mkdir build; cd build
$ cmake .. \
$   -DEigen3_DIR=/some/path
$ ctest
```

This will also output

# Testing

```console
$ mkdir build; cd build
$ cmake ..
$ make
$ ctest
```

Or you may optionally specify paths to folders with data for the tests:
```console
$ cmake .. -DDATA_PATH=/some/path
$ make
$ ctest
```

# Installation

You may install the python bindings like so (preferably using a virtual environment):
```console
$ python3 -m venv my_venv
$ source my_venv/bin/activate
$ python -m pip install --editable .
```
