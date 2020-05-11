# Parallel Least Squares

Optimizes multiple inputs in parallel using the Levenberg-Marquadt algorithm.

# Requirements

- CMake >= 3.12
- Eigen3
- python3
- numpy

# Testing

```console
$ mkdir build; cd build
$ cmake ..
$ ctest
```

You may also run `TestDriver` directly. There are preset directories in which the test
driver will look for matrix/vector files, however you may provide additional paths
as arguments.

# Installation

Preferably using a virtual environment:
```console
$ python3 -m venv my_venv
$ source my_venv/bin/activate
$ python -m pip install --user --editable .
```
