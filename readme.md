# Parallel Least Squares

Optimizes multiple inputs in parallel using the Levenberg-Marquadt algorithm.

# Requirements

- CMake >= 3.12
- C++17-capable C++ compiler 

# Building

You may have to specify directories of dependencies, like so:

# Installation

You may install the python bindings like so (preferably using a virtual environment):
```console
$ python3 -m venv venv
$ source venv/bin/activate
$ python setup.py build
$ python setup.py install
$ python
>>> import levmar_eigen
>>> 
```
