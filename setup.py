import sys
import os
import glob

from setuptools import find_packages, setup
from tools.setup_helpers.cmake_setuptools import (
        CMakeBuild,
        CMakeInstall,
        CMakeExtension,
        )
from tools.setup_helpers.extra_setuptools_commands import CleanBuild

libs = [ i for i in glob.glob(os.path.join('levmar_eigen', '*'))
        if not i.endswith('.py') and not '__pycache__' in i ]
print(libs)

extra_cmake_args = dict()

args = []
for i, arg in enumerate(sys.argv):
    if arg.startswith('-D') and '=' in arg:
        k, v = arg.split('=')
        k = k.replace('-D', '')
        extra_cmake_args[k] = v
    else:
        args.append(arg)

sys.argv = args

setup(
    name='levmar-eigen',
    version='0.1.2',
    author='Asher Mancinelli',
    author_email='asher.mancinelli@pnnl.gov',
    description='Dispatcher for multiple parallel runs of an LM least-squares solver',
    long_description='',
    ext_modules=[CMakeExtension('ParLeastSquares', **extra_cmake_args)],
    packages=find_packages(),
    package_data=dict(levmar_eigen=libs),
    include_package_data=True,
    cmdclass=dict(
        build_ext=CMakeBuild,
        install=CMakeInstall,
        clean=CleanBuild,
        ),
    zip_safe=False,
)
