from setuptools_cmake import *
from tools.setup_helpers.cmake_setuptools import CMakeBuild, CMakeInstall

setup(
    name='levmar-eigen',
    version='0.1',
    author='Asher Mancinelli',
    author_email='asher.mancinelli@pnnl.gov',
    description='Dispatcher for multiple parallel runs of an LM least-squares solver',
    long_description='',
    ext_modules=[CMakeExtension('ParLeastSquares')],
    cmdclass=dict(
        build_ext=CMakeBuild,
        install=CMakeInstall,
        ),
    zip_safe=False,
)
