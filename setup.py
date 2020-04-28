from scripts import *

setup(
    name='ParLeastSquares',
    version='0.1',
    author='Asher Mancinelli',
    author_email='asher.mancinelli@pnnl.gov',
    description='Dispatcher for multiple parallel runs of an LM least-squares solver',
    long_description='',
    ext_modules=[CMakeExtension('ParLeastSquares')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
)
