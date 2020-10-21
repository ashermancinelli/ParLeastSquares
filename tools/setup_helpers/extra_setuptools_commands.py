
from distutils.command.sdist import sdist
from setuptools import Command
from .utils import *

class SDistChecked(sdist):
    """ check submodules on sdist to prevent incomplete tarballs """
    def run(self):
        init_submodules()
        sdist.run(self)


class CleanBuild(Command):
    '''Clean build directory'''

    user_options = []

    def initialize_options(self):
        ...

    def finalize_options(self):
        ...

    def run(self):
        build_dir = get_build_dir()
        rmtree(build_dir)


