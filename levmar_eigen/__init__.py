import sys, os, shutil, glob
import os.path as op

cwd = op.dirname(op.abspath(__file__))

sobjs = glob.glob(op.join(cwd, '..', '*.so*')) + \
    glob.glob(op.join(cwd, '..', '*.dylib'))

if len(sobjs) > 0:
    for sobj in sobjs:
        shutil.copy(sobj, cwd)

sys.path.insert(0, op.dirname(op.abspath(__file__)))

from __levmar_eigen_C import *
