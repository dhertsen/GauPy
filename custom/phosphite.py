import gaupy.log
# import gaupy.molecules
# from gaupy.patterns import h, c, o, s, me
# import molmod.molecular_graphs as g
# import molmod.graphs as gr
# from collections import OrderedDict
# import numpy as np
# import math
import warnings
# import logging

# Dirty trick to silence warnings of confliciting
# modules on the HPC cluster and to silence NumPy
# warnings.
warnings.filterwarnings('ignore')

__all__ = ['SixFour']


class Phosphite(gaupy.log.LOGFile):

    @classmethod
    def partition(cls, logs, add_patterns=[], patterns=[],
                  group_by_number=False):
        return super(Phosphite, cls).partition(
            logs, add_patterns=['ntms', 'otms', 'prc', 'adduct', 'opt',
                                'adductirc', 'prcirc', 'noirc'])
