import gaupy.log
import gaupy.molecules
from gaupy.patterns import h, c, n, p
import molmod.graphs as gr
import molmod.molecular_graphs as g
from collections import OrderedDict
import warnings
import logging
from gaupy.utils import cached
import numpy as np

# Dirty trick to silence warnings of confliciting
# modules on the HPC cluster and to silence NumPy
# warnings.
warnings.filterwarnings('ignore')

__all__ = ['Phosphite']


class Phosphite(gaupy.log.LOGFile):

    info = '''
Phosphite additions to mono- and diunsatured systems.
Atoms:      N, C2/3/4/5/6, P2/4/6 (closest to respective C)
Distances:  PC2, PC4, PC6
Methods:    charges()
    '''

    def __init__(self, filename):
        logging.debug('Phosphite.__init__(): %s' % filename)
        super(Phosphite, self).__init__(filename)
        try:
            self._classify()
        except:
            logging.error('Failed to classify %s' % self.file)

    @classmethod
    def partition(cls, logs, add_patterns=[], patterns=[],
                  group_by_number=False):
        return super(Phosphite, cls).partition(
            logs, add_patterns=['ntms', 'otms', 'prc', 'adduct', 'opt',
                                'adductirc', 'prcirc', 'noirc', 'esma',
                                'iprperpendicular', 'iprflat'])

    @cached
    def pc2(self):
        if self.geometry.p2:
            d = self.geometry.dist(self.geometry.c2.n, self.geometry.p2.n)
            if d < 3.0:
                return d

    @cached
    def pc4(self):
        if self.geometry.p4:
            d = self.geometry.dist(self.geometry.c4.n, self.geometry.p4.n)
            if d < 3.0: 
                return d

    @cached
    def pc6(self):
        if self.geometry.p6:
            d = self.geometry.dist(self.geometry.c6.n, self.geometry.p6.n)
            if d < 3.0: 
                return d

    def _classify(self):
        od = OrderedDict()
        self.geometry.initiate_match()
        od['n'] = n
        od['c2'] = gr.CritOr(g.HasNeighbors(n, c, h),
                            g.HasNeighbors(n, c, h, p))
        self.geometry.set_matches(od)
        self.geometry.set_match('c3', self.geometry.closest(
            6, self.geometry.c2.n, only=self.geometry.unparsed))
        self.geometry.set_match('c4', self.geometry.closest(
            6, self.geometry.c3.n, only=self.geometry.unparsed))
        if 'di' in self.file:
            self.geometry.set_match('c5', self.geometry.closest(
                6, self.geometry.c4.n, only=self.geometry.unparsed))
            self.geometry.set_match('c6', self.geometry.closest(
                6, self.geometry.c5.n, only=self.geometry.unparsed))
        # classify P atoms as p2, p4, p6 based on whether they're
        # closer to C2, C4, C6 respectively
        for pat in self.geometry.atoms(15):
            reactives = [self.geometry.c2.n, self.geometry.c4.n]
            if 'di' in self.file:
                reactives.append(self.geometry.c6.n)
            closest_reactive = self.geometry.closest(6, pat, only=reactives)
            # dirty trick to convert a list index to p2/4/6
            name = 'p%i' % ((reactives.index(closest_reactive) + 1) * 2)
            self.geometry.set_match(name, pat)

    def charges(self, hi=True, npa=True, npa_suffix=None, all_carbons=True,
                string=False):
        if all_carbons:
            atoms = (['c2', 'c3', 'c4'] if 'mono' in self.file
                     else ['c2', 'c3', 'c4', 'c5', 'c6'])
        else:
            atoms = (['c2', 'c4'] if 'mono' in self.file
                     else ['c2', 'c4', 'c6'])
        indices = [getattr(self.geometry, x).n for x in atoms]
        try:
            if hi:
                hi = [self.hi_charges[i] for i in indices]
        except:
            hi = []
            logging.debug('No HI charges found.')
        try:
            if npa:
                npa_log = self
                if npa_suffix:
                    npa_log = gaupy.log.LOGFile(self.file.add(npa_suffix))
                npa = [npa_log.npa_charges[i] for i in indices]
        except:
            npa = []
            logging.debug('No NPA charges found.')
        if string:
            return '%-25s' % self.file + (
                '%-8.3f' * (len(npa) + len(hi))) % tuple(npa + hi)
        else:
            return npa + hi
