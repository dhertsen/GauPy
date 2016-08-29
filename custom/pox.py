import gaupy.log
import gaupy.molecules
import gaupy.utils as utils
import warnings
# import numpy as np
# from molmod.units import angstrom
from gaupy.patterns import h, c, o, s, me, ch2
from collections import OrderedDict
import logging


# Dirty trick to silence warnings of confliciting
# modules on the HPC cluster and to silence NumPy
# warnings.
warnings.filterwarnings('ignore')

__all__ = ['Pox']


class Pox(gaupy.log.LOGFile):

    def __init__(self, filename):
        super(Pox, self).__init__(filename)

        self.residues = self.stoichiometry['N']
        '''
        Number of residues in the chain. Equals the number of nitrogen
        atoms in the molecule, since no side chains with nitrogen atoms
        are considered.
        '''

        self._classify()


    @utils.cached
    def _hi(self):
        return self.get_hi_charges()

    @utils.cached
    def n(self):
        self.n = self._hi[self.geometry.monomer_nitrogen.n]

    @utils.cached
    def c5(self):
        self.c5 = self._hi[self.geometry.cation_c5.c5]
        
    @utils.cached
    def electro(self):
        if self.cation_c5:
            return self._hi[self.cation_c5]
        else:
            return None

    @utils.cached
    def nucleo(self):
        if self.monomer_n:
            return self._hi[self.monomer_n]
        else:
            return None

    def _classify(self):

        logging.debug('start _classify()')
        self.geometry.initiate_match()
        # self.geometry.set_match('four_c1', c2_bonded[0])

        # find all oxazoline rings
        self.rings = self.geometry.nrings(5)
        for r in self.rings:
            # Classify ring atoms, works both for cationic and neutral rings
            nitrogen = filter(lambda i: self.geometry.numbers[i] == 7, r)[0]
            oxygen = filter(lambda i: self.geometry.numbers[i] == 8, r)[0]
            cs_close_to_n = set(
                self.geometry.closest(6, nitrogen, n=2, only=r))
            cs_close_to_o = set(self.geometry.closest(6, oxygen, n=2, only=r))
            c2 = (cs_close_to_n & cs_close_to_o).pop()
            c4 = (cs_close_to_n - set([c2])).pop()
            c5 = (cs_close_to_o - set([c2])).pop()
            # Find out what type of ring it is
            # distance of N to third nearest C
            dist = self.geometry.dist(
                nitrogen, self.geometry.closest(6, nitrogen, n=3)[2])
            # If this distance is smaller than 1.75 A, it is a real bond.
            if dist <= 1.75:
                species = 'cation'
            else:
                species = 'monomer'
            # Add 'cation_' or 'monomer_'.
            for atom in ['nitrogen', 'oxygen', 'c2', 'c4', 'c5']:
                name = '%s_%s' % (species, atom)
                self.geometry.set_match(name, eval(atom))
            try:
                self.geometry.set_match('butylc1', self.geometry.closest(6, c2))
                attached_to_c1 = self.geometry.closest(6, self.geometry.butylc1.n, n=3)
                self.geometry.set_match('butylc2', ch2(6, 6), only=attached_to_c1)
                self.geometry.set_match('butylc3', self.geometry.closest(6, self.geometry.butylc2.n))
                self.geometry.set_match('butylc4', self.geometry.closest(6, self.geometry.butylc3.n))
            except:
                logging.debug('No butyl chain.')

        logging.debug('end _classify()')
