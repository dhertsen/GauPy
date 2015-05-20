import gaupy.log
import gaupy.molecules
import gaupy.utils as utils
import warnings
import numpy as np
from molmod.units import angstrom

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

        # find all oxazoline rings
        self.rings = self.geometry.get_nrings(5)
        for r in self.rings:
            # Classify ring atoms, works both for cationic and neutral rings
            n = filter(lambda i: self.geometry.numbers[i] == 7, r)[0]
            o = filter(lambda i: self.geometry.numbers[i] == 8, r)[0]
            cs_close_to_n = set(self.geometry.closest(6, n, n=2, only=r))
            cs_close_to_o = set(self.geometry.closest(6, o, n=2, only=r))
            c2 = (cs_close_to_n & cs_close_to_o).pop()
            c4 = (cs_close_to_n - set([self.c2])).pop()
            c5 = (cs_close_to_o - set([self.c2])).pop()
            # Find out what type of ring it is
            # distance of N to third nearest C
            dist = self.geometry.dist(n, self.geometry.closest(6, n, n=3)[2])
            # If this distance is smaller than 1.75 A, it is a real bond.
            if dist <= 1.75:
                species = 'cation'
            else:
                species = 'monomer'
            # Add 'cation_' or 'monomer_'.
            for atom in ['n', 'o', 'c2', 'c4', 'c5']:
                setattr(self, '%s_%s' % (species, atom), eval(atom))

    @utils.cached
    def _hi(self):
        return self.get_hi_charges()

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
