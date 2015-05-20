import gaupy.log
import gaupy.molecules
from gaupy.patterns import h, c, o, s, me
import molmod.molecular_graphs as g
import molmod.graphs as gr
from collections import OrderedDict
import numpy as np
import math
import warnings
import logging

# Dirty trick to silence warnings of confliciting
# modules on the HPC cluster and to silence NumPy
# warnings.
warnings.filterwarnings('ignore')

__all__ = ['SixFour']


class SixFour(gaupy.log.LOGFile):

    def __init__(self, filename):
        super(SixFour, self).__init__(filename)

        # adjusted scaling necessary for thiophene molecular graphs
        self.scalings = [1.0, 1.6]

        # stoichiometry determines which study is being studied
        if self.stoichiometry == 'C14H19O':
            self.system = 'furan'
        elif self.stoichiometry == 'C14H19S':
            self.system = 'thiophene'
        elif self.stoichiometry == 'C18H27O':
            self.system = 'ipr-me'
        elif self.stoichiometry == 'C17H25O':
            self.system = 'ipr-h'

        self._species()
        if 'reactant' not in self.species:
            self._exoendo()
            self._cistrans()

    def _classify(self):

        # custom patterns
        os = gr.CritOr(o, s)
        c_bonded_with_cch = g.HasNeighborNumbers(1, 6, 6)

        # patterns in specific order, since the groups recognized
        # will be removed from self.geometry.unparsed sequentially
        p = OrderedDict()
        p['six_hetero'] = os
        p['six_c5'] = g.HasNeighbors(os, me(6), c_bonded_with_cch)
        p['six_me'] = g.HasNeighbors(p['six_c5'], h, h, h)
        p['six_c4'] = g.HasNeighbors(p['six_c5'], h, c)
        p['six_c2'] = g.HasNeighbors(os, c, c)
        # In the first case, two methyl groups on the cationic position.
        # In the second case, a single one.
        p['six_cation'] = gr.CritOr(
            g.HasNeighbors(p['six_c2'], c, me(6), h),
            g.HasNeighbors(p['six_c2'], c, me(6), me(6)))
        # First case: before the second bond is formed.
        # Second case: after the second bond is formed.
        p['six_c3'] = gr.CritOr(
            g.HasNeighbors(p['six_c4'], p['six_c2'], h),
            g.HasNeighbors(p['six_c4'], p['six_c2'], h, c))
        # first bond has been formed when _classify() is called
        p['four_head2'] = g.HasNeighbors(p['six_cation'], c, c, h)
        p['four_c2'] = g.HasNeighbors(p['four_head2'], c, h, h)
        p['four_c1'] = gr.CritOr(
            g.HasNeighbors(p['four_c2'], c, h, h),
            g.HasNeighbors(p['four_c2'], c, c, h))
        p['four_c4'] = gr.CritOr(
            g.HasNeighbors(p['four_head2'], c, h),
            g.HasNeighbors(p['four_head2'], c, c))
        p['four_c3'] = g.HasNeighbors(p['four_c4'], c, h)
        # First case: before the second bond is formed.
        # Second case: after the second bond is formed.
        p['four_head1'] = gr.CritOr(
            g.HasNeighbors(p['four_c3'], c, h),
            g.HasNeighbors(p['four_c3'], c, c, h))

        # TODO
        '''
        - six and four instead of furan and diene
        - hetero instead of o
        - only diene_ipr further parsed, because it's necessay
          for cis/trans
        '''
        if 'ipr' in self.system:
            p['four_ipr'] = g.HasNeighbors(c, c, c, h)

        self.geometry.initiate_match()
        logging.debug('call set_matches from sixfour.SixFour._classify()')
        self.geometry.set_matches(p)

    def _exoendo(self):
        four = ((self.geometry.four_head1.xyz - self.geometry.four_c3.xyz)
                + (self.geometry.four_head2.xyz
                    - self.geometry.four_c4.xyz))
        four /= np.linalg.norm(four)
        six_subs = (self.geometry.six_c2.xyz
                    - self.geometry.six_hetero.xyz)
        six_subs /= np.linalg.norm(six_subs)
        angle = math.acos(np.dot(four, six_subs))
        # stereo attributes is stored both is self and self.geometry
        # (handy for making tables)
        if 0 <= angle <= (math.pi / 2):
            self.exoendo = self.geometry.stereo = 'endo'
        else:
            self.exoendo = self.geometry.stereo = 'exo'

    def _cistrans(self):

        if 'ipr' in self.system:
            # check on which side of the ring the iPr group is placed
            # unit vector perpendicular on ring, centered on six_c1
            n1 = np.cross(self.geometry.four_c2.xyz
                          - self.geometry.four_c1.xyz,
                          self.geometry.four_head1.xyz
                          - self.geometry.four_c1.xyz)
            n1 /= np.linalg.norm(n1)
            # unit vector from that six_c1 center to iPr
            ipr = self.geometry.four_ipr.xyz - self.geometry.four_c1.xyz
            ipr /= np.linalg.norm(ipr)
            # how are these two vector aligned with respect to eachother
            iprside = math.acos(np.dot(n1, ipr)) < (math.pi / 2)

            # do the same thing around the second bridge head
            # unit vector perpendicular on ring, centered on six_head2
            n2 = np.cross(self.geometry.four_c2.xyz
                          - self.geometry.four_head2.xyz,
                          self.geometry.four_c4.xyz
                          - self.geometry.four_head2.xyz)
            n2 /= np.linalg.norm(n2)
            # unit vector from that six_head center to cation
            head = self.geometry.six_cation.xyz - self.geometry.four_head2.xyz
            head /= np.linalg.norm(head)
            # alignment
            headside = math.acos(np.dot(n2, head)) < (math.pi / 2)

            # Define cis and trans. Don't bother with the signs, just check
            # whether == or != belongs to cis or trans.
            if iprside != headside:
                self.cistrans = self.geometry.cistrans = 'cis'
            else:
                self.cistrans = self.geometry.cistrans = 'trans'
        else:
            self.cistrans = None

    def _species(self):
        '''classify system: ts1, ts2, product, reactant, int
        product (imag), reactant (imag), int(imag), ts1(multi), ts2(multi)'''

        # reactant system consist of two seperate molecules
        if len(self.geometry.molecules()) == 2:
            if self.nimag == 0:
                self.species = 'reactant'
            else:
                self.species = 'reactant (imag)'
        # in all other cases, atoms have to be classified
        else:
            self._classify()
            # second bond formed?
            second = self.geometry.bonded(self.geometry.four_head1.n,
                                          self.geometry.six_c3.n)

            # let's classify single-molecule species

            # single imaginatry frequencies can be TSs or annoyances
            # in intermediates and products
            if self.nimag >= 1:

                # second bridge had been formed
                if second:
                    # still reacting?
                    print self.vibrating_atoms
                    if (set(self.vibrating_atoms)
                        == set([self.geometry.four_head1.n,
                                self.geometry.six_c3.n])):
                        if self.nimag == 1:
                            self.species = 'ts2'
                        else:
                            self.species = 'ts2(multi)'
                    # both bridges have been formed, but still a neg freq
                    else:
                        self.species = 'product(imag)'

                # only first bridge had been formed
                else:
                    # still reacting?
                    if (set(self.vibrating_atoms)
                        == set([self.geometry.four_head2.n,
                                self.geometry.six_cation.n])):
                        if self.nimag == 1:
                            self.species = 'ts1'
                        else:
                            self.species = 'ts1(multi)'
                    # already reacted, but still an annoying neg freq
                    else:
                        self.species = 'int(imag)'
            else:
                if second:
                    self.species = 'product'
                else:
                    self.species = 'intermediate'
