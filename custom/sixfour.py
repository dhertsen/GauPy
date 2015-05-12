import gaupy.log
import gaupy.molecules
from gaupy.patterns import h, c, o, me, ch2, ipr
import molmod.molecular_graphs as g
import molmod.graphs as gr
from collections import OrderedDict
import numpy as np
import math
import warnings

# Dirty trick to silence warnings of confliciting
# modules on the HPC cluster and to silence NumPy
# warnings.
warnings.filterwarnings('ignore')

__all__ = ['IsoPropyl']


class IsoPropyl(gaupy.log.LOGFile):

    def __init__(self, filename):
        super(IsoPropyl, self).__init__(filename)
        self._classify()
        self._exoendo()
        self._cistrans()

    def _classify(self):
        # Patterns
        p = OrderedDict()
        c_bonded_with_cch = g.HasNeighborNumbers(1, 6, 6)
        p['furan_o'] = o
        p['furan_c5'] = g.HasNeighbors(o, me(6), c_bonded_with_cch)
        p['furan_me'] = g.HasNeighbors(p['furan_c5'], h, h, h)
        p['furan_c4'] = g.HasNeighbors(p['furan_c5'], h, c)
        p['furan_c2'] = g.HasNeighbors(o, c, c)
        p['diene_c4'] = g.HasNeighbors(c_bonded_with_cch, c, me(6))
        p['ipr_tert'] = g.HasNeighbors(me(6), me(6), h, c)
        p['ipr_me1'] = me(ipr(6))
        p['ipr_me2'] = me(ipr(6))
        p['diene_me'] = g.HasNeighbors(p['diene_c4'], h, h, h)
        p['diene_c1'] = g.HasNeighbors(p['ipr_tert'], c, c, h)
        p['diene_c2'] = ch2(p['diene_c1'], 6)
        p['diene_head1'] = gr.CritOr(g.HasNeighbors(p['diene_c1'], c, h),
                                     g.HasNeighbors(p['diene_c1'], c, c, h))
        p['diene_c3'] = g.HasNeighbors(p['diene_head1'], p['diene_c4'], h)
        p['diene_head2'] = gr.CritOr(g.HasNeighbors(p['diene_c4'],
                                                    p['diene_c2'], h),
                                     g.HasNeighbors(p['diene_c4'],
                                                    p['diene_c2'], h, c))
        p['furan_c3'] = gr.CritOr(g.HasNeighbors(p['furan_c4'], c, h),
                                  g.HasNeighbors(p['furan_c4'], c, c, h))
        p['furan_c2'] = g.HasNeighbors(p['furan_c3'], o, c)
        p['furan_cation_me1'] = me(6)
        p['furan_cation_me2'] = me(6)

        # Pattern matching
        self.geometry.initiate_match()
        self.geometry.set_matches(p)
        if len(self.geometry.unparsed) == 1:
            self.geometry.furan_cation = gaupy.molecules.Match(
                self.geometry.unparsed[0], self.geometry)

    def _exoendo(self):
        diene = ((self.geometry.diene_head1.xyz
                  - self.geometry.diene_c3.xyz)
                 + (self.geometry.diene_head2.xyz
                    - self.geometry.diene_c4.xyz))
        diene /= np.linalg.norm(diene)
        dienophile_subs = (self.geometry.furan_o.xyz
                           - self.geometry.furan_c2.xyz)
        dienophile_subs /= np.linalg.norm(dienophile_subs)
        angle = math.acos(np.dot(diene, dienophile_subs))
        # stereo attributes is stored both is self and self.geometry
        # (handy for making tables)
        if 0 <= angle <= (math.pi / 2):
            self.exoendo = self.geometry.stereo = 'endo'
        else:
            self.exoendo = self.geometry.stereo = 'exo'

    def _cistrans(self):
        ipr = self.geometry.ipr_tert.xyz - self.geometry.diene_c1.xyz
        ipr /= np.linalg.norm(ipr)
        cis = self.geometry.furan_cation.xyz - self.geometry.diene_head2.xyz
        cis /= np.linalg.norm(cis)
        trans = self.geometry.diene_c4.xyz - self.geometry.diene_head2.xyz
        trans /= np.linalg.norm(trans)
        cisangle = math.acos(np.dot(cis, ipr))
        transangle = math.acos(np.dot(trans, ipr))
        if cisangle < transangle:
            self.cistrans = self.geometry.cistrans = 'cis'
        else:
            self.cistrans = self.geometry.cistrans = 'trans'
