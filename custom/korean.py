import gaupy.log
import gaupy.molecules
from gaupy.patterns import h, c, n
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

__all__ = ['SixThreeBicyclic']


class SixThreeBicyclic(gaupy.log.LOGFile):

    def __init__(self, filename):
        super(SixThreeBicyclic, self).__init__(filename)

        self.scalings = [1.0, 1.6]

        self._system()
        self._classify()
        self._pucker()
        self._species()

    def _system(self):
        # nucleophile
        if self.stoichiometry['Cl']:
            self.nucleophile = 'Cl'
            self.attacking_atom = 17
        else:
            C = self.stoichiometry['C']
            N = self.stoichiometry['N']
            if (C - N) == 10:
                self.nucleophile = 'N3'
                self.attacking_atom = 7
            elif (C - N) == 13:
                self.nucleophile = 'CN'
                self.attacking_atom = 6

    def _classify(self):
        p = OrderedDict()
        self.geometry.initiate_match()
        p['ph'] = gr.CritAnd(c, g.HasNeighbors(c, c, c))
        self.geometry.set_matches(p)

        # remove phenyl group to speed up pattern searching
        phgroup = self.geometry.closest(6, self.geometry.ph.n, n=2)
        phgroup += [self.geometry.closest(6, i, exclude=[self.geometry.ph.n])
                    for i in phgroup]
        phgroup += [self.geometry.closest(6, phgroup[-1], exclude=phgroup)]
        logging.debug(
            'korean.SixThreeBicyclic._classify(): phenyl group %s' % (
                phgroup + [self.geometry.ph.n]))
        for i in phgroup:
            self.geometry.unparsed.remove(i)

        self.geometry.set_match('c7', self.geometry.closest(
            6, self.geometry.ph.n, only=self.geometry.unparsed))
        self.geometry.set_match('me', self.geometry.closest(
            6, self.geometry.c7.n, only=self.geometry.unparsed))
        self.geometry.set_match('nplus', self.geometry.closest(
            7, self.geometry.c7.n, only=self.geometry.unparsed))

        nplus = gr.CritAnd(n, gr.CritOr(g.HasNeighbors(c, c, c, c),
                                        g.HasNeighbors(c, c, c)))
        c_nextto_nplus = g.HasNeighbors(nplus, h, h, c)
        p['c2'] = g.HasNeighbors(c_nextto_nplus, h, h, c)
        p['c3'] = g.HasNeighbors(p['c2'], h, h, c)
        p['c4'] = g.HasNeighbors(p['c3'], h, h, c)
        p['c1'] = g.HasNeighbors(nplus, h, h, p['c2'])
        self.geometry.set_matches(p)
        self.geometry.set_match(
            'c5', self.geometry.closest(6, self.geometry.c4.n,
                                        only=self.geometry.unparsed))
        self.geometry.set_match(
            'c6', self.geometry.closest(
                6, self.geometry.nplus.n, only=self.geometry.unparsed))
        self.geometry.set_match(
            'nu', self.geometry.closest(
                self.attacking_atom, self.geometry.c5.n,
                only=self.geometry.unparsed))

    def _pucker(self):
        middle = (self.geometry.c2.xyz + self.geometry.c4.xyz) / 2
        pucker = self.geometry.c3.xyz - middle
        pucker /= np.linalg.norm(pucker)

        plane1 = self.geometry.c2.xyz - self.geometry.c1.xyz
        plane2 = self.geometry.c3.xyz - self.geometry.c1.xyz
        cross = np.cross(plane1, plane2)
        cross /= np.linalg.norm(cross)

        dot = np.dot(pucker, cross)
        direction = math.acos(dot) / math.pi * 180
        if direction < 60:
            self.pucker = 'chair'
        else:
            self.pucker = 'boat'

    def _species(self):
        c5 = self.geometry.c5.n
        c6 = self.geometry.c6.n
        nplus = self.geometry.nplus.n
        if 'irc' in self.route_section:
            first = self.irc.geometries[0]
            last = self.irc.geometries[-1]
            deltac5n = last.dist(c5, nplus) - first.dist(c5, nplus)
            deltac6n = last.dist(c6, nplus) - first.dist(c6, nplus)
            if deltac5n > 0.05 and abs(deltac6n) < 0.05:
                self.species = 'irc3'
            elif deltac6n > 0.05 and abs(deltac5n) < 0.05:
                self.species = 'irc2'
            elif deltac5n + deltac6n < -0.05:
                self.species = 'irc1'
        else:
            if self.lowest_frequency < -200:
                c5n_length = self.geometry.dist(c5, nplus)
                c6n_length = self.geometry.dist(c6, nplus)
                ratio = c5n_length / c6n_length
                if ratio < 1:
                    self.species = 'tsi'
                else:
                    self.species = 'tsii'
            elif not self.nimag:
                c5n = self.geometry.bonded(c5, nplus)
                c6n = self.geometry.bonded(c6, nplus)
                if c5n and not c6n:
                    self.species = '2'
                elif c6n and not c5n:
                    self.species = '3'
                elif c5n and c6n:
                    self.species = '1'

    @classmethod
    def partition(cls, logs):
        return super(SixThreeBicyclic, cls).partition(
            logs, add_patterns=['tsi', 'tsii', '1', '2', '3', 'irc1', 'irc3',
                                'irc2'])
