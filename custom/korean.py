import gaupy.log
import gaupy.molecules
from gaupy.patterns import h, c, n, o
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
        logging.debug('SixThreeBicyclic.__init__(): %s' % filename)
        super(SixThreeBicyclic, self).__init__(filename)
        self.scalings = [1.0, 1.6]
        try:
            self._nucleophile()
        except:
            logging.error('%s SixThreeBicyclic.__init__(): nucleophile could not be parsed.' % self.file)
        try:
            self._classify()
        except:
            logging.error('%s SixThreeBicyclic.__init__(): system could not be classified.' % self.file)
        try:
            self._pucker()
        except:
            logging.error('%s SixThreeBicyclic.__init__(): pucker could not be determined.' % self.file)
        try:
            self._species()
        except:
            logging.error('%s SixThreeBicyclic.__init__(): species could not be determined.' % self.file)

    def _nucleophile(self):

        # determination of nucleophile based on stereochemistry,
        # regardless of the number of explicit CH3CN molecules
        hcn = self.stoichiometry['H'] - self.stoichiometry['C'] - self.stoichiometry['N'] 
        os = self.stoichiometry['O'] + self.stoichiometry['S']
        hcrat = self.stoichiometry['H'] / self.stoichiometry['C']
        logging.debug('SixThreeBicyclic._nucleophile(): hcn = %s' % hcn)
        logging.debug('SixThreeBicyclic._nucleophile(): os = %s' % os)
        logging.debug('SixThreeBicyclic._nucleophile(): hcrat = %s' % hcrat)
        
        if self.stoichiometry['Cl']:
            self.nucleophile = 'Cl'
            self.attacking_atom = 17
        elif self.stoichiometry['Br']:
            self.nucleophile = 'Br'
            self.attacking_atom = 35
        elif hcn == 2:
            if os == 0:
                self.nucleophile = 'N3'
                self.attacking_atom = 7
        elif hcn == 3:
            if os == 0:
                self.nucleophile = 'CN'
                self.attacking_atom = 6
            elif os == 1:
                self.nucleophile = 'SCN'
                self.attacking_atom = 16
        elif hcn == 5:
            if os == 2:
                self.nucleophile = 'acryl'
                self.attacking_atom = 8
        elif hcn == 6:
            if os == 0:
                self.nucleophile = 'NH2'
                self.attacking_atom = 7
            elif os == 1:
                self.nucleophile = 'OH'
                self.attacking_atom = 8
        elif hcn == 7:
            if os == 1:
                # holds at least up to 6 CH3CN
                if hcrat > 1.5:
                    self.nucleophile = 'OMe'
                    self.attacking_atom = 8
                elif hcrat < 1.5:
                    self.nucleophile = 'OPh'
                    self.attacking_atom = 8
        logging.debug('SixThreeBicyclic._system(): nucleophile: %s' % self.nucleophile)

    def _classify(self):
        p = OrderedDict()
        self.geometry.initiate_match()
        p['ph'] = gr.CritAnd(c, g.HasNeighbors(c, c, c))
        self.geometry.set_matches(p)

        if self.nucleophile == 'OPh':
            # find the O just to check which ph was parsed
            self.geometry.set_match('o', o)
            # wrong phenyl group selected, select the other one
            if self.geometry.dist(self.geometry.ph.n, self.geometry.o.n) < 3.0:
                self.geometry.set_match('ph', gr.CritAnd(c, g.HasNeighbors(c, c, c)))

        # remove phenyl group to speed up pattern searching
        phgroup = self.geometry.closest(6, self.geometry.ph.n, n=2)
        phgroup += [self.geometry.closest(6, i, exclude=[self.geometry.ph.n])
                    for i in phgroup]
        phgroup += [self.geometry.closest(6, phgroup[-1], exclude=phgroup)]
        logging.debug(
            'SixThreeBicyclic._classify(): phenyl group %s' % (
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
        if self.attacking_atom:
            self.geometry.set_match(
                'nu', self.geometry.closest(
                    self.attacking_atom, self.geometry.c5.n,
                    only=self.geometry.unparsed))
        else:
            logging.debug('SixThreeBicyclic._classify(): nu not parsed, since no nucleophile was found.')

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
