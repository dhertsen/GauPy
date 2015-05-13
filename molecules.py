import copy
import warnings
import subprocess
import os
import utils
import logging

# Dirty trick to silence warnings of confliciting
# modules on the HPC cluster and to silence NumPy
# warnings.
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    import molmod.graphs as grp
    import molmod.molecular_graphs as mgrp
    import numpy as np
    import molmod.molecules as mol
    import molmod.units as units
    import molmod.ic as ic


class Match(object):
    '''
    used in the classification of a molecules
    logfile.geometry.somecenter.xyz
    '''

    def __init__(self, number, supermolecule):
        self.n = number
        self.xyz = supermolecule.coordinates[number]
        self.z = supermolecule.numbers[number]


class SuperMolecule(mol.Molecule):
    ''' mention from_file, write_to_file '''

    def __init__(self, *pargs, **kargs):
        super(SuperMolecule, self).__init__(*pargs, **kargs)

    @classmethod
    def from_Molecule(cls, molecule):
        '''
        Make SuperMolecule from molmod.molecules.Molmod instance.
        '''
        newmolecule = copy.deepcopy(molecule)
        newmolecule.__class__ = SuperMolecule
        return newmolecule

    @classmethod
    def from_string(cls, string, ignore_header=False):
        ''' handles both symbols and numbers,
        take three last cols and first col'''
        atoms = []
        coordinates = []
        for i, l in enumerate(string.strip().split('\n')):
            if ignore_header and i < 2:
                continue
            split = l.split()
            atoms.append(utils.anum(split[0]))
            coordinates.append(np.array(map(float, split[-3:]))
                               * units.angstrom)
        return cls.from_Molecule(mol.Molecule(atoms, coordinates))

    def to_string(self, header=False, comment='', symbols=False):
        lines = []
        if header:
            lines += [str(self.size), str(comment)]
        for a, c in zip(self.numbers, self.coordinates):
            a = utils.asym(a) if symbols else a
            x, y, z = c / units.angstrom
            lines.append('%-6s %15.6f %15.6f %15.6f' % (a, x, y, z))
        return '\n'.join(lines)

    def dist(self, atom1, atom2):
        '''
        Calculate distance between atom1 and atom2 in molecule geom
        (molod.molecules.Molecule).

        Arguments:
            atom1/2:        number of atoms in xyz matrix
        Return:
            float
        '''
        return self.distance_matrix[atom1][atom2] / units.angstrom

    def angle(self, *atoms):
        return ic.bend_angle([self.coordinates[i]
                              for i in atoms])[0] / units.deg

    def dihedral(self, *atoms):
        ''' IUPAC sign convention'''
        return ic.dihed_angle([self.coordinates[i]
                               for i in atoms])[0] / units.deg

    def opbend_angle(self, *atoms):
        '''angle plane 012 and vector 03'''
        return ic.opbend_angle([self.coordinates[i]
                                for i in atoms])[0] / units.deg

    def opbend_dist(self, *atoms):
        ''' distance plane 012 and vector 03'''
        return ic.opbend_dist([self.coordinates[i]
                               for i in atoms])[0] / units.angstrom

    def closest(self, atom_type, reference, n=1, exclude=False, only=False):
        '''
        If n=1: return the number in the xyz matrix of the atom of atom_type
        which is the closest to a reference atom in a geometry.
        If n>1, return a list of the n closest atoms of atom_type, which are
        closest to the reference atom.

        Arguments:
            atom_type:      atom number or symbol of the atoms whose number
                            in the xyz matrix
                            must be returned
            reference:      number of the reference atom in the xyz matrix
            n:              number of atoms to return
            exclude:        exclude these atoms
            only:           limit search to these atoms
        '''
        # If search is not restricted, look for any atom.
        if not only:
            only = range(self.size)
        atom_type = utils.anum(atom_type)
        distances = self.distance_matrix[reference]
        out = [np.where(distances == d)[0][0] for d in sorted(distances)
               if self.numbers[np.where(distances == d)[0][0]] == atom_type
               and np.where(distances == d)[0][0] != reference
               and np.where(distances == d)[0][0] in only]
        if exclude:
            out = [o for o in out if o not in exclude]
        if n > 1:
            return out[:n]
        elif n == 1:
            return out[0]

    def nonhydrogens(self):
        return [i for i, at in enumerate(self.numbers) if at != 1]

    def avogadro(self):
        ''' writes a file'''
        self.write_to_file('.avogadro.xyz')
        process = subprocess.Popen(['avogadro', '.avogadro.xyz'],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        process.wait()
        os.remove('.avogadro.xyz')

    def part(self, atoms):
        '''
        make a new SuperMolecule
        '''
        newnumbers = []
        newcoord = np.empty([len(atoms), 3])
        for i, a in enumerate(atoms):
            newnumbers.append(self.numbers[a])
            newcoord[i] = (self.coordinates[a])
        return type(self)(newnumbers, newcoord)

    def _graph(self):
        if self.graph:
            return
        else:
            try:
                self.graph = mgrp.MolecularGraph.from_geometry(
                    self, scaling=self.scaling)
            except:
                self.set_default_graph()

    def get_nrings(self, n):
        '''
        Look for all n-membered rings
        '''
        self._graph()
        rings = []
        match_generator = grp.GraphSearch(mgrp.NRingPattern(n),)
        for ring_match in match_generator(self.graph):
            ring = ([b for a, b in ring_match.forward.items()])
            if all([not utils.is_translation(ring, x) for x in rings]):
                rings.append(ring)
        return rings

    def initiate_match(self, ignore_hydrogens=True):
        self._graph()
        if ignore_hydrogens:
            self.unparsed = self.nonhydrogens()
        else:
            self.unparsed = range(self.size)

    def set_match(self, name, pattern):
        logging.debug('SuperMolecule.set_match(): begin pattern %s' % name)
        for u in self.unparsed:
            if pattern(u, self.graph):
                setattr(self, name, Match(u, self))
                self.unparsed.remove(u)
                logging.debug(
                    'SuperMolecule.set_match(): end pattern %s' % name)
                return

    def set_matches(self, patterns):
        '''
        patterns is OrderedDict {attributes: patterns}
        returns the unparsed ones
        '''
        for name, pattern in patterns.iteritems():
            self.set_match(name, pattern)

    def molecules(self):
        '''
        Return the molecules as seperate SuperMolecules
        '''
        self._graph()
        molecules = []
        for molecule in self.graph.independent_vertices:
            molecules.append(self.part(molecule))
        return molecules

    def break_bond(self, atom1, atom2):
        '''
        Return two SuperMolecules by breaking bond
        '''
        self._graph()
        mol1, mol2 = self.graph.get_halfs(atom1, atom2)
        return self.part(mol1), self.part(mol2)

    def neighbors(self, atom1, atom2):
        '''
        Is there a bond between atom1 and atom2? Returns a boolean.
        '''
        self._graph()
        return atom2 in self.graph.neighbors[atom1]
