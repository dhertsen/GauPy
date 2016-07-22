import copy
import warnings
import subprocess
import os
import utils
import logging
import numpy
from itertools import chain
from molmod.transformations import Translation, Rotation

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


class Stoichiometry(object):

    def __init__(self, molecule):
        self._array = np.bincount(molecule.numbers)
        '''
        index is the atomic number, item is the number
        of atoms with that atom number in the molecule
        '''
        self._string = ''
        '''string representation'''
        # print these atom numbers first in the string representation
        first = [z for z in [6, 1, 8, 7] if z < len(self._array)]
        rest = [z for z in range(len(self._array)) if z not in first]
        # for z in first + rest:
        for i, z in enumerate(first + rest):
            self._string += self._atomstr(z)

    def _atomstr(self, z):
        '''print things like 'C6' based on the atom number (z) or even
        the atom symbol'''
        try:
            n = self._array[utils.anum(z)]
        except:
            n = 0
        if n:
            if n > 1:
                return '%s%i' % (utils.asym(z).upper(), n)
            elif n == 1:
                return '%s' % utils.asym(z).upper()
        else:
            return ''

    def __str__(self):
        '''returns the string representation'''
        return self._string

    def __getitem__(self, i):
        '''both Stoichio['C'] and Stoichio[6] will yield the number
        of carbon atoms'''
        try:
            return self._array[utils.anum(i)]
        # What happens if the number of Hf atoms is requested,
        # but only atoms up to Cl are included?
        except IndexError:
            return 0

    def __eq__(self, other):
        '''numpy arrays, other stoichiometries or string representations
        can be compared'''
        if isinstance(other, Stoichiometry):
            return np.array_equal(self._array, other._array)
        elif isinstance(other, str):
            return self._string.strip().lower() == other.strip().lower()
        else:
            try:
                return np.array_equal(self._array, other)
            except:
                return False

    def __ne__(self, other):
        return not self.__eq__(other)


class SuperMolecule(mol.Molecule):
    ''' mention from_file, write_to_file '''

    def __init__(self, *pargs, **kargs):
        super(SuperMolecule, self).__init__(*pargs, **kargs)
        self.scaling = 1.0
        self.stoichiometry = Stoichiometry(self)

    def scale(self, nimag, factor=[1.0, 1.5]):
        if nimag == 0:
            self.scaling = factor[0]
        elif nimag == 1:
            self.scaling = factor[1]

    @classmethod
    def from_Molecule(cls, molecule):
        '''
        Make SuperMolecule from molmod.molecules.Molmod instance.
        '''
        newmolecule = copy.deepcopy(molecule)
        newmolecule.__class__ = SuperMolecule
        newmolecule.stoichiometry = Stoichiometry(newmolecule)
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

    @classmethod
    def from_xyz(cls, filename):
        with open(filename) as f:
            return cls.from_string(f.read(), ignore_header=True)

    def to_string(self, header=False, comment='', symbols=True, width=15):
        lines = []
        if header:
            lines += [str(self.size), str(comment)]
        for a, c in zip(self.numbers, self.coordinates):
            a = utils.asym(a) if symbols else a
            x, y, z = c / units.angstrom
            #lines.append('%-6s %15.6f %15.6f %15.6f' % (a, x, y, z))
            lines.append('%-6s %*.6f %*.6f %*.6f' % (a, width, x, width, y,
                                                     width, z))
        return '\n'.join(lines)

    def to_xyz(self, filename):
        with open(filename, 'w') as f:
            string = self.to_string(header=True)
            f.write(string)
            print('%s written' % filename)

    @classmethod
    def join(cls, *molecules):
        numbers = list(chain(*[m.numbers for m in molecules]))
        coordinates = list(chain(*[m.coordinates for m in molecules]))
        return cls(numbers, coordinates)

    def dist(self, *atoms):
        '''
        Calculate distance between atom1 and atom2 in molecule geom
        (molod.molecules.Molecule). In angstrom!

        Arguments:
            atom1/2:        number of atoms in xyz matrix
        Return:
            float
        '''
        return self.distance_matrix[atoms[0]][atoms[1]] / units.angstrom

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

    def nrings(self, n):
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
        '''no hydrogens in self.unparsed!'''
        self._graph()
        if ignore_hydrogens:
            self.unparsed = self.nonhydrogens()
        else:
            self.unparsed = range(self.size)

    def set_match(self, name, pattern):
        '''pattern can be both a molmod pattern and a the number of the
        atom in the xyz matrix'''
        logging.debug('SuperMolecule.set_match(): begin pattern %s' % name)
        # try a molmod pattern
        try:
            for u in self.unparsed:
                if pattern(u, self.graph):
                    setattr(self, name, Match(u, self))
                    self.unparsed.remove(u)
                    logging.debug(
                        'SuperMolecule.set_match(): %s at atom %i' % (name, u))
                    return
            else:
                logging.debug(
                    'SupperMolecule.set_match(): pattern %s not found' % name)
        # if not, a number?
        except:
            try:
                if pattern in self.unparsed:
                    setattr(self, name, Match(pattern, self))
                    self.unparsed.remove(pattern)
                    logging.debug(
                        'SuperMolecule.set_match(): %s at atom %i' % (
                            name, pattern))
                else:
                    logging.warning(
                        'SuperMolecule.set_match(): '
                        + '%s: atom %i not in unparsed list' % (name, pattern))
            # neither a molmod pattern, nor a number
            except:
                logging.warning(
                    'SuperMolecule.set_match(): '
                    + '%s: invalid formatting of %s pattern' % (name, pattern))

    def set_matches(self, patterns):
        '''
        patterns is OrderedDict {attributes: patterns}
        If pattern with name is already set, ignore. Want to redo it?
        Use initiate_match() again. This way you can flush an OrderedDict
        in between, but don't have to empty it every time. You can still
        access previous patterns (not emptied) and already access certain
        parsed atoms (flushed).
        '''
        for name, pattern in patterns.iteritems():
            if not hasattr(self, name):
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

    def bonded(self, atom1, atom2):
        '''
        Is there a bond between atom1 and atom2? Returns a boolean.
        '''
        self._graph()
        return atom2 in self.graph.neighbors[atom1]

    def non_hydrogen_neighbors(self, atom):
        ''' return non-hydrogen neighbors of an atom'''
        self._graph()
        return [n for n in self.graph.neighbors[atom] if self.numbers[n] != 1]

    def transform(self, transformation):
        return self.__init__(self.numbers,
                             transformation.apply_to(self.coordinates))

    def translate(self, vector):
        self.transform(Translation(vector))

    def rotate(self, center, axis, angle):
        self.transform(Translation(-np.array(center)))
        self.transform(Rotation.from_properties(angle, axis, invert=False))
        self.transform(Translation(center))

    def randomize(self, center):
        self.transform(Translation(-np.array(center)))
        self.transform(Rotation.random())
        self.transform(Translation(center))

    def mirror(self, plane='xy'):
        if plane == 'xy':
            x = y = 1
            z = -1
        elif plane == 'xz':
            x = z = 1
            y = -1
        elif plane == 'yz':
            y = z = 1
            x = -1
        reflection = np.transpose(np.array([x * self.coordinates[:, 0],
                                            y * self.coordinates[:, 1],
                                            z * self.coordinates[:, 2]]))
        self.coordinates = reflection

    def atoms(self, at):
        '''return the indices of the atoms with the given atom number or symbol'''
        atom_type = utils.anum(at)
        return numpy.where(self.numbers == atom_type)[0]
