'''Parsing Gaussian output files'''

import re
import warnings
import filenames
import molecules
import fnmatch
from gaupy.utils import find_all, cached
import numpy as np
import molmod.molecules as mol
import molmod.units as units
import molmod.bonds as bonds
import logging
from math import log10

# Dirty trick to silence warnings of confliciting
# modules on the HPC cluster and to silence NumPy
# warnings.
warnings.filterwarnings('ignore')


def energy(value, type='gibbs'):
    '''
    Convert a energy or a filename to an energy. If a filename was provided,
    energy of type ``type`` will be parsed. ``float('nan')`` values are
    allowed and will even be returned if an error occurs. Passing LOGFile
    objects also works.
    '''
    try:
        return float(value)
    except:
        try:
            return float(getattr(LOGFile(value), type))
        except:
            try:
                return float(getattr(value, type))
            except:
                return float('nan')


def energies(values, type='gibbs'):
    '''
    Convert a list or NumPy array of energies and filenames to a NumPy array
    of energies. Energy of type ``type`` will be taken from log files.
    ``float('nan')`` values are allowed.
    '''
    if values:
        return np.array([energy(e, type) for e in values])
    else:
        return np.array([])


def relative_energies(values, reference='min', absolute=0, relative=0,
                      conversion='kjmol', type='gibbs'):
    '''
    Convert a list of energies and filanemes to a list of relative energies.

        First of all, the internal reference is applied. The original
        ``energies`` array is transformed to **internal reference**:
            ------------------------    -------------------------------------
            ``reference=``              reference
            ------------------------    -------------------------------------
            ``'min'``                   smallest energy in array
            ``'max'``                   largest energy in array
            ``'none'``                  use absolute values
            ``i``                       i-th energy of array
            ------------------------    -------------------------------------

        If any filename is present, the energy will be parsed from the
        corresponding log file first. The type of energy can be specified
        (see Gaussian white pages for more info):
            ------------------------    -------------------------------------
            ``type=``            type
            ------------------------    -------------------------------------
            ``'electronic'``            electronic
            ``'gibbs'``                 Gibbs free
            ``'enthalpy'``              enthalpy
            ``'zpe'``                   ZPE
            ``'zpesum'``                electronic energy + ZPE
            ``'thermal'``               thermal energy
            ``'gibbscorrection'``       Gibbs - electronic
            ``'enthalpycorrection'``    enthalpy - electronic
            ``'thermalcorrection'``     thermal - electronic
            ------------------------    -------------------------------------

        ``float('nan')`` values are allowed in the array. They will result
        in float('nan') values in the output array.

        The **internal reference** is then further converted as such:
            (**internal reference** + ``absolute``) * ``conversion``
            + ``relative``

        ``absolute`` may be a number or a filename, ``relative``
        is expected to be a number.
    '''

    # Necessary if [] or np.array([]) is passed
    if values:
        # Convert to numerical energies
        e = energies(values, type=type)
        # Internal reference
        if reference == 'min':
            e = e - np.nanmin(e)
        elif reference == 'max':
            e = e - np.nanmax(e)
        else:
            try:
                e = e - e[reference]
            # Failure of internal reference will result in the
            # parsed numerical energies
            except:
                pass
        # Absolute and relative references, conversion
        try:
            return ((e + energy(absolute)) / units.parse_unit(conversion)
                    + relative)
        # If the calculation fails, the numerical energies will be returned
        except:
            return e
    else:
        return np.array([])


class LOGFile(object):
    '''
    Parsed Gaussian output file from ``filename``. ``'.log'`` will be added
    to filenames without extension or will replace any other extension
    before parsing. Only this true log file will be parsed.
    '''

    # Default info message for special parser classes.
    info = 'No special info about this type.'

    def __init__(self, filename):
        logging.debug('log.LOGFile.__init__(): %s' % filename)
        self.files = filenames.GaussianFile(filename)
        '''
        :class:`gaussian.filenames.GaussianFile` object constructed
        from ``filename``. This attribute can be used to quickly access
        other extensions without the need for replace commands:
            =======================   ===========================
            Attribute                 Filename
            =======================   ===========================
            ``foo.files.root``        ``'calculation'``
            ``foo.files.log``         ``'calculation.log'``
            ``foo.files.chk``         ``'calculation.chk'``
            ``foo.files.arbitrary``   ``'calculation.arbitrary'``
            =======================   ===========================
        '''
        self.file = self.files.root
        '''
        Root of filename (filename without extension, with directory),
        equivalent with and shortcut for ``self.files.root``.
        '''
        self.scalings = [1.0, 1.5]
        '''
        Default scaling factors used for molecular graphics: first element
        for self.nimag==0, second element for self.nimag==1.
        '''

    # Methods

    def _get_geometry(self, pos, reverse=False, orientation='orientation'):
        '''
        Return the first geometry in the Gaussian output file after
        (``Reverse=False``) or before (``Reverse=True``) position ``pos``
        as a :class:`gaussian.molecules.SuperMolecule`.

        A geometry is marked by 'orientation' in the output file. Both Standard
        and Input Orientations can be found in most output files. The standard
        implementation of this function parses all of them. Reset orientation
        (``Standard orientation``, ``Input orientation`` if only a single
        type is wanted.
        '''

        # A regex to parse the coordinates from an 'orientation' section
        reg = re.compile((r'%s:.*$\n^ -*$\n^.*$\n^.*$\n^ -*'
                          + r'$\n([-\s\.0-9\n]*[0-9]\n)') % orientation,
                         re.MULTILINE)

        # Perform the reverse search if applicable. Necessary because reverse
        # regex searches are not possible.
        if reverse:
            npos = self._full.rfind('orientation', 0, pos) - 5
        else:
            npos = pos

        # Parse the coordinate lines with the regex
        lines = reg.search(self._full, npos).group(1)[:-1].split('\n')

        # Seperate these lines into atoms and coordinates. Remember:
        # molmod.molecules.Molecule (and gaussian.molecules.SuperMolecule)
        # have coordinates in atomic units, while Gaussian prints its
        # coordinates in angstrom.
        atoms = [l.split()[1] for l in lines]
        coordinates = [np.array(map(float, l.split()[-3:])) * units.angstrom
                       for l in lines]

        # Create a SuperMolecule
        geometry = mol.Molecule(atoms, coordinates)
        geometry = molecules.SuperMolecule.from_Molecule(geometry)
        geometry.scale(self.nimag, factor=self.scalings)
        return geometry

    def _get_geometries(self, orientation='orientation'):
        '''
        Return all geometries in the Gaussian output file as a list of
        :class:`gaussian.molecules.SuperMolecule`.

            Two types of coordinate systems are used by Gaussian:
                (1)  coordinate systems that respond to those used in the
                     input
                       (a)  ``orientation='Z-Matrix orientation'`` -- parse
                            all coordinates printed in internal coordinates,
                            only (?) if input is provided in internal
                            cordinates
                       (b)  ``orientation='Input orientation'`` -- parse all
                            coordinates printed with respect to the same
                            cartesian coordinate system as used in the input,
                            only (?) if input is provided as cartesian
                            coordinates
                (2)   coordinates systems that are used in the calculations
                      and are transformed for the sake of symmetry
                            ``orientation='Standard orientation'`` -- parse
                            all coordinates printed in these coordinates
        '''
        try:
            return [self._get_geometry(p, orientation=orientation)
                    for p in find_all(self._full, orientation)]
        except:
            return None

    def _get_energy(self, pos, reverse=False):
        '''
        Get the first electronic energy (Ha) before (``reverse=True``)  or
        after (``reverse=False``) a given position ``pos`` in the Gaussian
        output file.

            Electronic energies correspond to ``SCF Done`` entries in the
            Gaussian log file.

            Should work for force fields as well, although not
            thoroughly tested. Looks for ``Energy=``
            in that case.
        '''

        if 'SCF Done' in self._full:
            reg = re.compile(r'SCF Done: *E.* = *([-.0-9]*)')
            start = 'SCF Done'
        else:
            reg = re.compile(r'Energy= *([-.0-9]*)')
            start = 'NIter='
        # Perform an initial non-regex search because (1) reverse regex
        # searches are not possible and (2) a normal search is more
        # efficient than a regex one.
        if reverse:
            npos = self._full.rfind(start, 0, pos) - 50
        else:
            npos = self._full.find(start, pos) - 50
        return float(reg.search(self._full, npos).group(1))

    def _get_energy_correction(self, string):
        '''
        Search for a thermal correction to the electronic energy
        in the thermochemistry section of the log file.

            The ``string`` argument is the part before ``=``
            (e.g. ``'Thermal correction to Enthalpy'``,
            ``'Sum of electronic and thermal Energies'``).
        '''
        try:
            return float(re.search(string + r'= *([-.0-9]*)',
                                   self._thermochem_block).group(1))
        except:
            return float('nan')

    def to_input(self, filename=None):
        '''
        Return a Gaussian input file using the same methods and computational
        resources, starting from the end point of this calculation.

            Structure of the input file with (attributes)
            ::
                %nproc=(nproc)
                %mem=(mem)
                %chk=(chk)
                #(route_section)

                (comment)

                (charge) (multiplicity)
                (geometry)

                (modredundant)

                (scrf_nonstandard_input)

            If ``filename=None``, ``files.base`` will be used for ``(chk)``,
            else ``GaussianFile(filename).base`` will be used. In other words,
            the base (path without extension and directories) of the
            the ``filename`` argument or the current ``LOGFile``
            is used.

            Note that it is allowed to change attributes such as the route
            section, before you use the function. In this way, you can
            adapt calculation parameters before a subsequent input file
            is constructed.
        '''
        fn = filenames.GaussianFile(filename) if filename else self.files
        # This way, self.geometry = None cn be set, and no cartesian
        # coordinates are written. Can be handy if calculations are restarted
        # from a checkpoint file.
        try:
            geometry = self.geometry.to_string()
        except:
            geometry = ''
        inp = ''
        # Do not include %nproc and %memory if they were not provided in the
        # input file.
        if not np.isnan(self.nproc):
            inp += '%%nproc=%s\n' % self.nproc
        if self.memory:
            inp += '%%mem=%s\n' % self.memory
        # Do include %chk in all cases, even if it was not provided in the
        # input file.
        inp += ('%%chk=%(chk)s\n'
                '#%(route)s\n\n'
                '%(comment)s\n\n'
                '%(charge)i %(multiplicity)i\n'
                '%(geometry_text)s\n\n'
                '%(modredundant)s\n\n'
                '%(scrf)s\n\n') % {
                    'chk': fn.base,
                    'route': ' '.join(self.keywords),
                    'comment': self.comment,
                    'charge': self.charge,
                    'multiplicity': self.multiplicity,
                    'geometry_text': geometry,
                    'modredundant': self.modredundant,
                    'scrf': self.scrf_nonstandard_input}
        return inp.rstrip() + '\n\n'

    def write_input(self, filename=None):
        '''
        Write a Gaussian input file with the data parsed from the output
        file. See :meth:`gaussian.log.LOGFile.to_input` for more info.

            The base (path without extension and dirs) of the ``filename``
            argument or the current filename is used for the input file
            and the checkpoint line in the input file.
        '''
        fn = filenames.GaussianFile(filename) if filename else self.files
        f = open(fn.com, 'w')
        f.write(self.to_input(fn.base))
        f.close()
        print('%s written' % fn.com)

    def get_nbo_energy(self, atom1, atom2, debug=[]):
        '''
        Parse the NBO interaction energy between ``atom1`` and ``atom2``.
        Debug info will be stored in the ``debug`` list.
        '''

        # Block in which all interaction energies can be found.
        pertstr = ('Second Order Perturbation Theory Analysis '
                   + 'of Fock Matrix in NBO Basis')
        start = self._full.find('1.',
                                self._full.find('==\n',
                                                self._full.rfind(pertstr)))
        end = self._full.find('Natural Bond Orbitals (Summary)', start) - 4

        interaction = 0
        for line in self._full[start:end].split('\n'):
            if line and 'unit' not in line:
                line_groups = re.search('\)(.*)/.*\)(.*)\s+([\d.]+)\s+([\d.]+)'
                                        + '\s+([\d.]+)\s*', line).groups()
                energy = float(line_groups[2])
                atoms1 = [int(x) - 1 for x in re.findall('\d+',
                                                         line_groups[0])]
                atoms2 = [int(x) - 1 for x in re.findall('\d+',
                                                         line_groups[1])]
                if ((atom1 in atoms1 and atom2 in atoms2)
                        or (atom1 in atoms2 and atom2 in atoms1)):
                    interaction += energy
                    debug.append(line)

        # Return the total interaction energy.
        return interaction

    def get_hi_charges(self, csv_file=None):
        '''
        Hirshfeld I charges will be parsed from a csv file written by Horton
        and stored in the ``hi_charges`` attribute of the current
        instance as a list of atomic charges.

            If ``csv_file=None``, the csv file is assumed to have the same
            root as the current file. Else, the provided csv_file is used.
        '''
        if not csv_file:
            csv_file = self.files.csv
        csv = open(csv_file).read()
        begin = csv.find('Dataset,charges')
        # This will also work if the csv file is cut below the Dataset,charges
        # block (file compression). In that case end=-1, so that csv[begin:-1].
        end = csv.find('Dataset,history_charges')
        return map(float, csv[begin:end].split('\n')[2:-2])

    @cached
    def hi_charges(self):
        return self.get_hi_charges()

    @classmethod
    def parse_all(cls, *pargs):
        '''
        Return a list of :class:`gaussian.log.LOGFile` instances for the
        provided files, which may be provided with any or without extensions.
        '''
        return [cls(parg) for parg in pargs]

    def movie(self):
        '''
        Return an xyz movie from the Gaussian output file. This can be a movie
        of scan points, IRC points or geometry optimization steps. The latter
        will be created if no scan/IRC points are found. The comments contain
        the step number.
        '''
        base = self.irc or self.scan or self
        geometries = base.geometries
        return ''.join([g.to_string(header=True, comment=i) + '\n'
                        for i, g in enumerate(geometries)])

    @classmethod
    def select(cls, name, logfiles):
        '''
        From a list of :class:`gaussian.log.LOGFile`,
        select those where :attr:`gaussian.log.LOGFile.files.root`
        matches ``name`` (unix wildcards allowed).
        In case of a single find, return a single object instead of a list.
        '''
        selected = []
        for l in logfiles:
            if fnmatch.fnmatch(l.files.root, name):
                selected.append(l)
        if len(selected) == 1:
            selected = selected[0]
        return selected

    @classmethod
    def partition(cls, logs, add_patterns=[], patterns=[],
                  group_by_number=False):
        '''
        Partition a list of ``LOGFile``s according to their filename. Return
        an iterator of sets of ``LOGFile``s.
        See ``gaussian.filenames.GaussianFiles`` for more information.
        '''
        partition = filenames.GaussianFile.partition([
            l.files for l in logs], add_patterns=add_patterns,
            patterns=patterns, group_by_number=group_by_number)
        for p in partition:
            yield set([cls.select(n.root, logs) for n in p])

    @classmethod
    def set_relative_energies(cls, logs, reference='min', absolute=0,
                              relative=0, conversion='kjmol'):
        '''
        Set the relative electronic energy (``relenergy``),
        the relative enthalpy (``relenthalpy``) and the relative
        Gibbs free energy (``relgibbs``) attributes
        for each :class:`gaussian.log.LOGFile` in ``logs`` list.

        By default, the most stable system is chosen as a reference and
        the results are saved as kJ/mol. For more settings see
        :method:`gaussian.relative_energies`.
        '''
        reles = relative_energies(values=logs, reference=reference,
                                  absolute=absolute, relative=relative,
                                  conversion=conversion,
                                  type='energy')
        relhs = relative_energies(values=logs, reference=reference,
                                  absolute=absolute, relative=relative,
                                  conversion=conversion,
                                  type='enthalpy')
        relgs = relative_energies(values=logs, reference=reference,
                                  absolute=absolute, relative=relative,
                                  conversion=conversion,
                                  type='gibbs')
        for i, l in enumerate(logs):
            l.relenergy = reles[i]
            l.relenthalpy = relhs[i]
            l.relgibbs = relgs[i]

    @classmethod
    def set_rmsds(cls, logs):
        '''
        Set RMSDs (``rmsd`` attribute of :class:`gaussian.log.LOGFile`
        instance) for a number
        of Gaussian output files. First one is used as a reference.
        '''
        for l in logs:
            # Since set_rmsds() is calles in the construcor of GaussianTable
            # a try except is necessary. If molecules don't have corresponding
            # geometries, obtaining an RMSD will cause an error.
            try:
                l.rmsd = list(logs)[0].geometry.rmsd(l.geometry)[2]
            except:
                l.rmsd = float('nan')

    def __getattr__(self, name):
        '''
        Return ``None`` for missing attributes.
        '''
        return None

    # Cached attributes

    @cached
    def _full(self):
        '''
        The full text of the Gaussian output file.

            This is a large file, so all attributes have the cached
            decorator. This causes the attributes to be parsed only after
            their first request, thereby eliminating unnecessary searches
            of the long file. Note that the full text and the summary block
            (:attr:`gaussian.log.LOGFile._summary_block`) are also only
            constructed after their first request.
        '''
        try:
            return self.files.output()
        except:
            return

    @cached
    def exists(self):
        if self._full:
            return True

    @cached
    def _summary_block(self):
        '''
        Plain text of the summary block that appears at the end
        of each calculation.

            Only the first summary block of the log file is fetched. In a
            combined calculation (e.g. ``opt freq``), this block will
            contain the original route section, but evidently won't contain
            results from the second part.

            It is faster to retrieve information from this small block
            than from the full text. For several attributes (final
            energy, final geometry), the summary
            block is the first source. If this fails (abrupted calculation
            without a summary) or if more info is needed (all energies
            and geometries), the full text is searched.
        '''
        end = self._full.find('@')
        begin = self._full.rfind('GINC', 0, end)
        return self._full[begin:end].replace('\n ', '')

    @cached
    def _unspaced_summary(self):
        return self._summary_block.replace(' ', '')

    @cached
    def route_section(self):
        '''
        Route section of the calculation (after # in the input file), fully
        lowercase, # included
        '''
        if self._summary_block:
            return self._summary_block.split(r'\\')[1]
        else:
            try:
                # Find the first '#' the file.
                first_char_route_section = self._full.find('#')
                # Find the '----' before that point. This should mark the route
                # section.
                after_route_section = self._full.find('-------',
                                                      first_char_route_section)
                route_section = self._full[
                    first_char_route_section:after_route_section]
                # Join the lines in the route_section. Whitespace at the end
                # of the line is stripped. This may lead to errors when two
                # keywords are concatenated. The whole route section is
                # lowercased.
                return ''.join([x.strip()
                                for x in route_section.split('\n')]).lower()
            except:
                return None

    @cached
    def keywords(self):
        '''
        List of keywords in the route section, a simple (lowercase) split
        of the routesection, without #
        '''
        return self.route_section.replace('#', '', 1).split()

    def _options(self, keyword):
        '''
        Parse keyword(a,b=c,d=e,f) to {'a':True, 'b':'c', 'd':'e', 'f':True}
        '''
        for kw in self.keywords:
            if kw.startswith(keyword):
                options = dict()
                for option in re.split(',|\(|\)', kw)[1:]:
                    if '=' in option:
                        param, value = option.split('=')
                        options[param] = value
                    # necessary to avoid '':True entry in dictionary
                    elif option:
                        options[option] = True
                return options

    @cached
    def _charge_and_multiplicity(self):
        '''
        ``(charge, multiplicity)`` tuple used during parsing,
        since charge and multiplicity occur next to eachother
        '''
        try:
            if self._summary_block:
                return map(int, self._unspaced_summary.split(r'\\')[3]
                           .split('\\')[0].split(','))
            else:
                reg_chargemulti = re.compile(
                    r'Charge = *([-0-9]*) *Multiplicity = *(\S*)')
                chargemulti = reg_chargemulti.search(self._full,
                                                     self._full.find('Charge'))
                return map(float, (chargemulti.group(1), chargemulti.group(2)))
        except:
            return float('nan'), float('nan')

    @cached
    def charge(self):
        '''
        Total charge (e)
        '''
        return self._charge_and_multiplicity[0]

    @cached
    def multiplicity(self):
        '''
        Multiplicity (number of paired electrons + 1)
        '''
        return self._charge_and_multiplicity[1]

    @cached
    def modredundant(self):
        '''
        ModRedundant input section
        '''
        try:
            if self.scan:
                a = self._full.find('The following ModRedundant input')
                b = self._full.find('\n', a) + 2
                # Look for modredundant input in the first 1000 characters
                # after the start of the modredundant section. This ought to be
                # enough.
                modinput = []
                for line in self._full[b:b+1000].split('\n'):
                    split = line.split()
                    # Looks for lines that start correctly
                    if len(split) > 2:
                        if (split[0] in ['X', 'B', 'A', 'D', 'L']
                                and split[1].isdigit()):
                            modinput.append(line)
                return '\n'.join(modinput)
            else:
                return ''
        except:
            return ''

    @cached
    def geometries(self):
        '''
        All geometries in the log file (list of
        :class:`gaussian.molecules.SuperMolecule`)

            Geometries are printed after ``Input orientation`` (same
            cartesian coordinate system as input geometry, only if input
            in cartesian coordinates),
            ``Z-Matrix orientation`` (internal coordinates, only if input
            is provided as such) or ``Standard orientation``
            (new coordinate system) in the log file. By default, the
            input orientations are stored in
            :attr:`gaussian.log.LOGFile.geometries` and
            the standard orientations in
            :attr:`gaussian.log.LOGFile.geometries_standard`.
        '''
        # try standard orientations first
        g = self.standard_orientations
        if not g:
            # obtain input orientations if no standards were found
            g = self.input_orientations
        return g

    @cached
    def input_orientations(self):
        '''
        Return all geometries in input orientation, see
        :attr:`gaussian.log.LOGFile.geometries`.
        '''
        return self._get_geometries(orientation='Input orientation')

    @cached
    def standard_orientations(self):
        '''
        Return all geometries in standard orientation, see
        :attr:`gaussian.log.LOGFile.geometries`.
        '''
        return self._get_geometries(orientation='Standard orientation')

    @cached
    def geometry(self):
        '''
        Last geometry of the file (:class:`gaussian.molecules.SuperMolecule`).

            May be parsed as the geometry in the first summary block
            (Standard orientation?, not sure) or as the last element of
            :attr:`gaussian.log.LOGFile.geometries` (Input orientation).
            See :attr:`gaussian.log.LOGFile.geometries` for more info about
            different orientations.
        '''
        # If there is a summary block, parse last geometry from this block
        # for performance. In this way, not all geometries have to be
        # parsed, since self.geometries is only constructed when it is first
        # called.
        if self._summary_block:
            lines = self._unspaced_summary.split(r'\\')[3].split('\\')[1:]
            newblock = '\n'.join([l.replace(',', '\t') for l in lines])
            geom = molecules.SuperMolecule.from_string(newblock)
            geom.scale(self.nimag, factor=self.scalings)
            return geom
        # If not, fetch last geometry from self.geometries.
        else:
            logging.debug('log.geometry: no summary block found')
            return self.geometries[-1]

    @cached
    def energies(self):
        '''
        List of electronic energies (Ha), which can be found
        after each occurence of 'SCF Done'
        '''
        try:
            if 'SCF Done' in self._full:
                return [self._get_energy(p)
                        for p in find_all(self._full, 'SCF Done')]
            else:
                return [self._get_energy(p)
                        for p in find_all(self._full, 'Energy=')]
        except:
            return None

    @cached
    def steps(self):
        '''
        Number of geometrical steps
        '''
        try:
            return len(self.energies)
        except:
            return float('nan')

    @cached
    def energy(self):
        '''
        Final electronic energy (Ha)
        '''
        # If there is a summary block, parse last electronic energy
        # from this block for performance. In this way,
        # not all electronic energies have to be parsed, since
        # self.energies is only constructed when it is first
        # called.
        if self._summary_block:
            try:
                start = self._unspaced_summary.find('HF') + 3
                end = self._unspaced_summary.find('\\', start)
                return float(self._unspaced_summary[start:end])
            # If there is a summary block, but this approach fails,
            # try to fetch all electronic energies and pick last one.
            except:
                return self.energies[-1]
        # Do the same thing if no summary block is present.
        else:
            return self.energies[-1]

    @cached
    def memory(self):
        '''
        Requested memory (e.g. 1MB)
        '''
        # The link0 section is repeated literally at the top
        # of the log file.

        # Do not use sfind. If '%mem' is not found, '=' will still be found
        # after position -1. This will lead to nonsense if no memory limit was
        # provided in the input file.

        mem = self._full.find('%mem')
        if mem != -1:
            start = self._full.find('=', mem) + 1
            end = self._full.find('\n', start)
            return self._full[start:end]

    @cached
    def _thermochem_block(self):
        '''
        Plain text block containing the thermal corrections to the
        electronic energy
        '''
        if 'freq' in self.route_section:
            start = self._full.find('Vibrational temperatures')
            end = self._full.find('E (Thermal)')
            return self._full[start:end]

    @cached
    def zpe(self):
        '''
        Zero point correction (Ha)
        '''
        return self._get_energy_correction('Zero-point correction')

    @cached
    def thermalcorrection(self):
        '''
        Thermal correction to Energy (Ha)
            =  Zero-point correction
                    +  translational energy
                    +  vibrational energy
                    +  rotational energy
        '''
        return self._get_energy_correction('Thermal correction to Energy')

    @cached
    def enthalpycorrection(self):
        '''
        Thermal correction to Enthalpy (Ha)
            =  Thermal correction to Energy
                    +   PV
            =   Zero-point correction
                    +   translational energy
                    +   vibrational energy
                    +   rotational energy
                    +   PV
        '''
        return self._get_energy_correction('Thermal correction to Enthalpy')

    @cached
    def gibbscorrection(self):
        '''
        Thermal correction to Gibbs Free Energy (Ha)
            =  Thermal correction to Enthalpy
                    -   TS
            =   Zero-point correction
                    +   translational energy
                    +   vibrational energy
                    +   rotational energy
                    +   PV
                    +   (-TS)
        '''
        return self._get_energy_correction(
            'Thermal correction to Gibbs Free Energy')

    @cached
    def zpesum(self):
        '''
        Electronic energy + ZPE (Ha)
        '''
        return self._get_energy_correction(
            'Sum of electronic and zero-point Energies')

    @cached
    def thermal(self):
        '''
        Electronic energy + Thermal correction to Energy
        '''
        return self._get_energy_correction(
            'Sum of electronic and thermal Energies')

    @cached
    def enthalpy(self):
        '''
        Enthalpy (Ha)
            = Electronic energy + Thermal correction to Enthalpy
        '''
        return self._get_energy_correction(
            'Sum of electronic and thermal Enthalpies')

    @cached
    def gibbs(self):
        '''
        Gibbs free energy (Ha)
            = Electronic energy + Thermal correction to Free Energy
        '''
        return self._get_energy_correction(
            'Sum of electronic and thermal Free Energies')

    @cached
    def entropy(self):
        '''
        Entropy (Ha), calculated as (1/T) * (H - G)
        '''
        return (1.0 / self.temperature) * (self.enthalpy - self.gibbs)

    @cached
    def _temperature_and_pressure(self):
        '''
        (temperature (K), pressure (atm)) tuple, used in parsing
        since they appear in the same place
        '''
        if 'freq' in self.keywords:
            try:
                start = self._full.find('Temperature',
                                        self._full.find('Thermochemistry'))
                end = self._full.find('\n', start)
                split = self._full[start:end].split()
                temperature = float(split[1])
                pressure = float(split[4])
            except:
                temperature = float('nan')
                pressure = float('nan')
        else:
            temperature = float('nan')
            pressure = float('nan')
        return temperature, pressure

    @cached
    def temperature(self):
        '''
        Temperature (K) used during frequency calculations
        '''
        return self._temperature_and_pressure[0]

    @cached
    def pressure(self):
        '''
        Pressure (atm) used during frequency calculations
        '''
        return self._temperature_and_pressure[1]

    @cached
    def irc(self):
        '''
        :class:`gaussian.log.IRCPath` instance that contains relevant
        IRC energies and geometries. ``None`` if this isn't a IRC
        calculation.
        '''
        if 'irc' in self.route_section:
            ircpath = IRCPath(self)
            ircpath.maxcycle = self._options('irc').get('maxcycle', 20)
            ircpath.stepsize = self._options('irc').get('stepsize', 10)
            ircpath.maxpoints = self._options('irc').get('maxpoints', 10)
            return ircpath

    @cached
    def scan(self):
        '''
        :class:`gaussian.log.ScanPath` instance that contains relevant
        scan energies and geometries. ``None`` if this isn't a scan.
        '''
        if 'modredundant' in self.route_section:
            return ScanPath(self)

    @cached
    def frequencies(self):
        '''
        Vibrational frequencies (cm-1).
        '''
        frequencies = None

        if 'freq' in self.route_section:
            try:
                # List of frequencies
                freq_lines_start = list(find_all(
                    self._full, 'Frequencies --'))
                freq_lines_end = [self._full.find('\n', i)
                                  for i in freq_lines_start]
                freq_lines = [self._full[a:freq_lines_end[i]]
                              for i, a in enumerate(freq_lines_start)]
                frequencies = []
                for freq_line in freq_lines:
                    frequencies.extend((freq_line.split()[2:]))
                frequencies = map(float, frequencies)
                # If opt(calcall) is specified, frequencies are reported in the
                # optimisation section. If another, superfluous frequency
                # calculation was performed afterwards, discard it.
                if ('calcall' in self.route_section
                        and 'freq' in self.route_section):
                    frequencies = frequencies[:len(frequencies) / 2]
            except:
                pass
        return frequencies

    @cached
    def nimag(self):
        '''
        Number of imaginary frequencies
        '''
        try:
            return len([x for x in self.frequencies if x < 0])
        except:
            return float('nan')

    @cached
    def lowest_frequency(self):
        '''
        Lowest (least positive, most negative) vibrational frequency (cm-1)
        '''
        try:
            return self.frequencies[0]
        except:
            return float('nan')

    @cached
    def chk_consistent(self):
        '''
        Return ``True`` if the checkpoint file specified in the input file
        corresponds to the name of the output file (e.g. `calculation.log` and
        `calculation.chk`).
        '''
        return self.chk == self.files.chk

    @cached
    def npa_charges(self):
        '''
        NPA charges (list of atomic charges in, e)
        '''

        try:
            # The block in which the NPA charges can be found.
            start = self._full.find('---\n',
                                    self._full.find(
                                        'Total', self._full.rfind(
                                            'Summary of Natural Popul'
                                            + 'ation Analysis'))) + 4
            end = self._full.find('==', start) - 2

            # The third column in this block.
            return np.array([float(l.split()[2]) for l
                             in self._full[start:end].split('\n')])
        except:
            return None

    @cached
    def mulliken_charges(self):
        '''
        Mulliken charges (list of atomic charges, e)
        '''

        try:
            # Look for instances of 'Mulliken', start at the bottom. These
            # are clustered together. If a large jump in position is
            # encountered, i.e. jumping between clusters, stop searching.
            # The top two instances of 'Mulliken' in the cluster at the
            # bottom of the file mark the atomic Mulliken charges.
            finder = find_all(self._full, 'Mulliken', reverse=True)
            finds = [next(finder)]
            diff = 0
            while abs(diff) < 5000:
                finds.append(next(finder))
                diff = finds[-1] - finds[-2]
            lines = self._full[finds[-2]:finds[-3]].split('\n')[2:-1]
            return [float(l.split()[2]) for l in lines]
        except:
            return None

    @cached
    def error(self):
        '''
        Termination status of the calculation

            A normal termination will result in ``None``. An error termination
            will result in either one of the keywords below, depending on
            which error message can be found in the log file.

            ==================  =============================================
            Keyword             Error message
            ==================  =============================================
            ``convergence``     Convergence failure
            ``corrector``       Maximum number of corrector steps
            ``formbx``          FormBX had a problem
            ``intcoord``        Error in internal coordinate system
            ``multi/charge``    The combination of multiplicity
            ``ntrerr``          NtrErr Called from FileIO
            ``steps``           Number of steps exceeded
            ``svd``             SVD failed
            ``syntax``          A syntax error was detected in the input line
            ``unknown``         The error could not be assigned
            ==================  =============================================
        '''
        possible_errors = {'The combination of multiplicity': 'multi/charge',
                           'Convergence failure': 'convergence',
                           'Maximum number of corrector steps': 'corrector',
                           'Number of steps exceeded': 'steps',
                           'Error in internal coordinate system': 'intcoord',
                           'SVD failed': 'svd',
                           'FormBX had a problem': 'formbx',
                           'NtrErr Called from FileIO': 'ntrerr',
                           'A syntax error was detected in the input line':
                           'syntax',
                           'FileIO operation on non-existent file.':
                           'nochk'}

        # even a Gaussian calculation?
        if not self._full.startswith(' Entering Gaussian System'):
            return 'notgaussian'

        # bottom of the log file contains termination information
        bottom = self._full[-len(self._full)/2:]

        # terminated, with or without errors
        if 'termination' in bottom:

            # test for all possible errors
            for possible_error in possible_errors:
                if possible_error in bottom:
                    return possible_errors[possible_error]

            # All okay?
            # Last 200 characters should include 'Normal termination'.
            if 'Normal termination' in bottom[-200:]:
                return

            # terminated, but an unknown error
            return 'unknown'

        # did not terminate
        else:

            # 'Argument list too long' error will stop the calculation
            # abruptly. This will only show in the standard error file.
            stderr = self.files.get_stderr()
            if stderr:
                if 'Argument list too long' in open(stderr).read():
                    return 'argumentlist'
            return 'abrupted'

    @cached
    def termination_count(self):
        '''
        Number of 'termination's
        '''
        return len(list(find_all(self._full, 'termination')))

    @cached
    def stoichiometry(self):
        '''
        Stoichiometry
        '''
        try:
            return self.geometry.stoichiometry
        except:
            return None

    @cached
    def nproc(self):
        '''
        Number of processors used during a Gaussian calculation
        '''
        try:
            nproc_index = self._full.find('%nproc')
            return int(re.search('%nproc=([0-9]*)', self._full[
                nproc_index:nproc_index + 100]).group(1))
        except:
            return float('nan')

    @cached
    def predicted(self):
        '''
        Predicted energy changes in Hartree
        '''
        try:
            predicted = []
            for i in find_all(self._full, 'Predicted'):
                begin = self._full.find('=', i) + 1
                end = self._full.find('\n', begin)
                predicted.append(
                    float(self._full[begin:end].replace('D', 'E')))
            return predicted
        except:
            return []

    @cached
    def predicted_last(self):
        '''
        Logarithm (base 10) of the absolute
        value of the last predicted energy change
        '''
        try:
            return round(log10(abs(self.predicted[-1])), 1)
        except:
            return float('nan')

    @cached
    def comment(self):
        '''
        Comment line provided in input
        '''

        try:
            first_character_route_section = self._full.find('#')
            after_route_section = self._full.find(
                '-------', first_character_route_section)
            just_after_route_block = self._full.find(
                '------\n', after_route_section) + 10
            begin_comment = self._full.find('--\n', just_after_route_block) + 3
            end_comment = self._full.find('--', begin_comment + 1)
            return self._full[begin_comment:end_comment].strip()
        except:
            return None

    @cached
    def chk(self):
        '''
        Checkpoint filename from Gaussian output file.
        '''
        try:
            begin_checkpoint = self._full.find('chk=') + 4
            end_checkpoint = self._full.find('\n', begin_checkpoint)
            checkpoint = self._full[begin_checkpoint:end_checkpoint].strip()
            return filenames.GaussianFile(checkpoint).root
        except:
            return None

    @cached
    def _method_and_basis(self):
        '''
        (method, basis) tuple, used internally for parsing purposes because
        the two appear in the same place.
        '''
        method = False
        basis = False

        # Possible aliases for basis sets. The first element of every tuple
        # will be substituted by the second. List multiple tuples in order
        # of replacement, since certain aliases may overlap, e.g. g* and g**.
        basis_aliases = {'g**': 'g(d,p)', 'g*': 'g(d)'}
        try:
            method, basis = [x for x in self.keywords
                             if re.search('^[^\(]+/', x)
                             ][0].split('/')
            for b in basis_aliases:
                basis = basis.replace(b, basis_aliases[b])
            lot = '%s/%s' % (method, basis)
        except:
            # If no 'method/basis' was found, assume that the first word
            # of the route section is the method. This is the casef for
            # semiempirical methods and forcefields.
            basis = None
            method = lot = self.keywords[0]
        return method, basis, lot

    @cached
    def method(self):
        '''
        Method (Gaussian keyword: DFT functional, semiempirical method, etc.)
        '''
        return self._method_and_basis[0]

    @cached
    def basis(self):
        '''
        Basis set (Gaussian keyword). ``None`` for semiempirical methods.
        '''
        return self._method_and_basis[1]

    @cached
    def lot(self):
        '''
        Method/basis for ab initio, method for semiempirics
        '''
        return self._method_and_basis[2]

    @cached
    def dipole(self):
        '''
        Total dipole moment (field-independent basis, Debye)
        '''
        try:
            start = self._full.find('Tot=',
                                    self._full.find('Dipole moment')) + 4
            end = self._full.find('\n', start)
            return float(self._full[start:end])
        except:
            return float('nan')

    @cached
    def scrf(self):
        '''
        SCRF part (implicit solvation) in the route section or empty string
        '''
        scrf = re.search(r'(scrf\S*)', self.route_section)
        if not scrf:
            scrf = ''
        else:
            scrf = scrf.group(1)
        return scrf

    @cached
    def scrf_nonstandard_input(self):
        '''
        Nonstandard SCRF block in input file or empty string
        '''
        try:
            start = self._full.find(
                'Using the following non-standard input for PCM:') + 48
            end = self._full.find('--- end of non-standard input.', start)
            if not start == 47 and not end == -1:
                return self._full[start:end].strip()
            else:
                return ''
        except:
            return ''

    @cached
    def cputime(self):
        '''
        CPU time
        '''
        try:
            begin = self._full.find('Job cpu time') + 13
            end = self._full.find('\n', begin)
            formatted_time = map(float,
                                 self._full[begin:end].strip().split()[::2])
            return sum(formatted_time * np.array([24*3600, 3600, 60, 1]))
        except:
            return None

    @cached
    def normal_coordinates(self):
        '''Returns normal coordinates for three lowest
        vibrational frequencies'''
        # Do not check for 'freq' in the route section. Some people use fake
        # Gaussian output files, with frequencies, but without freq.
        try:
            start = self._full.find('Frequencies')
            end = self._full.find('Frequencies', start+1)
            block = self._full[start:end].split('\n')[5:-3]
            return [np.array([l.split()[2:5] for l in block], np.float64),
                    np.array([l.split()[5:8] for l in block], np.float64),
                    np.array([l.split()[8:11] for l in block], np.float64)]
        except:
            return None

    @cached
    def normal_lengths(self):
        '''Norms of the normal coordinates for lowest vibrational freq'''
        return np.linalg.norm(self.normal_coordinates, axis=1)

    @cached
    def reacting(self, treshold=[1.2, 2.0], only=None):
        '''
        return a list of bonds (i.e. sets of atom pairs) that could be
        considered to be reacting in the lowest frequency. In this context
        two atoms react if:
            - force vectors for this frequency are almost parallel
            - distance between the atoms is between treshold[0] and
              treshold[1] times the VDW radius for these atoms
        search can be limited to a subset of atoms via only
        '''
        logging.debug('log.LOGFile.reacting')
        # limit search to these atoms
        if not only:
            only = range(self.geometry.size)
        # normalize all force vectors
        normal = []
        for vec in self.normal_coordinates[0]:
            denom = np.linalg.norm(vec) or 1
            normal.append(vec / denom)
        reacting = []
        # for all atom pairs
        for i, ni in enumerate(normal):
            for j, nj in enumerate(normal):
                if i < j and i in only and j in only:
                    # if force the vectors are parallel
                    if abs(np.dot(ni, nj)) > 0.90:
                        # VDW distance between atom i and atom j
                        vdw = bonds.bonds.get_length(
                            self.geometry.numbers[i],
                            self.geometry.numbers[j]) / units.angstrom
                        # real distance
                        real = self.geometry.dist(i, j)
                        # real distance should be quite close to VDW distance
                        if treshold[0] < real / vdw < treshold[1]:
                            reacting.append(set([i, j]))
        logging.debug('reacting atoms in %s are %s' % (self.file, reacting))
        return reacting


class Path(object):
    '''
    Generic path class for scan and IRC paths.

        Needs ``_point_positions``
        and ``_logfile`` (:class:`gaussian.log.LOGFile`)
        attributes to be functional. At each of the point positions, the
        electronic energy and geometry above that position will be selected
        as a point. The parsing of ``_point_positions`` determines the type
        of path.
    '''

    @cached
    def geometries(self):
        '''
        Geometries of the path points
        (list of :class:`gaussian.molecules.SuperMolecule`,
        ``Input orientations``)
        '''
        return [self._logfile._get_geometry(x, reverse=True)
                for x in self._point_positions]

    @cached
    def energies(self):
        '''
        Electronic energies of the path points (Ha)
        '''
        return [self._logfile._get_energy(x, reverse=True)
                for x in self._point_positions]

    @cached
    def relativeenergies(self):
        '''
        Relative electronic energies of the path points (kJ/mol)
        '''
        return relative_energies(self.energies)

    @cached
    def length(self):
        '''
        Number of points
        '''
        return len(self.energies)


class IRCPath(Path):
    '''
    IRC path derived from :class:`gaussian.log.Path`
    '''

    def __init__(self, logfile):
        self._logfile = logfile

    @cached
    def _point_positions(self):
        '''
        Fetch _point_positions (see parent class :class:`gaussian.log.Path`)
        '''
        # Find all 'Point Number: x    Path Number: 1'.
        path1 = list(find_all(self._logfile._full, 'Path Number:   1'))
        path2 = list(find_all(self._logfile._full, 'Path Number:   2'))
        return path2[::-1] + path1

    @cached
    def direction(self):
        '''
        Direction of the IRC path (``'both'``, ``'forward'`` or ``'reverse'``)
        '''
        if 'both directions' in self._logfile._full:
            return 'both'
        elif 'path following direction' in self._logfile._full:
            if ('fwd' in self._logfile.route_section
                    or 'forward' in self._logfile.route_section):
                return 'forward'
            elif ('rev' in self._logfile.route_section
                    or 'reverse' in self._logfile.route_section):
                return 'reverse'


class ScanPath(Path):
    '''
    Scan path derived from :class:`gaussian.log.Path`
    '''

    def __init__(self, logfile):
        self._logfile = logfile

    @cached
    def _point_positions(self):
        '''
        Fetch _point_positions (see parent class :class:`gaussian.log.Path`)
        '''
        return list(find_all(self._logfile._full, 'Stationary point found'))
