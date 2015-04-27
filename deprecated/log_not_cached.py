'''Parsing Gaussian output files'''
import re
import warnings
import filenames
import molecules
import fnmatch
import utils
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    import numpy as np
    import molmod.units as units
    import molmod.molecules as mol
    # TODO make cached attributes
    # from molmod.utils import cached
# TODO IRC and scan: make a path class with geometries, energies, relative
# energies


class LOGFile(object):
    '''Parsed Gaussian output file from ``filename``'''

    def __init__(self, filename):
        self.files = filenames.GaussianFile(filename)
        '''
        ``GaussianFile`` object constructed from ``filename``
        '''
        self.file = self.files.root
        '''
        Root of filename (filename without extension, with directory),
        equivalent with ``self.files.root``.
        '''
        self._full = open(filename).read()
        '''Full text of the Gaussian output file'''
        self.error = self._get_error()
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


            Typically, only normally terminated files are processed. For
            example::
                l = LOGFile(f)
                if not l.error:
                    # do something

        '''
        route_section, keywords = self._get_route_section()
        self.route_section = route_section
        '''
        Route section of the calculation (after # in the input file), fully
        lowercase

            The following code prints the filenames for which a frequency
            calculation was requested::
                for f in filenames:
                    l = LOGFile(f)
                    if 'freq' in l.route_section:
                        print l

        '''
        self.keywords = keywords
        '''
        Keywords in route section, lowercase, a simple split of
        :attr:`gaussian.log.LOGFile.route_section`
        '''
        self.geometries = self._get_geometries(orientation='Input orientation')
        '''
        All geometries in the log file (list of
        :class:`gaussian.molecules.SuperMolecule`)

            Geometries are printed after ``Input orientation`` (same
            orientation as input geometry) or ``Standard orientation``
            (new coordinate system) in the log file. By default, the
            input orientations are stored in ``self.geometries``.

            This code checks whether the number of molecules is the same
            in the first and the last geometry of the calculation::
                l = LOGFile(f)
                nbegin = len(l.geometries[0].molecules())
                nend = len(l.geometries[-1].molecules())
                assert nbegin == nend

        '''
        try:
            self.geometry = self.geometries[-1]
            '''
            Last geometry (:class:`gaussian.molecules.SuperMolecule`), i.e.
            the last element of :attr:`gaussian.log.LOGFile.geometries`
            '''
            self.size = self.geometry.size
            '''
            Number of atoms in the systems, equivalent to
            ``self.geometry.size``
            '''
        except:
            self.geometry = None
            self.size = float('nan')
        charge, multi = self._get_charge_and_multiplicity()
        self.charge = charge
        '''Charge of system (e)'''
        self.multiplicity = multi
        '''Multiplicity of the systems, as defined in input file'''
        self.stoichiometry = self._get_stoichiometry()
        '''Stoichiometry, e.g. ``'C18H27O(1+)'``'''
        self.nproc = self._get_nproc()
        '''Number of processors used for the calculation'''
        self.comment = self._get_comment()
        '''Comment line from Gaussian input file'''
        self.chk = self._get_chk()
        '''
        Name of the checkpoint file (without .chk extension) which
        was given in the Gaussian
        input file and used during the calculation.

            Generally, the checkpoint which was used (``self.chk + '.chk'``)
            will be equal to the name of the output file in which the log
            extension is replaced by chk (``self.files.chk``).
            In some cases, it may be wise to check whether this was actually
            the case, before the further use of (formatted) checkpoint files::
                l = LOGFile(f)
                assert l.chk == l.files.chk
        '''
        method, basis, lot = self._get_method_basis()
        self.method = method
        '''Method (DFT functional, MP2, etc.)'''
        self.basis = basis
        '''Basis set'''
        self.lot = lot
        '''
        Level of theory, combination of method and basis, e.g. B3LYP/6-31+G**,
        MP2/aug-cc-pVDZ, PM6
        '''
        self.dipole = self._get_dipole()
        '''Dipole moment (Debye)'''
        scrf, scrf_nonstandard_input = self._get_scrf()
        self.scrf = scrf
        '''SCRF part in route section'''
        self.scrf_nonstandard_input = scrf_nonstandard_input
        '''Nonstandard SCRF information provided in Gaussian input file'''
        self.electronics = self._get_energies()
        '''Electronic energies after each 'SCF Done' (Ha)'''
        try:
            self.electronic = self.electronics[-1]
            '''Last electronic energy (Ha), i.e. last element of
            :attr:`gaussian.log.LOGFile.electronics`'''
            self.steps = len(self.electronics)
            '''Number of electronic steps'''
        except:
            self.electronic = float('nan')
            self.steps = float('nan')
        self.gibbs = float('nan')
        '''Gibbs free energy = electronic energy + Gibbs correction (Ha)'''
        self.enthalpy = float('nan')
        '''Enthalpy = electronic energy + enthalpy correction (Ha)'''
        self.zpe = float('nan')
        '''Zero-point energy (Ha)'''
        self.thermalcorrection = float('nan')
        '''Thermal correction (Ha)'''
        self.enthalpycorrection = float('nan')
        '''Enthalpy correction (Ha)'''
        self.gibbscorrection = float('nan')
        ''' Gibbs correction (Ha)'''
        self.zpesum = float('nan')
        '''Electronic energy + zero-point energy (Ha)'''
        self.thermal = float('nan')
        '''Electronic energy + thermal correction (Ha)'''
        self._set_energy_corrections()
        temp, pres = self._get_temperature_and_pressure()
        self.temperature = temp
        '''Temperature of energy corrections (K)'''
        self.pressure = pres
        '''Pressure of energy corrections (atm)'''
        freqs, lowfreq, nimag = self._get_frequencies()
        self.frequencies = freqs
        '''List of vibrational frequencies (cm-1)'''
        self.lowest_frequency = lowfreq
        '''Smallest or most negative vibrational frequency (cm-1)'''
        self.nimag = nimag
        '''Number of imaginary vibrational frequencies'''
        irc, ircg, irce, ircrele, irclen, ircdir = self._get_irc()
        self.irc = irc
        '''Is this an IRC calculation (``bool``)?'''
        self.ircgeometries = ircg
        '''IRC geometries (list of ``SuperMolecule``)'''
        self.ircenergies = irce
        '''IRC electronic energies (Ha)'''
        self.ircrelativeenergies = ircrele
        '''IRC relative electronic energies (kJ/mol)'''
        self.irclength = irclen
        '''Length of IRC path'''
        self.ircdirection = ircdir
        '''Direction of IRC path (``forward``, ``reverse``, ``both``)'''
        scan, scang, scane = self._get_scan()
        self.scan = scan
        '''Is this a scan (``bool``)?'''
        self.scangeometries = scang
        '''Scan geometries (list of ``SuperMolecule``)'''
        self.scanenergies = scane
        '''Scan energies (Ha)'''
        self.cputime = self._get_cputime()
        '''CPU time (s)'''
        self.memory = self._get_memory()
        '''Requested memory (``str``, e.g. 1MB)'''
        # Default scaling behaviour for molecular graphs of TSs
        if self.nimag == 1:
            self.geometry.scaling = 1.5
            # TODO default scaling behaviour help in SuperMolecule

    def to_input(self, filename=None):
        '''
        Convert to Gaussian input file. ``filename`` is used for the checkpoint
        file. If ``None``, use the current filename.
        '''
        # TODO make a Gaussian input file class, so that you don't have to care
        # about the change in self.files.

        fn = filenames.GaussianFile(filename) if filename else self.files
        # This way, self.geometry = None cn be set, and no cartesian
        # coordinates are written. Can be handy if calculations are restarted
        # from a checkpoint file.
        try:
            geometry = self.geometry.to_string()
        except:
            geometry = ''
        return '''%%nproc=%(nproc)s
%%mem=%(mem)s
%%chk=%(chk)s
#%(route)s

%(comment)s

%(charge)s %(multiplicity)s
%(geometry_text)s

%(scrf)s
''' % {'nproc': self.nproc,
            'mem': self.memory,
            'chk': fn.base,
            'route': ' '.join(self.keywords),
            'comment': self.comment,
            'charge': self.charge,
            'multiplicity': self.multiplicity,
            'geometry_text': geometry,
            'scrf': self.scrf_nonstandard_input}

    def write_input(self, filename=None):
        '''
        Write a Gaussian input file with the data parsed from the output
        file. This is most likely a restart file. Use ``filename`` for
        the checkpoint file and the input file (.com). If ``None``, use the
        current filename.
        '''
        fn = filenames.GaussianFile(filename) if filename else self.files
        f = open(fn.com, 'w')
        f.write(self.to_input(fn.base))
        f.close()

    def __getattr__(self, name):
        '''
        Return ``None`` for missing attributes.
        '''
        return None

    def _get_geometry(self, pos, reverse=False, orientation='orientation'):
        '''
        Return the first geometry in the Gaussian output file after
        (``Reverse=False``) or before (``Reverse=True``) position ``pos``
        as a ``SuperMolecule``.

        A geometry is marked by 'orientation' in the output file. Both Standard
        and Input Orientations can be found in most output files. The standard
        implementation of this function parses all of them. Reset orienation
        (``Standard orientation``, ``Input orientation`` if only a single
        type is wanted.
        '''

        # The regex used for obtaining an unprocessed geometry
        reg = re.compile((r'%s:.*$\n^ -*$\n^.*$\n^.*$\n^ -*'
                          + r'$\n([-\s\.0-9\n]*[0-9]\n)') % orientation,
                         re.MULTILINE)

        # Obtain unprocessed geometry. Stops searching if found.
        if reverse:
            actual_pos = self._full.rfind('orientation', 0, pos) - 5
        else:
            actual_pos = pos
        unprocessed_geometry = reg.search(self._full, actual_pos).group(1)[:-1]

        lines = unprocessed_geometry.split('\n')
        atoms = [l.split()[1] for l in lines]
        coordinates = [np.array(map(float, l.split()[-3:])) * units.angstrom
                       for l in lines]
        geometry = mol.Molecule(atoms, coordinates)

        return molecules.SuperMolecule.from_Molecule(geometry)

    def _get_geometries(self, orientation='orientation'):
        '''
        Parse all geometries from the Gaussian output file as a list of
        ``SuperMolecule``.

        See ``_get_geometry`` for the ``orientation`` argument.
        '''
        try:
            return [self._get_geometry(p, orientation=orientation)
                    for p in utils.find_all(self._full, orientation)]
        except:
            return None

    def _get_energy(self, pos, reverse=False):
        '''
        Get the first electronic energy (Ha) before (``reverse=True``)  or
        after (``reverse=False``) a given position ``pos`` in the Gaussian
        output file.
        '''

        if 'SCF Done' in self._full:
            reg_energy = re.compile(r'SCF Done: *E.* = *([-.0-9]*)')
            start_string = 'SCF Done'
        else:
            reg_energy = re.compile(r'Energy= *([-.0-9]*)')
            start_string = 'NIter='
        if reverse:
            actual_pos = self._full.rfind(start_string, 0, pos) - 50
        else:
            actual_pos = self._full.find(start_string, pos) - 50
        found_energy = reg_energy.search(self._full, actual_pos)
        energy = found_energy.group(1)
        return float(energy)

    def _get_energies(self):
        '''
        Get all electronic energies (Ha) in the Gaussian output file.
        '''
        try:
            if 'SCF Done' in self._full:
                return [self._get_energy(p)
                        for p in utils.find_all(self._full, 'SCF Done')]
            else:
                return [self._get_energy(p)
                        for p in utils.find_all(self._full, 'Energy=')]
        except:
            return None

    def _get_route_section(self):
        '''
        Parse the route section from the Gaussian output file and the keywords
        therein.

        A tuple of the full route section (``str``) and a splitted one
        (list of ``str``) is returned.
        '''
        try:
            # Find the first '#' the file.
            first_char_route_section = self._full.find('#')
            # Find the '----' before that point. This should mark the route
            # section.
            after_route_section = self._full.find('-------',
                                                  first_char_route_section)
            route_section = self._full[first_char_route_section
                                       + 1:after_route_section]
            # Join the lines in the route_section. Whitespace at the end of the
            # line is stripped. This may lead to errors when two keywords are
            # concatenated. The whole route section is lowercased.
            route_section = ''.join([x.strip()
                                    for x in route_section.split(
                                        '\n')]).lower()
            # A list of keywords by splitting the route section.
            keywords = route_section.split()
            return route_section, keywords
        except:
            return None, None

    def _get_memory(self):
        '''
        Return the memory requested for the calculation.
        '''
        a = self._full.find('%mem')
        a = self._full.find('=', a) + 1
        b = self._full.find('\n', a)
        return self._full[a:b]

    def _set_energy_corrections(self):
        '''
        Parse all thermal energy corrections (Ha) and reset the attributes
        of the ``LOGFile`` accordingly.

        The following attributes will be set: zpe, thermalcorrection,
        enthalpycorrection, gibbscorrection, zpesum, thermal, enthalpy, gibbs.
        '''
        names = ['zpe', 'thermalcorrection', 'enthalpycorrection',
                 'gibbscorrection', 'zpesum', 'thermal', 'enthalpy', 'gibbs']
        strings = ['Zero-point correction', 'Thermal correction to Energy',
                   'Thermal correction to Enthalpy',
                   'Thermal correction to Gibbs Free Energy',
                   'Sum of electronic and zero-point Energies',
                   'Sum of electronic and thermal Energies',
                   'Sum of electronic and thermal Enthalpies',
                   'Sum of electronic and thermal Free Energies']
        if 'freq' in self.route_section:
            start = self._full.find('Vibrational temperatures')
            end = self._full.find('E (Thermal)')
            ecorr_section = self._full[start:end]
            for name, string in zip(names, strings):
                try:
                    e = float(re.search(string + r'= *([-.0-9]*)',
                                        ecorr_section).group(1))
                    setattr(self, name, e)
                except:
                    setattr(self, name, float('nan'))

    def _get_temperature_and_pressure(self):
        '''
        Return the temperature (K) and the pressure (atm) used in the
        frequency calculations.

        ``temperature, pressure`` is returned.
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

    def _get_irc_path(self, path):
        '''
        Obtain the geometries and the energies of an IRC path.

        Arguments:
            logfile:        filename or full text
            path:           path number (1,2)
        Return:
            {'geometries': list of molmod.molecules.Molecule,
            'energies': list of floats}
        '''

        geometries = []
        energies = []

        # Find all 'Point Number: x    Path Number: 1'.
        # Store in MatchObject foundPointNumbersPath.
        foundPointNumbersPath = re.finditer('Point Number.*Path Number:   '
                                            + str(path) + '', self._full)

        # Iterate over these findings.
        for foundPointNumberPath in foundPointNumbersPath:

            # Find the 'Input orientation' right before 'Point Number' for the
            # points in path 1. This is the geometry of the point. Store index
            # at which this orientation is found in InputOrientation.
            foundOrientation = self._full.rfind('orientation', 0,
                                                foundPointNumberPath.start())

            # Obtain that very 'Input orientation' (i.e. the corresponding
            # geometry)
            ProcessedGeometry = self._get_geometry(foundOrientation)

            # Obtain the corresponding electronic energy, i.e. the SCF Done
            # after 'Input orientation'
            SCFDoneEnergy = self._get_energy(foundOrientation)

            # Add geometry and energies to results
            geometries.append(ProcessedGeometry)
            energies.append(SCFDoneEnergy)

        # Return geometries and energies as a dictionary
        return {'geometries': geometries, 'energies': energies}

    def _get_irc(self):
        '''
        Get all IRC data.

        Arguments:
            logfile:        filename or full text
            route:          provide parsed route option (optionally),
                            increases efficiency (default: None)
        Return:
            {'ircenergies': list of floats, 'ircrelativeenergies':
                list of floats, 'ircgeometries': list of molmod.molecules.
                Molecule, 'irclength': int, 'ircdirection': 'forward',
                'reverse, 'both'}
        '''

        irc = False
        ircenergies = None
        ircgeometries = None
        ircrelativeenergies = None
        irclength = float('nan')
        ircdirection = None

        if 'irc' in self.route_section:
            irc = True
            try:
                # Remark: IRC analysis should work if freq is erroneously used
                # after the irc
                # Important remark: only 'Input orientation's are stored.
                # Some geometries and electronic energies should be present.

                # Obtain energies and geometries from IRC paths. If only one
                # (reverse, forward) path was calculated, path2 will be empty.
                path1 = self._get_irc_path(1)
                path2 = self._get_irc_path(2)

                # Merge the energies and geometries of both paths.
                ircenergies = path2['energies'][::-1] + path1['energies']
                ircgeometries = (path2['geometries'][::-1]
                                 + path1['geometries'])

                # Relative IRC energies
                ircrelativeenergies = utils.relative_energies(
                    ircenergies)
                # Number of IRC path points
                irclength = len(ircenergies)

                # Combined forward-reverse or single path
                if 'both directions' in self._full:
                    ircdirection = 'both'
                elif 'path following direction' in self._full:
                    if ('fwd' in self.route_section
                            or 'forward' in self.route_section):
                        ircdirection = 'forward'
                    elif ('rev' in self.route_section
                            or 'reverse' in self.route_section):
                        ircdirection = 'reverse'
            except:
                pass

        return (irc, ircgeometries, ircenergies, ircrelativeenergies,
                irclength, ircdirection)

    def _get_scan(self):
        '''
        Get all geometries and energies from a PES scan (modredundant).

        Arguments:
            logfile:        filename or full text
            route_section:  provide parsed route option (optionally),
                            increases efficiency (default: None)
        Return:
            return {'geometries': list of molmod.molecules.Molecule,
            'energies': list of floats}
        '''

        if 'modredundant' in self.route_section:
            scan = True
            try:
                # Different way of handling output if a FF was used.
                if 'SCF Done' in self._full:
                    force_field = False
                else:
                    force_field = True
                points = [i for i in utils.find_all(self._full,
                                                    'Stationary point found')]
                scangeometries = [self._get_geometry(i, reverse=True)
                                  for i in points]
                scanenergies = [self._get_energy(
                    i, reverse=True, force_field=force_field)
                    for i in points]
            except:
                scangeometries = None
                scanenergies = None
        else:
            scan = False
            scangeometries = None
            scanenergies = None
        return scan, scangeometries, scanenergies

    def _get_frequencies(self):
        '''
        Obtain vibrational frequencies.

        Arguments:
            logfile:        filename or full text
            route_section:  provide parsed route option (optionally),
                            increases efficiency (default: None)
        Return:
            {'frequencies': None or list of floats, 'nimag': None or int,
            'lowest_frequency': None or float}
        '''
        frequencies = None
        nimag = None
        lowest_frequency = None

        if 'freq' in self.route_section:
            try:
                # List of frequencies
                freq_lines_start = list(utils.find_all(
                    self._full, 'Frequencies --'))
                freq_lines_end = [self._full.find('\n', i)
                                  for i in freq_lines_start]
                freq_lines = [self._full[a:freq_lines_end[i]]
                              for i, a in enumerate(freq_lines_start)]
                frequencies = []
                for freq_line in freq_lines:
                    frequencies.extend((freq_line.split()[2:]))
                frequencies = map(float, frequencies)
                lowest_frequency = frequencies[0]
                # Number of imaginary frequencies
                nimag = len(filter(lambda x: x < 0, frequencies))
            except:
                pass
        return frequencies, lowest_frequency, nimag

    def get_npa_charges(self):
        '''
        Extract NPA charges from Gaussian output file.

        Arguments:
            logfile:        Gaussian output filename or full text
        Return:
            list of floats
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
            self.npa_charges = np.array([float(l.split()[2]) for l
                                        in self._full[start:end].split('\n')])
        except:
            self.npa_charges = None

    def _get_mulliken_charges(self):
        '''
        Extract Mulliken charges from Gaussian output file.

        Arguments:
            logfile:        Gaussian output filename or full text
        Return:
            list of floats
        '''

        try:
            start = self._full.find('Mulliken charges')
            end = self._full.find('Sum of Mulliken charges')
            self.mulliken_charges = [float(l.split()[2]) for l in
                                     self._full[start:end].split('\n')[2:-1]]
        except:
            self.mulliken_charges = None

    def get_nbo_energy(self, atom1, atom2, debug=[]):
        '''
        Parse NBO interaction energies between two atoms from Gaussian output
        file.

        Arguments:
            logfile:        filename or full text
            atom1:          first atom
            atom2:          second atom
            debug           list that will be used to store all parsed lines
        Return:
            float
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
        default: self.files.csv
        '''
        if not csv_file:
            csv_file = self.files.csv
        csv = open(csv_file).read()
        begin = csv.find('Dataset,charges')
        # This will also work if the csv file is cut below the Dataset,charges
        # block (file compression). In that case end=-1, so that csv[begin:-1].
        end = csv.find('Dataset,history_charges')
        self.hi_charges = map(float, csv[begin:end].split('\n')[2:-2])

    def _get_first_summary_block(self):
        end = self._full.find('@')
        begin = self._full.rfind('GINC')
        return self._full[begin:end].replace('\n', '')

    def _get_error(self):
        '''
        Check whether a Gaussian calculation finished normally.

        Arguments:
            logfile:        filename or full text
        Return:
            normal termination (None) or error in ['The combination of
            multiplicity', 'Convergence failure', 'Maximum number of corrector
            steps', 'Number of steps exceeded', 'Error in internal coordinate
            system', 'SVD failed', 'Not terminated']
        '''
        error = None
        possible_errors = {'The combination of multiplicity': 'multi/charge',
                           'Convergence failure': 'convergence',
                           'Maximum number of corrector steps': 'corrector',
                           'Number of steps exceeded': 'steps',
                           'Error in internal coordinate system': 'intcoord',
                           'SVD failed': 'svd',
                           'FormBX had a problem': 'formbx',
                           'NtrErr Called from FileIO': 'ntrerr',
                           'A syntax error was detected in the input line':
                           'syntax'}

        # Only search the bottom of the log file for termination information.
        bottom = self._full[-len(self._full)/2:]

        # If the calculation terminated at all.
        if 'termination' in bottom:

            # Test for all possible errors, but only if no error was found yet,
            # i.e. elif in a for loop.
            for possible_error in possible_errors:
                if not error and possible_error in bottom:
                    error = possible_errors[possible_error]

            # The last 200 characters should include 'Normal termination',
            # since otherwise a combined calculation (e.g. opt freq) could not
            # be evaluated.
            if not error:
                if 'Normal termination' in bottom[-200:]:
                    error = None

                # Other cases are unknown errors.
                else:
                    error = 'unknown'

        # No 'termination in ongoing calculations.
        else:
            error = 'abrupted'

        return error

    def _get_charge_and_multiplicity(self):
        '''
        Parse charge and multiplicity from Gaussian output file.

        Arguments:
            logfile:        filename or full text
        Return:
            {'charge': int, 'multiplicity': int}
        '''
        reg_chargemulti = re.compile(
            r'Charge = *([-0-9]*) *Multiplicity = *(\S*)')
        chargemulti = reg_chargemulti.search(self._full,
                                             self._full.find('Charge'))
        try:
            return chargemulti.group(1), chargemulti.group(2)
        except:
            return float('nan'), float('nan')

    def _get_stoichiometry(self):
        '''
        Parse stoichiometry from Gaussian output file.

        Arguments:
            logfile:        filename or full text
        Return:
            str
        '''
        try:
            begin = self._full.find('Stoichiometry') + 13
            end = self._full.find('\n', begin)
            return self._full[begin:end].strip()
        except:
            return None

    def _get_nproc(self):
        '''
        Get the number of processors used during a Gaussian calculation
        from a Gaussian output file.

        Arguments:
            logfile:        filename or full text
        Return:
            int
        '''
        try:
            nproc_index = self._full.find('%nproc')
            return re.search('%nproc=([0-9]*)', self._full[
                nproc_index:nproc_index + 100]).group(1)
        except:
            return float('nan')

    def _get_comment(self):
        '''
        Get comment line from Gaussian output file.

        Arguments:
            logfile:        filename or full text
        Return:
            str
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

    def _get_chk(self):
        '''
        Get checkpoint filename from Gaussian output file.

        Arguments:
            logfile:        filename or full text
        Return:
            str
        '''
        try:
            begin_checkpoint = self._full.find('chk=') + 4
            end_checkpoint = self._full.find('\n', begin_checkpoint)
            checkpoint = self._full[begin_checkpoint:end_checkpoint].strip()
            return checkpoint.replace('.chk', '')
        except:
            return None

    def _get_method_basis(self):
        '''
        Retrieve method and basis from Gaussian output file or from
        previously parsed keywords (see gaussian.get_route_section).
        At least one of these needs to be provided.

        Arguments:
            logfile:        filename or full text (optional, default: None)
            keywords:       previously parsed keywords (optional,
                            default: None)
        Return:
            {'basis': str, 'method': str}
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
                basis.replace(b, basis_aliases[b])
            lot = '%s/%s' % (method, basis)
        except:
            # If no 'method/basis' was found, assume that the first word
            # of the route section is the method. This is the casef for
            # semiempirical methods and forcefields.
            basis = None
            method = lot = self.keywords[0]
        return method, basis, lot

    def _get_dipole(self):
        '''
        Fetch the total dipole moment (field-independent basis, Debye) from
        Gaussian output file.

        Arguments:
            logfile:        filename or full text

        Return:
            float or None
        '''
        try:
            start = self._full.find('Tot=',
                                    self._full.find('Dipole moment')) + 4
            end = self._full.find('\n', start)
            return float(self._full[start:end])
        except:
            return float('nan')

    def _get_scrf(self):
        '''
        Retrieve SCRF keywords and non-standard input from Gaussian
        output file.

        Arguments:
            logfile:        filename or full text
            route_section:  provide parsed route option (optionally),
                            increases efficiency (default: None)
        Return:
            {'scrf': str, 'scrf_nonstandard_input': str}
        '''
        scrf = re.search(r'(scrf\S*)', self.route_section)
        if not scrf:
            scrf = ''
        else:
            scrf = scrf.group(1)
        try:
            start = self._full.find(
                'Using the following non-standard input for PCM:') + 48
            end = self._full.find('--- end of non-standard input.', start)
            if not start == 47 and not end == -1:
                scrf_nonstandard_input = self._full[start:end].strip()
            else:
                scrf_nonstandard_input = ''
        except:
            scrf_nonstandard_input = ''
        return scrf, scrf_nonstandard_input

    def _get_cputime(self):
        '''
        Parse CPU time from Gaussian output file.

        Arguments:
            log:            full text or filename
        Return:
            number of seconds (float)
        '''
        try:
            begin = self._full.find('Job cpu time') + 13
            end = self._full.find('\n', begin)
            formatted_time = map(float,
                                 self._full[begin:end].strip().split()[::2])
            return sum(formatted_time * np.array([24*3600, 3600, 60, 1]))
        except:
            return None

    @classmethod
    def parse_all(cls, *pargs):
        '''without extension'''
        '''wildcards mogelijk op pythonniveau'''
        '''make generator of it?'''
        '''often multiple operatons on same list ==> better not a generator?'''
        return [cls(parg) for parg in pargs]

    def movie(self):
        '''
        Return an xyz movie from the Gaussian output file. This can be a movie
        of scan point, IRC points or geometry optimization steps. The latter
        will be created if no scan/IRC points are found. The comments contain
        the step number.

        Arguments:
            logfile:            filename or full text of Gaussian output file
        Return:
            str
        '''
        geometries = (self.ircgeometries or self.scangeometries
                      or self.geometries)
        return ''.join([g.to_string(header=True, comment=i) + '\n'
                        for i, g in enumerate(geometries)])

    @classmethod
    def select(cls, name, logfiles):
        '''
        From a list of parsed logfiles, select those where logfile.files.root
        matches name (unix wildcards). In case of a single find, return one
        lgofile.
        '''
        selected = []
        for l in logfiles:
            if fnmatch.fnmatch(l.files.root, name):
                selected.append(l)
        if len(selected) == 1:
            selected = selected[0]
        return selected

    @classmethod
    def partition(cls, logs):
        '''
        Partition a list of ``LOGFile``s according to their filename. Return
        an iterator of sets of ``LOGFile``s. See ``GaussianFiles`` for
        more information.
        '''
        partition = filenames.GaussianFile.partition(
            [l.files for l in logs])
        for p in partition:
            yield set([cls.select(n.root, logs) for n in p])

    @classmethod
    def set_relative_energies(cls, logs):
        '''
        Set relative electronic energies (``relelectronic`` attribute of
        ``LOGFile`` instance), enthalpies (``relenthalpy``) and Gibbs
        free energies (``relgibbs``) for a number of Gaussian ouput
        files.
        '''
        reles = utils.relative_energies([k.electronic for k in logs])
        relhs = utils.relative_energies([k.enthalpy for k in logs])
        relgs = utils.relative_energies([k.gibbs for k in logs])
        for i, l in enumerate(logs):
            l.relelectronic = reles[i]
            l.relenthalpy = relhs[i]
            l.relgibbs = relgs[i]

    @classmethod
    def set_rmsds(cls, logs):
        '''
        Set RMSDs (``rmsd`` attribute of ``LOGFile`` instance) for a number
        of Gaussian output files.
        '''
        for l in logs:
            # Since set_rmsds() is calles in the construcor of GaussianTable
            # a try except is necessary. If molecules don't have corresponding
            # geometries, obtaining an RMSD will cause an error.
            try:
                l.rmsd = list(logs)[0].geometry.rmsd(l.geometry)[2]
            except:
                l.rmsd = float('nan')
