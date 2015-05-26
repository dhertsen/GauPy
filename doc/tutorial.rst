Using the Python modules
==========================

Load output files
-----------------

:class:`gaupy.log.LOGFile` is the main class for analyzing output files. Fetching attributes from the output files is straightforward::

    >>> from gaupy.log import LOGFile
    >>> lf = LOGFile('output.log')
    >>> lf.multiplicity
    1

For the sake of efficiency, attributes are parsed and cached when they are first requested. Depending on the attribute, a different default value will be returned if it could not be retrieved from the output file. If a single numeric value is expected, ``nan`` is the default value. In all other cases, a default with a ``False`` boolean value is returned (``None``, ``False``, ``[]``). The latter allows for pythonic constructs::

    lf = LOGFile('output.log')
    if lf.energies and not lf.geometries:
        print('What is this? Energies found, but no geometries')

This is probably a good time to talk about filenames and :class:`gaupy.filenames.GaussianFile`. It does not matter which file type is passed to the ``LOGFile`` constructor: ``LOGFile('output.log')``, ``LOGFile('output.chk')``, ``LOGFile('output')``, etc. are equivalent and will look for an output file with a log extension. For each calculation, a ``GaussianFile`` instance is created and stored as the ``files`` attribute of the ``LOGFile``. This object allows for a number of common filename operations and quick access to different extensions. As the root of the filename ``lf.files.root`` is used frequently, ``lf.file`` is provided as a quick alias.

::

    >>> lf = LOGFile('folder/output.extension')
    >>> lf.files.root
    folder/output
    >>> lf.files.base
    output
    >>> lf.otherextension
    folder/output.otherextension



Geometries
-----------

Several types of geometries can be fetched:

.. include:: geometry-table.rst

Since requesting all geometries can be slow, it is recommended to load only the last one via ``geometry`` unless more are needed. All these geometries are :class:`gaupy.molecules.SuperMolecule` instances, a subclass of `molmod.molecules.Molecule <http://molmod.github.io/molmod/reference/basic.html#molmod-molecules-molecular-systems>`_. Several inherited attributes:

================================================================================                ============================================================
Attribute                                                                                       Description
================================================================================                ============================================================
:attr:`coordinates <gaupy.molecules.SuperMolecule.coordinates>`                                 cartesian coordinates (atomic units)
:attr:`mass <gaupy.molecules.SuperMolecule.mass>`                                                molecular mass
:attr:`numbers <gaupy.molecules.SuperMolecule.numbers>`                                          atomic numbers
:attr:`size <gaupy.molecules.SuperMolecule.size>`                                                number of atoms
================================================================================                ============================================================

A lot of functionality has been added in :class:`SuperMolecule <gaupy.molecules.SuperMolecule>`, both by extending the core functionality of `Molecule <http://molmod.github.io/molmod/reference/basic.html#molmod-molecules-molecular-systems>`_ and by creating some aliases. For most atomic properties, the rows of two-dimensional and the elements of one-dimensional `NumPy <http://www.numpy.org>`_ arrays correspond to atoms in the XYZ matrix in the log files in the same order. Several geometric parameters can be calculated from a :class:`SuperMolecule <gaupy.molecules.SuperMolecule>`:

================================================================================                ============================================================
Function                                                                                        Description
================================================================================                ============================================================
:func:`dist(*atoms) <gaupy.molecules.SuperMolecule.dist>`                                       distance between in angstrom
:func:`angle(*atoms) <gaupy.molecules.SuperMolecule.angle>`                                     angle between three atoms in degrees
:func:`dihedral(*atoms) <gaupy.molecules.SuperMolecule.dihedral>`                               dihedral angle between four atoms in degrees
:func:`opbend_angle(*atoms) <gaupy.molecules.SuperMolecule.opbend_angle>`                       angle between the plane defined by the first three atoms
                                                                                                and the vector defined by the first and the fourth atom
                                                                                                in degrees
:func:`opbend_dist(*atoms) <gaupy.molecules.SuperMolecule.opbend_dist>`                         distance between the plane defined by the first three atoms
                                                                                                and the vector defined by the first and the fourth atom
                                                                                                in angstroms
================================================================================                ============================================================

More complex molecular operations can be performed as well:

=================================================================================================================   ========================================================================================================================================================================================================================
Function                                                                                                                Description
=================================================================================================================   ========================================================================================================================================================================================================================
:func:`closest(atom_type, reference, n=1, exclude=False, only=False) <gaupy.molecules.SuperMolecule.closest>`           Return a list of the ``n`` atoms of atom type ``atom_type`` that are closest to atom with number ``reference``. This search can be limited to a list of atoms (``only``) and/or exclude a list of atoms (``exclude``).
:func:`nonhydrogens() <gaupy.molecules.SuperMolecule.nonhydrogens>`                                                      return a list of all heavy atoms
:func:`part(atoms) <gaupy.molecules.SuperMolecule.part>`                                                                return a new :class:`SuperMolecule <gaupy.molecules.SuperMolecule>` that consist of the atoms in ``atoms``
:func:`nrings(n) <gaupy.molecules.SuperMolecule.nrings>`                                                            returns a list of all ``n`` membered rings in the molecule
:func:`molecules() <gaupy.molecules.SuperMolecule.molecules>`                                                           returns all molecules in the system as a list of :class:`SuperMolecule <gaupy.molecules.SuperMolecule>`
:func:`break_bond(atom1, atom2) <gaupy.molecules.SuperMolecule.break_bond>`                                              breaks the bond between ``atom1`` and ``atom2`` and return both fragments as new <gaupy.molecules.SuperMolecule>`
:func:`bonded(atom1, atom2) <gaupy.molecules.SuperMolecule.bonded>`                                                     checks whether there is a bond between ``atom1`` and ``atom2``, return a boolean
:func:`non_hydrogen_neighbors(atom) <gaupy.molecules.SuperMolecule.non_hydrogen_neighbors>`                                        returns a list of heavy atom neighbors of atom ``atom``
=================================================================================================================   ========================================================================================================================================================================================================================

Energies
--------

Multiple energies are available:

.. include:: energy-table.rst

::

    for s in ['system1', 'system2', 'system3']:
        l = LOGFile(s)
        lspe = LOGFile(s + '-spe')
        gibbs = l['energy'] + lspe['gibbscorrection']
