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

Several types of geometries can be fetched. It is recommended to use ``geometry`` for the last geometry of the file, since this is generally a lot faster.

.. include:: geometry-table.rst

All these geometries are ``gaupy.molecules.SuperMolecule`` instances. 


Energies
---------

Multiple energies are available:

.. include:: energy-table.rst
