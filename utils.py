'''Common functions'''
import warnings
import log

# Dirty trick to silence warnings of confliciting
# modules on the HPC cluster and to silence NumPy
# warnings.
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    import numpy as np
    import molmod.units as units
    import molmod.periodic as periodic


def anum(atom):
    '''
    convert atomic symbol or number to atomic number
    '''
    return periodic.periodic[atom].number


def asym(atom):
    '''
    covert atomic symbol or number to atomic symbol
    '''
    return periodic.periodic[atom].symbol


def get_full_attr(obj, attr):
    '''
    Make constructs like _getattr(somobject, 'a.b.c.d') possible. Whenever
    there is an error fetching the terminal attribute, None is returned.
    '''
    try:
        if '.' not in attr:
            return getattr(obj, attr)
        else:
            superattr, subattr = attr.split('.', 1)
            return get_full_attr(getattr(obj, superattr), subattr)
    except:
        return None


def translate(l, shift):
    '''
    list to list
    '''
    return [l[(i - shift) % len(l)] for i, x in enumerate(l)]


def all_translations(l):
    '''
    list to generator
    '''
    for x in range(0, len(l)):
        yield translate(l, x)
        yield translate(l[::-1], x)


def is_translation(l1, l2):
    for t in all_translations(l1):
        if l2 == t:
            return True
    else:
        return False


def find_all(string, sub):
    '''
    Generator to find all indices at which a substring occurs in a string.

    Arguments:
        string:         parent string
        sub:            substring
    Return:
        generator, yields int
    '''
    start = 0
    while True:
        start = string.find(sub, start)
        if start == -1:
            return
        yield start
        start += len(sub)


def energy(value, energy_type='gibbs'):
    '''
    Convert a energy or a filename to an energy. If a filename was provided,
    energy of type ``energy_type`` will be parsed. ``float('nan')`` values are
    allowed and will even be returned if an error occurs.
    '''
    try:
        return float(energy)
    except:
        try:
            return getattr(log.LOGFile(value), energy_type)
        except:
            return float('nan')


def energies(values, energy_type='gibbs'):
    '''
    Convert a list or NumPy array of energies and filenames to a NumPy array
    of energies. Energy of type ``energy_type`` will be taken from log files.
    ``float('nan')`` values are allowed.
    '''
    if values:
        return np.array([energy(e, energy_type) for e in values])
    else:
        return np.array([])


def relative_energies(values, reference='min', absolute=0, relative=0,
                      conversion='kjmol', energy_type='gibbs'):
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
            ``energy_type=``            type
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
    if energies:
        # Convert to numerical energies
        e = energies(values, energy_type=energy_type)
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


class cached(object):
    """
    Descriptor (non-data) for building an attribute on-demand on first use.
    """
    def __init__(self, factory):
        """
        <factory> is called such: factory(instance) to build the attribute.
        """
        self._attr_name = factory.__name__
        self._factory = factory

    def __get__(self, instance, owner):
        # Build the attribute.
        attr = self._factory(instance)

        # Cache the value; hide ourselves.
        setattr(instance, self._attr_name, attr)

        return attr
