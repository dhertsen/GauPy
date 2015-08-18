'''Common functions'''
import warnings
from ast import literal_eval

# TODO import gaupy.utils in interactive session throws an exception

# Dirty trick to silence warnings of confliciting
# modules on the HPC cluster and to silence NumPy
# warnings.
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
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
    Should also work for unnested lists:
    get_full_attr(someobject, 'energies[3]')
    '''
    try:
        if '.' not in attr:
            if '[' in attr and ']' in attr:
                listname, index = attr.split('[', 1)
                index = index.split(']', 1)[0]
                return getattr(obj, listname)[index]
            else:
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


def find_all(string, sub, index=None, reverse=False):
    '''
    Generator to find all indices at which a substring occurs in a string.

    Arguments:
        string:         parent string
        sub:            substring
    Return:
        generator, yields int
    '''
    if not reverse:
        start = index or 0
        while True:
            start = string.find(sub, start)
            if start == -1:
                return
            yield start
            start += len(sub)
    else:
        start = index or len(string)
        while True:
            start = string.rfind(sub, 0, start)
            if start == -1:
                return
            yield start
            start -= len(sub)


def sfind(string, *pargs):
    '''sequential search function'''
    index = 0
    for p in pargs:
        index = string.find(p, index + 1)
    return index


def liteval(string):
    try:
        return literal_eval(string)
    except:
        return string


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
