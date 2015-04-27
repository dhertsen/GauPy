import molmod.graphs as g
import molmod.molecular_graphs as mg
import utils
import warnings

# Dirty trick to silence warnings of confliciting
# modules on the HPC cluster and to silence NumPy
# warnings.
warnings.filterwarnings('ignore')

h = mg.HasAtomNumber(1)
c = mg.HasAtomNumber(6)
n = mg.HasAtomNumber(7)
o = mg.HasAtomNumber(8)

csp3 = g.CritAnd(c, mg.HasNumNeighbors(4))
csp2 = g.CritAnd(c, mg.HasNumNeighbors(3))
csp = g.CritAnd(c, mg.HasNumNeighbors(2))


def _convert(atom):
    try:
        return mg.HasAtomNumber(utils.anum(atom))
    except:
        return atom


def me(atom):
    '''
    at can be 'C', 'c', 6 or a pattern
    '''
    at = _convert(atom)
    return g.CritAnd(mg.HasNeighbors(h, h, h, at), c)


def ipr(atom):
    at = _convert(atom)
    return g.CritAnd(c, mg.HasNeighbors(me(c), me(c), h, at))


def ch2(x, y):
    x = _convert(x)
    y = _convert(y)
    return g.CritAnd(c, mg.HasNeighbors(h, h, x, y))
