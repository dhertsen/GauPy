import os
import re
import glob
import shutil
import logging


class GaussianFile(object):
    '''
    split [a-zA-Z0-9]
    possible to override reaction_patterns
    all extensions possible
    filenames in different folders can be compared
    split is without path, root with
    '''

    complex_patterns = ['fwd', 'rev', 'react', 'int', 'prc', 'prod']
    irc_patterns = ['irc' + x for x in complex_patterns]
    ts_patterns = ['ts', 'ts1', 'ts2', 'tsopt']
    reaction_patterns = (complex_patterns + irc_patterns + ts_patterns
                         + ['spe'])

    def __init__(self, name):
        '''
        should work both with a GaussianFile instance and and strings
        '''
        try:
            # constructor with GaussianFile argument
            new = type(self)(name.root)
            self.root = new.root
            self.base = new.base
            self.directory = new.directory
            self.split = new.split
        except:
            # constructor from filename (str)
            self.root = os.path.splitext(name)[0]
            self.base = os.path.basename(self.root)
            self.directory = os.path.dirname(self.root)
            self.split = re.findall(r"[a-zA-Z0-9]+", self.base)

    def __getattr__(self, extension):
        if extension not in ['root', 'base', 'directory', 'split']:
            return '%s.%s' % (self.root, extension)

    def get_all(self):
        return glob.iglob('%s.*' % self.root)

    def _reset_root(self, dir):
        self.root = os.path.join(dir, os.path.basename(self.root))

    def replace(self, find, replace):
        newsplit = []
        for i, element in enumerate(self.split):
            newsplit.append(element if not element == find else replace)
        return GaussianFile(self.directory + '-'.join(newsplit))

    def add(self, element):
        return GaussianFile('%s-%s' % (self.root, element))

    def remove(self, *elements):
        newsplit = self.split
        for el in elements:
            newsplit = [s for s in self.newsplit if s != el]
        return GaussianFile(self.directory + '-'.join(newsplit))

    def ts_to_irc(self, replace=False):
        '''
        yields two values, cause a single TS yields two IRCs
        dictionary with keys 'fwd', 'rev'
        replace: try ts_patterns -> irc_patterns first, then add
        else: just add
        '''
        ircs = dict()
        for i in ['fwd', 'rev']:
            irc = 'irc' + i
            ircs[i] = GaussianFile('%s-%s' % (self.root, irc))
            if replace:
                for t in self.ts_patterns:
                    if t in self.split:
                        ircs[i] = self.replace(t, irc)
        return ircs

    def irc_to_end(self):
        '''
        Returns a single value, since one IRC corresponds to a single complex.
        '''
        for irc, cplx in zip(self.irc_patterns, self.complex_patterns):
            if irc in self.split:
                return self.replace(irc, cplx)
        else:
            return self.add('opt')

    def move_all(self, dir):
        for f in self.get_all():
            to = os.path.join(dir, os.path.basename(f))
            shutil.move(f, to)
            print('Moved %s to %s' % (f, to))
        self._reset_root(dir)

    def move_ext(self, dir, *extensions):
        ''' does not trigger reset of root'''
        for ext in extensions:
            ext = getattr(self, ext)
            to = os.path.join(dir, os.path.basename(ext))
            shutil.move(ext, to)
            print('Moved %s to %s' % (ext, to))

    def copy_all(self, dir):
        for f in self.get_all():
            to = os.path.join(dir, os.path.basename(f))
            shutil.copy(f, to)
            print('Copied %s to %s' % (f, to))
        self._reset_root(dir)

    def copy_ext(self, dir, *extensions):
        ''' does not trigger reset of root'''
        for ext in extensions:
            ext = getattr(self, ext)
            to = os.path.join(dir, os.path.basename(ext))
            shutil.copy(ext, to)
            print('Copied %s to %s' % (ext, to))

    def remove_all(self):
        for f in self.get_all():
            os.remove(f)
            print('Removed %s' % f)

    def remove_ext(self, *extensions):
        for ext in extensions:
            ext = getattr(self, ext)
            os.remove(ext)
            print('Removed %s' % ext)

    def same_reaction(self, other):
        '''other must be GaussianFile'''
        unique_self = set(self.split) - set(self.reaction_patterns)
        unique_other = set(other.split) - set(self.reaction_patterns)
        return unique_self == unique_other

    @classmethod
    def partition(cls, files):
        # haalt ook dubbels eruit
        # is een generator
        l = set(files)
        while l:
            a = l.pop()
            b = filter(lambda x: a.same_reaction(x), l)
            c = set([a] + b)
            l = l - c
            yield c

    @classmethod
    def max_width(cls, files):
        ''' works on string filenames as well '''
        try:
            return max([len(f.root) for f in files])
        except:
            return max([len(cls(f).root) for f in files])
