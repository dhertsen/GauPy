import log
import warnings
import utils
import numpy as np
import hpc

# Dirty trick to silence warnings of confliciting
# modules on the HPC cluster and to silence NumPy
# warnings.
warnings.filterwarnings('ignore')


def _ignore(value):
    '''
    Check whether a data point can be ignored.
    '''
    if value in [None, 'na', 'nan', '']:
        return True
    else:
        try:
            if np.isnan(value):
                return True
            else:
                return False
        except:
            return False


class ColumnLayout(object):
    def __init__(self, layout):
        self.layout = layout
        self.width = len(layout % 0)
        self.string = (layout[:layout.find('.')] + 's'
                       if '.' in layout else layout[:-1] + 's')

    def format(self, data):
        try:
            if _ignore(data):
                return ' ' * self.width
            else:
                return self.layout % data
        except:
            return self.string % data


class GaussianTable(object):
    alias = {'gibbs': 'G', 'electronic': 'E', 'irc.length': '#IRC',
             'gibbs_relative': 'Grel', 'electronic_relative': 'Erel',
             'nimag': '#IM', 'lowest_frequency': 'FREQ', 'relgibbs': 'Grel',
             'frequencies[1]': 'FREQ2', 'frequencies[2]': 'FREQ3',
             'frequencies[3]': 'FREQ4', 'frequencies[4]': 'FREQ5',
             'relenthalpy': 'Hrel', 'temperature': 'T', 'pressure': 'P',
             'multiplicity': 'MULTI', 'termination_count': 'NTERM',
             'predicted_last': 'PRED', 'scan.length': '#SCAN'}

    def __init__(self, logs, columns, layout, partition=False,
                 partition_numbers=False,
                 partition_energies_rmsds=False,
                 sort_gibbs=False, type=log.LOGFile):
        self.type = type
        self.logs = [l for l in self.type.parse_all(*logs) if l.exists]
        if sort_gibbs:
            self.logs.sort(key=lambda x: x.gibbs)
        self.type.set_relative_energies(self.logs)
        self.type.set_rmsds(self.logs)
        self.columns = columns
        if 'hpc' in self.columns:
            hpc.set_status(*self.logs)
        self._layout(layout)
        self._empty_columns()
        self.width = sum([l.width for l in self.layout])
        self.caption = ''.join([l.string % self.alias.get(c, c.upper())
                                for c, l in zip(self.columns, self.layout)])
        self.line = self.width * '-'
        self._rows(partition=partition, partition_numbers=partition_numbers,
                   partition_energies_rmsds=partition_energies_rmsds)

    def _layout(self, layout):
        '''
        set self.layout from a list of strings,
        and fit columns
        fit+i is allowed
        '''
        self.layout = []
        for i, c in enumerate(layout):
            if 'fit' in c:
                try:
                    extra = int(c.replace('fit', ''))
                except:
                    extra = 0
                # maximum width of content
                contentwidth = max([len('%s' % utils.get_full_attr(
                    l, self.columns[i])) for l in self.logs])
                # make sure column header is not wider,
                # also look at aliases
                width = max(contentwidth, len(self.alias.get(
                    self.columns[i], self.columns[i])))
                fitted = '%%-%is' % (width + extra)
                self.layout.append(ColumnLayout(fitted))
            else:
                self.layout.append(ColumnLayout(c))

    def _empty_columns(self):
        '''
        delete empty and invalid columns
        '''
        for i, c in enumerate(self.columns):
            empty = [all([_ignore(utils.get_full_attr(lg, col))
                          for lg in self.logs])
                     for col in self.columns]
            self.layout = [lay for i, lay in enumerate(self.layout)
                           if not empty[i]]
            self.columns = [col for i, col in enumerate(self.columns)
                            if not empty[i]]

    def _row(self, logfile):
        return(''.join([l.format(utils.get_full_attr(logfile, c))
                        for c, l in zip(self.columns, self.layout)]))

    def _rows(self, partition, partition_numbers, partition_energies_rmsds):
        self.rows = [self.line, self.caption, self.line]
        if partition:
            for p in self.type.partition(self.logs,
                                         group_by_number=partition_numbers):
                if partition_energies_rmsds:
                    self.type.set_relative_energies(p)
                    self.type.set_rmsds(p)
                for l in p:
                    self.rows.append(self._row(l))
                self.rows.append(self.line)
        else:
            for l in self.logs:
                row = self._row(l)
                if row:
                    self.rows.append(row)
            self.rows.append(self.line)

    def __str__(self):
        return '\n'.join(self.rows)
