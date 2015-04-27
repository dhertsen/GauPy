import gaussian


def ordered_rings(pox):
    '''
    Return the ring atoms of a pox in the (O, C2, N, C3, C4) way. This will
    probably not work for transition states.

    Arguments:
        pox:    filename with arbitrary extension

    Return:
        List of dictionaries {'N':i, 'C2':j, 'N':k, 'C4':l, 'C5':m}
    '''
    ordered_rings = []
    f = gaussian.filenames(pox)
    geom = gaussian.analyze_log(f['log'])['geometry']
    for xyznumbers in gaussian.get_nrings(geom, 5):
        atomnumbers = [geom.numbers[n] for n in xyznumbers]
        try:
            o = atomnumbers.index(8)
            n = atomnumbers.index(7)
            c_closest_o = set(gaussian.closest(6, xyznumbers[o], geom, n=2,
                                               only=xyznumbers))
            c_closest_n = set(gaussian.closest(6, xyznumbers[n], geom, n=2,
                                               only=xyznumbers))
            o = xyznumbers[o]
            n = xyznumbers[n]
            c2 = list(c_closest_o & c_closest_n)[0]
            c4 = list(c_closest_n - c_closest_o)[0]
            c5 = list(c_closest_o - c_closest_n)[0]
            ordered_rings.append({'O': o, 'C2': c2, 'N': n, 'C4': c4,
                                  'C5':  c5})
        except:
            pass
    return ordered_rings


def ring_charges(pox):
    '''
    Return the Hirshfeld I charges in all oxazoline rings. This will probably
    not work for transition states.

    Arguments:
        pox:    filename with arbitrary extension

    Return:
        List of dictionaries [{'N': charge_of_N_in_this ring,..}, ...]
    '''
    f = gaussian.filenames(pox)
    rings = ordered_rings(pox)
    for ring in rings:
        for atom in ring:
            ring[atom] = gaussian.hi_charges(f['csv'])[ring[atom]]
    return rings


def cation_monomer_ring_charges(pox):
    '''
    Return the Hirshfeld I charges in all oxazoline rings. The rings will
    be classified as monomer or cation rings. This will probably not work for
    transition states.

    Arguments:
        pox:    filename with arbitrary extension

    Return:
        List of dictionaries ['cation': {'N': charge_of_N_in_cation_ring,..}
        , 'monomer'; {}]
    '''
    identified_rings = dict()
    charges = ring_charges(pox)
    for ring in charges:
        if ring['N'] > -0.4:
            identified_rings['cation'] = ring
        else:
            identified_rings['monomer'] = ring
    return identified_rings
