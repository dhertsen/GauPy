import gaussian
import numpy as np


def classify(geometry):
    '''
    0   X
    1   reactive C
    2   reactive N
    3   C attached to reactive C, closest to an N atom
    4   other Cs in ring
    5       ''
    6       ''
    7       ''
    '''
    out = []
    at = geometry.numbers
    co = geometry.coordinates

    # Find halogen.
    for i, a in enumerate(at):
        if a > 8:
            out.append(i)
    # Carbon closest to halogen.
    out.append(gaussian.closest(6, out[0], geometry))
    # Nitrogen closest to that carbon.
    out.append(gaussian.closest(7, out[1], geometry))
    # Next carbon, attached to the critical carbon, and closest to a N atom.
    carbons = gaussian.closest(6, out[1], geometry, n=2)
    closest_n = [np.linalg.norm(co[c] - co[gaussian.closest(7, c, geometry)])
                 for c in carbons]
    if min(closest_n) < 1.8:
        dinitro = True
    else:
        dinitro = False
    out.append(gaussian.sort_by(carbons, closest_n))
    # The other ring carbons.
    for i in range(4):
        out.append(gaussian.closest(6, out[-1], geometry, exclude=out))
    # First nitro group
    out.append(gaussian.closest(7, out[5], geometry))
    out.extend(gaussian.closest(8, out[-1], geometry, n=2))
    # Second nitro group
    if dinitro:
        out.append(gaussian.closest(7, out[3], geometry))
        out.extend(gaussian.closest(8, out[-1], geometry, n=2))

    return out


def nitro(geometry):
    dc = 1.3
    n = range(len(geometry.numbers))
    nitros = []
    for i in n:
        if geometry.numbers[i] == 7:
            for j in n:
                if (geometry.numbers[j] == 8
                        and gaussian.dist(i, j, geometry) < dc):
                    for k in n:
                        if (geometry.numbers[k] == 8
                                and gaussian.dist(i, k, geometry) < dc
                                and k > j):
                            nitros.append((i, j, k))
    return nitros


def ammonia(geometry):
    n = range(len(geometry.numbers))
    ammonias = []
    dc = 1.15
    for i in n:
        if geometry.numbers[i] == 7:
            for j in n:
                if (geometry.numbers[j] == 1
                        and gaussian.dist(i, j, geometry) < dc):
                    for k in n:
                        if (geometry.numbers[k] == 1
                                and gaussian.dist(i, k, geometry) < dc
                                and k > j):
                            for l in n:
                                if (geometry.numbers[l] == 1
                                        and gaussian.dist(i, l, geometry) < dc
                                        and l > k):
                                    ammonias.append((i, j, k, l))
    return ammonias


def reactive_nc(geometry):
    '''
    Returns the nitrogen and carbon atom in a C-NH2 or C-NH3 motif. Can be
    used to find the reactive carbon and nitrogen atom in transition states,
    Meiseinheimer complexes and products.
    '''
    n = range(len(geometry.numbers))
    nh = 1.15
    nc = 2.5
    for i in n:
        if geometry.numbers[i] == 7:
            for j in n:
                if (geometry.numbers[j] == 1
                        and gaussian.dist(i, j, geometry) < nh):
                    for k in n:
                        if (geometry.numbers[k] == 1
                                and gaussian.dist(i, k, geometry) < nh
                                and k > j):
                            for l in n:
                                if (geometry.numbers[l] == 6
                                        and gaussian.dist(i,
                                                          l, geometry) < nc):
                                    return {'N': i, 'C': l}


def planarity(geometry):
    cl = classify(geometry)
    p1 = geometry.coordinates[cl[4]]
    p2 = geometry.coordinates[cl[5]]
    p3 = geometry.coordinates[cl[6]]
    p4 = geometry.coordinates[cl[1]]
    return point_to_plane_distance(p1, p2, p3, p4)


def point_to_plane_distance(p1, p2, p3, p4):
    '''
    Return distance between plane defined by p[1-3] and point p4.
    '''
    return abs(np.dot(np.cross(p2-p1, p3-p1) / np.linalg.norm(
        np.cross(p2-p1, p3-p1)), p4-p1))
