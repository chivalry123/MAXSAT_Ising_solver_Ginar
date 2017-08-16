#!/usr/bin/env python

"""
Convert Cluster Expansion parameters in CASM's format to the format used
by Wenxuan's ground-state finder code.

Restrictions:
- The PRIM file has to contain fractional (direct) coordinates.
- The components of all atomic coordinates have to be in [0,1).

"""

from __future__ import print_function, division

__author__ = "Alexander Urban, Wenxuan Huang, Ziqin Rong"
__email__ = "alexurba@mit.edu, key01027@mit.edu, rongzq08@mit.edu"
__date__ = "2015-04-02"
__version__ = "0.1.1"

import argparse
import numpy as np
import sys


def read_PRIM(prim_file):
    with open(prim_file, 'r') as fp:
        fp.readline()
        a = float(fp.readline())
        avec = np.zeros((3, 3))
        avec[0] = map(float, fp.readline().split())
        avec[1] = map(float, fp.readline().split())
        avec[2] = map(float, fp.readline().split())
        avec *= a
        nsites = map(int, fp.readline().split())
        cooformat = fp.readline()
        if cooformat.strip() not in ["Direct", "direct"]:
            print("Error: This script can only handle direct coordinates!")
            print(cooformat)
            sys.exit()
        coords = []
        species = []
        for i in range(np.sum(nsites)):
            line = fp.readline()
            coords.append(map(float, line.split()[:3]))
            # FIXME: make sure that all coordinates are in [0:1]
            species.append(line.split()[3:])
        coords = np.array(coords)
    return (avec, coords, species)


def read_ECIs(eci_out_file):
    with open(eci_out_file, 'r') as fp:
        for i in range(7):
            line = fp.readline()
        eci = []
        cluster = []
        line = fp.readline()
        while line:
            if abs(float(line.split()[1]))>1e-20:
                eci.append(float(line.split()[1]))
                cluster.append(int(line.split()[2]))
            line = fp.readline()
    return (eci, cluster)


def get_CE_sites(coords, species):
    sites = []
    species_index = []
    for i in range(len(coords)):
        if species[i] > 1:
            sites.append(i)
            if species[i] not in species_index:
                species_index.append(species[i])
    return (sites, species_index)


def parse_orbit(fp, J, nclus, coords, sites, species_index, eci):
    # print ("\n in parse_orbit ")

    fp.readline()
    (size, num) = map(int, fp.readline().split()[1:3])
    for i in range(num):
        # print ("\n i is {}".format(i))
        fp.readline()
        T = []
        pos = []
        spc = []
        order=[]
        # read and convert a single cluster
        for j in range(size):
            # print ("\n j is {}".format(j))
            line = fp.readline()
            coo = np.array(map(float, line.split()[:3]))
            coo = np.where(np.abs(coo)<1.0e-10, 0, coo)
            site_species = line.split()[3:-1]
            spc.append(species_index.index(site_species))
            # NOTE: only works if original coordinates are in [0:1]
            T.append(map(int, np.floor(coo)))
            order.append(line.split()[-1])
            for s in sites:
                if (np.linalg.norm(coords[s] + T[j] - coo) < 0.0001):
                    pos.append(s)
        # print ("\n non array T is ")
        # print (T)


        T = np.array(T)
        T_min = np.min(T, axis=0)

        # print ("\n array T is ")
        # print (T)

        T = T + (np.array([1, 1, 1]) - T_min)

        # print ("\n T_min is")
        # print (T_min)
        # print ("\n np.array([1, 1, 1]) is")
        # print (np.array([1, 1, 1]))


        nclus += 1
        print("Cluster {}".format(nclus))
        for c in range(size):
            print("{},{},{},{},{},{}   ".format(
                # T[c, 0], T[c, 1], T[c, 2], pos[c]+1, spc[c]+1), end="")
                T[c, 0], T[c, 1], T[c, 2], pos[c]+1, 1 ,order[c]), end="")
        print("\nJ={}\n".format(J))
    return nclus


def skip_orbit(fp):
    fp.readline()
    (size, num) = map(int, fp.readline().split()[1:3])
    for i in range(num):
        fp.readline()
        for j in range(size):
            fp.readline()


def parse_FCLUST(fbclust_file, coords, sites, species_index, eci, cluster):
    if cluster[0] != 0:
        print("Constant {}\n".format(0))
    else:
        print("Constant {}\n".format(eci[0]))
    nclus = 0
    with open(fbclust_file, 'r') as fp:
        norbit = int(fp.readline())
        for iorb in range(1, norbit+1):
            if iorb in cluster:
                J = eci[cluster.index(iorb)]
                # print ("\n outside parse_orbit with iorb={} ".format(iorb))
                nclus = parse_orbit(fp, J, nclus, coords, sites,
                                    species_index, eci)
            else:
                skip_orbit(fp)


def convert(prim_file, fbclust_file, eci_out_file):

    (avec, coords, species) = read_PRIM(prim_file)
    # print ("\navec is")
    # print (avec)
    # print ("\ncoords is")
    # print (coords)
    # print ("\nspecies is")
    # print (species)

    (eci, cluster) = read_ECIs(eci_out_file)

    # print ("\neci is")
    # print (eci)
    # print ("\ncluster is")
    # print (cluster)

    # print ("\nbefore get_CE_sites(coords, species)\n")
    (sites, species_index) = get_CE_sites(coords, species)

    # print ("\nsites is")
    # print (sites)
    # print ("\n species_index is ")
    # print (species_index)


    parse_FCLUST(fbclust_file, coords, sites, species_index, eci, cluster)


if (__name__ == "__main__"):

    parser = argparse.ArgumentParser(
        description=__doc__+"\n{} {}".format(__date__, __author__),
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument(
        "--PRIM", "-p",
        help="Path to a CASM PRIM file (default: PRIM).",
        default="PRIM")

    parser.add_argument(
        "--FBCLUST", "-f",
        help="Path to a CASM FCLUST file (default: FBCLUST).",
        default="FBCLUST")

    parser.add_argument(
        "--eci-out", "-e",
        help="Path to a CASM eci.out file (default: eci.out).",
        default="eci.out")

    args = parser.parse_args()

    convert(args.PRIM, args.FBCLUST, args.eci_out)

    #note output to J_in_tern_casm.in

