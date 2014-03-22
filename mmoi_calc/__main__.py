# -*- coding: utf-8 -*-
"""
Command line utility to calculate moments of inertia of a molecule described
by a "connection table" (found in .mol and .sdf files).
"""

from __future__ import print_function
import argparse

import numpy

from . import ctab


def main():
    """
    Parses the provided file and calculates inertia.
    """
    parser = argparse.ArgumentParser(
        description="Calculate msoi from a ctab.")
    parser.add_argument('file', help="the molecule file to process")
    parser.add_argument('--principal_axes', action='store_true',
                        help="also show the molecule's principal axes")
    args = parser.parse_args()

    with open(args.file, 'r') as molfile:
        molecule = ctab.Parser().molfile(molfile)

    numpy.set_printoptions(precision=2, suppress=True)
    print("Molecule: {}".format(molecule))
    if args.principal_axes:
        print("Principal moments of inertia: {}\nPrincipal axes:\n{}".format(
            *molecule.inertia(moments_only=False)))
    else:
        print("Principal moments of inertia: {}".format(molecule.inertia()))


if __name__ == "__main__":
    main()
