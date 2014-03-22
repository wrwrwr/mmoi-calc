# -*- coding: utf-8 -*-
import argparse

import numpy

import ctab


parser = argparse.ArgumentParser(description="Calculate msoi from a ctab.")
parser.add_argument('file', help="the molecule file to process")
parser.add_argument('--principal_axes', action='store_true',
                    help="also show the molecule's principal axes")
args = parser.parse_args()

with open(args.file, 'r') as file:
    molecule = ctab.Parser().molfile(file)
   
numpy.set_printoptions(precision=2, suppress=True)
print "Molecule: {}".format(molecule)
if args.principal_axes:
    print "Principal moments of inertia: {}\nPrincipal axes:\n{}".format(
        *molecule.inertia(moments_only=False))
else:
    print "Principal moments of inertia: {}".format(molecule.inertia())
