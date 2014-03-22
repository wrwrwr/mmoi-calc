# -*- coding: utf-8 -*-
import argparse

import numpy

import ctab


parser = argparse.ArgumentParser(description='Calculate msoi from a ctab.')
parser.add_argument('file', help='the molecule file to process')
args = parser.parse_args()

with open(args.file, 'r') as file:
    molecule = ctab.Parser().molfile(file)
   
numpy.set_printoptions(precision=2, suppress=True)
print """
Molecule: {}
Principal moments of inertia: {}
""".format(molecule, molecule.inertia())
