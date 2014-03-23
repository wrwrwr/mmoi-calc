# -*- coding: utf-8 -*-
"""
Rudimentary parser and data structures that can represent a "connection table"
(contained in some .mol, .sdf and other files).

TODO: Currently only processes the atoms block.
"""

import logging
import re

import numpy
from .elements import ELEMENTS


# Connection table (for the legacy V2000 format).
# Reference: http://download.accelrys.com/freeware/ctfile-formats/.
# Note: The v2 format seems to be designed with fixed position splitting
#       in mind, admittedly using regular expressions is an overshot.
V2_TYPES = {
    'int2': r'(?=\ *-?\d*\ *)[-\d ]{2}',
    'int3': r'(?=\ *-?\d*\ *)[-\d ]{3}',
    'real': r'(?=\ *-?\d*\.?\d*\ *)[-.\d ]{10}',
    'bool': r'(?=\ *[01 ]\ *)[01 ]{3}'}

# Connection table header, specifies counts of lines in following blocks.
# aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
V2_COUNTS = re.compile(r'''
    (?P<atoms>{int3})             # Number of atoms.
    (?P<bonds>{int3})             # Number of bond lines.
    (?P<lists>{int3})             # Number of atom lists.
    .{{3}}                        # Obsolete.
    (?P<chiral>{bool})            # Chirality flag (1 = chiral).
    (?P<stexts>{int3})            # Structural text lines.
    .{{12}}                       # Obsolete x4.
    (?P<props>{int3})             # Number of additional property lines.
    \ ?[vV]2000\s*                # Version identifier.
'''.format(**V2_TYPES), re.VERBOSE)

# A single atom line of the legacy connection table.
# xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
V2_ATOM = re.compile(r'''
    (?P<x>{real})                 # Coordinates (resulting bond lengths look
    (?P<y>{real})                 # like in Angstroms, but that doesn't seem
    (?P<z>{real})                 # to be documented anywhere). Also note the
    \                             # space separating further entries.
    (?P<symbol>[\w\d* ]{{3}})     # Atom symbol, group label or
                                  # query -- L, A, Q, LP or *.
    (?P<mass_diff>{int2})         # Difference from mass in the periodic
                                  # table; should be within the -3..4 range.
    (?P<charge>{int3})            # Charge (1 = +3, 2 = +2, 3 = +1, 4 = double
                                  # radical, 5 = -1, 6 = -2, 7 = -3).
    (?P<stereo>{int3})            # Stereo parity (1 = odd, 2 = even,
                                  # 3 = either or unmarked center).
    (?P<hydrogen>{int3})          # Hydrogen count (in excess of explicit;
                                  # queries only).
    (?P<stereo_care>{bool})       # Whether to consider stereo for double
                                  # bonds (1 = match; queries only).
    (?P<valence>{int3})           # Number of bonds including implied
                                  # hydrogens (0 = default, 15 = 0).
    (?P<no_hydrogen>{bool})       # Legacy zero additional hydrogen designator
                                  # (1 = no hydrogen allowed).
    .{{6}}                        # Unused x2.
    (?P<atom_atom>{int3})         # Atom-atom mapping number (reactions only).
    (?P<conf>{int3})              # Inversion / retention flag (1 = inverted,
                                  # 2 = retained; reactions only).
    (?P<exact>{bool})\s*          # Exact change flag (1 = changes must be
                                  # exact).
'''.format(**V2_TYPES), re.VERBOSE)


class Atom:
    """
    A line in the "atoms" block."
    """
    def __init__(self, match):
        """
        Atom object from a regex match.

        Drops in element data under ``element`` property and
        sets ``mass`` taking into account mass difference entry.

        TODO: Additional properties may also include mass modifications.
        """
        self.symbol = match.group('symbol').strip()
        self.element = ELEMENTS[self.symbol]
        self.coords = numpy.array(match.group('x', 'y', 'z'), dtype=float)
        mass_difference = match.group('mass_diff').strip()
        if mass_difference == '':
            mass_difference = 0
        self.mass = self.element.mass + int(mass_difference)

    def __str__(self):
        return self.symbol

    def __repr__(self):
        return "{} {}".format(self.symbol, self.coords)


class Molecule:
    """
    Collects all information from a "connection table" (which does not need to
    actually describe a molecule, but some more general atoms collection).
    """
    def __init__(self, atoms):
        """
        Molecule object from atoms and bonds data.
        """
        self.atoms = atoms
        self.mass = sum(a.mass for a in atoms)

    def center_of_mass(self):
        """
        Calculates the molecule's center of mass (in data coordinates).
        """
        return sum(a.mass * a.coords for a in self.atoms) / self.mass

    def inertia(self, principal=True, moments_only=True):
        """
        Calculates the moments of inertia of the molecule (Dalton Angstrom^2).

        With ``principal`` set to false returns inertia matrix in the data
        coordinates (with non-zero products of inertia), with ``moments_only``
        set to false, returns principal moments and axes (as matrix columns).

        See: https://en.wikipedia.org/wiki/Moment_of_inertia
            #The_inertia_matrix_for_spatial_movement_of_a_rigid_body
        """
        center_of_mass = self.center_of_mass()
        inertia = numpy.zeros((3, 3))
        for atom in self.atoms:
            cmx, cmy, cmz = atom.coords - center_of_mass
            inertia += - atom.mass * numpy.matrix(
                [[ 0,   -cmz,  cmy],
                 [ cmz,  0,   -cmx],
                 [-cmy,  cmx,  0  ]])**2
        if principal:
            # TODO: SciPy has an ``eigvals_only`` argument to eigh.
            if moments_only:
                return numpy.linalg.eigvalsh(inertia)
            else:
                return numpy.linalg.eigh(inertia)
        else:
            return inertia

    def __str__(self):
        return ' '.join(str(a) for a in self.atoms)

    def __repr__(self):
        return repr(self.atoms)


class ParseError(Exception):
    """
    Signals a critical parsing failure.
    """
    pass


class Parser:
    """
    Can produce a ``Molecule`` from a file containing a "connection table".
    """
    def __init__(self):
        # TODO: Some parser configuration?
        pass

    def molfile(self, lines):
        """
        Takes a file containing a V2000 ctab and returns a ``Molecule``.

        The ``lines`` argument may be a list, an open ``file`` object or
        anything else that can be iterated in a line-wise fashion (like a
        ``StringIO``).
        """
        # Find the "counts" line.
        for line in lines:
            counts = V2_COUNTS.match(line)
            if counts:
                break
        else:
            raise ParseError("Couldn't find the ctab header; "
                             "is this a V2000 file?")

        # Parse atoms data (assumed to follow).
        atoms_count = int(counts.group('atoms'))
        atoms = []
        for line in lines:
            atom = V2_ATOM.match(line)
            if not atom:
                logging.warn("Couldn't parse atom line: {}.".format(line))
            else:
                atoms.append(Atom(atom))
            if len(atoms) == atoms_count:
                break
        else:
            logging.warn("More atoms declared than could be parsed (counts: "
                         "{}, found: {}).".format(atoms_count, len(atoms)))
        logging.info("Atoms: {}.".format(atoms))

        return Molecule(atoms)
