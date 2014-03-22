# -*- coding: utf-8 -*-
import logging
import re

from elements import ELEMENTS
import numpy


# Connection table header (for the legacy V2000 format).
# Reference: http://download.accelrys.com/freeware/ctfile-formats/.
v2_counts = re.compile(r'''\s*
    (?P<atoms>\d+)\s+             # Number of atoms.
    (?P<bonds>\d+)\s+             # Number of bond lines.
    (?P<lists>\d+)\s+             # Number of atom lists.
    \d+\s+                        # Obsolete.
    (?P<chiral>[01])\s+           # Chirality flag (1 = chiral).
    (?P<stexts>\d+)\s+            # Structural text lines.
    (?:\d+\s+){,4}                # Obsolete x4.
    (?P<props>\d+)\s+             # Number of additional property lines.
    V2000\s*                      # Version identifier.
''', re.VERBOSE)

# Atom line of the legacy connection table.
v2_atom = re.compile(r'''\s*
    (?P<x>-?\d+(?:\.\d+)?)\s+     # X coordinate.
    (?P<y>-?\d+(?:\.\d+)?)\s+     # Y coordinate.
    (?P<z>-?\d+(?:\.\d+)?)\s+     # Z coordinate (resulting bond lengths look
                                  # like in Angstroms, but that doesn't seem
                                  # to be documented anywhere).
    (?P<symbol>\w+|R\#\w+|\*)\s+  # Atom symbol, group label or
                                  # query -- L, A, Q, LP or *.
    (?P<mass>-?\d+)\s+            # Difference from mass in the periodic
                                  # table; should be within the -3..4 range.
    (?P<charge>[0..7])\s+         # Charge (1 = +3, 2 = +2, 3 = +1, 4 = double
                                  # radical, 5 = -1, 6 = -2, 7 = -3).
    (?P<stereo>[0..3])\s+         # Stereo parity (1 = odd, 2 = even,
                                  # 3 = either or unmarked center).
    (?P<hydrogen>[0..5])\s+       # Hydrogen count (in excess of explicit;
                                  # queries only).
    (?P<stereo_care>[01])\s+      # Whether to consider stereo for double
                                  # bonds (1 = match; queries only).
    (?P<valence>\d+)\s+           # Number of bonds including implied
                                  # hydrogens (0 = default, 15 = 0).
    (?P<no_hydrogen>[01])\s+      # Legacy zero additional hydrogen designator
                                  # (1 = no hydrogen allowed).
    \d+\s+\d+\s+                  # Unused x2.
    (?P<atom_atom>\d+)\s+         # Atom-atom mapping number (reactions only).
    (?P<conf>[0..2])\s+           # Inversion / retention flag (1 = inverted,
                                  # 2 = retained; reactions only).
    (?P<exact>[01])\s*            # Exact change flag (1 = changes must be
                                  # exact).
''', re.VERBOSE)


class Atom:
    def __init__(self, match):
        """
        Atom object from a regex match.
        
        Drops in element data under ``element`` property and
        sets ``mass`` taking into account mass difference data.
        
        TODO: Additional properties may also include mass modifications.
        """
        self.symbol = match.group('symbol')
        self.element = ELEMENTS[self.symbol]
        self.coords = numpy.array(match.group('x', 'y', 'z'), dtype=float)
        self.mass = self.element.mass + int(match.group('mass'))
        
    def __str__(self):
        return self.symbol

    def __repr__(self):
        return "{} {}".format(self.symbol, self.coords)
        
        
class Molecule:
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
        
        With ``pricipal`` set to false returns inertia matrix in the data
        coordinates (with non-zero products of inertia), with ``moments_only``
        set to false, returns principal moments and axes (as matrix columns).
        
        See: https://en.wikipedia.org/wiki/Moment_of_inertia
            #The_inertia_matrix_for_spatial_movement_of_a_rigid_body
        """
        com = self.center_of_mass()
        inertia = numpy.zeros((3, 3))
        for atom in self.atoms:
            r = atom.coords - com
            inertia += - atom.mass * numpy.matrix(
                [[ 0,   -r[2], r[1]],
                 [ r[2], 0,   -r[0]],
                 [-r[1], r[0], 0   ]])**2
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
    """Signals a critical parsing failure."""
    pass

        
class Parser:
    def molfile(self, file):
        """
        Takes a file containing a V2000 ctab and returns a ``Molecule``.
        
        Argument should be an open ``File`` object (or anything implementing 
        the same interface like a ``StringIO``).
        """
        # Find the "counts" line.
        for line in file:
            counts = v2_counts.match(line)
            if counts:
                break
        else:
            raise ParseError("Couldn't find the ctab header; "
                             "is this a V2000 file?")
            
        # Parse atoms data (assumed to follow).
        atoms_count = int(counts.group('atoms'))
        atoms = []
        for line in file:
            atom = v2_atom.match(line)
            if not atom:
                logging.warn("Couldn't parse atom line: {}.", line)
            else:
                atoms.append(Atom(atom))
            if len(atoms) == atoms_count:
                break
        else:
            logging.warn("More atoms declared than could be parsed (counts: "
                         "{}, found: {}).".format(atoms_count, len(atoms)))
        logging.info("Atoms: {}.".format(atoms))
        
        # TODO: Parse bonds etc.
        
        return Molecule(atoms)
