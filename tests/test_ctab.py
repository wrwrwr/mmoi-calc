from __future__ import unicode_literals
import io
import unittest

from mmoi_calc import ctab


CO2_CTAB = '''
  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C    0  0  0  0  0  0  0  0  0  0  0  0
    1.1600    0.0000    0.0000 O    0  0  0  0  0  0  0  0  0  0  0  0
   -1.1600    0.0000    0.0000 O    0  0  0  0  0  0  0  0  0  0  0  0
'''


class TestParser(unittest.TestCase):
    def test_simple_ctab(self):
        molecule = ctab.Parser().molfile(io.StringIO(CO2_CTAB))
        self.assertIsInstance(molecule, ctab.Molecule)
        self.assertEqual(
            sorted(a.symbol for a in molecule.atoms), sorted(('C', 'O', 'O')))


class TestMolecule(unittest.TestCase):
    def setUp(self):
        self.co2 = ctab.Parser().molfile(io.StringIO(CO2_CTAB))

    def test_center_of_mass(self):
        self.assertEqual(list(self.co2.center_of_mass()), [0, 0, 0])

    def test_inertia(self):
        # Calculated using mass = 15.9994 u and bonds lenghts = 116 pm.
        self.assertEqual(
            list(self.co2.inertia()), [0, 43.05758528, 43.05758528])


if __name__ == '__main__':
    unittest.main()
