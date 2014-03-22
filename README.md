# Molecule moments of inertia calculator

A small utility that can take a .mol (or .sdf) file and calculate
principal moments of inertia of the described molecule.


### Requirements

[Python](https://www.python.org/) (tested with 2.7 and 3.2) and
[NumPy](http://www.numpy.org/) (tested with 1.8).

Windows users may prefer to install one of
[distributions with NumPy](http://www.scipy.org/install.html).


### Command line usage

    $ python -m mmoi_calc <ctab-file>


### Running tests

    $ python -m unittest discover
