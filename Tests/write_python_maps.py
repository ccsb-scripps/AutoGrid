from PyAutoDock.AutoGrid import AutoGrid4
from MolKit import Read

mol = Read("hsg1_sm.pdbqt")[0]
mol.buildBondsByDistance()

ag = AutoGrid4(mol, atom_types = ['A','C','HD','NA', 'N','OA' ], npts =
[5,5,5], center= [2.5, 6.5, -7.5])
ag.write_maps()

