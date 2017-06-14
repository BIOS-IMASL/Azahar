"""
Functions to compute and change torsionals angles for glycans
"""
import pymol
from pymol import cmd
import numpy as np


def get_phi(obj, bond, state=1):
    """
    Get the dihedral angle phi (OR-C1-O'x-C'x)
    """
    resi_i, resi_j, atom_i, atom_j = bond[0], bond[2], bond[4], bond[5]

    if atom_i > atom_j:
        atom_j = atom_i
        resi_i, resi_j = resi_j, resi_i
    atom1 = '%s and resi %s and name OR' % (obj, resi_i)
    atom2 = '%s and resi %s and name C1' % (obj, resi_i)
    atom3 = '%s and resi %s and name O%s' % (obj, resi_j, atom_j)
    atom4 = '%s and resi %s and name C%s' % (obj, resi_j, atom_j)
    value = cmd.get_dihedral(atom1, atom2, atom3, atom4, state=state)
    return value


def get_psi(obj, bond, state=1):
    """
    Get the dihedral angle psi (C1-O'x-C'x-C'x-1)
    """
    resi_i, resi_j, atom_i, atom_j = bond[0], bond[2], bond[4], bond[5]
    if atom_i > atom_j:
        atom_j = atom_i
        resi_i, resi_j = resi_j, resi_i
    atom1 = '%s and resi %s and name C1' % (obj, resi_i)
    atom2 = '%s and resi %s and name O%s' % (obj, resi_j, atom_j)
    atom3 = '%s and resi %s and name C%s' % (obj, resi_j, atom_j)
    atom4 = '%s and resi %s and name C%s' % (obj, resi_j, (atom_j - 1))
    value = cmd.get_dihedral(atom1, atom2, atom3, atom4, state=state)
    return value


def get_omega(obj, bond, state=1):
    """
    Get the dihedral angle omega (O6-C6-C5-C4)
    """
    resi_i, resi_j, atom_i, atom_j = bond[0], bond[2], bond[4], bond[5]
    atom1 = '%s and resi %s and name O6' % (obj, resi_j)
    atom2 = '%s and resi %s and name C6' % (obj, resi_j)
    atom3 = '%s and resi %s and name C5' % (obj, resi_j)
    atom4 = '%s and resi %s and name C4' % (obj, resi_j)
    value = cmd.get_dihedral(atom1, atom2, atom3, atom4, state=state)
    return value


def get_chi():
    pass


def set_phi(obj, bond, angle, state=1):
    """
    Set the dihedral angle phi (OR-C1-O'x-C'x)
    """
    resi_i, resi_j, atom_i, atom_j = bond[0], bond[2], bond[4], bond[5]
    if atom_i > atom_j:
        atom_j = atom_i
        resi_i, resi_j = resi_j, resi_i
    atom1 = '%s and resi %s and name OR' % (obj, resi_i)
    atom2 = '%s and resi %s and name C1' % (obj, resi_i)
    atom3 = '%s and resi %s and name O%s' % (obj, resi_j, atom_j)
    atom4 = '%s and resi %s and name C%s' % (obj, resi_j, atom_j)
    cmd.set_dihedral(atom1, atom2, atom3, atom4, angle, state=state)


def set_psi(obj, bond, angle, state=1):
    """
    Set the dihedral angle psi (C1-O'x-C'x-C'x-1)
    """
    resi_i, resi_j, atom_i, atom_j = bond[0], bond[2], bond[4], bond[5]
    if atom_i > atom_j:
        atom_j = atom_i
        resi_i, resi_j = resi_j, resi_i
    atom1 = '%s and resi %s and name C1' % (obj, resi_i)
    atom2 = '%s and resi %s and name O%s' % (obj, resi_j, atom_j)
    atom3 = '%s and resi %s and name C%s' % (obj, resi_j, atom_j)
    atom4 = '%s and resi %s and name C%s' % (obj, resi_j, (atom_j - 1))
    cmd.set_dihedral(atom1, atom2, atom3, atom4, angle, state=state)


def set_omega(obj, bond, angle, state=1):
    """
    Get the dihedral angle omega (O6-C6-C5-C4)
    """
    resi_i, resi_j, atom_i, atom_j = bond[0], bond[2], bond[4], bond[5]
    atom1 = '%s and resi %s and name O6' % (obj, resi_j)
    atom2 = '%s and resi %s and name C6' % (obj, resi_j)
    atom3 = '%s and resi %s and name C5' % (obj, resi_j)
    atom4 = '%s and resi %s and name C4' % (obj, resi_j)
    cmd.set_dihedral(atom1, atom2, atom3, atom4, angle, state=state)


def set_chi(obj, bond, state=1):
    """
    Set the dihedral angles chi (Cx-1-Cx-O'x-H'x) to one of tree predefinied
    angles (-60, 60, 180). The perturbed chi angle is selected at random.
    This function considers the omega angle as a chi angle.
    """
    resi_i, resi_j = bond[0], bond[2]
    chises = [1, 2, 3, 4, 5, 6]
    np.random.shuffle(chises)
    angle = np.random.choice([-60, 60, 180])
    resi_i = np.random.choice([resi_i, resi_j])
    for chi in chises:
        if chi == 1:
            try:
                cmd.set_dihedral(
                    '%s and resi %s and name OR' %
                    (obj, resi_i), '%s and resi %s and name C%S' %
                    (obj, resi_i, chi), '%s and resi %s and name O%s' %
                    (obj, resi_i, chi), '%s and resi %s and name H%so' %
                    (obj, resi_i, chi), angle, state=state)
                break
            except BaseException:
                pass
        elif chi <= 5:
            try:
                cmd.set_dihedral(
                    '%s and resi %s and name C%s' %
                    (obj, resi_i, chi - 1), '%s and resi %s and name C%S' %
                    (obj, resi_i, chi), '%s and resi %s and name O%s' %
                    (obj, resi_i, chi), '%s and resi %s and name H%so' %
                    (obj, resi_i, chi), angle, state=state)
                break
            except BaseException:
                pass
        else:
            try:
                set_omega(obj, bond, angle)
                break
            except BaseException:
                pass
