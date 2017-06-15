"""
Functions to compute and optimize the energy of molecules
"""
from __future__ import division
from pymol import cmd
import openbabel as ob
import numpy as np
from scipy.spatial.distance import cdist


def minimize(selection='tmp', forcefield='MMFF94', method='steepest descent',
             nsteps=2000, conv=1E-6, cutoff=False, cut_vdw=6.0, cut_elec=8.0,
             rigid_geometry=True):
    """
    Use openbabel to minimize the energy of a molecule.
    """
    pdb_string = cmd.get_pdbstr(selection)
    name = cmd.get_legal_name(selection)
    obconversion = ob.OBConversion()
    obconversion.SetInAndOutFormats('pdb', 'pdb')
    mol = ob.OBMol()
    obconversion.ReadString(mol, pdb_string)
    if rigid_geometry:
        constraints = ob.OBFFConstraints()
        for angle in ob.OBMolAngleIter(mol):
            b, a, c = [mol.GetAtom(x + 1) for x in angle]
            value = mol.GetAngle(a, b, c)
            b, a, c = [x + 1 for x in angle]
            constraints.AddAngleConstraint(a, b, c, value)
        for i in ob.OBMolBondIter(mol):
            a, b = (i.GetBeginAtomIdx(), i.GetEndAtomIdx())
            value = i.GetLength()
            constraints.AddDistanceConstraint(a, b, value)
        ff = ob.OBForceField.FindForceField(forcefield)
        ff.Setup(mol, constraints)
        ff.SetConstraints(constraints)
    else:
        ff = ob.OBForceField.FindForceField(forcefield)
        ff.Setup(mol)
    if cutoff:
        ff.EnableCutOff(True)
        ff.SetVDWCutOff(cut_vdw)
        ff.SetElectrostaticCutOff(cut_elec)
    if method == 'conjugate gradients':
        ff.ConjugateGradients(nsteps, conv)
    else:
        ff.SteepestDescent(nsteps, conv)
    ff.GetCoordinates(mol)
    nrg = ff.Energy()
    pdb_string = obconversion.WriteString(mol)
    cmd.delete(name)
    if name == 'all':
        name = 'all_'
    cmd.delete(selection)
    cmd.read_pdbstr(pdb_string, selection)
    return nrg


def generate_sphere_points(n):
    """
    Returns an array of evenly distributed points on the
    surface of a the unit sphere using the Golden-Section
    Spiral algorithm.

    Parameters
    ----------
    n : int. Number of points
    """

    golden_angle = np.pi * (3 - 5**0.5)

    theta = golden_angle * np.arange(n)
    z = np.linspace(1 - 1.0 / n, 1.0 / n - 1, n)

    radius = np.sqrt(1 - z * z)
    points = np.zeros((n, 3))
    # change from polar to cartesian coordinates
    points[:, 0] = radius * np.cos(theta)
    points[:, 1] = radius * np.sin(theta)
    points[:, 2] = z

    const = 4.0 * np.pi / n

    return points, const


def assign_atom_types(selection='all'):
    """
    TODO document me!

    read http://openbabel.org/dev-api/classOpenBabel_1_1OBAtom.shtml#ae09ed28481ac044dab3f31c8605b44a9
    for available functions provided by openbabel to extract atom properties. There is a GetType() function
    but I found that function very limited.
    """
    atom_types = []
    pdb_string = cmd.get_pdbstr(selection)
    mol = ob.OBMol()
    obconversion = ob.OBConversion()
    obconversion.SetInAndOutFormats('pdb', 'pdb')
    obconversion.ReadString(mol, pdb_string)
    rings = mol.GetSSSR()
    for at in ob.OBMolAtomIter(mol):
        ring_member = [ring.IsMember(at) for ring in rings]
        neighbors = [neighbor.GetAtomicNum()
                     for neighbor in ob.OBAtomAtomIter(at)]
        atom_types.append(
            (at.GetIndex(),
             at.GetAtomicNum(),
             at.GetHvyValence(),
             any(ring_member),
                at.IsAromatic(),
                at.MemberOfRingCount(),
                (neighbors)))
    return atom_types


def assign_params(atom_types):
    """
    Assign solvation parameters for each chemical group. See: DOI: 10.1002/prot.340140112
    """
    SRA_model = {'hydroxyl_carboxyl H': (2.85, -0.0487), 
                 'amine_amide H': (2.85, -0.0487),
                 'thiol H': (2.85, 0.0690),
                 'aliphatic CH3': (0.946, 0.0676),
                 'aliphatic CH2': (0.946, 0.0193),
                 'aliphatic CH': (0.946, 0.00438),
                 'aliphatic_alicyclic C': (0.946, 0.01691),
                 'alicyclic CH2': (0.946, 0.0190),
                 'alicyclic CH': (0.946, 0.0190),
                 'aromatic CH': (0.946, -0.0110),
                 'aromatic C': (0.946, -0.137),
                 'aromatic C of fused ring': (0.946, -0.103),
                 'aromatic CH of fused ring': (0.946, 0.2342),
                 'carbonyl_carboxylic C': (0.946, -0.104),
                 'N primary amine': (4.10, -0.0225),
                 'N secondary amine': (4.10, -0.0213),
                 'N cyclic amine': (4.10, -0.0254),
                 'aromatic N': (4.10, -0.0234),
                 'N amide': (4.10, -0.0337),
                 'ether_hydroxyl O': (2.83, -0.0282),
                 'carboxylic O': (2.83, 0.0394),
                 'carbonyl O': (2.83, -0.0489),
                 'amide carbonyl O': (2.83, -0.0920),
                 'thiol or sulfide S': (7.37, -0.0020)}
    params = {}
    for atom in atom_types:
        at_index, elem, heavy_nb, in_ring, is_arom, ring_membership, neighbors = atom
        if elem == 1:
            if neighbors[0] == 8:
                params[at_index] = SRA_model['hydroxyl_carboxyl H']
            elif neighbors[0] == 7:
                params[at_index] = SRA_model['amine_amide H']
            elif neighbors[0] == 16:
                params[at_index] = SRA_model['thiol H']
        elif elem == 6:
            if neighbors.count(1) == 3:
                params[at_index] = SRA_model['aliphatic CH3']
            elif neighbors.count(1) == 2:
                if in_ring:
                    params[at_index] = SRA_model['alicyclic CH2']
                else:
                    params[at_index] = SRA_model['aliphatic CH2']
            elif neighbors.count(1) == 1:
                if in_ring:
                    params[at_index] = SRA_model['alicyclic CH']
                else:
                    params[at_index] = SRA_model['aliphatic CH']
                if is_arom:
                    if ring_membership > 1:
                        params[at_index] = SRA_model['aromatic CH of fused ring']
                    else:
                        params[at_index] = SRA_model['aromatic CH']
            elif neighbors.count(1) == 0:
                if is_arom:
                    if ring_membership > 1:
                        params[at_index] = SRA_model['aromatic C of fused ring']
                    else:
                        params[at_index] = SRA_model['aromatic C']
                else:
                    params[at_index] = SRA_model['aliphatic_alicyclic C']
            elif neighbors.count(8) > 1:
                params[at_index] = SRA_model['carbonyl_carboxylic C']
        elif elem == 7:
            if neighbors.count(1) == 2:
                if is_arom:
                    params[at_index] = SRA_model['N cyclic amine']
                else:
                    params[at_index] = SRA_model['N primary amine']
            elif neighbors.count(1) == 1:
                params[at_index] = SRA_model['N secondary amine']
            elif is_arom:
                params[at_index] = SRA_model['aromatic N']
            # elif # N de laamide
                #params.append(SRA_model['aromatic N'])
        elif elem == 8:
            if len(neighbors) == 2:
                params[at_index] = SRA_model['ether_hydroxyl O']
            else:
                params[at_index] = SRA_model['carboxylic O']
            # C=O
            # C=O of ester
            # O de la amide
        if elem == 16:
            params[at_index] = SRA_model['thiol or sulfide S']
    return params


def get_neighbors(coord, probe, k, params):
    """
    Returns list of indices of neighbors.

    Parameters
    ----------
    coord : Array (n,3). Cartesian coordinates
    probe : float. Radius of the solvent
    k: int. Residue from which to compute neighbors
    """
    dist = cdist(coord, coord, metric='euclidean')
    neighbor_indices = []
    radius = params[k][0] + probe * 2
    for key, values in params.items():
        if dist[key, k] < radius + values[0]:
            neighbor_indices.append(key)
    return neighbor_indices


def get_sasa(params, points, const, selection='all', probe=1.4):
    """
    Returns the solvent-accessible-surface area and empirical solvatation term

    Parameters
    ----------
    probe : float. Radius of the solvent
    """
    # compute the area each point represents
    areas = []
    energies = []
    model = cmd.get_model(selection)
    coord = np.array(model.get_coord_list())
    for key, values in params.items():
        # scale the points to the correct radius
        radius = values[0] + probe
        points_scaled = coord[key] + points * radius
        # get all the neighbors of the i residue
        neighbors = get_neighbors(coord, probe, key, params)
        # compute the distance between points and neighbors
        d_matrix = cdist(points_scaled, coord[neighbors],
                         metric='euclidean')
        # create a matrix and store the vdW radii for each neighbor
        nb_matrix = np.zeros((len(points), len(neighbors)))
        for nb_i, nb in enumerate(neighbors):
            nb_matrix[:, nb_i] = values[0]
        # compute the number of buried points, we have to be carefull
        # because we have counted how many times a point is buried
        # and we only need how many points are buried
        buried = np.sum(np.sum(d_matrix < nb_matrix + probe, axis=1) > 0)
        exposed = len(points) - buried
        area_per_atom = const * exposed * radius**2
        energy_per_atom = area_per_atom * values[1]
        areas.append(area_per_atom)
        energies.append(energy_per_atom)
    return sum(energies), sum(areas)


def set_sasa(n=1000):
    """
    Sets solvent-accessible-surface area
    """
    #cmd.alter("all", "type='ATOM'")
    atom_types = assign_atom_types()
    params = assign_params(atom_types)
    points, const = generate_sphere_points(n)
    return params, points, const
