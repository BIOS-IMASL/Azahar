from __future__ import division
import pymol
from pymol import cmd, stored
from pymol.cgo import *
import numpy as np
import time
from utils import get_glyco_bonds_within_chain_and_model, writer

def elapsed(starting, s):
    e = time.time() - starting
    print(s, e)
    return e

def find_rings(resn_list, chain, model):
    """determine wich atoms define the sugar rings"""
    matrix_rings = []
    start = time.time()
    for resi in resn_list:
        ring = []
        # identify the oxygens that belong to the ring
        o_atm = cmd.get_model(f'model {model} and chain {chain} and resi {resi} and name C1 extend 1 and (resi {resi} and name O* and not name O1*)')
        ring.append(o_atm.atom[0].name)
        for carbon in range(1, 10):
            if cmd.select(
                'tmp', 'not hydrogen and (neighbor (model %s and chain %s and resi %s and name c%s ))' %
                    (model, chain, resi, carbon)) > 2:
                ring.append('C%s' % carbon)
        while True:
            if cmd.select('tmp', 'model %s and chain %s and resi %s and name %s extend 1 and name %s' % (
                    model, chain, resi, ring[0], ring[-1])) == 0:
                ring.pop()
            else:
                break
        matrix_rings.append(ring)
    return matrix_rings


def get_ring_coords(resn_list, matrix, chain, model):
    """obtain coordinates of sugar rings"""
    matrix_coords = []
    for state in range(1, cmd.count_states() + 1):
        coords = []
        for i, resi in enumerate(resn_list):
            stored.coords = []
            cmd.iterate_state(state, 'model %s and chain %s and resi %s and name %s' % (
                model, chain, resi,  '+'.join(matrix[i])), 'stored.coords.append([x,y,z])')
            coords.append(stored.coords)
        matrix_coords.append(coords)
    return matrix_coords


def get_bonds_coords(resn_list, matrix, chain, model):
    """obtain coordinates of the atoms in the glycosidic bond"""
    matrix_coords = []
    for state in range(1, cmd.count_states() + 1):
        coords = []
        for bond in matrix:
            stored.pos = []
            #Ex bond:
            # (3101, 'NAG', 3100, 'NAG', 1, 4)
            # (475, 'MAN', 471, 'MAN', 1, 6) Exocyclic
            if bond[4] == 6:
                cmd.iterate_state(
                    state, 'model %s and chain %s and resi %s and name C%s or model %s and chain %s and resi %s and name C%s' %
                    (model, chain, bond[0], 5, model, chain, bond[2], bond[5]), 'stored.pos.append((x,y,z))')
            elif bond[5] == 6:
                cmd.iterate_state(
                    state, 'model %s and chain %s and resi %s and name C%s or model %s and chain %s and resi %s and name C%s' %
                    (model, chain, bond[0], bond[4], model, chain, bond[2], 5), 'stored.pos.append((x,y,z))')
            else:
                cmd.iterate_state(
                    state, 'model %s and chain %s and resi %s and name C%s or model %s and chain %s and resi %s and name C%s' %
                    (model, chain, bond[0], bond[4], model, chain, bond[2], bond[5]), 'stored.pos.append((x,y,z))')
            x1, y1, z1 = stored.pos[0]
            x2, y2, z2 = stored.pos[1]
            coords.append((x1, y1, z1, x2, y2, z2))
        matrix_coords.append(coords)
    return matrix_coords


def get_colors_c1(resn_list, color, chain, model):
    """obtain colors of c1 atoms"""
    matrix_colors = []
    if color == 'auto':
        for state in range(1, cmd.count_states() + 1):
            colors = []
            for i, resi in enumerate(resn_list):
                stored.colors = []
                cmd.iterate_state(
                    state, 'model %s and chain %s and resi %s and name C1' %
                           (model, chain, resi), 'stored.colors.append(color)')
                colors.extend(stored.colors)
            matrix_colors.append(colors)
    else:
        for state in range(1, cmd.count_states() + 1):
            matrix_colors.append([color] * len(resn_list))
    return matrix_colors


def get_bonds_colors(resn_list, matrix, color, chain, model):
    """obtain colors for the bonds"""
    matrix_colors = []
    if color == 'auto':
        for state in range(1, cmd.count_states() + 1):
            colors = []
            for bond in matrix:
                stored.colors = []
                cmd.iterate_state(
                    state, 'model %s and chain %s and resi %s and name C1 or model %s and chain %s and resi %s  and name C1' %
                    (model, chain, bond[0], model, chain, bond[2]), 'stored.colors.append(color)')
                colors.append((stored.colors[0], stored.colors[1]))
            matrix_colors.append(colors)
    else:
        for state in range(1, cmd.count_states() + 1):
            matrix_colors.append([(color, color)] * (len(resn_list) - 1))
    return matrix_colors


def hexagon(obj, coords, colors, rep, radius):
    """draw the rings"""
    for idx, coord in enumerate(coords):
        r, g, b = cmd.get_color_tuple(colors[idx])
        coord.append(coord[0])
        coord = np.array(coord)
        mean = np.mean(coord[:-1], axis=0)
        x3, y3, z3 = mean

        # average the normals of two neighbouring triangles
        normals = [0.5 * (
            np.cross((coord[i] - coord[i - 1]), (mean - coord[i])) +
            np.cross((coord[i + 1] - coord[i]), (mean - coord[i])))
            for i in range(len(coord) - 1)]

        centernormal = np.mean(normals, axis=0).tolist()
        # add first value to be able to cycle trought the triangles
        normals.append(normals[0])

        for i in range(len(coord) - 1):
            x1, y1, z1 = coord[i]
            x2, y2, z2 = coord[i + 1]
            # Triangles
            if rep == 'cartoon':
                tri = [
                    BEGIN, TRIANGLES,
                    COLOR, r, g, b,
                    NORMAL, normals[i][0], normals[i][1], normals[i][2],
                    NORMAL, normals[i + 1][0], normals[i + 1][1], normals[i + 1][2],
                    NORMAL, centernormal[0], centernormal[1], centernormal[2],
                    VERTEX, x1, y1, z1,
                    VERTEX, x2, y2, z2,
                    VERTEX, x3, y3, z3,
                    END
                ]
                obj.extend(tri)
            obj.extend([CYLINDER, x1, y1, z1, x2, y2,
                        z2, radius, r, g, b, r, g, b])
            obj.extend([COLOR, r, g, b, SPHERE, x1, y1, z1, radius])
    return obj


def beads(obj, coords, colors, radius):
    """draw beads"""
    for idx, coord in enumerate(coords):
        r, g, b = cmd.get_color_tuple(colors[idx])
        x1, y1, z1 = np.mean(coord, axis=0)
        obj.extend([COLOR, r, g, b, SPHERE, x1, y1, z1, radius])
    return obj


def cylinder(obj, coords, colors, radius):
    """draw the bonds between rings"""
    for idx, coord in enumerate(coords):
        r1, g1, b1 = cmd.get_color_tuple(colors[idx][0])
        r2, g2, b2 = cmd.get_color_tuple(colors[idx][1])
        x1, y1, z1, x2, y2, z2 = coord
        obj.extend([CYLINDER, x1, y1, z1, x2, y2, z2,
                    radius, r1, g1, b1, r2, g2, b2])
    return obj

@cmd.extend
def cartoonize(chains, color='auto', rep='cartoon', show_bonds=False, model = "all"):
    """draw a cartoon representation of glycans, iterate over models and all chains."""

    obj_names = []
    stored.all_models = []

    cmd.iterate('(all)', 'stored.all_models.append((model))')
    models = set(stored.all_models)

    if model.lower() == "all":
        for model in models:
            for chain in chains:
                obj_names.extend(cartoonize_model_and_chain(model, chain, color, rep, show_bonds))
    elif model in models:
        for chain in chains:
            obj_names.extend(cartoonize_model_and_chain(model, chain, color, rep, show_bonds))
    else:
        print("model", model, "not in list of models:",models)
        return
    obj_names.sort()
    for obj_name in obj_names:
        cmd.group("azahar", obj_name)

@cmd.extend
def cartoonize_model_and_chain(model, chain, color='auto', rep='cartoon', show_bonds=False):
    """draw a cartoon representation of glycans"""
    reps = ['cartoon', 'wire', 'beads']
    colors = ['auto', 'green', 'red', 'blue', 'yellow', 'cyan', 'magenta',
                   'orange', 'grey', 'black', 'white']

    if rep not in reps:
        print("Rep not understood. Available: ", reps)
    if color not in colors:
        print("color not understood. Available: ", colors)


    start = time.time()
    stored.ResiduesNumber = []
    atoms = cmd.get_model(f'model {model} and chain {chain} and name C1')
    resn_list = [int(at.resi) for at in atoms.atom]

    #con_matrix = writer2(bonds)
    rings = find_rings(resn_list, chain, model)
    rings_coords = get_ring_coords(resn_list, rings, chain, model)
    colors = get_colors_c1(resn_list, color, chain, model)

    cmd.set('suspend_updates', 'on')

    obj = []

    #Set the radius
    if rep == 'beads':
        radius = 0.18
    elif rep == "cartoon":
        radius = 0.075
    else:
        radius = 0.035

    #Show bonds is currently the slowest step.
    ## It doesn't work well with branching (1-6 is too short)
    ## It doesn't work with Glycan-protein connections.
    obj = []
    names = []
    if show_bonds:
        bonds = get_glyco_bonds_within_chain_and_model(resn_list[0], resn_list[-1] + 1, chain, model)
        con_matrix = writer(bonds)

        bonds_colors = get_bonds_colors(resn_list, con_matrix, color, chain, model)
        bonds_coords = get_bonds_coords(resn_list, con_matrix, chain, model)

        for state, coords in enumerate(rings_coords):

            obj = cylinder(
                obj,
                bonds_coords[state],
                bonds_colors[state],
                radius)
            name = 'cbb' + chain + "_" + model
            names.append(name)
            cmd.load_cgo(obj,  name ,state + 1)

    #Independant objects to toggle either one on/off.
    obj = []
    for state, coords in enumerate(rings_coords):

        if rep == 'beads':
            obj = beads(obj, coords, colors[state], 1.8)
        else:
            obj = hexagon(obj, coords, colors[state], rep, radius)


        name = 'cbr' + chain + "_" + model
        names.append(name)
        cmd.load_cgo(obj, name, state + 1)

    cmd.select('glycan', 'byres name C1')
    cmd.delete('glycan')
    cmd.delete('tmp')
    cmd.set('two_sided_lighting', 1)
    cmd.set('suspend_updates', 'off')
    return names
