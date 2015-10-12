import pymol
from pymol import cmd, stored
from pymol.cgo import *
import numpy as np


def get_glyco_bonds(resn_list):
    """
    Obtain glycosidic bonds from a pymol object
    """       
    stored.nb = []
    for res_i in resn_list:
        # TODO In the future we should be able to deal with glyco-conjugates!
        cmd.iterate( "not polymer and (neighbor resi %s)" % res_i,  'stored.nb.append((%s, int(resi), name[1], resn))' % res_i)
    return stored.nb


def writer(bonds):
    """
    Write conectivity matrix
    """
    con_matrix = []
    for i, ibond in enumerate(bonds):
        for jbond in bonds[i+1:]:
            if ibond[0] == jbond[1] and ibond[1] == jbond[0]:
                con_matrix.append((jbond[0], jbond[3], ibond[0],
                 ibond[3], ibond[2], jbond[2]))        
    return con_matrix


def find_rings(resn_list):
    """determine wich atoms define the sugar rings"""
    matrix_rings = []
    for resi in resn_list:
        ring = []
        stored.oxy = []
        #identify the oxygens that belong to the ring
        cmd.iterate('resi %s and name c1 extend 1 and (resi %s and name O* and not name O1*)' % (resi, resi), 'stored.oxy.append(name)')
        ring.append(stored.oxy[0])
        for carbon in range(1, 10):
            if cmd.select('tmp', 'not hydrogen and (neighbor (resi %s and name c%s))' % (resi, carbon)) > 2:
                ring.append('C%s' % carbon)
        while True:
            if cmd.select('tmp', 'resi %s and name %s extend 1 and name %s' % (resi, ring[0], ring[-1])) == 0:
                ring.pop()
            else:
                break
        matrix_rings.append(ring)
    return matrix_rings
    

def get_ring_coords(resn_list, matrix):
    """obtain coordinates of sugar rings"""
    matrix_coords = []
    for state in range(1, cmd.count_states()+1):
        coords = []
        for i, resi in enumerate(resn_list):
            stored.coords = []
            cmd.iterate_state(state, 'resi %s and name %s' % (resi, '+'.join(matrix[i])), 'stored.coords.append([x,y,z])')
            coords.append(stored.coords)
        matrix_coords.append(coords)
    return matrix_coords

 
def get_bonds_coords(resn_list, matrix):
    """obtain coordinates of the atoms in the glycosidic bond"""
    matrix_coords = []
    for state in range(1, cmd.count_states()+1):
        coords = []
        for bond in matrix:
            stored.pos = []
            if bond[4] == '6':
                cmd.iterate_state(state, 'resi %s and name c%s or resi %s and name c%s' % (bond[0], 5, bond[2], bond[5]), 'stored.pos.append((x,y,z))')
            elif bond[5] == '6':
                cmd.iterate_state(state, 'resi %s and name c%s or resi %s and name c%s' % (bond[0], bond[4], bond[2], 5), 'stored.pos.append((x,y,z))')
            else:
                cmd.iterate_state(state, 'resi %s and name c%s or resi %s and name c%s' % (bond[0], bond[4], bond[2], bond[5]), 'stored.pos.append((x,y,z))')
            x1, y1, z1 = stored.pos[0]
            x2, y2, z2 = stored.pos[1]
            coords.append((x1, y1, z1, x2, y2, z2))
            print bond[0], bond[4], bond[2], bond[5]
            
        matrix_coords.append(coords)
    return matrix_coords        


def hexagon(coords):
    """draw the rings"""
    radius = 0.075
    r, g, b = 0, 1, 0    
    for coord in coords:
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

        for i in range(len(coord)-1):
            x1, y1, z1 = coord[i]
            x2, y2, z2 = coord[i+1]

            # Triangles
            tri = [
    		    BEGIN, TRIANGLES,
    		    COLOR, r, g, b,
    		    NORMAL, normals[i][0], normals[i][1], normals[i][2],
    		    NORMAL, normals[i+1][0], normals[i+1][1], normals[i+1][2],
    		    NORMAL, centernormal[0], centernormal[1], centernormal[2],
    		    VERTEX, x1, y1, z1,
    		    VERTEX, x2, y2, z2,
    		    VERTEX, x3, y3, z3,
                END
                ]   
                         
            obj.extend(tri)      
            obj.extend([CYLINDER, x1, y1, z1, x2, y2, z2, radius, r, g, b, r, g, b])
            obj.extend([SPHERE, x1, y1, z1, radius, COLOR, r, g, b])
    return obj


def cylinder(coords):
    """draw the bonds between rings"""
    r, g, b = 0, 1, 0
    radius = 0.075
    for coord in coords:
        x1, y1, z1, x2, y2, z2 = coord
        obj.extend([CYLINDER, x1, y1, z1, x2, y2, z2, radius, r, g, b, r, g, b])
    return obj


def cartoonize():
    """draw a cartoon representation of glycans"""
    global obj

    stored.ResiduesNumber = []
    cmd.iterate('name c1', 'stored.ResiduesNumber.append((resi))')
    resn_list = [int(i) for i in stored.ResiduesNumber]
    bonds = get_glyco_bonds(resn_list)
    con_matrix = writer(bonds)
    rings = find_rings(resn_list)
    rings_coords = get_ring_coords(resn_list, rings)
    bonds_coords = get_bonds_coords(resn_list, con_matrix)

    cmd.set('suspend_updates', 'on')
    for state, coords in enumerate(rings_coords):
        obj = []
        obj = hexagon(coords)
        obj = cylinder(bonds_coords[state])
        cmd.load_cgo(obj,'cgo01', state+1)

    cmd.select('glycan', 'byres name C1')
    cmd.hide('everything', 'glycan')
    cmd.delete('glycan')
    cmd.delete('tmp')
    cmd.set('two_sided_lighting', 1)
    cmd.set('suspend_updates', 'off')
