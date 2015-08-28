"""
Functions to compute several glycan properties
The currently available functions are:
Radius of gyration
Ramachandran plot
"""
from __future__ import division  
import pymol
from pymol import cmd, stored
import numpy as np
import matplotlib.pyplot as plt
import os, sys
path = os.path.dirname(__file__)
sys.path.append(path)
from torsionals import get_phi, get_psi


def pose_from_pdb(pdb_file):
    """
    Obtain residue indexes 
    """
    stored.ResiduesNumber = []
    cmd.iterate('name c1', 'stored.ResiduesNumber.append((resi))')
    first = int(stored.ResiduesNumber[0])
    last = first+len(stored.ResiduesNumber)
    return first, last
    
    
def get_glyco_bonds(first, last):
    """
    Obtain glycosidic bonds from a pymol object
    """       
    stored.nb = []
    for res_i in range(first, last):
        # TODO In the future we should be able to deal with glyco-conjugates!
        cmd.iterate( "not polymer and (neighbor resi %s)" % res_i,  'stored.nb.append((%s, int(resi), name[-1], resn))' % res_i)
    return stored.nb
    

def writer(bonds):
    """
    Write conectivity matrix
    """
    con_matrix = []
    for i, ibond in enumerate(bonds):    
        for jbond in bonds[i+1:]:
            if ibond[0] == jbond[1] and ibond[1] == jbond[0]:
                con_matrix.append((jbond[0], ibond[3], ibond[0],
                jbond[3], int(ibond[2]), int(jbond[2])))        
    return con_matrix


def r_gyration(selection='all', from_state=1, to_state=1,  visual=True, by_state=True):
    """
    Calculates radius of gyration for a PyMOL object
    
    Arguments:
    -------------------------------------------------
    selection: key word. Selects and object 
    from_state: int. First state to calculate RG
    to_state: int. Last state to calculate RG
    visual: boolean. 
    by_state: boolean.
    """
    fd = open('Rg.dat', 'w')
    radii = []
    centers = []
    for state in range(from_state, to_state+1):
        model = cmd.get_model(selection)
        xyz = np.array(model.get_coord_list())
        center = np.average(xyz, 0)
        centers.append(center)
        rg = np.sqrt(np.sum((xyz - center)**2)/len(xyz))
        fd.write('%9d%8.2f\n'% (state, rg))
        radii.append(rg)
    rg_mean = sum(radii)/len(radii)
    centers_mean = sum(centers)/len(centers)
    fd.write('Rg_mean = %8.2f\n'% rg_mean)   
    fd.close()
    print 'Rg_mean = %8.2f\n' % rg_mean

    if visual:
        if by_state:
            cmd.set('defer_updates', 'on')
            for state in range(from_state, to_state+1):
                cmd.pseudoatom('sphere_rg', pos=tuple(center), vdw=radii[state-1],
                 color='blue', state=state)
            cmd.set('defer_updates', 'off')
        else:
            cmd.pseudoatom('sphere_rg', pos=tuple(centers_mean), vdw=rg_mean,
             color='blue')
        cmd.show("spheres", 'sphere_rg')


def rama_plot(selection='all', from_state=1, to_state=1):
    """ 
    Makes a scatter plot with the phi and psi angle pairs
    """  
    first, last = pose_from_pdb(selection)
    bonds = get_glyco_bonds(first, last)
    
    con_matrix = writer(bonds)
    
    phi = []
    psi = []
    for state in range(from_state, to_state+1):
        for element in con_matrix:
            phi.append(get_phi(selection, element))
            psi.append(get_psi(selection, element))
    
    plt.figure(figsize=(6,6))
    plt.scatter(phi, psi)
    plt.xlabel('$\phi$', fontsize=16)
    plt.ylabel('$\psi$', fontsize=16, rotation=0)
    plt.xlim(-180, 180) 
    plt.ylim(-180, 180)
    plt.show()

