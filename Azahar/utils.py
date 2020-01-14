"""
Functions to compute several glycan properties
The currently available functions are:
Radius of gyration
Ramachandran plot
"""
from __future__ import division
try:
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
except ImportError:
    print('matplotlib not found')
import sys
if sys.version_info[0] < 3:
    import Tkinter
else:
    import tkinter as Tkinter
from pymol.cgo import *
import pymol
from pymol import cmd, stored
import numpy as np
import os
path = os.path.dirname(__file__)
sys.path.append(path)
from torsionals import get_phi, get_psi


def analyse(type_analysis, selection, from_state, to_state, step, visual,
            by_state,
            cut_off):
    if type_analysis == 'Rama scatter':
        rama_plot(selection, from_state, to_state, step, scatter=True)
    elif type_analysis == 'Rama hex':
        rama_plot(selection, from_state, to_state, step, scatter=False)
    elif type_analysis == ' Rg':
        r_gyration(selection, from_state, to_state, step, visual, by_state)
    elif type_analysis == 'Hydrogen_bonds':
        hydro_pairs(selection, cut_off)


def pose_from_pdb(pdb_file):
    """
    Obtain residue indexes
    """
    stored.ResiduesNumber = []
    cmd.iterate('name C1', 'stored.ResiduesNumber.append((resi))')
    if stored.ResiduesNumber:
        first = int(stored.ResiduesNumber[0])
        last = first + len(stored.ResiduesNumber)
        return first, last
    else:
        print("GlycanNotFound", "There is no glycan molecule \
or the atoms in your molecule have a non-standard nomenclature")
        return None, None

def get_glyco_bonds(first, last):
    """
    Obtain glycosidic bonds from a pymol object
    """
    stored.nb = []
    for res_i in range(first, last):
        # TODO In the future we should be able to deal with glyco-conjugates!
        cmd.iterate("not polymer and (neighbor resi %s and chain %s)" %
                    res_i,
                    'stored.nb.append((%s, int(resi), name[-1], resn))' %
                    res_i)
    return stored.nb

def get_glyco_bonds_within_chain_and_model(first, last, chain, model):
    """
    Obtain glycosidic bonds from a pymol object
    """
    stored.nb = []
    for res_i in range(first, last):
        # TODO In the future we should be able to deal with glyco-conjugates!
        cmd.iterate("not polymer and chain %s and (neighbor resi %s and chain %s)" %
                    (chain, res_i, chain),
                    'stored.nb.append((%s, int(resi), name[-1], resn))' %
                    res_i)
    return stored.nb


def writer(bonds):
    """
    Write conectivity matrix
    """
    con_matrix = []
    for i, ibond in enumerate(bonds):
        for jbond in bonds[i + 1:]:
            if ibond[0] == jbond[1] and ibond[1] == jbond[0]:
                con_matrix.append((jbond[0], jbond[3], ibond[0],
                                   ibond[3], int(ibond[2]), int(jbond[2])))
    return con_matrix


def r_gyration(selection='all', from_state=1, to_state=1, step=1, visual=True,
               by_state=True):
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
    states = cmd.count_states(selection)
    if states:
        fd = open('Rg.dat', 'w')
        radii = []
        centers = []
        for state in range(from_state, to_state + 1, step):
            model = cmd.get_model(selection, state)
            xyz = np.array(model.get_coord_list())
            center = np.average(xyz, 0)
            centers.append(center)
            rg = np.mean(np.sum((xyz - center)**2, 1))** 0.5
            fd.write('%9d%8.2f\n' % (state, rg))
            radii.append(rg)
        try:
            rg_mean = sum(radii) / len(radii)
            centers_mean = sum(centers) / len(centers)
        except ZeroDivisionError:
            rg_mean = np.nan
            centers_mean = np.nan
        fd.write('Rg_mean = %8.2f\n' % rg_mean)
        fd.close()
        print('Rg_mean = %8.2f\n' % rg_mean)
    
        if visual:
            cmd.delete('sphere_rg')
            r, g, b = 0, 0, 1
            if by_state:
                cmd.set('defer_updates', 'on')
                count = 0
                for state in range(from_state, to_state + 1, step):
                    x1, y1, z1 = tuple(centers[count])
                    radius = radii[count]
                    obj = [COLOR, r, g, b, SPHERE, x1, y1, z1, radius]
                    cmd.load_cgo(obj, 'sphere_rg', state)
                    count += 1
                cmd.set('defer_updates', 'off')
                # workaround. find a better way to fix the cgo persistent for
                # the last states
                for i in range(state + 1, to_state + 1):
                    cmd.load_cgo([], 'sphere_rg', i)
            else:
                x1, y1, z1 = tuple(centers_mean)
                radius = rg_mean
                obj = [COLOR, r, g, b, SPHERE, x1, y1, z1, radius]
                cmd.load_cgo(obj, 'sphere_rg')


def rama_plot(selection='all', from_state=1, to_state=1, step=1, scatter=True):
    """
    Makes a scatter plot with the phi and psi angle pairs
    """
    first, last = pose_from_pdb(selection)
    if first or last:
        bonds = get_glyco_bonds(first, last)

        con_matrix = writer(bonds)

        phi = []
        psi = []
        for state in range(from_state, to_state + 1, step):
            for element in con_matrix:
                phi.append(get_phi(selection, element, state))
                psi.append(get_psi(selection, element, state))

        if scatter:
            plt.scatter(phi, psi)
        else:
            gridsize = 100
            #gridsize = int(2*len(phi)**(1/3))
            # if gridsize < 36:
            #    gridsize = 36
            plt.hexbin(phi,
                       psi,
                       gridsize=gridsize,
                       cmap=plt.cm.summer,
                       mincnt=1)

        plt.xlabel('$\phi$', fontsize=16)
        plt.ylabel('$\psi$', fontsize=16, rotation=0)
        plt.xlim(-180, 180)
        plt.ylim(-180, 180)
        plt.show()


def hydro_pairs(selection, cut_off):
    """
    Find hydrogen bonds for a given selection
    """
    states = cmd.count_states(selection)
    if states:
        hb = []
        for i in range(1, states + 1):
            p = cmd.find_pairs("%s and (name O* and neighbor elem H) or"
                               "(name N* and neighbor elem H)"  % selection, 
                               "%s and name N* or name O*" % selection,
                               cutoff=3.5,
                               mode=1,
                               angle=70,
                               state1=i,
                               state2=i)
            hb.append(p)

        if any([j for i in hb for j in i]):
            seen = {}
            for state in hb:
                for pairs in state:
                    if pairs in seen:
                        seen[pairs] = seen[pairs] + 1
                    else:
                        seen[pairs] = 1
            occurrence = seen.items()
            occurrence = sorted(occurrence, key=lambda x: x[1], reverse=True)
        
            fd = open("Hydrogen_bonds.dat", "w")
        
            fd.write("--------------------------------------------------------"
                     "------------\n")
            fd.write("-----------------------------Results--------------------"
                     "------------\n")
            fd.write("--------------------------------------------------------"
                     "------------\n")
            fd.write("            Donor             |"
                     "            Aceptor           |\n")
            fd.write("   Object  Residue  Atom Index|"
                     "   Object  Residue  Atom Index| % occurrence\n")
        
            stored.donors = []
            stored.aceptors = []
            sub = []
        
            occurrence_list = []
            for i in range(len(occurrence)):
                cmd.iterate("index %s" % (occurrence[i][0][0][1]),
                            "stored.donors.append((resn, elem))")
                cmd.iterate("index %s" % (occurrence[i][0][1][1]),
                            "stored.aceptors.append((resn, elem))")
                if (occurrence[i][1] * 100 / states) >= cut_off:
                    fd.write("%8s%8s%6s%8s|%8s%8s%6s%8s|%5s\n" %
                             (occurrence[i][0][0][0],
                              stored.donors[i][0],
                              stored.donors[i][1],
                              occurrence[i][0][0][1],
                              occurrence[i][0][1][0],
                              stored.aceptors[i][0],
                              stored.aceptors[i][1],
                              occurrence[i][0][1][1],
                              "%.2f" % (occurrence[i][1] * 100 /states)
                              ))
                    occurrence_list.append(occurrence[i][0])
            fd.close()
        
            cmd.set('suspend_updates', 'on')
            for state, bonds in enumerate(hb):
                for bond in bonds:
                    if bond in occurrence_list:
                        cmd.distance('HB',
                                     'index %s' % bond[0][1],
                                     'index %s' % bond[1][1],
                                     state=state,
                                     cutoff=3.5)
        
            cmd.hide("labels", "HB")
            cmd.set('suspend_updates', 'off')
            print("Check working directory for hydrogen bonds text file report")
        else:
            print("No hydrogen bonds detected for this selection")
