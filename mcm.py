from __future__ import division
import pymol
from pymol import cmd
import openbabel as ob
import numpy as np
import glob, os, sys
path = os.path.dirname(__file__)
sys.path.append(path)
from torsionals import *
from energy import minimize, set_sasa, get_sasa
from utils import pose_from_pdb, get_glyco_bonds, writer



def mcm_run(pose, mc_steps, SASA):
    #cmd.set('defer_updates', 'on')
    cmd.feedback('disable', 'all', 'everything')   ##uncomment for debugging
    cmd.set('pdb_conect_all', 1)
    ################################# MCM Parameters ##########################
    T = 300. # Temperature 
    k = 0.0019872041 # Boltzmann constant
    NRG_old = 1E6 # Ridicously large initial energy, always accep first step
    angles_prob = [1/3, 1/3, 1/3] # probability to sample phi, psi or chi
    ############################################################################

    first, last = pose_from_pdb(pose)
    glyco_bonds = get_glyco_bonds(first, last)
    con_matrix = writer(glyco_bonds)
    
    if SASA:
        params, points, const = set_sasa(n=1000)


    ## randomize initial conformation
    for i in range(len(con_matrix)-1):
        bond = con_matrix[i]
        angle_values = np.random.uniform(-180, 180, size=2)
        set_psi(pose, bond, angle_values[0])
        set_phi(pose, bond, angle_values[1])
        for i in range(6):
            set_chi(pose, bond)
    cmd.save('init.pdb') # XXX
        

    # Remove previous pdb files
    prev_files = glob.glob('test_*.pdb')
    try:
        os.remove('result.pdb')
    except:
        pass
    for prev_file in prev_files:
        os.remove(prev_file)


    ## MCM
    fd = open("log.txt", "w")
    count = 0
    for i in range(0, mc_steps):
        random_angle = np.random.choice(['phi', 'psi', 'chi'], p=angles_prob)
        random_res = np.random.random_integers(0, len(con_matrix)-1)
        bond = con_matrix[random_res]
        cmd.copy('tmp', pose)
        if random_angle == "phi":
            phi = get_phi('tmp', bond)
            angle_value = np.random.normal(phi, 30)
            set_phi('tmp', bond, angle_value)
        elif random_angle == "psi":
            psi = get_psi('tmp', bond)
            angle_value = np.random.normal(psi, 30)
            set_psi('tmp', bond, angle_value)
        else:
            set_chi('tmp', bond)
            set_chi('tmp', bond)
        NRG_new = minimize('tmp', nsteps=100, rigid_geometry=False)
        if SASA:
            solvatation_nrg = get_sasa(params, points, const, selection='all', probe=0)[0]
            NRG_new = NRG_new + solvatation_nrg
        if NRG_new < NRG_old:
            NRG_old = NRG_new 
            fd.write('%6d%10.2f%6s\n' % (count, NRG_new, random_angle))
            cmd.copy(pose, 'tmp')
            cmd.delete('tmp')
            cmd.save('test_%05d.pdb' % count)
            count += 1
        else:
            delta = np.exp(-(NRG_new-NRG_old)/(T*k))
            if delta > np.random.uniform(0, 1):
                NRG_old = NRG_new
                fd.write('%6d%10.2f%6s\n' % (count, NRG_new, random_angle))

                cmd.copy(pose, 'tmp')
                cmd.delete('tmp')
                cmd.save('test_%05d.pdb' % count)
                count += 1
        cmd.delete('tmp')
    fd.close()

    cmd.delete('all')
    # Save all conformations into a single file
    for i in range(0, count):
        cmd.load('test_%05d.pdb' % i)
        cmd.create('trace', 'test_%05d' % i, 1, i+1)
        cmd.delete('test_%05d.pdb' % i)
    cmd.save('result.pdb', 'trace', state=0)
    cmd.delete('all')
    cmd.load('result.pdb')
    #cmd.set('defer_updates', 'off')
