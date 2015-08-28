import __main__
__main__.pymol_argv = ['pymol','-qc']
import pymol
from pymol import cmd 
import os, sys
pymol.finish_launching()

cmd.load('1.pdb')
offset = 100
cmd.alter('1', 'ID = str(int(ID)+100)')
cmd.save('1', '1')


