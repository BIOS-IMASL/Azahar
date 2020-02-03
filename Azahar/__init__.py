"""
BioModelBuilder
Described at PyMOL wiki: http://www.pymolwiki.org/index.php/Azahar
Authors: Agustina Arroyuelo and Osvaldo Martin
email: agustinaarroyuelo@gmail.com, aloctavodia@gmail.com
Date: august 2015
License: GNU General Public License
Version 0.8
"""


import sys
if sys.version_info[0] < 3:
    import Tkinter
else:
    import tkinter as Tkinter
    
import Pmw
from pymol import cmd, stored
import os
import numpy as np
path = os.path.dirname(__file__)
sys.path.append(path)
db_path = os.path.join(path, 'db_glycans')
from BuildOligo import read_input, builder
from cartoonize import cartoonize
from utils import analyse
from mcm import mcm_run


def __init__(self):
    """ Adds this plugin to the Pymol menu """
    self.menuBar.addmenuitem('Plugin', 'command',
                             'Azahar',
                             label='Azahar',
                             command=lambda: mainDialog(self.root))


def mainDialog(root=None):
    """ Creates the GUI """

    ################# Define some necessary functions #########################
    def add(mol_name):
        """
        Takes the elements selected by the user with the GUI and uses
        them to write a matrix, which specifies the way residues are
        bonded in the final molecule. It assigns an index and a name
        to each residue and a bond to each pair.
        """
        residue0 = selected_res0.get()
        residue1 = selected_res1.get()
        bond = selected_bond.get()
        linear_residues = int(n_residues.get())
        first = int(first_res.get())
        total = int(total_res.get())

        if total == 0 and os.path.isfile('%s_matrix.dat' % mol_name):
            os.remove('%s_matrix.dat' % mol_name)
        conectivity_matrix = open('%s_matrix.dat' % mol_name, 'a')

        for i in range(linear_residues):
            text = '%5s%30s%5s%30s%2s%2s\n' % (first + i,
                                               residue0,
                                               total + i + 1,
                                               residue1,
                                               bond[1], bond[3])
            print('%s' % text),
            conectivity_matrix.write(text)
        total = total + i + 1
        first_res.set(total)
        total_res.set(total)
        conectivity_matrix.close()

    def create(mol_name):
        """ create the defined molecule """
        if os.path.isfile('%s_matrix.dat' % mol_name):
            residues, bonds = read_input('%s_matrix.dat' % mol_name)
            builder(residues, bonds, mol_name)
            to_state.set(cmd.count_states(sel0_value.get()))
            cmd.zoom()
            cmd.util.chainbow(mol_name)
        else:
            print("FileNotFound", "You should add residues\
            before creating a molecule.")

    def reset():
        n_residues.set(1)
        first_res.set(0)
        total_res.set(0)

    def enable_disable_rg():
        """
        enables the checkbutton by_state if rg visualization is enabled and
        if there is more than one state to visualize
        """
        if vis_rg_value.get() and (to_state.get() - from_state.get()) > 1:
            entry_by_state.configure(state='normal')
        else:
            entry_by_state.configure(state='disabled')
        entry_by_state.update()

    def enable_disable_cutoff(items):
        """
        enable the cut_off entry if the hydrogen_bonds computation is selected
        """
        if items == 'Hydrogen_bonds':
            entry_cut_off.configure(state='normal')
        else:
            entry_cut_off.configure(state='disabled')
        entry_by_state.update()

    ############################ Create the GUI ##############################
    master = Tkinter.Toplevel(root)
    master.title("Azahar")
    w = Tkinter.Label(master, text="\nLife is Sweet.\n",
                      background='black',
                      foreground='dark orange',
                      height=1)
    w.pack(expand=1, fill='both', padx=4, pady=4)
    ############################ NoteBook ####################################
    Pmw.initialise()
    nb = Pmw.NoteBook(master, hull_width=400, hull_height=260)
    p1 = nb.add('Creation')
    p2 = nb.add('Visualization')
    p3 = nb.add('Calculations')
    p4 = nb.add('Monte Carlo')
    p10 = nb.add('    About   ')
    nb.pack(padx=5, pady=5, fill=Tkinter.BOTH, expand=1)
    ############################ Creation TAB ################################
    group = Pmw.Group(p1, tag_text='options')
    group.pack(fill='both', expand=1, padx=5, pady=5)
    # Select first residue
    selected_res0 = Tkinter.StringVar(master=group.interior())
    selected_res0.set('a-D-glucose')

    residues_templates = sorted([os.path.splitext(x)[0]
                                 for x in os.listdir(db_path)])
    Pmw.OptionMenu(group.interior(),
                   labelpos='w',
                   label_text='Residue 1',
                   menubutton_textvariable=selected_res0,
                   items=residues_templates,
                   menubutton_width=18,
                   ).grid(row=0, columnspan=2)
    # Select Bond
    selected_bond = Tkinter.StringVar(master=group.interior())
    selected_bond.set("(1,4)")
    Pmw.OptionMenu(group.interior(),
                   labelpos='w',
                   label_text='   Bond    ',
                   menubutton_textvariable=selected_bond,
                   items=["(1,4)", "(1,6)", "(1,3)", "(1,2)", "(1,5)", "(1,1)",
                          "(6,1)", "(4,1)", "(3,1)", "(2,1)", "(5,1)"],
                   menubutton_width=18,
                   ).grid(row=1, columnspan=2)
    # Select Second residue
    selected_res1 = Tkinter.StringVar(master=group.interior())
    selected_res1.set("a-D-glucose")
    Pmw.OptionMenu(group.interior(),
                   labelpos='w',
                   label_text='Residue 2',
                   menubutton_textvariable=selected_res1,
                   items=residues_templates,
                   menubutton_width=18,
                   ).grid(row=2, columnspan=2)
    # Select Repetitions
    Tkinter.Label(group.interior(), text='Repetitions').grid(row=3, column=0)
    n_residues = Tkinter.StringVar(master=group.interior())
    n_residues.set(1)
    entry_n_residues = Tkinter.Entry(group.interior(),
                                     textvariable=n_residues,
                                     width=5)
    entry_n_residues.grid(row=3, column=1)
    # entry_n_residues.update()
    # Select position to insert ramification
    Tkinter.Label(group.interior(), text='Position').grid(row=4, column=0)
    first_res = Tkinter.StringVar(master=group.interior())
    first_res.set(0)
    entry_first_res = Tkinter.Entry(group.interior(),
                                    textvariable=first_res,
                                    width=5)
    entry_first_res.grid(row=4, column=1)
    # entry_first_res.update()
    # Total number of residues, hidden
    total_res = Tkinter.StringVar(master=group.interior())
    total_res.set(0)
    Tkinter.Label(group.interior(), text='New molecule').grid(row=5, column=0)
    name = Tkinter.StringVar(master=group.interior())
    name.set('carb')
    mol_name = Tkinter.Entry(group.interior(), textvariable=name, width=5)
    mol_name.grid(row=5, column=1)
    mol_name.update()

    Tkinter.Button(p1,
                   text="  add  ",
                   command=lambda: add(mol_name.get())
                   ).pack(side=Tkinter.LEFT)
    Tkinter.Button(p1,
                   text="create",
                   command=lambda: create(mol_name.get())
                   ).pack(side=Tkinter.LEFT)
    Tkinter.Button(p1, text="reset", command=reset).pack(side=Tkinter.RIGHT)
    ############################ Visualization TAB ###########################
    group = Pmw.Group(p2, tag_text='options')
    group.pack(fill='both', expand=1, padx=5, pady=5)

    colors_list = ['auto', 'green', 'red', 'blue', 'yellow', 'cyan', 'magenta',
                   'orange', 'grey', 'black', 'white']
    Tkinter.Label(group.interior(), text='Colors').grid(row=0, column=0)
    colors = Tkinter.StringVar(master=group.interior())
    colors.set(colors_list[0])
    Pmw.OptionMenu(group.interior(),
                   labelpos='w',
                   menubutton_textvariable=colors,
                   items=colors_list,
                   menubutton_width=12,
                   ).grid(row=0, column=1)

    rep_list = ['cartoon', 'wire', 'beads']
    Tkinter.Label(group.interior(), text='Representations').grid(row=1, column=0)
    rep = Tkinter.StringVar(master=group.interior())
    rep.set(rep_list[0])
    Pmw.OptionMenu(group.interior(),
                   labelpos='w',
                   menubutton_textvariable=rep,
                   items=rep_list,
                   menubutton_width=12,
                   ).grid(row=1, column=1)

    stored.all_models = []
    cmd.iterate('(all)', 'stored.all_models.append((model))')

    models = ["all"]
    models.extend(np.unique(stored.all_models))
    Tkinter.Label(group.interior(), text='Models').grid(row=2, column=0)
    model = Tkinter.StringVar(master=group.interior())
    model.set("all")
    Pmw.OptionMenu(group.interior(),
                   labelpos='w',
                   menubutton_textvariable=model,
                   items=models,
                   menubutton_width=12,
                   ).grid(row=2, column=1)


    Tkinter.Label(group.interior(), text='Chains').grid(row=3, column=0)
    chains = Tkinter.StringVar(master=group.interior())
    stored.iterchains = []
    cmd.iterate('(all)', 'stored.iterchains.append((chain))')
    all_chains = "".join(set(stored.iterchains))
    entry_chain = Tkinter.Entry(group.interior(),
                                  textvariable=chains,
                                  width=12)
    entry_chain.insert(0, 'A')
    chains.set(all_chains)
    entry_chain.grid(row=3, column=1)
    entry_chain.configure(state='normal')
    entry_chain.update()

    Tkinter.Label(group.interior(), text='Draw Bonds?').grid(row=4, column=0)
    bond_val = Tkinter.BooleanVar(master=group.interior())
    bond_val.set(False)
    check_bond = Tkinter.Checkbutton(
        group.interior(),
        variable=bond_val)
    check_bond.grid(row=4, column=1)
    check_bond.configure(state='normal')
    check_bond.update()

    Tkinter.Button(p2,
                   text="visualize",
                   command=lambda: cartoonize(
                                              chains.get(),
                                              colors.get(),
                                              rep.get(),
                                              bond_val.get(),
                                              model.get())
                   ).pack()
    ############################Calculations TAB##############################
    group = Pmw.Group(p3, tag_text='options')
    group.pack(fill='both', expand=1, padx=5, pady=5)

    # List available computations
    Tkinter.Label(group.interior(), text=' Analysis ').grid(row=0, column=0)
    type_analysis = Tkinter.StringVar(master=group.interior())
    type_analysis.set('Rama scatter')
    Pmw.OptionMenu(group.interior(),
                   labelpos='w',
                   menubutton_textvariable=type_analysis,
                   items=['Rama scatter', 'Rama hex', ' Rg', 'Hydrogen_bonds'],
                   menubutton_width=12,
                   command=enable_disable_cutoff,
                   ).grid(row=0, column=1)

    Tkinter.Label(group.interior(), text='Selection').grid(row=1, column=0)
    sel0_value = Tkinter.StringVar(master=group.interior())
    sel0_value.set('all')
    entry_sel0_value = Tkinter.Entry(
        group.interior(),
        textvariable=sel0_value,
        width=12)
    entry_sel0_value.grid(row=1, column=1)
    entry_sel0_value.configure(state='normal')
    # entry_sel0_value.update()
    entry_sel0_value.bind(
        "<Return>", lambda e: to_state.set(
            cmd.count_states(
                sel0_value.get())))

    Tkinter.Label(group.interior(), text='From state').grid(row=2, column=0)
    from_state = Tkinter.IntVar(master=group.interior())
    from_state.set(1)
    entry_from_state = Tkinter.Entry(
        group.interior(),
        textvariable=from_state,
        width=12)
    entry_from_state.grid(row=2, column=1)

    Tkinter.Label(group.interior(), text='To state').grid(row=3, column=0)
    to_state = Tkinter.IntVar(master=group.interior())
    to_state.set(cmd.count_states(sel0_value.get()))
    entry_to_state = Tkinter.Entry(group.interior(),
                                   textvariable=to_state,
                                   width=12)
    entry_to_state.grid(row=3, column=1)

    Tkinter.Label(group.interior(), text='step').grid(row=4, column=0)
    step = Tkinter.IntVar(master=group.interior())
    step.set(1)
    entry_step = Tkinter.Entry(group.interior(), textvariable=step, width=12)
    entry_step.grid(row=4, column=1)

    Tkinter.Label(group.interior(), text='Cut off').grid(row=4, column=2)
    cut_off = Tkinter.IntVar(master=group.interior())
    cut_off.set(0)
    entry_cut_off = Tkinter.Entry(group.interior(),
                                  textvariable=cut_off,
                                  width=12)
    entry_cut_off.grid(row=4, column=3)
    entry_cut_off.configure(state='disabled')
    entry_cut_off.update()

    Tkinter.Label(group.interior(), text='Rg Visualization').grid(row=5, column=0)
    vis_rg_value = Tkinter.BooleanVar(master=group.interior())
    vis_rg_value.set(False)
    entry_rg = Tkinter.Checkbutton(
        group.interior(),
        variable=vis_rg_value,
        command=enable_disable_rg)
    entry_rg.grid(row=5, column=1)
    entry_rg.update()
    by_state_value = Tkinter.BooleanVar(master=group.interior())
    by_state_value.set(False)
    entry_by_state = Tkinter.Checkbutton(
        group.interior(),
        text='By state',
        variable=by_state_value)
    entry_by_state.grid(row=6, columnspan=1)
    entry_by_state.configure(state='disabled')
    entry_by_state.update()

    Tkinter.Button(p3,
                   text="run analysis",
                   command=lambda: analyse(type_analysis.get(),
                                           sel0_value.get(),
                                           int(entry_from_state.get()),
                                           int(entry_to_state.get()),
                                           int(entry_step.get()),
                                           vis_rg_value.get(),
                                           by_state_value.get(),
                                           int(entry_cut_off.get()))
                   ).pack()
    ############################## MCM TAB ####################################
    group = Pmw.Group(p4, tag_text='options')
    group.pack(fill='both', expand=1, padx=5, pady=5)
    # Selection
    Tkinter.Label(group.interior(), text='Molecule').grid(row=0, column=0)
    molecule = Tkinter.StringVar(master=group.interior())
    molecule.set('carb')
    entry_molecule = Tkinter.Entry(group.interior(), textvariable=molecule, width=5)
    entry_molecule.grid(row=0, column=1)
    # Number of iteration
    Tkinter.Label(group.interior(), text='Iterations').grid(row=1, column=0)
    iterations = Tkinter.IntVar(master=group.interior())
    iterations.set(10000)
    entry_iterations = Tkinter.Entry(
        group.interior(),
        textvariable=iterations,
        width=5)
    entry_iterations.grid(row=1, column=1)
    # use solvent or not
    Tkinter.Label(group.interior(), text='SASA').grid(row=2, column=0)
    sasa = Tkinter.BooleanVar(master=group.interior())
    sasa.set(True)
    entry_sasa = Tkinter.Checkbutton(group.interior(), variable=sasa)
    entry_sasa.grid(row=2, column=1)
    entry_sasa.configure(state='active')
    # randomize initial structure
    Tkinter.Label(group.interior(), text='Randomize').grid(row=3, column=0)
    randomize = Tkinter.BooleanVar(master=group.interior())
    randomize.set(True)
    entry_rnd = Tkinter.Checkbutton(group.interior(), variable=randomize)
    entry_rnd.grid(row=3, column=1)
    entry_rnd.configure(state='active')
    Tkinter.Button(p4,
                   text="run MCM",
                   command=lambda: mcm_run(entry_molecule.get(),
                                           int(entry_iterations.get()),
                                           bool(sasa.get()),
                                           bool(randomize.get()))
                  ).pack()
    ############################  About TAB   ################################
    group = Pmw.Group(p10, tag_text='About')
    group.pack(fill='both', expand=1, padx=5, pady=5)
    text = """This plugin was developed by Agustina
Arroyuelo and Osvaldo Martin as part of the
Warren L. DeLano Memorial PyMOL
Open-Source Fellowship.
For instructions on how to use this plugin, please
read:
http://www.pymolwiki.org/index.php/Azahar
"""
    interior_frame = Tkinter.Frame(group.interior())
    bar = Tkinter.Scrollbar(interior_frame)
    text_holder = Tkinter.Text(interior_frame,
                               yscrollcommand=bar.set,
                               foreground="#cecece",
                               background="#000000",
                               font="Times 12")
    bar.config(command=text_holder.yview)
    text_holder.insert(Tkinter.END, text)
    text_holder.pack(side=Tkinter.LEFT, expand="yes", fill="both")
    bar.pack(side=Tkinter.LEFT, expand="yes", fill="y")
    interior_frame.pack(expand="yes", fill="both")
