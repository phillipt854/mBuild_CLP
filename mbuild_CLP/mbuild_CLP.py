#!/usr/bin/env python
# coding: utf-8

# In[7]:


from __future__ import print_function
import mbuild as mb
import builtins as __builtin__
import numpy as np
from foyer import Forcefield
#All the old import from compound
import collections
from collections import OrderedDict, defaultdict
from copy import deepcopy
import itertools
import os
import sys
import tempfile
from warnings import warn
##
from mbuild.utils.io import run_from_ipython, import_
r0BB = 0.5 # BB-BB bond length of 0.5 sigma for CLP model 
r0BBHB = 0.37


class BB(mb.Compound):
    def __init__(self):
        super(BB, self).__init__(pos=[0,0,0],name = 'bb')#Initizlize an instance of abstract class
        
        bead = mb.Particle(pos=[0.0, 0.0, 0.0], name='BB')
        self.add(bead)
        
        port = mb.Port(anchor=list(self.particles_by_name('BB'))[0], orientation=[0, 1, 0], separation=r0BB/2)
        self.add(port, 'up')
        port = mb.Port(anchor=list(self.particles_by_name('BB'))[0], orientation=[1, 0, 0], separation=r0BBHB/2)
        self.add(port,'toHB')
        port = mb.Port(anchor=list(self.particles_by_name('BB'))[0], orientation=[0, -1, 0], separation=r0BB/2)
        self.add(port, 'down')
        

class BBP(BB):
    def __init__(self):
        super(BBP,self).__init__()
        for par in self.particles():
            par.name = '_BBP'
            #print(par.name)
        
class BBO(BB):
    def __init__(self):
        super(BBO,self).__init__()
        for par in self.particles():
            par.name = '_BBO'

class BBG(BB):
    def __init__(self):
        super(BBG,self).__init__()
        for par in self.particles():
            par.name = '_BBG'
class BBK(BB):
    def __init__(self):
        super(BBK,self).__init__()
        for par in self.particles():
            par.name = '_BBK'

class BBD(BB):
    def __init__(self):
        super(BBD,self).__init__()
        for par in self.particles():
            par.name = '_BBD'        

# Define hydrogen bonding (HB) class        
class HB(mb.Compound):
    def __init__(self):
        super(HB, self).__init__(pos=[0.0, 0.0, 0.0], name='hb')
        
        bead = mb.Particle(pos=[0.0, 0.0, 0.0], name='HB')
        self.add(bead)
        
        port = mb.Port(anchor=list(self.particles_by_name('HB'))[0], orientation=[-1, 0, 0], separation=r0BBHB/2)
        self.add(port, 'toBB')      

class HBP(HB):
    def __init__(self):
        super(HBP,self).__init__()
        for par in self.particles():
            par.name = '_HBP'
        
class HBG(HB):
    def __init__(self):
        super(HBG,self).__init__()
        for par in self.particles():
            par.name = '_HBG'
        
# Define amino acid (AA) class - combine backbone and hydrogen bonding beads
class AAP(mb.Compound):
    def __init__(self):
        super(AAP,self).__init__()
        bb = BBP()
        hb = HBP()
        self.add((bb,hb))
        #Move
        mb.force_overlap(move_this=hb, from_positions=hb['toBB'],to_positions=bb['toHB'])

class AAO(mb.Compound):
    def __init__(self):
        super(AAO,self).__init__()
        bb = BBO()
        self.add(bb)

class AAG(mb.Compound):
    def __init__(self):
        super(AAG,self).__init__()
        bb = BBG()
        hb = HBG()
        self.add((bb,hb))
        #Move
        mb.force_overlap(move_this=hb, from_positions=hb['toBB'],to_positions=bb['toHB'])
        
class AAK(mb.Compound):
    def __init__(self):
        super(AAK,self).__init__()
        bb = BBK()
        self.add(bb)

class AAD(mb.Compound):
    def __init__(self):
        super(AAD,self).__init__()
        bb = BBD()
        hb = HBP() # P and D share the same HB bead type 
        self.add((bb,hb))
        #Move
        mb.force_overlap(move_this=hb, from_positions=hb['toBB'],to_positions=bb['toHB'])


def get_AA(type='P'):
    if type == 'P':
        return AAP()
    if type == 'O':
        return AAO()
    if type == 'G':
        return AAG()
    if type == 'K':
        return AAK()
    if type == 'D':
        return AAD
    return None

    




class CLP(mb.Compound):
    sequence = None
    def __init__(self,seq=None):
        self.sequence = seq
        super(CLP,self).__init__()
        seq_len = len(seq)
        if seq_len == 0:
            pass
        last_AA = get_AA(seq[0])
        self.add(last_AA)
        
        for letter in seq[1:]:
            new_AA = get_AA(letter)
            self.add(new_AA)
            #Here we need to pay some attention to the sequencing
            #Define an convention of the down port of "earlier beads" connects to the up port of "later beads"
            #Have to move new_AA to last_AA because the positions are screwed up otherwise. This is a bug and should be expected to be fixed soon
            #Have to always refer to the last available port since first last_AA has two available ports
            mb.force_overlap(move_this = new_AA, from_positions=(new_AA.all_ports())[0],to_positions = (last_AA.all_ports())[-1])
            last_AA = new_AA
    def print_seq(self):
        __builtin__.print(self.sequence)
            #Questions:
            #move last AA or last_AA[0]?
            #better way to reference children?
            #Do I need to translate? How to translate?
    def CLP_visualize_py3dmol(self, show_ports=False):
        """Visualize the Compound using py3Dmol.
        Allows for visualization of a Compound within a Jupyter Notebook.
        Parameters
        ----------
        show_ports : bool, optional, default=False
            Visualize Ports in addition to Particles
        color_scheme : dict, optional
            Specify coloring for non-elemental particles
            keys are strings of the particle names
            values are strings of the colors
            i.e. {'_CGBEAD': 'blue'}
        Returns
        ------
        view : py3Dmol.view
        """
        py3Dmol = import_('py3Dmol')
        remove_digits = lambda x: ''.join(i for i in x if not i.isdigit()
                                              or i == '_')

        #modified_color_scheme = {}
        

        for particle in self.particles():
            #particle.name = remove_digits(particle.name).upper()
            if not particle.name:
                particle.name = 'UNK'
                
        view = py3Dmol.view()
        rad = {'_BBP':.5,'_BBK':.5,'_BBG':.5,'_BBO':.5,'_BBD':.5,'_HBP':0.22,'_HBG':0.22}
        col = {'_BBP':'#0000FF','_BBO':'#FF8000','_BBG':'#00FF00','_BBK':'#000000','_BBD':'#FF0000','_HBP':'#FFFF00'
                 ,'_HBG':'#FFFF00'}
        
        for p in self.particles(include_ports=False):
            view.addSphere({
                'center': {'x':p.pos[0], 'y':p.pos[1], 'z':p.pos[2]},
                'radius' :rad[p.name],
                'color': col[p.name],
                'alpha': 0.9})
        view.zoomTo()
        view.show()

        return view


# In[38]:

class CLP_helix(mb.Compound):
    def __init__(self,sequences=[]):
        super(CLP_helix,self).__init__()
        seq_len = len(sequences)
        if seq_len == 0:
            pass
        # Create 3 CLP strands
        new_CLP = []
        new_CLP.append(CLP(sequences))
        new_CLP.append(CLP(sequences))
        new_CLP.append(CLP(sequences))
        
        # Create triple helix
        spacingPos = 1.04/(np.sqrt(3))
        bBLength = 0.5
        positions = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])
        spin_angle = [-np.pi/6.,-2*np.pi/3.,2.*np.pi/3.]
        for l in range(3):
            th = l*2.0*np.pi/3.0
            offset = l*bBLength
            positions[l] = np.add(positions[l],[-spacingPos,offset,0.])
            positions[l] = np.dot(positions[l],[[np.cos(th),0.,np.sin(th)],[0.,1.0,0.],[-np.sin(th),0.,np.cos(th)]])
            new_CLP[l].translate_to(positions[l])
            new_CLP[l].spin(spin_angle[l],[0.,1.,0,])
            self.add(new_CLP[l])
        
    def visualize(self, show_ports=False):
        """Visualize the Compound using py3Dmol.
        Allows for visualization of a Compound within a Jupyter Notebook.
        Parameters
        ----------
        show_ports : bool, optional, default=False
            Visualize Ports in addition to Particles
        color_scheme : dict, optional
            Specify coloring for non-elemental particles
            keys are strings of the particle names
            values are strings of the colors
            i.e. {'_CGBEAD': 'blue'}
        Returns
        ------
        view : py3Dmol.view
        """
        py3Dmol = import_('py3Dmol')
        remove_digits = lambda x: ''.join(i for i in x if not i.isdigit()
                                              or i == '_')

        #modified_color_scheme = {}
        

        for particle in self.particles():
            #particle.name = remove_digits(particle.name).upper()
            if not particle.name:
                particle.name = 'UNK'
                
        view = py3Dmol.view()
        rad = {'_BBP':.5,'_BBK':.5,'_BBG':.5,'_BBO':.5,'_BBD':.5,'_HBP':0.22,'_HBG':0.22}
        col = {'_BBP':'#0000FF','_BBO':'#FF8000','_BBG':'#00FF00','_BBK':'#000000','_BBD':'#FF0000','_HBP':'#FFFF00'
                 ,'_HBG':'#FFFF00'}
        
        for p in self.particles(include_ports=False):
            view.addSphere({
                'center': {'x':p.pos[0], 'y':p.pos[1], 'z':p.pos[2]},
                'radius' :rad[p.name],
                'color': col[p.name],
                'alpha': 0.9})
        view.zoomTo()
        view.show()

        return view
        
         

class CLP_box(mb.Compound):
    def __init__(self,sequences = [],dim = [0,0,0]):
        super(CLP_box,self).__init__()
        if len(sequences) >0 and dim[0]*dim[1]*dim[2] != len(sequences):
            dim = [len(sequences),1,1]
        seq_num = 0
        for i in range(dim[0]):
            for j in range(dim[1]):
                for k in range(dim[2]): 
                    # Create 3 CLP strands
                    new_CLP_helix = CLP_helix(sequences[seq_num])
                    new_CLP_helix.translate_to([i*3.,j*3.,k*3.])
                    self.add(new_CLP_helix)
                   
                        
             
                    seq_num += 1

    def visualize(self):
        py3Dmol = import_('py3Dmol')
        view = py3Dmol.view()
        rad = {'_BBP':.5,'_BBK':.5,'_BBG':.5,'_BBO':.5,'_BBD':.5,'_HBP':0.22,'_HBG':0.22}
        col = {'_BBP':'#0000FF','_BBO':'#FF8000','_BBG':'#00FF00','_BBK':'#000000','_BBD':'#FF0000','_HBP':'#FFFF00'
                 ,'_HBG':'#FFFF00'}
        
        remove_digits = lambda x: ''.join(i for i in x if not i.isdigit()
                                              or i == '_')

        #modified_color_scheme = {}
        
        for chain in self.children:
            for particle in chain.particles():
                #particle.name = remove_digits(particle.name).upper()
                if not particle.name:
                    particle.name = 'UNK'
                
        
        
                for p in chain.particles(include_ports=False):
                    view.addSphere({
                        'center': {'x':p.pos[0], 'y':p.pos[1], 'z':p.pos[2]},
                        'radius' :rad[p.name],
                        'color': col[p.name],
                        'alpha': 0.9})
        view.zoomTo()
        view.show()
        return view

    def write_lammps(self,filename):
        
        cgff = Forcefield(forcefield_files='mbuild_CLP/CLP.xml')
        #apply it to the mbuild structure
        test_box_typed =cgff.apply(self,assert_dihedral_params=False)
        #Output a LAMMPS data file
        mb.formats.lammpsdata.write_lammpsdata(test_box_typed,filename,atom_style='full')

    
    def create_lammps_input_script(self,T=4.6,script_num=1,sim_name='test.lammps',time_steps=2e6,atom_style='full',seed_1='42',seed_2='42'):
        re = import_('re')
	    # Calculate dielectric constant (using equation from Ahmad's soft matter paper)
        T_real = T*(0.1/8.314)*4184
        eps_r = 249.4 - 0.788*T_real + 0.00072*(T_real**2)
        eps_r = np.round(eps_r,decimals=5)

	# Specify the name of the sample input script 
        sample_name = 'mbuild_CLP/sample.in'

	# Loop over the total number of scripts and create files 
        for i in range(script_num):
        # Specify input file name 
            filename = 'test_sample' + str(i+1) + '.in'

		# Treat the first input script separately
            if (i == 0):
			# Define substitutions
                subs = [
        				('T_ref',str(T)),
        				('eps_r',str(eps_r)),
        				('sample.dcd',sim_name + '.' + str(i+1) + '.dcd'),
        				('sample.xyz',sim_name + '.' + str(i+1) + '.xyz'),
        				('time_steps_a', str(int(time_steps + 1))),
        				('time_steps_b', str(int(time_steps))),
        				('info.dat', 'info' + str(i+1) + '.dat'),
        				('atom_style full', 'atom_style ' + atom_style),
        				('seed_1',seed_1),
        				('seed_2',seed_2)
                                        ]
            else:
           # Replace relevant fields (Don't forget to replace read_data and data.lammps)
                subs = [
                    ('T_ref',str(T)),
                    ('eps_r',str(eps_r)),
                    ('read_data '+ sim_name,'read_restart restart.* remap'),
                    ('sample.dcd',sim_name + '.' + str(i+1) + '.dcd'),
                    ('sample.xyz',sim_name + '.' + str(i+1) + '.xyz'),
                    ('time_steps_a', str(int(time_steps + 1))),
                    ('time_steps_b', str(int(time_steps))),
                    ('info.dat', 'info' + str(i+1) + '.dat'),
                    ('atom_style full', 'atom_style ' + atom_style),
                    ('seed_1',seed_1),
                    ('seed_2',seed_2)
                    ]
    
		# Read sample input script 
            with open(sample_name) as fp:
                text = fp.read()
    
            lines = text.splitlines()
            for line_no,line in enumerate(lines):
                for pattern,replace in subs:
                # Do replacements
                    lines[line_no] = re.sub(pattern,replace,lines[line_no])

                # Don't re-assign velocities for i > 0 (Also, don't energy minimize)
                if (i > 0):
                    del lines[105:109]
	
                with open(filename,'w') as fp: 
                    fp.write('\n'.join(lines))
		








