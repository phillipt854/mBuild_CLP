import mbuild as mb
import py3Dmol
import mbuild_CLP
from mbuild_CLP.mbuild_CLP import CLP_box

test_box = CLP_box(['POG'],[1,1,2])
test_box.visualize()
#test_box.write_lammps('test.lammps')
#test_box.create_lammps_input_script(sim_name='test.lammps')

