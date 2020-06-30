import mbuild as mb
import numpy as np
import pytest
from mbuild_CLP import *


def test_sc_sequence():
    chain = CLP('POG')
    name_list = []
    for aa in chain.children:
        name_list.append(aa.name)
    assert name_list == ['AAP', 'AAO', 'AAG'] and chain.sequence == 'POG'


def test_sc_particles():
    chain = CLP('POG')
    name_list = []
    for aa in chain.children:
        for p in aa.children:
            name_list.append(p.name)
    assert name_list == ['_bbp', '_hbp', '_bbo',
                         '_bbg', '_hbg']


def test_HBBB_bond():
    # Check there is only 1 bond that has the only 2 particles in the system 
    # connected
    chain = CLP('P')
    particle_set = set([p for p in chain.particles()])
    bonds = [b for b in chain.bonds()]
    firstBond = set(bonds[0])
    assert len(bonds) == 1 and firstBond == particle_set


def test_triplehelix():
    ds = CLP_helix('POG')
    chain_list = []
    for chain in ds.children:
        chain_list.append(chain)
    assert len(chain_list) == 3

