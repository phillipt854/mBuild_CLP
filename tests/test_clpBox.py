import mbuild as mb
import numpy as np
import pytest
from mbuild_CLP import *

def test_boxChains():
    testBox = CLP_box(['POG'], dim=[1, 1, 1])
    chains = [c.sequence for c in testBox.children]
    testBox.save('clpbox.gsd')
    assert chains == ['POG']
