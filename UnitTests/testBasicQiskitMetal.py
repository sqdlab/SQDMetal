import os
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"
os.environ["PMIX_MCA_gds"]="hash"
# Import useful packages
import qiskit_metal as metal
from qiskit_metal import designs, draw
from qiskit_metal import MetalGUI, Dict, open_docs
from qiskit_metal.toolbox_metal import math_and_overrides
from qiskit_metal.qlibrary.core import QComponent
from collections import OrderedDict
# To create plots after geting solution data.
import numpy as np
# Packages for the simple design
from SQDMetal.Comps.Xmon import Xmon
from SQDMetal.Comps.Junctions import JunctionDolanPinStretch

# from SQDMetal.Utilities.QUtilities import QUtilities
# from qiskit_metal.toolbox_python.attr_dict import Dict
# #
# from SQDMetal.Comps.Capacitors import CapacitorInterdigitalPinStretch
# from qiskit_metal.qlibrary.terminations.launchpad_wb import LaunchpadWirebond

import shutil

import unittest

class TestBasic(unittest.TestCase):
    ERR_TOL = 5e-13
    
    def test_XmonDesign(self):
        # Set up chip design as planar, multiplanar also available
        design = designs.DesignPlanar({}, overwrite_enabled=True)
        # Set up chip dimensions 
        design.chips.main.size.size_x = '500um'
        design.chips.main.size.size_y = '500um'
        design.chips.main.size.size_z = '500um'
        design.chips.main.size.center_x = '0mm'
        design.chips.main.size.center_y = '0mm'
        # Create the x-mon
        xmon = Xmon(design, 'x-mon', options=Dict(pos_x=0, pos_y=0,
                                            vBar_width='24um', hBar_width='24um', vBar_gap=f'{16}um', hBar_gap=f'{16}um',
                                            cross_width=f'{60*2+24}um', cross_height=f'{60*2+24}um',
                                            gap_up='24um', gap_left='24um', gap_right='24um', gap_down='24um'))

        # Create the Josephson junction
        JunctionDolanPinStretch(design, 'junction', options=Dict(pin_inputs=Dict(start_pin=Dict(component=f'x-mon',pin='right')),
                                                                dist_extend='24um',
                                                                layer=2,
                                                                finger_width='0.4um', t_pad_size='0.385um',
                                                                squid_width='5.4um', prong_width='0.9um'));

    # def test_Gallery(self):
    #     design = designs.DesignPlanar({}, overwrite_enabled=True)
    #     design.delete_all_components()
    #     LaunchpadWirebond(design, 'LP1', options = dict(chip='main', orientation='45', lead_length='20um', pad_height='20um', 
    #                         pad_width='40um', pad_gap='20um'))
    #     ##########################
    #     CapacitorInterdigitalPinStretch(design, 'leCap', options=Dict(pin_inputs=Dict(start_pin=Dict(component=f'LP1',pin='tie')), dist_extend='100um',
    #                                                         cpw_width='15um', len_flat='10um', fing_len='10um', fing_len_gap='1um', fing_wid='2um',
    #                                                         len_diag='7um', fing_wid_gap='1um', N_total=9, larger_first=True))
    #     ##########################
    #     QUtilities.plot_highlight_component('leCap', design)


if __name__ == '__main__':
    TestBasic().test_XmonDesign()
    unittest.main()
