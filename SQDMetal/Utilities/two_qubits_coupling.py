import qiskit_metal as metal
from qiskit_metal import Dict, QComponent
from qiskit_metal.qlibrary.tlines.meandered import RouteMeander
from qiskit_metal.qlibrary.couplers.coupled_line_tee import CoupledLineTee
from SQDMetal.Utilities.QUtilities import QUtilities
from qiskit_metal.analyses import cpw_calculations
from SQDMetal.Utilities.Materials import Material
from SQDMetal.Comps.Junctions import JunctionSingleDolanPinStretch
from SQDMetal.Utilities.CpwParams import CpwParams
from SQDMetal.Comps.Bandages import BandageRectPin
from SQDMetal.Utilities.QubitDesigner import ResonatorHalfWave
from qiskit_metal.qlibrary.tlines.straight_path import RouteStraight
from qiskit_metal.qlibrary.terminations.launchpad_wb import LaunchpadWirebond
from SQDMetal.Comps.Qubits import TransmonTaperedInsets
from SQDMetal.Comps.Joints import Joint, JointExtend
from SQDMetal.Comps.Markers import MarkerSquare4, MarkerDicingCross
from SQDMetal.Comps.Wires import WirePins
from SQDMetal.Comps.Labels import LabelText
from SQDMetal.Comps import Junctions
import numpy as np
from qiskit_metal import designs, draw
from qiskit_metal import MetalGUI, open_docs
import matplotlib.pyplot as plt
from SQDMetal.Comps.Capacitors import Smooth_rectangle, Smooth_Capacitor_Semicircle, Smooth_synapse
from SQDMetal.Comps.Resonators import ResonatorMeanderPinPin
from SQDMetal.Utilities.MakeGDS import MakeGDS

## This file contains a function that gives a 2 qubit design given a coupling strength J12 and the capacitor choice. 

def two_qubits_coupling(
    J12_target,
    capacitor_type='rode',
    design=None,
    position=[0.0, 0.0],
    freq=[4.5e9, 4.5e9],
    print_results=False,
    rebuild_design=True,
    readout_pad=False,
    ensemble_orientation_deg=0.0,
    flux_line_left=False,
    flux_line_right=False,
):

    ''' This function gives a 2 qubit design given a coupling strength J12 and the capacitor choice. The coupling strength is determined by the capacitor geometry and the frequencies of the qubits.
        This function is limited to the frequency 4.5 GHz, 5 GHz and coupling strength J12 is confined between different values for different capacitor types and frequencies difference.

        rode capacitor: same frequency 4.5 GHz: J12 can be tuned between 
                        diff frequency 4.5 GHz and 5 GHz: J12 can be tuned between  5 to 30 MHz
        
        semi-circle capacitor: same frequency 4.5 GHz: J12 can be tuned between 
                                diff frequency 4.5 GHz and 5 GHz: J12 can be tuned between

        Synapse capacitor: same frequency 4.5 GHz: J12 can be tuned between 
                            diff frequency 4.5 GHz and 5 GHz: J12 can be tuned between

        The function takes in the target coupling strength J12_target, the capacitor type, the design object (if None, a new design will be created), the position of the qubits, the frequencies of the qubits, 
        and whether to print the results of the capacitor parameters.

        *J12_target in MHz
        *capacitor_type: 'rode', 'semicircle' or 'synapse'
        *design: Qiskit Metal design object
        *position: List of x and y coordinates for the qubit positions, the position will be at the point between the two qubits, default is [0.0, 0.0] mm
        *freq: List of frequencies for the two qubits
        *print_results: Boolean indicating whether to print the capacitor parameters


        The function returns the design object, which can be used to visualize and export the design.'''
    #################### Target feasibility check ###################################################################################
    
    diff_freq = False
    if freq[0] != freq[1]:
        diff_freq = True

    if capacitor_type == 'rode' and not diff_freq:
        if J12_target < 5 or J12_target > 50.5:
            raise ValueError("For rode capacitor with same frequencies, J12 must be between 5 and 50.5 MHz.")
    elif capacitor_type == 'rode' and diff_freq:
        if J12_target < 5 or J12_target > 35.5:
            raise ValueError("For rode capacitor with different frequencies, J12 must be between 5 and 35.5 MHz.")
    elif capacitor_type == 'semicircle' and not diff_freq:
        if J12_target < 8 or J12_target > 30.5:
            raise ValueError("For semi-circle capacitor with same frequencies, J12 must be between 8 and 30.5 MHz.")
    elif capacitor_type == 'semicircle' and diff_freq:
        if J12_target < 8 or J12_target > 20.5:
            raise ValueError("For semi-circle capacitor with different frequencies, J12 must be between 8 and 20.5 MHz.")
    elif capacitor_type == 'synapse' and not diff_freq:
        if J12_target < 0 or J12_target > 25.5:
            raise ValueError("For synapse capacitor with same frequencies, J12 must be between 0 and 25.5 MHz.")
    elif capacitor_type == 'synapse' and diff_freq:
        if J12_target < 0.5 or J12_target > 17.5:
            raise ValueError("For synapse capacitor with different frequencies, J12 must be between 0.5 and 17.5 MHz.")
    else:
        raise ValueError("Invalid capacitor type. Must be 'rode', 'semicircle' or 'synapse'.")


    #################### Design initialization ###################################################################################
    
    qubit1_exist = 'Qubit1' in design.components

    if design is None:
        design = designs.DesignPlanar({}, overwrite_enabled=True)
        # Set up chip dimensions 
        design.chips.main.size.size_x = '5.0mm'
        design.chips.main.size.size_y = '3.5mm'
        design.chips.main.size.size_z = '-280um'
        design.chips.main.size.center_x = '0.0mm'
        design.chips.main.size.center_y = '0mm'

       # Resonator and feedline gap width (W) and center conductor width (S) 
        design.variables['cpw_width'] = '10 um' #S
        design.variables['cpw_gap'] = '6 um' #W 

    # universal capacitor parameters
    gap_side = 25
    gap_front = 25
    gap_back = 25
    fillet_radius = 10
    chrgln_pin_x_offset = 0
    chrgln_pin_y_offset =  0
    taper_width_base = 200
    semi_radius = 61
    res_freq_GHZ = freq[1] + 1.3 *1e9
    #if qubit1_exist: 
    #    res_freq_GHZ = freq[3] + 1.3 *1e9
    resonator = ResonatorHalfWave(res_freq_GHZ)
    # Rotate the pair around the ensemble center (position) by ensemble_orientation_deg.
    theta = np.deg2rad(float(ensemble_orientation_deg))
    c, s = np.cos(theta), np.sin(theta)
    x_c, y_c = position[0], position[1]
    if flux_line_left:
        taper_width_base = 100
        junction_centered_left = 'left'
    else: 
        junction_centered_left = True
    if flux_line_right:
        taper_width_base = 100
        junction_centered_right = False
    else:
        junction_centered_right = True

    dx1, dy1 = -1.5, 0.0
    dx2, dy2 = 1.5, 0.0
    pos_x1 = x_c + dx1 * c - dy1 * s
    pos_y1 = y_c + dx1 * s + dy1 * c
    pos_x2 = x_c + dx2 * c - dy2 * s
    pos_y2 = y_c + dx2 * s + dy2 * c

    q_orient = str(int(round(ensemble_orientation_deg)) % 360)
    cap_orient_0 = int(round(ensemble_orientation_deg)) % 360
    cap_orient_180 = (180 + int(round(ensemble_orientation_deg))) % 360


    m = Material("siliconcryo")

    # CPW parameters
    cpw_width = 10e-6
    cpw = CpwParams(m.permittivity, 500e-6)
    cpw_gap = cpw.get_gap_from_width(cpw_width)

    l_fullwave, _, _ = cpw_calculations.guided_wavelength(freq=resonator.f0,line_width=cpw_width,
                line_gap=cpw_gap,substrate_thickness=cpw.dielectric_thickness,film_thickness=200e-9,)
    l_halfwave = l_fullwave / 2

    readout_pad_options = Dict(
        connection_pads=Dict(
            readout=Dict(
                cpw_extend="5um",
                pad_gap="90um",
                cpw_width=f"{cpw_width*1e6}um",
                cpw_gap=f"{cpw_gap*1e6}um",
                loc_W="0",
                loc_H="-1",
            )
        )
    )
     
    surname = '' 
    if qubit1_exist:
        surname = '_2nd'

    ################# Transmons defintion ###################################################################################

    TransmonTaperedInsets(design, f'Qubit1{surname}', options = Dict(pos_x = f'{pos_x1}mm', pos_y = f'{pos_y1}mm', orientation = q_orient,
                                                       pocket_lower_tighten = '-100.0um',pocket_height = '500um', pocket_width = '1000um', chrgln_pin_x_offset=f"{chrgln_pin_x_offset}um", chrgln_pin_y_offset=f"{chrgln_pin_y_offset}um",
                                                       taper_width_base = f"{taper_width_base}um", junction_centered = junction_centered_left,
                                                       **readout_pad_options))


    TransmonTaperedInsets(design, f'Qubit2{surname}', options = Dict(pos_x = f'{pos_x2}mm', pos_y = f'{pos_y2}mm', orientation = q_orient,
                                                      pocket_lower_tighten = '-100.0um',pocket_height = '500um', pocket_width = '1000um', chrgln_pin_x_offset=f"{chrgln_pin_x_offset}um", chrgln_pin_y_offset=f"{chrgln_pin_y_offset}um",
                                                      taper_width_base = f"{taper_width_base}um", junction_centered = junction_centered_right,
                                                       **readout_pad_options))


    ################ Capacitor parameters ###################################################################################
    if capacitor_type == 'rode':
        rect_width = 80
        if not diff_freq:
           rect_length = (J12_target + 21.207)/0.337
           
        elif diff_freq:
           rect_length = (-0.041 + np.sqrt(0.041**2 - 4*5.66e-4*(-0.251-J12_target)))/(2*5.66e-4)


        Smooth_rectangle(design, f'Capacitor1{surname}',options=Dict(pin_inputs=Dict(start_pin=Dict(component=f'Qubit1{surname}',pin='top_right_pin')), orientation= 180,
                                                              rect_width = f'{rect_width}um',rect_length = f'{rect_length}um', gap_side = f'{gap_side}um', gap_back = f'{gap_back}um', gap_front = f'{gap_front}um',
                                                              fillet_radius = f'{fillet_radius}um',fillet_resolution = 20, gap_fillet_radius = f'{fillet_radius}um'))

        Smooth_rectangle(design, f'Capacitor2{surname}',options=Dict(pin_inputs=Dict(start_pin=Dict(component=f'Qubit2{surname}',pin='top_left_pin')), orientation= 0,
                                                              rect_width = f'{rect_width}um',rect_length = f'{rect_length}um',gap_side = f'{gap_side}um', gap_back = f'{gap_back}um', gap_front = f'{gap_front}um',
                                                              fillet_radius = f'{fillet_radius}um',fillet_resolution = 20, gap_fillet_radius = f'{fillet_radius}um'))


    if capacitor_type == 'semicircle': 
        rect_length = 170
        rect_width = 190
        if not diff_freq:
            circle_offset = -(J12_target - 25.850)/0.480 # (0.389 - np.sqrt(0.389**2 - 4*(-4.516e-3)*(26.242-J12_target)))/(2*(-4.516e-3)) #-(J12_target - 25.850)/0.480
            if circle_offset < -10:
               circle_offset = -10

        elif diff_freq:
            circle_offset =  (0.232 - np.sqrt(0.232**2 - 4*1.137e-3*(15.482-J12_target)))/(2*1.137e-3)
            if circle_offset < -10:
               circle_offset = -10

        Smooth_Capacitor_Semicircle(design, f'Capacitor1{surname}',options=Dict(pin_inputs=Dict(start_pin=Dict(component=f'Qubit1{surname}',pin='top_right_pin')), orientation=cap_orient_180,
                                                              rect_width = f'{rect_width}um',rect_length = f'{rect_length}um', semi_radius = f'{semi_radius}um', circle_offset = f'{circle_offset}um',
                                                              fillet_radius = f'{fillet_radius}um',fillet_resolution = 20, gap_side = f'{gap_side}um', gap_back = f'{gap_back}um'))

        Smooth_Capacitor_Semicircle(design, f'Capacitor2{surname}',options=Dict(pin_inputs=Dict(start_pin=Dict(component=f'Qubit2{surname}',pin='top_left_pin')), orientation=cap_orient_0,
                                                              rect_width = f'{rect_width}um',rect_length = f'{rect_length}um', semi_radius = f'{semi_radius}um', circle_offset = f'{circle_offset}um',
                                                              fillet_radius = f'{fillet_radius}um',fillet_resolution = 20, gap_side = f'{gap_side}um', gap_back = f'{gap_back}um'))

    if capacitor_type == 'synapse':

        bulb_radius = 100
        rect_length = 40
        rect_width = 10
        circle_offset = -40
        big_fillet_radius = 20
        gap_side = 10
        gap_front = 10

        if not diff_freq:
           bulb_scale_y = 1.1

        elif diff_freq:
           bulb_scale_y = 1.2

        Smooth_synapse(design, f'Capacitor1{surname}', options = Dict(pin_inputs=Dict(start_pin=Dict(component=f'Qubit1{surname}',pin='top_right_pin')),
                               pos_x = '0.0mm', pos_y = '0.0mm', orientation = str(cap_orient_0), fillet_resolution = 16, bulb_radius = f"{bulb_radius}um",
                                                   rect_length = f"{rect_length}um", rect_width = f"{rect_width}um", bulb_scale_y = f"{bulb_scale_y}", circle_offset = f"{circle_offset}um",
                                                   gap_side = f"{gap_side}um", gap_front = f"{gap_front}um", fillet_radius = f"{fillet_radius}um", semi_radius = f"{semi_radius}um", big_fillet_radius = f"{big_fillet_radius}um"))

        Smooth_synapse(design, f'Capacitor2{surname}', options = Dict(pin_inputs=Dict(start_pin=Dict(component=f'Qubit2{surname}',pin='top_left_pin')),
                               pos_x = '0.0mm', pos_y = '0.0mm', orientation = str(cap_orient_180), fillet_resolution = 16, bulb_radius = f"{bulb_radius}um",
                                                   rect_length = f"{rect_length}um", rect_width = f"{rect_width}um", bulb_scale_y = f"{bulb_scale_y}", circle_offset = f"{circle_offset}um",
                                                   gap_side = f"{gap_side}um", gap_front = f"{gap_front}um", fillet_radius = f"{fillet_radius}um", semi_radius = f"{semi_radius}um", big_fillet_radius = f"{big_fillet_radius}um"))


    if capacitor_type not in ['rode', 'semicircle', 'synapse']:
        raise ValueError("Invalid capacitor type. Must be 'rode', 'semicircle' or 'synapse'.")


    ############### Resonator definition ###################################################################################

    ResonatorMeanderPinPin(design, f'Resonator1{surname}', options = Dict(pin_inputs=Dict(start_pin=Dict(component=f'Capacitor1{surname}', pin='a'),end_pin=Dict(component=f'Capacitor2{surname}', pin='a')),
                                            total_length = f"{l_halfwave*1e3}mm", start_left = True, const_radius = '10um', const_width_max = '200um', fillet_padding = '50 um',
                                            trace_width = f"{cpw_width*1e6}um", gap_width = f"{cpw_gap*1e6}um"))
    

    if print_results:
        if capacitor_type == 'rode':
            print("Rode Capacitor parameters:")
            print("===========================")
            print(f"Rode capacitor parameters for J12 = {J12_target} MHz:")
            print(f"Rectangular length: {rect_length} um")
            print(f"Rectangular width: {rect_width} um")
        elif capacitor_type == 'semicircle':
            print("Semi-circle Capacitor parameters:")
            print("===========================")
            print(f"Semi-circle capacitor parameters for J12 = {J12_target} MHz:")
            print(f"Rectangular length: {rect_length} um")
            print(f"Rectangular width: {rect_width} um")
            print(f"Circle offset: {circle_offset} um")
        elif capacitor_type == 'synapse':
            print("Synapse Capacitor parameters:")
            print("===========================")
            print(f"Synapse capacitor parameters for J12 = {J12_target} MHz:")
            print(f"Bulb radius: {bulb_radius} um")
            print(f"Rectangular length: {rect_length} um")
            print(f"Rectangular width: {rect_width} um")
            print(f"Circle offset: {circle_offset} um")
            print(f"Bulb scale Y: {bulb_scale_y}")


    if rebuild_design:
        design.rebuild()

    return design