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
from SQDMetal.Utilities.two_qubits_coupling import two_qubits_coupling

class two_Qubits_chip(QComponent):
    """A class to generate a two qubits chip with a coupling capacitor between them. The coupling capacitor is designed using a Capacitor component. 
    The design of the coupling capacitor is based on the parameters defined in the two_qubits_coupling function, which calculates the parameters based on the desired coupling strength and qubit frequencies. The class also includes options for adding bandages to the qubits and adjusting the position of the coupling capacitor. The design can be easily modified by changing the parameters in the two_qubits_coupling function and the options in the class.
    """

    # default options
    default_options = Dict(
        lp_width_um=400,
        component_prefix="",
        sample_label="",
        offset_x_mm=0,  # origin offset
        offset_y_mm=0,  # origin offset
        chip_x_mm=10,
        chip_y_mm=10,
        lp_inset_mm=0.9,
        q_y_coord_mm=2.0,
        q_x_coord_mm=0.0,
        coupler_qubit_offset_mm=0.3,
        with_labels=False,
        with_junctions=True,
        with_markers=True,
        with_patches=True,
        with_sample_label=True,
        font="sans",
        font_size_sublabels="80um",
        font_size_samplelabel="120um",
        readout_res_fillet="40um",
        meander_spacing="150um",
        J_C_uA_um2=0.45, #critical current density
        L_JJ_nH_list=[15.8,12.9,14.5,14.5],  #list of junction inductances in nH
        J12_target_MHz=20, # target coupling strength in MHz
        cpw_width_um=10,
    )

    def make(self):
        """
        All components in the design are placed here
        """

        # CLASS OPTIONS
        p = self.p

        # coordinates
        offset_x_mm = p.offset_x_mm * 1e3
        offset_y_mm = p.offset_y_mm * 1e3
        chip_x_mm = p.chip_x_mm
        chip_y_mm = p.chip_y_mm
        q_y_coord_mm = p.q_y_coord_mm
        q_x_coord_mm = p.q_x_coord_mm

        # labels
        sample_label = p.sample_label
        if sample_label == "":
            sample_label = self.name
        prefix = p.component_prefix
        if prefix != "":
            prefix = f"{prefix}_"
        else:
            prefix = f"{self.name}_"

        #setup  2 coupled qubits coordinates

        qubit12_coords = [
            q_x_coord_mm, q_y_coord_mm]
        
        qubit34_coords = [
            -q_x_coord_mm, q_y_coord_mm]

        lp_abs_x = chip_x_mm / 2 - p.lp_inset_mm
        lp_coords_mm = [(-lp_abs_x, 0), (lp_abs_x, 0)]  # launchpads
        label_offset_y_mm = 0.8  # label offset (from qubit centre)
        coupler_qubit_offset_mm = p.coupler_qubit_offset_mm

        # labels
        font_size = p.font_size_sublabels

        # CPW parameters
        m = Material("siliconcryo")
        cpw_width = p.cpw_width_um * 1e-6
        cpw = CpwParams(m.permittivity, 500e-6)
        cpw_gap = cpw.get_gap_from_width(cpw_width)
        gap_ratio = cpw_gap / cpw_width  # gap = gap_ratio * width
        film_thickness_nm = 200

        # Junction inductances
        # define junction inductances in Henries for each qubit
        L_JJ_nH = np.array(p.L_JJ_nH_list) # single junction qubits

        # Junction critical current density
        J_C_uA_um2 = p.J_C_uA_um2

        # Calculate required junction areas corresponding to `L_JJ_H` array
        JJ_widths_um = Junctions.get_square_JJ_width(
            J_C_uA_um2=J_C_uA_um2, target_LJ_nH=L_JJ_nH, rounding=True
        )
        # print(f"Josephson junction widths (um): {JJ_widths_um}")

        #
        #
        ##OPTIONS
        #
        #
        if p.with_junctions == False:
            p.with_patches = False

        # component options
        #qubit_options = Dict(pad_gap="100um",connection_pads=Dict(readout=Dict(cpw_extend="5um", #   pad_gap="100um",
        #                                      pad_gap="90um", # changed in NQCT-02
        #                                      cpw_width=f"{cpw_width*1e6}um",cpw_gap=f"{cpw_gap*1e6}um")),)

        junction_options = Dict(
            dist_extend="60um",
            t_pad_size="0.385um",
            # squid_width="10um",
            prong_width="1um",
            prong_length="5um",
            stem_width="1um",
            # t_pad_length="3um",
            t_finger_length="3um",
            layer=2,
        )

        bandage_options = Dict(width="6um", height="6um", layer=3)

        coupler_options = Dict(
            coupling_space="5um",
            coupling_length="220um",
            down_length="200um",
            prime_width=f"{cpw_width*1e6}um",
            prime_gap=f"{cpw_gap*1e6}um",
            second_width=f"{cpw_width*1e6}um",
            second_gap=f"{cpw_gap*1e6}um",
            open_termination=True,
        )

        lp_width_um = 400
        lp_gap_um = gap_ratio * lp_width_um
        lp_options = dict(
            chip="main",
            lead_length="100um",
            pad_height="300um",
            pad_width=f"{lp_width_um}um",
            pad_gap=f"{lp_gap_um}um",
            taper_height="250um",
            trace_width=f"{cpw_width*1e6}um",
            trace_gap=f"{cpw_gap*1e6}um",
        )


        # Qubit coupling
        freq = [4.5e9, 5.0e9, 4.7e9, 4.7e9] # frequencies of qubits in Hz
        res_freq_GHZ = freq[1] + 1.3e9
        res2_freq_GHZ = freq[3] + 1.3e9
        J12_target = p.J12_target_MHz 
        # define resonator 1
        r = ResonatorHalfWave(res_freq_GHZ * 1e9)
        l_full, _, _ = cpw_calculations.guided_wavelength(
            freq=r.f0,
            line_width=cpw_width,
            line_gap=cpw_gap,
            substrate_thickness=cpw.dielectric_thickness,
            film_thickness=film_thickness_nm * 1e-9,
        )
        l_half = l_full / 2  # lambda/2 resonator length

        two_qubits_coupling(
            J12_target,
            capacitor_type='rode',
            design=self.design,
            position=qubit12_coords,
            print_results=False,
            freq=freq,
            rebuild_design=False,
            readout_pad=True,
        )
        
        # junctions
        if p.with_junctions:
            JunctionSingleDolanPinStretch(self.design,"Qubit1_JJ",
                options=Dict(pin_inputs=Dict(start_pin=Dict(component="Qubit1", pin="pin_island")),
                    t_finger_width=f"{JJ_widths_um[0]}um",finger_width=f"{JJ_widths_um[0]}um",
                    **junction_options,),
            )
            JunctionSingleDolanPinStretch(self.design,"Qubit2_JJ",
                options=Dict(pin_inputs=Dict(start_pin=Dict(component="Qubit2", pin="pin_island")),
                    t_finger_width=f"{JJ_widths_um[1]}um",finger_width=f"{JJ_widths_um[1]}um",
                    **junction_options,),
            )
        # patches
        if p.with_patches:
            BandageRectPin(
                self.design, "Qubit1_patch1",
                options=Dict( target_comp="Qubit1_JJ", target_pin="t", **bandage_options),
            )
            BandageRectPin(
                self.design,"Qubit1_patch2",
                options=Dict(target_comp="Qubit1_JJ", target_pin="f", **bandage_options ),
            )

        if p.with_patches:
            BandageRectPin(self.design,"Qubit2_patch1",
                options=Dict(target_comp="Qubit2_JJ", target_pin="t", **bandage_options),)
            BandageRectPin(self.design,"Qubit2_patch2",
                options=Dict(target_comp="Qubit2_JJ", target_pin="f", **bandage_options),
            )

        # coupler
        coupler = CoupledLineTee(
            self.design,"Coupler1",
            options=dict(
                pos_x=f"{offset_x_mm + q_x_coord_mm - coupler_qubit_offset_mm - 1.5}mm",
                pos_y=f"{offset_y_mm}mm",
                orientation="180",
                **coupler_options,
            ),
        )

        coupler2 = CoupledLineTee(
            self.design, "Coupler2",
            options=dict(
                pos_x=f"{offset_x_mm + q_x_coord_mm - coupler_qubit_offset_mm + 1.5}mm",
                pos_y=f"{offset_y_mm}mm",
                orientation="180",
                **coupler_options,
            ),
        )
        # calculate extra resonator length due to coupling regions
        res_coupler_l_mm = QUtilities.calc_length_of_path(
            self.design,
            component_name="Coupler1",
            trace_name="prime_cpw_sub",
        )

        q_coupler_l_mm = 0 #qubit.get_resonator_length_mm()
        extra_length_mm = res_coupler_l_mm + q_coupler_l_mm
        # readout resonator 1
        readout = RouteMeander(
            self.design,"Resonator1_Read_Q",
            options=Dict(
                hfss_wire_bonds=True,
                pin_inputs=Dict( start_pin=Dict(component="Qubit1", pin="readout"),
                    end_pin=Dict(component="Coupler1", pin="second_end"),),
                lead=Dict(
                    start_straight="50um",
                    end_straight="50um",
                ),
                meander=Dict(
                    asymmetry="0um",
                    spacing=p.meander_spacing,
                ),
                trace_width=f"{cpw_width*1e6}um",
                trace_gap=f"{cpw_gap*1e6}um",
                fillet=p.readout_res_fillet,
                total_length=f"{(l_half * 1e3) - extra_length_mm:.2f}mm",
            ),
        )
         # readout resonator2
        readout2 = RouteMeander(
            self.design,
            "Resonator2_Read_Q",
            options=Dict(
                hfss_wire_bonds=True,
                pin_inputs=Dict(
                    start_pin=Dict(component="Qubit2", pin="readout"),
                    end_pin=Dict(
                        component="Coupler2", pin="second_end"
                    ),
                ),
                lead=Dict(
                    start_straight="50um",
                    end_straight="50um",
                ),
                meander=Dict(
                    asymmetry="0um",
                    spacing=p.meander_spacing,
                ),
                trace_width=f"{cpw_width*1e6}um",
                trace_gap=f"{cpw_gap*1e6}um",
                fillet=p.readout_res_fillet,
                total_length=f"{(l_half * 1e3) - extra_length_mm:.2f}mm",
            ),
        )

         ######################## 2nd two qubits ensemble ########################

        surname = '_2nd' 
         # define resonator 2
        r2 = ResonatorHalfWave(res2_freq_GHZ * 1e9)
        l2_full, _, _ = cpw_calculations.guided_wavelength(
            freq=r2.f0,
            line_width=cpw_width,
            line_gap=cpw_gap,
            substrate_thickness=cpw.dielectric_thickness,
            film_thickness=film_thickness_nm * 1e-9,
            )
        l2_half = l2_full / 2  # lambda/2 resonator length
         
        two_qubits_coupling(
            J12_target,
            capacitor_type='rode',
            design=self.design,
            position=qubit34_coords,
            print_results=False,
            freq=freq,
            rebuild_design=False,
            readout_pad=True,
        )
        
        '''# junctions
        if p.with_junctions:
            JunctionSingleDolanPinStretch(self.design,f"Qubit1_JJ{surname}",
                options=Dict(pin_inputs=Dict(start_pin=Dict(component=f"Qubit1{surname}", pin="pin_island")),
                    t_finger_width=f"{JJ_widths_um[2]}um",finger_width=f"{JJ_widths_um[2]}um",
                    **junction_options,),
            )
            JunctionSingleDolanPinStretch(self.design,f"Qubit2_JJ{surname}",
                options=Dict(pin_inputs=Dict(start_pin=Dict(component=f"Qubit2{surname}", pin="pin_island")),
                    t_finger_width=f"{JJ_widths_um[3]}um",finger_width=f"{JJ_widths_um[3]}um",
                    **junction_options,),
            )
        # patches
        if p.with_patches:
            BandageRectPin(
                self.design, f"Qubit1_patch1{surname}",
                options=Dict( target_comp=f"Qubit1_JJ{surname}", target_pin="t", **bandage_options),
            )
            BandageRectPin(
                self.design,f"Qubit1_patch2{surname}",
                options=Dict(target_comp=f"Qubit1_JJ{surname}", target_pin="f", **bandage_options ),
            )

        if p.with_patches:
            BandageRectPin(self.design,f"Qubit2_patch2{surname}",
                options=Dict(target_comp=f"Qubit2_JJ{surname}", target_pin="f", **bandage_options),)
            BandageRectPin(self.design,f"Qubit2_patch1{surname}",
                options=Dict(target_comp=f"Qubit2_JJ{surname}", target_pin="t", **bandage_options),
            )

        # coupler
        coupler3 = CoupledLineTee(
            self.design,"Coupler3",
            options=dict(
                pos_x=f"{offset_x_mm + q_x_coord_mm - coupler_qubit_offset_mm - 1.5}mm",
                pos_y=f"{offset_y_mm}mm",
                orientation="180",
                **coupler_options,
            ),
        )

        coupler4 = CoupledLineTee(
            self.design, "Coupler4",
            options=dict(
                pos_x=f"{offset_x_mm + q_x_coord_mm - coupler_qubit_offset_mm + 1.5}mm",
                pos_y=f"{offset_y_mm}mm",
                orientation="180",
                **coupler_options,
            ),
        )
        # calculate extra resonator length due to coupling regions
        res_coupler_l_mm = QUtilities.calc_length_of_path(
            self.design,
            component_name="Coupler3",
            trace_name="prime_cpw_sub",
        )

        q_coupler_l_mm = 0 #qubit.get_resonator_length_mm()
        extra_length_mm = res_coupler_l_mm + q_coupler_l_mm
        # readout resonator 3
        readout = RouteMeander(
            self.design,"Resonator3_Read_Q",
            options=Dict(
                hfss_wire_bonds=True,
                pin_inputs=Dict( start_pin=Dict(component=f"Qubit1{surname}", pin="readout"),
                    end_pin=Dict(component="Coupler3", pin="second_end"),),
                lead=Dict(
                    start_straight="50um",
                    end_straight="50um",
                ),
                meander=Dict(
                    asymmetry="0um",
                    spacing=p.meander_spacing,
                ),
                trace_width=f"{cpw_width*1e6}um",
                trace_gap=f"{cpw_gap*1e6}um",
                fillet=p.readout_res_fillet,
                total_length=f"{(l2_half * 1e3) - extra_length_mm:.2f}mm",
            ),
        )
         # readout resonator4
        readout4 = RouteMeander(
            self.design,
            "Resonator4_Read_Q",
            options=Dict(
                hfss_wire_bonds=True,
                pin_inputs=Dict(
                    start_pin=Dict(component=f"Qubit2{surname}", pin="readout"),
                    end_pin=Dict(
                        component="Coupler4", pin="second_end"
                    ),
                ),
                lead=Dict(
                    start_straight="50um",
                    end_straight="50um",
                ),
                meander=Dict(
                    asymmetry="0um",
                    spacing=p.meander_spacing,
                ),
                trace_width=f"{cpw_width*1e6}um",
                trace_gap=f"{cpw_gap*1e6}um",
                fillet=p.readout_res_fillet,
                total_length=f"{(l2_half * 1e3) - extra_length_mm:.2f}mm",
            ),
        ) '''

        ############################# Lanchpad, feedline and Markers ################################ 
        # Launchpads
        lp_left = LaunchpadWirebond(
            self.design,
            "LP_left",
            options=dict(
                **lp_options,
                pos_x=offset_x_mm + lp_coords_mm[0][0],
                pos_y=offset_y_mm,
                orientation="0",
            ),
        )
        lp_right = LaunchpadWirebond(
            self.design,
            "LP_right",
            options=dict(
                **lp_options,
                pos_x=offset_x_mm - lp_coords_mm[0][0],
                pos_y=offset_y_mm,
                orientation="180",
            ),
        )
        #
        #
        #
        # Feedline
        feedline = RouteStraight(
            self.design,
            "feedline_main",
            options=dict(
                pin_inputs=dict(
                    start_pin=dict(component="LP_left", pin="tie"),
                    end_pin=dict(component="LP_right", pin="tie"),
                ),
                trace_width=f"{cpw_width*1e6}um",
                trace_gap=f"{cpw_gap*1e6}um",
            ),
        )

        ''' # Markers
        if p.with_markers:
            # place markers in corners
            corners_mm = [
                (-0.5 * chip_x_mm, 0.5 * chip_y_mm),
                (0.5 * chip_x_mm, 0.5 * chip_y_mm),
                (0.5 * chip_x_mm, -0.5 * chip_y_mm),
                (-0.5 * chip_x_mm, -0.5 * chip_y_mm),
            ]
            inset_mm = 0.5
            corners_inset = [
                (-0.5 * chip_x_mm + inset_mm, 0.5 * chip_y_mm),
                (0.5 * chip_x_mm - inset_mm, 0.5 * chip_y_mm),
                (0.5 * chip_x_mm - inset_mm, -0.5 * chip_y_mm),
                (-0.5 * chip_x_mm + inset_mm, -0.5 * chip_y_mm),
            ]

            cross_markers = []
            square_markers = []
            for i, coords in enumerate(corners_mm):
                square_markers.append(
                    MarkerSquare4(
                        self.design,
                        f"{prefix}marker_square_{i}",
                        options=Dict(
                            pos_x=f"{offset_x_mm + corners_inset[i][0]}mm",
                            pos_y=f"{offset_y_mm + corners_inset[i][1]}mm",
                            layer=4,
                        ),
                    )
                )
                cross_markers.append(
                    MarkerDicingCross(
                        self.design,
                        f"{prefix}marker_cross_{i}",
                        options=Dict(
                            pos_x=f"{offset_x_mm + coords[0]}mm",
                            pos_y=f"{offset_y_mm + coords[1]}mm",
                            layer=0,
                        ),
                    )
                )'''