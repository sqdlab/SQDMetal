from qiskit_metal import designs
from SQDMetal.Utilities.Materials import Material 
from SQDMetal.Utilities.MakeGDS import MakeGDS
from SQDMetal.Utilities.CpwParams import CpwParams
from SQDMetal.Utilities.QUtilities import QUtilities
import numpy as np
import os
from datetime import datetime

class MultiDieChip:
    @staticmethod
    def make_resonator_chip(export_filename, 
                            export_path="", 
                            export_type="all",
                            frequency_range=(6e9, 7e9), 
                            num_resonators=5, 
                            cpw_width="9um", 
                            coupling_gap="20um", 
                            tl_y="0um",
                            res_vertical="1500um",
                            lp_to_res="300um",
                            lp_inset="0um",
                            lp_dimension="600um",
                            lp_taper="300um",
                            substrate_material="silicon", 
                            substrate_thickness="0.5mm", 
                            film_thickness="100nm",
                            chip_dimension=("20mm", "20mm"), 
                            chip_border="500um",     
                            die_dimension=("7.0mm", "4.0mm"),    
                            die_num=[1,1],
                            fill_chip=True,
                            markers_on=True,
                            text_label="",
                            text_size=800,
                            text_position=None,
                            print_all_infos=True
                            ):
        
        t = datetime.now().strftime("%Y%m%d_%H%M")  # add timestamp to export

        print(f"{t}\nBuilding chip \"{export_filename}\" with the following options:\n")
        print(locals())
        print('\n\n')

        # assign other values form input parameters
        substrate = Material(substrate_material)
        er = substrate.permittivity

        # initialise and clear design
        design = designs.DesignPlanar(
            metadata={}, overwrite_enabled=True, enable_renderers=True
        )
        design.delete_all_components()

        # set chip dimension and material
        design.chips.main.size.size_x = chip_dimension[0]
        design.chips.main.size.size_y = chip_dimension[1]
        design.chips.main.size.size_z = '-' + substrate_thickness
        design.chips.main.material = [substrate_material]

        # enact
        design.chips.main

        # initialise CpwParams object for calculations
        c = CpwParams(rel_permittivity=er, dielectric_thickness=QUtilities.parse_value_length(substrate_thickness))
        cpw_params = c.fromQDesign(design=design)

        # calculate gap from specified width for 50 Ohm impedence
        gap = str(f"{c.get_gap_from_width(trace_width=QUtilities.parse_value_length(cpw_width)) * 1e6:.2f}um")

        # calculate gap width ratio
        width_to_gap_ratio = QUtilities.parse_value_length(cpw_width) / QUtilities.parse_value_length(gap)

        # calculate cpw impedance
        cpw_impedance = CpwParams.calc_impedance(tr_wid=QUtilities.parse_value_length(cpw_width), tr_gap=QUtilities.parse_value_length(gap), er=er, h=QUtilities.parse_value_length(substrate_thickness))

        # calculate number of die in x and y if fill_chip is True
        if fill_chip:
            die_num[0] = int(
                (QUtilities.parse_value_length(chip_dimension[0]) - QUtilities.parse_value_length(chip_border))
                // QUtilities.parse_value_length(die_dimension[0])
            )
            die_num[1] = int(
                (QUtilities.parse_value_length(chip_dimension[1]) - QUtilities.parse_value_length(chip_border))
                // QUtilities.parse_value_length(die_dimension[1])
            )
        total_die_num = die_num[0] * die_num[1] # total number of die

        # list of tuples containing die centre coordinates
        die_coords = QUtilities.calc_die_coords(chip_dimension, die_dimension, die_num)

        # calculate frequencies
        freq_list = np.linspace(frequency_range[0], frequency_range[1], num_resonators, endpoint=True)

        # calculate total launchpad width (from edge of die)
        lp_extent = (
            QUtilities.parse_value_length(lp_inset)
            + QUtilities.parse_value_length(lp_dimension) / width_to_gap_ratio
            + QUtilities.parse_value_length(lp_dimension)
            + QUtilities.parse_value_length(lp_taper)
            + (25 * 1e-6)
        )

        print('CPW properties:')
        print(f'\tCPW gap: {gap}')
        print(f'\tCPW width: {cpw_width}')
        print(f'\tCPW impedance: {cpw_impedance:.2f} Ohm')
        print('\nLaunchpad properties:')
        print(f'\tLaunchpad extent: {lp_extent * 1e6:.2f}um\n\n')

        print('~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~\n')

        # place and store launchpads, resonators, transmission lines
        launchpads = []
        resonators_list, resonator_vals, resonator_names, capacitances, inductances = [], [], [], [], []
        transmission_lines = []

        # loop through dies
        for i, origin in enumerate(die_coords):
            
            # draw launchpads
            launchpads.append(
                QUtilities.place_launchpads(
                    design=design,
                    cpw_gap=gap,
                    cpw_width=cpw_width,
                    die_origin=origin,
                    die_dimension=die_dimension,
                    die_number=(i + 1),
                    dimension=lp_dimension,
                    inset=lp_inset,
                    taper=lp_taper,
                )
            )

            # draw resonators
            resonators, resonator_val, resonator_name, capacitance, inductance = QUtilities.place_resonators_hanger(design=design,
                gap=gap,
                width=cpw_width,
                num_resonators=num_resonators,
                frequencies=freq_list,
                die_origin=origin,
                die_dimension=die_dimension,
                die_index=i,
                film_thickness=film_thickness,
                launchpad_extent=lp_extent,
                coupling_gap=coupling_gap,
                transmission_line_y=tl_y,
                launchpad_to_res="250um",
                min_res_gap="50um",
                LC_calculations=True,
                print_statements=True
                )
            
            resonators_list.append(resonators)
            resonator_vals.append(resonator_val)
            resonator_names.append(resonator_name)
            capacitances.append(capacitance)
            inductances.append(inductance)

            # draw transmission lines
            transmission_lines.append(QUtilities.place_transmission_line_from_launchpads(
                design=design,
                tl_y=tl_y,
                gap=gap,
                width=cpw_width,
                die_index=i))

            # draw markers
            QUtilities.place_markers(design=design, die_origin=origin, die_dim=die_dimension)

        # setup GDS export (positive)
        gds_export = MakeGDS(design, export_type='positive')

        # add text label
        gds_export.add_text(text_label="Full chip maker test", size=600, position=(0.05, 0.93))

        # setup export path based on user inputs
        if export_path=="":
            full_export_path = os.path.join(f"{export_filename}_{t}.gds")
        else:
            full_export_path = os.path.join(export_path, f"{export_filename}_{t}.gds")

        print('\n\n~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~\n')

        # do export
        gds_export.export(full_export_path)

        print(f"Exported at {full_export_path}")