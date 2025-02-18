from qiskit_metal import designs
from SQDMetal.Utilities.Materials import Material
from SQDMetal.Utilities.MakeGDS import MakeGDS
from SQDMetal.Utilities.CpwParams import CpwParams
from SQDMetal.Utilities.QUtilities import QUtilities
import numpy as np
import matplotlib.pyplot as plt
import os
from pprint import pprint
from datetime import datetime
from shapely.geometry import Polygon


class MultiDieChip:

    def __init__(self, export_filename=None):
        self.export_filename = export_filename
        self.design = None

    # @staticmethod
    def make_resonator_chip(
        self,
        export_filename=None,
        export_path="",
        export_type="all",
        export_threshold=1e-9,
        frequency_range=(6e9, 7e9),
        num_resonators=5,
        cpw_width="9um",
        feedline_upscale=1.0,
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
        die_dimension=("7.1mm", "4.4mm"),
        die_num=[1, 1],
        fill_chip=False,
        markers_on=True,
        text_label="",
        text_size=600,
        text_position=None,
        print_all_infos=True,
        date_stamp_on_export=True,
        single_circuit_for_simulation=False,
        plot_inline_mpl=True
    ):
        """
        Creates a `.gds` full-wafer layout file for a simple coplanar waveguide $\lambda/4$ resonator chip containing a number of resonators (usually 5) capacitively coupled to a transmission line.

        Inputs:
            - export_filename - Filename for gds export (e.g. "test")
            - export_path - Path for export (e.g. 'exports'); the file will then be output to /exports/test.gds
            - export_type - (Defaults to "all") Export type for lithography as per `MakeGDS` (options: "all", "positive", "negative")
            - export_threshold - (Defaults to 1e-9) The smallest feature width that can exist; anything smaller will get culled out. This is to help remove artefacts from floating-point inaccuracies. It is given in metre
            - frequency_range - (Defaults to (6e9, 7e9)) Tuple containing minimum and maximum resonator frequencies in Hz
            - num_resonators - (Defaults to 5) Number of resonators per die
            - cpw_width - (Defaults to "9um") Width of the central trace on the resonators. The gap will be automatically calculated for 50 Ohm impedance based on the `substrate_material`. If feedline_upscale==0, the feedline width will match the resonator width. You can pass a list of strings if you want to scale each resonator seperately.
            - feedline_upscale - (Defaults to 1.0) scale of the feedline gap and width as a multiple of the resonator dimensions
            - coupling_gap - (Defaults to "20um") Amount of ground plane in the coupling gap between the feedline and the resonator
            - tl_y - (Defaults to "0um") The die-relative y-value for the main straight of the feedline
            - res_vertical - (Defaults to "1500um") Vertical length of resonator meanders
            - lp_to_res - (Defaults to "300um") Minimum distance between the launchpad taper and the coupling length of the left-most resonator
            - lp_inset - (Defaults to "0um") Inset of the launchpads along the x-axis relative to the die boundaries
            - lp_dimension - (Defaults to "600um") Width of the launchpads' conductive centre pad (the launchpad gap scales accordingly)
            - lp_taper - (Defaults to "300um") Length of the taper from launchpad to feedline
            - substrate_material - (Defaults to "silicon") Substrate material (currently only "silicon" and "sapphire" are supported)
            - substrate_thickness - (Defaults to "0.5mm") Substrate thickness
            - film_thickness - (Defaults to "100nm") Film thickness
            - chip_dimension - (Defaults to ("20mm", "20mm")) Dimensions of the chip as an (x, y) Tuple
            - chip_border - (Defaults to "500um") Chip border to leave un-patterned
            - die_dimension - (Defaults to ("7.1mm", "4.4mm")) Dimensions of the die as an (x, y) Tuple
            - die_num - (Defaults to [1, 1])) Die layout in [x, y] as a length-2 list
            - fill_chip - (Defaults to True) Boolean to choose whether the full chip is automatically populated with dies (over-rides die_num if True)
            - markers_on - (Defaults to True) Print dicing markers on export
            - text_label - (Optional) Text label to print on chip. Prints a default message if text_label="". Set text_label=None for no text label.
            - text_size - (Defaults to 600) Text size
            - text_position - (Optional) Tuple of text label location as normalised (x, y) (e.g. (0.1, 0.9) will place the text label 1/10th of the way along the chip in the x-direction, and 9/10ths of the way up the chip in the y-direction)
            - print_all_infos - (Defaults to True) Choose whether to print info as the `.gds` is being generated
            - date_stamp_on_export - (Defaults to True) If True, will print date and time of export in exported .gds filename

        Outputs:
            - design - Qiskit Metal design object for the generated chip
        """

        assert self.export_filename != None or export_filename != None, "Please provide an export filename"
        if self.export_filename == None:
            self.export_filename = export_filename
        elif export_filename == None:
                export_filename = self.export_filename
        
        if single_circuit_for_simulation and fill_chip:
            print("\nWarning: single_circuit_for_simulation=True will override fill_chip=True. We will only print a single chip (die_num=[1, 1]).\n")
            fill_chip = False
            die_num = [1, 1]

        # check if chip geometry is constant or scaled per-resonator
        assert QUtilities.is_string_or_list_of_strings(cpw_width)
        scaled_geometry = isinstance(cpw_width, list)
        if scaled_geometry == True:
            assert (
                len(cpw_width) > 1
            ), "If cpw_width is a list, it must have more than one element"
            assert (
                len(cpw_width) == num_resonators
            ), "If cpw_width is a list, it must have the same length as num_resonators"

        # TODO: add automatic Palace sim
        # TODO: add per-die labels for easy ID during fabrication

        t = datetime.now().strftime("%Y%m%d_%H%M")  # add timestamp to export

        print(f'{t}\nBuilding chip "{export_filename}" with the following options:\n')
        if print_all_infos:
            pprint(locals(), sort_dicts=False)

        # assign other values form input parameters
        substrate = Material(substrate_material)
        er = substrate.permittivity

        # initialise and clear design
        design = designs.DesignPlanar(
            metadata={}, overwrite_enabled=True, enable_renderers=True
        )
        design.delete_all_components()

        # convert some values to floats for ease-of-use
        chip_dim_x = QUtilities.parse_value_length(chip_dimension[0])
        chip_dim_y = QUtilities.parse_value_length(chip_dimension[1])
        die_dim_x = QUtilities.parse_value_length(die_dimension[0])
        die_dim_y = QUtilities.parse_value_length(die_dimension[1])
        chip_border_xy = QUtilities.parse_value_length(chip_border)
        res_vertical_y = QUtilities.parse_value_length(res_vertical)

        # set chip dimension and material
        design.chips.main.size.size_x = chip_dimension[0]
        design.chips.main.size.size_y = chip_dimension[1]
        design.chips.main.size.size_z = "-" + substrate_thickness
        design.chips.main.material = [substrate_material]

        # enact
        design.chips.main

        # initialise CpwParams object for calculations
        c = CpwParams(
            rel_permittivity=er,
            dielectric_thickness=QUtilities.parse_value_length(substrate_thickness),
        )
        cpw_params = c.fromQDesign(design=design)

        # calculate gap from specified width for 50 Ohm impedence
        if scaled_geometry == True:
            gap = [
                str(
                    f"{c.get_gap_from_width(trace_width=QUtilities.parse_value_length(i)) * 1e6:.2f}um"
                )
                for i in cpw_width
            ]
            width_0 = cpw_width[0]
            gap_0 = gap[0]
        else:
            gap = str(
                f"{c.get_gap_from_width(trace_width=QUtilities.parse_value_length(cpw_width)) * 1e6:.2f}um"
            )
            gap_0 = gap
            width_0 = cpw_width

        # calculate gap width ratio
        width_to_gap_ratio = QUtilities.parse_value_length(
            width_0
        ) / QUtilities.parse_value_length(gap_0)

        # calculate cpw impedance
        cpw_impedance = CpwParams.calc_impedance(
            tr_wid=QUtilities.parse_value_length(width_0),
            tr_gap=QUtilities.parse_value_length(gap_0),
            er=er,
            h=QUtilities.parse_value_length(substrate_thickness),
        )

        # calculate number of die in x and y if fill_chip is True
        if fill_chip:
            die_num[0] = int((chip_dim_x - chip_border_xy) // die_dim_x)
            die_num[1] = int((chip_dim_y - chip_border_xy) // die_dim_y)
            assert (
                die_num[0] > 0 and die_num[1] > 0
            ), "Chip dimensions must be larger than twice the chip border if fill_chip=True"
            print(
                f"\nChip fill enabled - overriding die_num. Writing {die_num[0]} x {die_num[1]} dies to fill chip.\n"
            )

        # list of tuples containing die centre coordinates (floats in units of meters)
        die_coords = QUtilities.calc_die_coords(chip_dimension, die_dimension, die_num)

        # calculate frequencies
        freq_list = np.linspace(
            frequency_range[0], frequency_range[1], num_resonators, endpoint=True
        )

        # calculate total launchpad width (from edge of die)
        lp_extent = (
            QUtilities.parse_value_length(lp_inset)
            + QUtilities.parse_value_length(lp_dimension) / width_to_gap_ratio
            + QUtilities.parse_value_length(lp_dimension)
            + QUtilities.parse_value_length(lp_taper)
            + (25 * 1e-6)
        )

        # calculate feedline dimensions (if different from resonators)
        if feedline_upscale != 1.0:
            feedline_cpw_width = (
                f"{(QUtilities.parse_value_length(width_0) * feedline_upscale) * 1e6}um"
            )
            feedline_cpw_gap = (
                f"{(QUtilities.parse_value_length(gap_0) * feedline_upscale) * 1e6}um"
            )
        else:
            feedline_cpw_width = width_0
            feedline_cpw_gap = gap_0

        print(f"Chip Size: {chip_dimension[0]} x {chip_dimension[1]}")
        print(f"Die Size:  {die_dimension[0]} x {die_dimension[1]}\n")
        print(f"Resonator properties:")
        print(f" Frequencies    : {[f'{i * 1e-9} GHz' for i in freq_list]}")
        print(f" Gap            : {gap}")
        print(f" Width          : {cpw_width}")
        print(f" Impedance      : {cpw_impedance:.2f} Ω\n")
        print(f"Launchpad extent: {lp_extent * 1e6:.2f} µm\n")

        # place and store launchpads, resonators, transmission lines
        launchpads = []
        resonators_list, resonator_vals, resonator_names, capacitances, inductances = (
            [],
            [],
            [],
            [],
            [],
        )
        transmission_lines = []

        # loop through dies
        for i, origin in enumerate(die_coords):

            # draw launchpads
            launchpads.append(
                QUtilities.place_launchpads(
                    design=design,
                    cpw_gap=feedline_cpw_gap,
                    cpw_width=feedline_cpw_width,
                    die_origin=origin,
                    die_dimension=die_dimension,
                    die_number=(i + 1),
                    dimension=lp_dimension,
                    inset=lp_inset,
                    taper=lp_taper,
                    print_checks=print_all_infos,
                )
            )

            # TODO: update calculation of resonator start position (in y) to account for up-scaled feedline as given by feedline_upscale argument.

            # draw resonators
            resonators, resonator_val, resonator_name, capacitance, inductance = (
                QUtilities.place_resonators_hanger(
                    design=design,
                    gap=gap,
                    width=cpw_width,
                    num_resonators=num_resonators,
                    frequencies=freq_list,
                    die_origin=origin,
                    die_dimension=die_dimension,
                    die_index=i,
                    film_thickness=film_thickness,
                    launchpad_extent=lp_extent,
                    feedline_upscale=feedline_upscale,
                    coupling_gap=coupling_gap,
                    transmission_line_y=tl_y,
                    launchpad_to_res=lp_to_res,
                    res_vertical=res_vertical,
                    min_res_gap="50um",
                    LC_calculations=True,
                    print_statements=print_all_infos,
                )
            )

            resonators_list.append(resonators)
            resonator_vals.append(resonator_val)
            resonator_names.append(resonator_name)
            capacitances.append(capacitance)
            inductances.append(inductance)

            # draw transmission lines
            transmission_lines.append(
                QUtilities.place_transmission_line_from_launchpads(
                    design=design,
                    tl_y=tl_y,
                    gap=feedline_cpw_gap,
                    launchpads=launchpads,
                    width=feedline_cpw_width,
                    die_index=i,
                    die_origin=origin,
                )
            )

            # draw markers
            if markers_on:
                QUtilities.place_markers(
                    design=design, die_origin=origin, die_dim=die_dimension
                )

        # add text label
        if text_label != None:
            if text_position == None:
                gds_export.add_text(
                    text_label=text_label, size=text_size, position=(0.05, 0.9)
                )
            else:
                gds_export.add_text(
                    text_label=text_label, size=text_size, position=text_position
                )

        if single_circuit_for_simulation:
            # # Create the chip boundary (polygon)
            sim_border = 0.1  # mm
            minxy = [np.inf, np.inf]
            maxxy = [-1*np.inf, -1*np.inf]
            # Iterate over all components in the design to check smalles boundary
            for comp in design.components.values():
                # Get the shape of the component
                (minx, miny, maxx, maxy) = comp.qgeometry_bounds()
                if minx < minxy[0]:
                    minxy[0] = minx
                if miny < minxy[1]:
                    minxy[1] = miny
                if maxx > maxxy[0]:
                    maxxy[0] = maxx
                if maxy > maxxy[1]:
                    maxxy[1] = maxy
            # Add a border around the chip
            x_bounds = max(abs(maxxy[0] + sim_border), abs(minxy[0] - sim_border))
            y_bounds = max(abs(minxy[1] - sim_border), abs(maxxy[1] + sim_border))
            print(f"Resizing chip for simulation: {(x_bounds*2)}mm x {(y_bounds*2)}mm\n")
            design.chips.main.size.size_x = f"{(x_bounds*2)}mm"
            design.chips.main.size.size_y = f"{(y_bounds*2)}mm"
            design.chips.main
            design.rebuild()
        
        # setup GDS export (positive)
        gds_export = MakeGDS(
            design,
            export_type=export_type,
            print_statements=print_all_infos,
            threshold=export_threshold,
            precision=export_threshold,
        )

        # only export if a filename is passed
        if export_filename != None:
            # setup export path based on user inputs
            if export_path == "":
                if date_stamp_on_export == True:
                    full_export_path = os.path.join(f"{export_filename}_{t}.gds")
                else:
                    full_export_path = os.path.join(f"{export_filename}.gds")
            else:
                if date_stamp_on_export == True:
                    full_export_path = os.path.join(
                        export_path, f"{export_filename}_{t}.gds"
                    )
                else:
                    full_export_path = os.path.join(
                        export_path, f"{export_filename}.gds"
                    )
            # do export
            gds_export.export(full_export_path)
            print(f"Exported at {os.path.abspath(full_export_path)}")
            print(f"")

        # assign self.design
        self.design = design

        if plot_inline_mpl == True:
            fig, ax = plt.subplots(figsize=(chip_dim_x*3e3, chip_dim_y*3e3))
            for comp in design.components.values():
                comp.qgeometry_plot(ax)
            full_export_path = os.path.join(
                        export_path, f"{export_filename}.png"
                    )
            fig.savefig(full_export_path)
            print(f"Circuit diagram saved at {os.path.abspath(full_export_path)}\n")

        # return design regardless of whether the GDS was exported or not
        return design
    
    def run_palace_sims(self):
        assert self.design != None, "Please run make_resonator_chip() first."
        # add simulation here
        print(f"Not yet implemented.")
