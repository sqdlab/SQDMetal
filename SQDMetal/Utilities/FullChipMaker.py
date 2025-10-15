# Copyright 2025 Prasanna Pakkiam
# SPDX-License-Identifier: Apache-2.0

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
import psutil


class MultiDieChip:

    def __init__(self, export_filename=None):
        self.export_filename = export_filename
        self.design = None
        self.substrate_material = None
        self.film_material = None
        self.start_freq = None
        self.num_resonators = None
        self.eigen_sim = None

    def make_resonator_chip(
        self,
        export_filename=None,
        export_path="",
        export_type="all",
        export_threshold=1e-9,
        frequency_range=(6e9, 7e9),
        frequency_list=None,
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
        film_material="aluminium",
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
        plot_inline_mpl=True,
        fillet="85um",
        radius="100um",
    ):
        """
        Creates a `.gds` full-wafer layout file for a simple coplanar waveguide $\lambda/4$ resonator chip containing a number of resonators (usually 5) capacitively coupled to a transmission line.

        Inputs:
            - export_filename - Filename for gds export (e.g. "test")
            - export_path - Path for export (e.g. 'exports'); the file will then be output to /exports/test.gds
            - export_type - (Defaults to "all") Export type for lithography as per `MakeGDS` (options: "all", "positive", "negative")
            - export_threshold - (Defaults to 1e-9) The smallest feature width that can exist; anything smaller will get culled out. This is to help remove artefacts from floating-point inaccuracies. It is given in metre
            - frequency_range - (Defaults to (6e9, 7e9)) Tuple containing minimum and maximum resonator frequencies in Hz
            - frequency_list - (Defaults to None) List containing discrete frequencies (must be same length as num_resonators) in Hz
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
            - single_circuit_for_simulation - (Defaults to False) If true, the multi-chip layour will be cropped to a single chip that can be used for simulation

        Outputs:
            - design - Qiskit Metal design object for the generated chip
        """

        if frequency_list:
            assert len(frequency_list) == num_resonators, "frequency_list must be the same length as num resonators."

        assert film_material in [
            "aluminium",
            "tantalum",
            "TiN"
        ], "Only aluminium and tantalum are supported for the film material currently."

        assert (
            self.export_filename is not None or export_filename is not None
        ), "Please provide an export filename"
        if self.export_filename is None:
            self.export_filename = export_filename
        elif export_filename is None:
            export_filename = self.export_filename

        if single_circuit_for_simulation:
            print(
                "\nPreparing as single_circuit_for_simulation.\n"
            )
            fill_chip = False
            die_num = [1, 1]
            markers_on = False
            text_label = None
            text_position = None
            chip_border = "0um"
            export_type = "positive"

        # check if chip geometry is constant or scaled per-resonator
        assert QUtilities.is_string_or_list_of_strings(cpw_width)
        assert QUtilities.is_string_or_list_of_strings(coupling_gap), "If you want different coupling_gap per-resonator, pass as a list of strings (e.g. ['5um, '10um', ...])"
        if isinstance(coupling_gap, list):
            assert len(coupling_gap) == num_resonators
        if len(cpw_width) == 1:
            cpw_width = cpw_width[0]
        scaled_geometry = isinstance(cpw_width, list)
        if scaled_geometry:
            assert (
                len(cpw_width) == num_resonators
            ), "If cpw_width is a list, it must have the same length as num_resonators"

        # update class variables
        self.substrate_material = substrate_material
        self.film_material = film_material
        self.start_freq = frequency_range[0]
        self.num_resonators = num_resonators
        self.fillet = fillet
        self.radius = radius

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
        if scaled_geometry:
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
        if frequency_list is None:
            freq_list = np.linspace(
                frequency_range[0], frequency_range[1], num_resonators, endpoint=True
            )
        elif frequency_list:
            freq_list = frequency_list

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
        print("Resonator properties:")
        print(f" Frequencies    : {[f'{i * 1e-9:.3f} GHz' for i in freq_list]}")
        print(f" Gap            : {gap}")
        print(f" Width          : {cpw_width}")
        print(f" Impedance      : {cpw_impedance:.2f} Ω")
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

            if print_all_infos:
                print(f'\nPrinting die {i}')

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

            # draw resonators
            '''
            Add resonators to the design. In progress: can we route better to avoid any bugs (WIP)? Testing filleting options, radius (?), etc. 
            TODO: finish this function and ensure the results are same as older designs... 
            There should be an option for qiskit_metal or sqdmetal routing (hopefully). 
            '''
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
                    radius=self.radius,
                    fillet=self.fillet
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
                    design=design, die_origin=origin, die_dim=die_dimension, die_index=i
                )

        if single_circuit_for_simulation:
            # update chip boundary
            sim_border = 0.4  # mm
            minxy = [np.inf, np.inf]
            maxxy = [-1 * np.inf, -1 * np.inf]
            print("\nCropping chip for simulation:")
            # iterate over all components in the design to check smallest boundary
            for i, comp in enumerate(design.components.values()):
                (minx, miny, maxx, maxy) = comp.qgeometry_bounds()
                if minx < minxy[0]:
                    minxy[0] = minx
                if miny < minxy[1]:
                    minxy[1] = miny
                if maxx > maxxy[0]:
                    maxxy[0] = maxx
                if maxy > maxxy[1]:
                    maxxy[1] = maxy
                comp_lims = [minx, miny, maxx, maxy]
                print(f" Comp {i:02}: {str(comp.name):30} -->   " + " ".join(f"{v:8.4f}" for v in comp_lims) + " [mm]")
            # add a border around the minimum boundary and rebuild
            x_total = (maxxy[0] - minxy[0]) + (2 * sim_border)
            y_total = (maxxy[1] - minxy[1]) + (2 * sim_border)
            new_chip_centre = [(minxy[0] + maxxy[0]) / 2, (minxy[1] + maxxy[1]) / 2]
            print(
                f"\n Resized chip to {x_total:.2f} x {y_total:.2f} mm\n"
                f"  xmin = {minxy[0]:6.3f} mm\n"
                f"  xmax = {maxxy[0]:6.3f} mm\n"
                f"  ymin = {minxy[1]:6.3f} mm\n"
                f"  ymax = {maxxy[1]:6.3f} mm\n"
            )
            print(f"  New x-y chip centre : {new_chip_centre[0]:6.2f}, {new_chip_centre[1]:6.2f} [mm]")
            print(f"  New x-y dimensions  : {x_total:6.2f}, {y_total:6.2f} [mm]\n")
            design.chips.main.size.center_x = f"{new_chip_centre[0]}mm"
            design.chips.main.size.center_y = f"{new_chip_centre[1]}mm"
            design.chips.main.size.size_x = f"{x_total}mm"
            design.chips.main.size.size_y = f"{y_total}mm"
            design.chips.main
            design.rebuild()

        # setup GDS export
        gds_export = MakeGDS(
            design,
            export_type=export_type,
            print_statements=print_all_infos,
            threshold=export_threshold,
            precision=export_threshold,
        )

        # add text label
        if text_label is not None:
            if text_position is None:
                gds_export.add_text(
                    text_label=text_label, size=text_size, position=(0.05, 0.9)
                )
            else:
                gds_export.add_text(
                    text_label=text_label, size=text_size, position=text_position
                )

    
        # setup export path based on user inputs
        if export_path == "":
            if date_stamp_on_export:
                full_export_path = os.path.join(f"{export_filename}_{t}.gds")
            else:
                full_export_path = os.path.join(f"{export_filename}.gds")
        else:
            if date_stamp_on_export:
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
        print("")

        # assign self.design
        self.design = design

        if plot_inline_mpl:
            fig, ax = plt.subplots(figsize=(chip_dim_x * 3e3, chip_dim_y * 3e3))
            ax.set_aspect('equal', adjustable='datalim')
            for comp in design.components.values():
                comp.qgeometry_plot(ax)
            full_export_path = os.path.join(export_path, f"{export_filename}.png")
            fig.savefig(full_export_path)
            print(f"Circuit diagram saved at {os.path.abspath(full_export_path)}\n")

        # return design regardless of whether the GDS was exported or not
        return design

    def setup_palace_eigenmode_resonator(
        self,
        palace_binary: str,
        num_eigenmodes=None,
        run_locally=True,
        leave_free_cpu_num=None,
        num_cpus=None,
        fine_mesh_min_max=(14e-6, 250e-6),
        start_freq=None,
        num_saved_solns=None,
        solver_order=2,
        solver_tol=1.0e-8,
        **kwargs,
    ):
        """
        Setup and run Palace eigenmode simulations on a resonator chip (as per MultiDieChip.make_resonator_chip) Qiskit Metal design object. Automatically handles number of eigenmodes, start frequency, and fine-meshing around transmission line and resonators. Assumes that ports are setup on launchpads.
        """
        assert self.design is not None, "Please run make_resonator_chip() first."
        assert (
            run_locally
        ), "Only local simulations are supported at the moment."
        # CPU num
        assert (leave_free_cpu_num or num_cpus) != None, "You must supply either 'num_cpus' or 'leave_free_cpu_num'."
        if leave_free_cpu_num == None:
            num_cpus = num_cpus
        elif num_cpus == None:
            physical_cores = psutil.cpu_count(logical=False)
            num_cpus = physical_cores - leave_free_cpu_num, 
            # print(f"Number of physical CPU cores: {physical_cores}")
        else:
            num_cpus = num_cpus
        # import Palace
        try:
            from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation # type: ignore
            from SQDMetal.Utilities.Materials import MaterialInterface # type: ignore

            print(
                "SQDMetal Palace modules successfully imported. Continuing with simulation."
            )
        except ImportError:
            print("PALACE Eigenmode Simulation module not found. Exiting simulation")
            exit()
        if start_freq is None:
            start_freq = self.start_freq
        assert start_freq > 1e9, "Starting frequency should be given in Hz. It looks like you provided a value in GHz."
        # set number of eigenmodes
        if num_eigenmodes is None:
            num_eigenmodes = self.num_resonators + 1
        # set number of saved solutions
        if num_saved_solns is None:
            num_saved_solns = num_eigenmodes # by default, save every eigenmode
        # Fine meshing resonators, tranmission line and launchpads
        fine_mesh_components_1 = []
        print("\nResonator component names:")
        for name in self.design.components.keys():
            if name.endswith("GHz"):
                print(f" {name}")
                fine_mesh_components_1.append(name)
        print("\nTransmission line name:")
        for name in self.design.components.keys():
            if name.startswith("tl"):
                print(f" {name}")
                fine_mesh_components_1.append(name)
        fine_mesh_launchpads = []
        print("\nLaunchpad names:")
        for name in self.design.components.keys():
            if name.startswith("lp"):
                print(f" {name}")
                fine_mesh_launchpads.append(name)
        assert len(fine_mesh_launchpads) == 2, f"There should only be 2 launchpad components, currently there are {len(fine_mesh_launchpads)}."
        # Setup for local simulation
        if run_locally:
            sim_mode = "SimPC"  # noqa: F841 # abhishekchak52: unused variable sim_mode
        else:
            print("HPC setup not yet supported.")
            exit()
        # Eigenmode Simulation Options
        user_defined_options = {
            "mesh_refinement": 0,  # refines mesh in PALACE - essetially divides every mesh element in half
            "dielectric_material": self.substrate_material,  # choose dielectric material - 'silicon' or 'sapphire'
            "starting_freq": start_freq - 0.5e9,  # starting frequency in Hz
            "number_of_freqs": num_eigenmodes,  # number of eigenmodes to find
            "solns_to_save": num_saved_solns,  # number of electromagnetic field visualizations to save
            "solver_order": solver_order,  # increasing solver order increases accuracy of simulation, but significantly increases sim time
            "solver_tol": solver_tol,  # error residual tolerance foriterative solver
            "solver_maxits": 200,  # number of solver iterations
            "mesh_sampling": 130,  # number of points to mesh along a geometry
            "fillet_resolution": 12,  # Number of vertices per quarter turn on a filleted path
            "palace_dir": palace_binary,  # "PATH/TO/PALACE/BINARY",
            "num_cpus": num_cpus,  # number of CPUs to use for simulation
            "mesh_min": fine_mesh_min_max[0]*1e3,  # minimum mesh size in mm
            "mesh_max": fine_mesh_min_max[1]*1e3  # maximum mesh size in mm
        }
        user_defined_options.update(kwargs) # update any other options based on kwargs
        max_option_length = max(len(str(option)) for option in user_defined_options)
        print("\n\nRunning simulation with the following options:")
        for option, value in user_defined_options.items():
            print(f"{option:<{max_option_length}} : {value}")
        if self.export_filename.endswith(".gds"):
            sim_name = self.export_filename[:-4]
        else:
            sim_name = self.export_filename
        # Creat the Palace Eigenmode simulation
        self.eigen_sim = PALACE_Eigenmode_Simulation(
            name=sim_name,  # name of simulation
            metal_design=self.design,  # feed in qiskit metal design
            sim_parent_directory="",  # choose directory where mesh file, config file and HPC batch file will be saved
            mode="simPC",  # choose simulation mode 'HPC' or 'simPC'
            meshing="GMSH",  # choose meshing 'GMSH' or 'COMSOL'
            user_options=user_defined_options,  # provide options chosen above
            view_design_gmsh_gui=False,  # view design in GMSH gui
            create_files=True,
        )  # create mesh, config and HPC batch files
        self.eigen_sim.add_metallic(1)
        self.eigen_sim.add_ground_plane()
        # Add in the RF ports
        self.eigen_sim.create_port_CPW_on_Launcher("lp_L_die0", 30e-6)
        self.eigen_sim.create_port_CPW_on_Launcher("lp_R_die0", 30e-6)
        # Fine-mesh routed paths (resonators, transmission line)
        # TODO: add option for mesh_around_path (new option that David wrote for CPW meshing)
        self.eigen_sim.fine_mesh_around_comp_boundaries(
            fine_mesh_components_1, min_size=fine_mesh_min_max[0], max_size=fine_mesh_min_max[1]
        )
        # Fine-mesh launchpads (less fine)
        self.eigen_sim.fine_mesh_around_comp_boundaries(
            fine_mesh_launchpads, min_size=fine_mesh_min_max[1], max_size=250e-6
        )
        # setup interfaces
        if (self.substrate_material == "silicon") and (
            self.film_material == "aluminium"
        ):
            self.eigen_sim.setup_EPR_interfaces(
                metal_air=MaterialInterface("aluminiumair"),
                substrate_air=MaterialInterface("siliconair"),
                substrate_metal=MaterialInterface("aluminiumsilicon"),
            )
        elif (self.substrate_material == "silicon") and (
            self.film_material == "tantalum"
        ):
            self.eigen_sim.setup_EPR_interfaces(
                metal_air=MaterialInterface("tantalumair"),
                substrate_air=MaterialInterface("siliconair"),
                substrate_metal=MaterialInterface("tantalumsilicon"),
            )
        self.eigen_sim.prepare_simulation()
        print(f"\nGMSH mesh exported at {os.path.abspath(f'{sim_name}.msh')}\n")
    
    def run_palace_eigenmode_resonator(self):
        self.eigen_sim.run()
