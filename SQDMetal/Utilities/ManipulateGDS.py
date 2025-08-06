# import gdspy
import gdstk
from SQDMetal.Utilities.QUtilities import QUtilities
from matplotlib.font_manager import FontProperties
from time import gmtime, strftime
import warnings


class ManipulateGDS:
    def __init__(
        self, import_path_gds, export_path_gds, import_cells=None, origin=(0, 0)
    ):
        """
        Class that will load a GDS file and perform actions such as making an array of the passed design (or a subsection of)

        Inputs:
            import_path_gds - Path to GDS file for import
            export_path_gds - Path to export GDS
            import_cells - (Optional) Cells to import (must be cells present in the passed GDS), defaults to all
            origin - (Defaults to (0,0)) Tuple containing relative origin to import the design
        """

        # initialise ManipulateGDS
        self.infile = import_path_gds
        self.outfile = export_path_gds
        self.import_cells = import_cells
        self.origin = origin
        self.gds_units = 1e-6

        # initiate gds cell and temp cell to help with copying
        self.lib = gdstk.Library("GDS Export")
        self.cell = self.lib.new_cell("Main")
        self.lib.replace(self.cell)

        # self.gds_in = gdstk.GdsLibrary(infile=import_path_gds)
        self.gds_in = gdstk.read_gds(infile=import_path_gds)

        # list available cells and import cells
        print(f"Available cells:\n {[c.name for c in self.gds_in.cells]}\n")
        if import_cells is None:
            total_polys = sum([len(c.polygons) for c in self.gds_in.cells])
            print(f"Imported: All ({total_polys} polygons)\n")
        else:
            print(f"Imported:\n {import_cells}\n")


    def flatten_cells(self, cell_keys=None, export=False):
        """
        Takes `cell_keys` (or all cells by default) from a GdsLibrary object and flattens them as seperate layers to a single cell export object.

        Inputs:
            - cell_keys - (Optional) List of cell names to extract from the passed GDS file
            - export - (Defaults to False) Choose whether to export the generated file, or just to operate on the design (default)
        """

        self.lib.cells.clear()

        timestamp = strftime("%Y%m%d_%H%M", gmtime())
        cell_name = f"Flattened_{timestamp}"
        temp_cell_name = f"Temp_{timestamp}"

        self.cell_flattened = self.lib.new_cell(cell_name)
        self.cell_temp = self.lib.new_cell(temp_cell_name)

        if not isinstance(cell_keys, list):
            print(f"Single cell: {cell_keys}")
        # get list of import cells (default: all)
        if cell_keys is None:
            cell_keys = [c.name for c in self.gds_in.cells]
        else:
            cell_keys = list(cell_keys)
        # flatten all cells to single cell with different layers
        for layer_idx, key in enumerate(cell_keys):
            # copy polygons
            to_copy = self.gds_in[key]
            #array = gdstk.Reference(to_copy, origin=self.origin)
            # self.cell_flattened = self.cell_temp.flatten().copy(
            #     name=cell_name
            # )  # flatten
            self.cell_flattened.add(*to_copy.polygons)
            # reset temp cell
            self.cell_temp.remove(
                *[poly for poly in self.cell.polygons if poly.layer == layer_idx]
            )
            print(f"Copied {key:<13}-> layer {layer_idx} ({len(to_copy.polygons)} polygons)")
        # Remove temp cell
        target_cell = next((cell for cell in self.lib.cells if cell.name == temp_cell_name), None)
        if target_cell:
            self.lib.remove(target_cell)

        # export
        if export:
            print(
                f"\nOutputting cell {cell_name} ({len(self.cell_flattened.get_polygonsets())} polygons)\n--> {self.outfile}\n\n"
            )

            self.lib.write_gds(self.outfile, cells=[cell_name])
        else:
            print(
                f"\n{cell_keys} were flattened to a single cell called '{cell_name}' with {len(self.cell_flattened.polygons)} polygons.\n"
            )

    def make_array_onChip(
        self,
        columns,
        rows,
        spacing=("50um", "50um"),
        chip_dimension=("20mm", "20mm"),
        export=True,
        export_path=None,
        use_cells=None,
        add_labels=False,
        label_string=None,
        label_offset=None,
        label_size=50,
    ):
        """
        Makes an array of a design input and exports it as a new GDS. Labels each arrayed structure.

        Inputs:
            - columns - number of columns in the array
            - rows - number of rows in the array
            - spacing - (Defaults to ("50um", "50um")) Tuple containing x- and y-spacing of the array as strings with units
            - chip_dimension - (Defaults to ("20mm", "20mm")) Dimension of the chip on which to array the structure (WIP)
            - export - (Defaults to True) Choose whether to export the generated file (default), or just to operate on the design
            - export_path - (Optional) Export path (if different from self)
            - use_cells - (Optional) List of cells from the input GDS to include in the arrayed structure
            - add_label - (Defaults to False) Add a label beneath each arrayed structure (WIP)
            - label_string - (Optional) If None, defaults to R1C1, R1C2 etc. (R{row_index}C{column_index}). Otherwise, if a single string is passed, this will be printed with an incrementing number. Otherwise, a list corresponding in length to the total number of arrayed structures can be passed and will be printed along row 1 from left to right, then row 2 etc.
            - label_offset - (Optional) Describes how far beneath the origin for each array point the label is. Defaults to half of spacing in y. Can be passed as a string with units (e.g. 300um).
            - label_size - (Defaults to 50) Text size for the printed labels

        Output:
            - cell_array - gdspy cell array containing arrayed structures
        """

        # warning that labels do not work currently
        if add_labels:
            warnings.warn(
                "Labels are not working currently! Setting add_labels = False and continuing the array process.\n\n"
            )
            add_labels = False

        # clear curent cells
        self.lib.cells.clear()

        # import assertions
        assert isinstance(chip_dimension, tuple) and isinstance(
            chip_dimension[0], str
        ), r"Input argument 'chip_dimension' should be a tuple contining strings of the chip's (x, y) dimensions with units - e.g. 'chip_dimension=('20mm', '20mm')', as in qiskit and SQDMetal."
        if isinstance(label_string, list):
            assert len(label_string) == (
                columns * rows
            ), "If you are passing a list of labels, please ensure it is the same length as the total number of arrayed structures (i.e. rows * columns)."

        # parse values
        sp_x = QUtilities.parse_value_length(spacing[0])  # noqa: F841 # abhishekchak52: unused variable sp_x
        sp_y = QUtilities.parse_value_length(spacing[1])

        lbl_offset = sp_y / 2 if label_offset is None else QUtilities.parse_value_length(label_offset)  # noqa: F841 # abhishekchak52: unused variable lbl_offset

        # set font
        fp = FontProperties(family="serif", style="italic")  # noqa: F841 # abhishekchak52: unused variable fp

        # default: copy all cells. Else use function input (if any), finally use class init input (if any)
        if use_cells is not None:
            self.flatten_cells(cell_keys=use_cells)
        elif (self.import_cells is None) and (use_cells is None):
            self.flatten_cells(cell_keys=None)
        else:
            self.flatten_cells(cell_keys=self.import_cells)

        # initialise an array cell and do arraying
        self.cell_array = self.cell_flattened.copy(name="Array")
        array_ref = gdstk.Reference(self.cell_array,
                                    columns=columns,
                                    rows=rows,
                                    spacing=(
                                        QUtilities.parse_value_length(spacing[0]) / self.gds_units,
                                        QUtilities.parse_value_length(spacing[1]) / self.gds_units,
                                        )
                                    )
        self.cell_array.add(array_ref)
        self.cell_array.flatten()
        self.lib.add(self.cell_array)

        # print(f'Cells: {[c.name for c in self.lib.cells]}')
        # print(f'Polys: {[len(c.polygons) for c in self.lib.cells]}')

        # export if option is set
        if export:
            if export_path is None:
                array_export = self.outfile
            else:
                array_export = export_path

            # export array only
            if not add_labels:
                # Keep only the cell named "Array"
                self.lib.remove(
                    *[cell for cell in self.lib.cells if cell.name != "Array"]
                )
                self.lib.write_gds(array_export)
                print(
                    f"\nOutputting cell 'Array' ({len(self.cell_array.polygons)} polygons)\n --> {array_export}\n"
                )
            # export array and text
            if add_labels:
                self.lib.remove(
                    *[
                        cell
                        for cell in self.lib.cells
                        if cell.name not in ["Array", "Text"]
                    ]
                )
                self.lib.write_gds(array_export)
                print(
                    f"\nOutputting cell: 'Array' ({len(self.cell_array.polygons)} polygons)\nOutputting cell: 'Text' ({len(self.cell_text.polygons)} polygons)\n--> {array_export}"
                )

        print(
            f"\nFrom {self.infile}:\n A {columns} x {rows} array of the design was made\n Spacing: {spacing}\n\nThe array was exported to: {self.outfile}"
        )

        return self.cell_array

    # depreciated: no use Comps/Labels.py for text
    # # for rendering system fonts
    # def render_text(self, text,
    #                 size=None,
    #                 position=(0, 0),
    #                 font_prop=None,
    #                 tolerance=0.1):

    #     path = TextPath(position, text, size=size, prop=font_prop)
    #     polys = []
    #     xmax = position[0]
    #     for points, code in path.iter_segments():
    #         if code == path.MOVETO:
    #             c = gdspy.Curve(*points, tolerance=tolerance)
    #         elif code == path.LINETO:
    #             c.L(*points)
    #         elif code == path.CURVE3:
    #             c.Q(*points)
    #         elif code == path.CURVE4:
    #             c.C(*points)
    #         elif code == path.CLOSEPOLY:
    #             poly = c.get_points()
    #             if poly.size > 0:
    #                 if poly[:, 0].min() < xmax:
    #                     i = len(polys) - 1
    #                     while i >= 0:
    #                         if gdspy.inside(
    #                             poly[:1], [polys[i]], precision=0.1 * tolerance
    #                         )[0]:
    #                             p = polys.pop(i)
    #                             poly = gdspy.boolean(
    #                                 [p],
    #                                 [poly],
    #                                 "xor",
    #                                 precision=0.1 * tolerance,
    #                                 max_points=0,
    #                             ).polygons[0]
    #                             break
    #                         elif gdspy.inside(
    #                             polys[i][:1], [poly], precision=0.1 * tolerance
    #                         )[0]:
    #                             p = polys.pop(i)
    #                             poly = gdspy.boolean(
    #                                 [p],
    #                                 [poly],
    #                                 "xor",
    #                                 precision=0.1 * tolerance,
    #                                 max_points=0,
    #                             ).polygons[0]
    #                         i -= 1
    #                 xmax = max(xmax, poly[:, 0].max())
    #                 polys.append(poly)
    #     return polys
