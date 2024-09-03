import gdspy
from SQDMetal.Utilities.QUtilities import QUtilities
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath

class ManipulateGDS:
    def __init__(self, import_path_gds, export_path_gds, import_cells=None, origin=(0,0)):
        '''
        Class that will load a GDS file and perform actions such as making an array of the passed design (or a subsection of)

        Inputs:
            import_path_gds - Path to GDS file for import
            export_path_gds - Path to export GDS
            import_cells - (Optional) Cells to import (must be cells present in the passed GDS), defaults to all
            origin - (Defaults to (0,0)) Tuple containing relative origin to import the design
        '''
        # clear gds library
        gdspy.current_library.cells.clear()
        # initialise ManipulateGDS
        self.infile = import_path_gds
        self.outfile = export_path_gds
        self.import_cells = import_cells
        self.origin = origin
        self.gds_units = 1e-6
        # initiate gds cell and temp cell to help with copying
        self.lib = gdspy.GdsLibrary()
        # self.cell = self.lib.new_cell('Array')
        self.gds_in = gdspy.GdsLibrary(infile=import_path_gds)
        # list available cells and import cells
        print(f'\n\nAvailable cells:\n\t{list(self.gds_in.cells.keys())}\n')
        if import_cells == None:
            print(f'Imported: All\n')
        else:
            print(f'Imported:\n\t{import_cells}\n')


    def flatten_cells(self, cell_keys=None, export=False):
        """
        Takes `cell_keys` (or all cells by default) from a GdsLibrary object and flattens them as seperate layers to a single cell export object.

        Inputs:
            - cell_keys - (Optional) List of cell names to extract from the passed GDS file
            - export - (Defaults to False) Choose whether to export the generated file, or just to operate on the design (default)
        """
    
        self.cell_flattened = self.lib.new_cell('Flattened')
        self.cell_temp = self.lib.new_cell('Temp')

        if not isinstance(cell_keys, list):
            print(f'Single cell: {cell_keys}')
        # get list of import cells (default: all)
        if cell_keys==None:
            cell_keys = list(self.gds_in.cells.keys())
        else:
            cell_keys = list(cell_keys)
        # flatten all cells to single cell with different layers
        for layer_idx, key in enumerate(cell_keys):
            print(f"Copying cell '{key}' --> layer {layer_idx}")
            # copy polygons to temp cell (and shift origin if needed)
            self.cell_temp.add((gdspy.CellReference(self.gds_in.cells[key], origin=self.origin)).get_polygonsets())

            self.cell_temp.flatten(single_layer=layer_idx)
            # add to cell 'Flattened'
            self.cell_flattened.add(self.cell_temp.get_polygonsets())
            # reset temp cell
            self.cell_temp.remove_polygons(lambda pts, layer, datatype: layer==layer_idx)
        self.lib.remove(cell='Temp')
        if export:
            print(f"\nOutputting cell 'Flattened' ({len(self.cell_flattened.get_polygonsets())} polygons)\n\t--> {self.outfile}\n\n")  

            self.lib.write_gds(self.outfile, cells=['Flattened'])
        else:
            print(f"\nCells {cell_keys} were flattened to a single cell called 'Flattened.\n\n")
        


    def make_array_onChip(self, columns, rows, 
                          spacing=("50um", "50um"), chip_dimension=("20mm", "20mm"), export=True, export_path=None, use_cells=None, add_labels=True, label_string=None, label_offset=None, label_size=50):
        '''
        Makes an array of a design input and exports it as a new GDS. Labels each arrayed structure. 

        Inputs:
            - columns - number of columns in the array
            - rows - number of rows in the array
            - spacing - (Defaults to ("50um", "50um")) Tuple containing x- and y-spacing of the array as strings with units
            - chip_dimension - (Defaults to ("20mm", "20mm")) Dimension of the chip on which to array the structure (WIP)
            - export - (Defaults to True) Choose whether to export the generated file (default), or just to operate on the design
            - export_path - (Optional) Export path (if different from self)
            - use_cells - (Optional) List of cells from the input GDS to include in the arrayed structure
            - add_label - (Defaults to True) Add a label beneath each arrayed structure
            - label_string - (Optional) If None, defaults to R1C1, R1C2 etc. (R{row_index}C{column_index}). Otherwise, if a single string is passed, this will be printed with an incrementing number. Otherwise, a list corresponding in length to the total number of arrayed structures can be passed and will be printed along row 1 from left to right, then row 2 etc.
            - label_offset - (Optional) Describes how far beneath the origin for each array point the label is. Defaults to half of spacing in y. Can be passed as a string with units (e.g. 300um).
            - label_size - (Defaults to 50) Text size for the printed labels
        
        Output:
            - cell_array - gdspy cell array containing arrayed structures
        '''

        # TODO: add option to autofill chip (requires additional arguments on chip_border (number as a string with units), fill_chip (boolean))

        # import assertions
        assert (isinstance(chip_dimension, tuple) and isinstance(chip_dimension[0], str)), r"Input argument 'chip_dimension' should be a tuple contining strings of the chip's (x, y) dimensions with units - e.g. 'chip_dimension=('20mm', '20mm')', as in qiskit and SQDMetal."   
        if isinstance(label_string, list):
            assert len(label_string)==(columns * rows), "If you are passing a list of labels, please ensure it is the same length as the total number of arrayed structures (i.e. rows * columns)." 

        # parse values
        sp_x = QUtilities.parse_value_length(spacing[0])
        sp_y = QUtilities.parse_value_length(spacing[1])
        if label_offset != None: 
            lbl_offset = QUtilities.parse_value_length(label_offset)
        else:
            lbl_offset = sp_y/2
        
        # set font
        # TODO: add as argument
        fp = FontProperties(family="serif", style="italic")

        # default: copy all cells. Else use function input (if any), finally use class init input (if any)
        if use_cells != None:
            self.flatten_cells(cell_keys=use_cells) 
        elif (self.import_cells == None) and (use_cells == None):
            self.flatten_cells(cell_keys=None) 
        else:
            self.flatten_cells(cell_keys=self.import_cells)

        # initialise an array cell
        self.cell_array = self.lib.new_cell('Array')

        # do arraying
        arrays_of_polygonsets = gdspy.CellArray(self.cell_flattened, 
                                columns=columns, 
                                rows=rows,
                                spacing=(QUtilities.parse_value_length(spacing[0])/self.gds_units,QUtilities.parse_value_length(spacing[1])/self.gds_units)
                                ).get_polygonsets()
        self.cell_array.add(arrays_of_polygonsets)

        # do labelling (origin is in bottom left)
        if add_labels:
            self.cell_text = self.lib.new_cell("Text")
            index = 0
            text = []
            if isinstance(label_string, list): 
                labels = label_string
            if add_labels:
                for i in range(rows):
                    for j in range(columns):
                        # assign labels
                        if label_string==None:
                            label = f'R{i}C{j}'
                        elif isinstance(label_string, str):
                            label = f'{label_string}_{i}{j}' 
                        else:
                            label = label_string[index]

                        # calculate label position
                        current_origin_x = j * sp_x
                        current_origin_y = (rows - i) * sp_y
                        label_position = (current_origin_x, (current_origin_y-lbl_offset))

                        print(f'{label} --> {label_position}')

                        # add text to gds on layer 999
                        # text = gdspy.PolygonSet(self.render_text(text=label, size=10, position=label_position, font_prop=fp), layer=999)
                        text = gdspy.Text(label, label_size, label_position, layer=11)

                        self.cell_text.add(text)

                        index += 1

        # export if option is set
        if export:
            if export_path == None: 
                array_export = self.outfile
            else:
                array_export = export_path
            
            if not add_labels:
                self.lib.write_gds(array_export, cells=['Array'])
                print(f"\nOutputting cell 'Array' ({len(self.cell_array.get_polygonsets())} polygons)\n\t--> {array_export}\n")
            if add_labels: 
                self.lib.write_gds(array_export, cells=['Array', 'Text'])
                print(f"\nOutputting cell: 'Array' ({len(self.cell_array.get_polygonsets())} polygons)\nOutputting cell: 'Text' ({len(self.cell_text.get_polygonsets())} polygons)\n\t--> {array_export}\n")

        # TODO: add centering by default (i.e. so that the whole array is centred on (0,0))

        print(f"\nFrom {self.infile}:\n\tA {columns} x {rows} array of the design was made\n\tSpacing: {spacing}\n\nThe array was exported to:\n\t{self.outfile}")

        return self.cell_array
    

    # for rendering system fonts
    def render_text(self, text, 
                    size=None, 
                    position=(0, 0), 
                    font_prop=None, 
                    tolerance=0.1):
        
        path = TextPath(position, text, size=size, prop=font_prop)
        polys = []
        xmax = position[0]
        for points, code in path.iter_segments():
            if code == path.MOVETO:
                c = gdspy.Curve(*points, tolerance=tolerance)
            elif code == path.LINETO:
                c.L(*points)
            elif code == path.CURVE3:
                c.Q(*points)
            elif code == path.CURVE4:
                c.C(*points)
            elif code == path.CLOSEPOLY:
                poly = c.get_points()
                if poly.size > 0:
                    if poly[:, 0].min() < xmax:
                        i = len(polys) - 1
                        while i >= 0:
                            if gdspy.inside(
                                poly[:1], [polys[i]], precision=0.1 * tolerance
                            )[0]:
                                p = polys.pop(i)
                                poly = gdspy.boolean(
                                    [p],
                                    [poly],
                                    "xor",
                                    precision=0.1 * tolerance,
                                    max_points=0,
                                ).polygons[0]
                                break
                            elif gdspy.inside(
                                polys[i][:1], [poly], precision=0.1 * tolerance
                            )[0]:
                                p = polys.pop(i)
                                poly = gdspy.boolean(
                                    [p],
                                    [poly],
                                    "xor",
                                    precision=0.1 * tolerance,
                                    max_points=0,
                                ).polygons[0]
                            i -= 1
                    xmax = max(xmax, poly[:, 0].max())
                    polys.append(poly)
        return polys