import gdspy
import os
from SQDMetal.Utilities.MakeGDS import MakeGDS
from SQDMetal.Utilities.QUtilities import QUtilities

class ManipulateGDS:
    def __init__(self, import_path_gds, export_path_gds, import_cells=None, origin=(0,0)):
        '''
        Class that will load a GDS file and perform actions such as making an array of the passed design (or a subsection of)

        Inputs:
            import_path_gds - path to GDS file for import
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
        


    def make_array_onChip(self, columns, rows, 
                          spacing=("0um", "0um"), chip_dimension=("20mm", "20mm"), export=True, export_path=None, use_cells=None):
        # import assertions
        assert (isinstance(chip_dimension, tuple) and isinstance(chip_dimension[0], str)), r"Input argument 'chip_dimension' should be a tuple contining strings of the chip's (x, y) dimensions with units - e.g. 'chip_dimension=('20mm', '20mm')', as in qiskit and SQDMetal."    
        
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

        # export if option is set
        if export:
            if export_path == None: 
                array_export = self.outfile
            else:
                array_export = export_path
            self.lib.write_gds(array_export, cells=['Array'])
            print(f"\nOutputting cell 'Array' ({len(self.cell_array.get_polygonsets())} polygons)\n\t--> {array_export}\n\n")

        # TODO: add centering by default (i.e. so that the whole array is centred on (0,0))

        return self.cell_array


f = ManipulateGDS(
                  import_path_gds="/Users/uqzdegna/Documents/Uni/PhD/qDesignDesk/junctions/manhattan_one_geo_v1.gds",
                  export_path_gds="/Users/uqzdegna/Documents/Uni/PhD/qDesignDesk/junctions/ManipulateGDS_test.gds",
                  origin=(0, 4000)
                  )

m = f.make_array_onChip(columns=16,
                        rows=16,
                        export=True, 
                        spacing=("500um", "500um"),
                        use_cells=['TOP_main_2', 'TOP_main_1'],
                        export_path="/Users/uqzdegna/Documents/Uni/PhD/qDesignDesk/junctions/ManipulateGDS_array.gds")

exit()