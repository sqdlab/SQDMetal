import gdstk
import numpy as np
import shapely
from SQDMetal.Utilities.ShapelyEx import ShapelyEx
from SQDMetal.Utilities.QUtilities import QUtilities
from SQDMetal.Utilities.GeometryProcessors.GeomBase import GeomBase


class GeomGDS(GeomBase):
    def __init__(self, gds_file, **kwargs):
        self.gds_file = gds_file
        cell_ind = kwargs.get('gds_cell_index', 0)
        chip_thickness = kwargs.get('gds_chip_thickness', 500e-6)

        restrict_rect = kwargs.get('segment_rectangle', None)   #Given as [x1,y1,x2,y2]
        if isinstance(restrict_rect, list) or isinstance(restrict_rect, np.ndarray) or isinstance(restrict_rect, tuple):
            self._restrict_rect = [min([restrict_rect[0], restrict_rect[2]]), min([restrict_rect[1], restrict_rect[3]]),
                                   max([restrict_rect[0], restrict_rect[2]]), max([restrict_rect[1], restrict_rect[3]])]

            self._chip_size_x = self._restrict_rect[2]-self._restrict_rect[0]
            self._chip_size_y = self._restrict_rect[3]-self._restrict_rect[1]
            self._chip_centre = [(self._restrict_rect[0]+self._restrict_rect[2])*0.5,
                                 (self._restrict_rect[1]+self._restrict_rect[3])*0.5,
                                 0]
        else:
            self._restrict_rect = None
        self._chip_size_z = chip_thickness

    @property
    def ParserType(self):
        return "GDS"

    @property
    def chip_size_x(self):
        return self._chip_size_x

    @property
    def chip_size_y(self):
        return self._chip_size_y

    @property
    def chip_size_z(self):
        return self._chip_size_z

    @property
    def chip_centre(self):
        return self._chip_centre[:]

    def process_layers(self, metallic_layers, ground_plane, **kwargs):
        if not ground_plane['omit']:
            print("WARNING: Don't use add_ground_plane() when using GDS. It is being ignored here...")

        unit_conv = kwargs.get('unit_conv', 1)

        leGds = gdstk.read_gds(infile=self.gds_file)
        leUnits = leGds.unit / unit_conv
        metal_polys_all = []
        for cur_layer in metallic_layers:
            if cur_layer['type'] == 'design_layer':
                cell_ind = cur_layer['other'].get('cell_index', 0)
                leCell = leGds.cells[cell_ind]

                metal_polys = leCell.get_polygons(layer=cur_layer['layer_id'], datatype=0)
                metal_polys = [shapely.Polygon(cur_poly.points*leUnits) for cur_poly in metal_polys]
                if isinstance(self._restrict_rect, list):
                    clip_rect = [x/unit_conv for x in self._restrict_rect]
                    metal_polys = [shapely.clip_by_rect(cur_poly, *clip_rect) for cur_poly in metal_polys]
                metal_polys = ShapelyEx.fuse_polygons_threshold(metal_polys, cur_layer['fuse_threshold'])
                if cur_layer['threshold'] > 0:
                    metal_polys = metal_polys.simplify(cur_layer['threshold'])

                metal_polys_all += ShapelyEx.shapely_to_list(metal_polys)
            elif cur_layer['type'] == 'Uclip':
                pass
                # if cur_layer['clip_type'] == 'inplaneLauncher':
                #     Uclip = QUtilities.get_RFport_CPW_groundU_Launcher_inplane(self.design, cur_layer['qObjName'], cur_layer['thickness_side'], cur_layer['thickness_back'], cur_layer['separation_gap'])
                # elif cur_layer['clip_type'] == 'inplaneRoute':
                #     Uclip = QUtilities.get_RFport_CPW_groundU_Route_inplane(self.design, cur_layer['route_name'], cur_layer['pin_name'], cur_layer['thickness_side'], cur_layer['thickness_back'], cur_layer['separation_gap'])
                # metal_polys.append(shapely.Polygon(Uclip))

        #Given as min x, min y, max x, max y
        temp_metal_fuse = ShapelyEx.fuse_polygons_threshold(metal_polys_all, kwargs.get('fuse_threshold', 1e-12) / unit_conv)

        if not isinstance(self._restrict_rect, list):
            substrate_box_bnds = temp_metal_fuse.bounds
            self._chip_size_x = (substrate_box_bnds[2] - substrate_box_bnds[0])*unit_conv
            self._chip_size_y = (substrate_box_bnds[3] - substrate_box_bnds[1])*unit_conv
            self._chip_centre = [(substrate_box_bnds[2]+substrate_box_bnds[0])/2*unit_conv, (substrate_box_bnds[3]+substrate_box_bnds[1])/2*unit_conv,0]
        else:
            substrate_box_bnds = self._restrict_rect[:]
        substrate = ShapelyEx.rectangle(*substrate_box_bnds)

        dielectric_gaps = shapely.difference(substrate, temp_metal_fuse)
        dielectric_gaps = ShapelyEx.shapely_to_list(dielectric_gaps)

        return ShapelyEx.shapely_to_list(temp_metal_fuse), dielectric_gaps
