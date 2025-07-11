from qiskit_metal.qgeometries.qgeometries_handler import QGeometryTables
from qiskit_metal.renderers.renderer_mpl.mpl_renderer import QMplRenderer
from SQDMetal.Utilities.QiskitShapelyRenderer import QiskitShapelyRenderer
import numpy as np
import shapely
from SQDMetal.Utilities.ShapelyEx import ShapelyEx
from SQDMetal.Utilities.QUtilities import QUtilities
from SQDMetal.Utilities.GeometryProcessors.GeomBase import GeomBase


class GeomQiskitMetal(GeomBase):
    def __init__(self, design, **kwargs):
        self.design = design

        restrict_rect = kwargs.get('segment_rectangle', None)   #Given as [x1,y1,x2,y2]
        if isinstance(restrict_rect, list) or isinstance(restrict_rect, np.ndarray) or isinstance(restrict_rect, tuple):
            self._restrict_rect = [min([restrict_rect[0], restrict_rect[2]]), min([restrict_rect[1], restrict_rect[3]]),
                                  max([restrict_rect[0], restrict_rect[2]]), max([restrict_rect[1], restrict_rect[3]])]
            self._chip_size_x = self._restrict_rect[2]-self._restrict_rect[0]
            self._chip_size_y = self._restrict_rect[3]-self._restrict_rect[1]
            self._chip_centre = [(self._restrict_rect[0]+self._restrict_rect[2])*0.5,
                                (self._restrict_rect[1]+self._restrict_rect[3])*0.5,
                                QUtilities.parse_value_length(self.design.chips['main']['size']['center_z'])]
        else:
            self._chip_size_x = QUtilities.parse_value_length(self.design.chips['main']['size']['size_x'])
            self._chip_size_y = QUtilities.parse_value_length(self.design.chips['main']['size']['size_y'])
            self._chip_centre = [QUtilities.parse_value_length(self.design.chips['main']['size'][x]) for x in ['center_x', 'center_y', 'center_z']]
            self._restrict_rect = [self.chip_centre[0]-0.5*self._chip_size_x, self.chip_centre[1]-0.5*self._chip_size_y,
                                  self.chip_centre[0]+0.5*self._chip_size_x, self.chip_centre[1]+0.5*self._chip_size_y]
        self._chip_size_z = QUtilities.parse_value_length(self.design.chips['main']['size']['size_z'])

    @property
    def ParserType(self):
        return "Qiskit-Metal"

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
        unit_conv = kwargs.get('unit_conv', 1)
        fuse_threshold = kwargs.get('fuse_threshold', 1e-12)
        threshold = kwargs.get('threshold', 1e-9)
        fillet_resolution = kwargs.get('fillet_resolution', 12)

        qm_units = QUtilities.get_units(self.design)

        metal_polys = []
        #Handle the ground plane...
        if not ground_plane['omit']:
            qmpl = QiskitShapelyRenderer(None, self.design, None)
            gsdf = qmpl.get_net_coordinates(fillet_resolution)
            filt = gsdf.loc[gsdf['subtract'] == True]
            #TODO: It should only take stuff from the metallic layers specified when calculating spaces???
            space_polys = ShapelyEx.fuse_polygons_threshold(filt['geometry'].buffer(0), fuse_threshold / qm_units)
            space_polys = ShapelyEx.shapely_to_list(space_polys)
            space_polys = [shapely.affinity.scale(x, xfact=qm_units, yfact=qm_units, origin=(0,0)) for x in space_polys]
            #
            metal_surface = shapely.geometry.box(self._chip_centre[0] - 0.5*self._chip_size_x, self._chip_centre[1] - 0.5*self._chip_size_y,
                                                self._chip_centre[0] + 0.5*self._chip_size_x, self._chip_centre[1] + 0.5*self._chip_size_y)
            ground_plane_poly = shapely.difference(metal_surface, shapely.geometry.multipolygon.MultiPolygon(space_polys))
            if not ground_plane_poly.is_empty:
                metal_polys.append(ground_plane_poly)
        #Gather all polygons into contiguous groups
        for cur_layer in metallic_layers:
            if cur_layer['type'] == 'design_layer':
                cur_layer['resolution'] = fillet_resolution
                cur_layer['threshold'] = threshold
                # cur_layer['fuse_threshold'] = fuse_threshold  #Let the individual fuse_thresholds work here...
                cur_layer['restrict_rect'] = self._restrict_rect[:]
                metal_polys_all, metal_sel_ids = QUtilities.get_metals_in_layer(self.design, **cur_layer)
                metal_polys += metal_polys_all
            elif cur_layer['type'] == 'Uclip':
                if cur_layer['clip_type'] == 'inplaneLauncher':
                    Uclip = QUtilities.get_RFport_CPW_groundU_Launcher_inplane(self.design, cur_layer['qObjName'], cur_layer['thickness_side'], cur_layer['thickness_back'], cur_layer['separation_gap'])
                elif cur_layer['clip_type'] == 'inplaneRoute':
                    Uclip = QUtilities.get_RFport_CPW_groundU_Route_inplane(self.design, cur_layer['route_name'], cur_layer['pin_name'], cur_layer['thickness_side'], cur_layer['thickness_back'], cur_layer['separation_gap'])
                metal_polys.append(shapely.Polygon(Uclip))
        #Try to fuse any contiguous polygons...
        metal_polys = ShapelyEx.fuse_polygons_threshold(metal_polys, fuse_threshold)
        metal_polys = shapely.affinity.scale(metal_polys, xfact=1/unit_conv, yfact=1/unit_conv, origin=(0,0))
        simplify_edge_min_angle_deg = kwargs.get("simplify_edge_min_angle_deg",-1)
        if simplify_edge_min_angle_deg >= 0:
            metal_polys = ShapelyEx.simplify_poly_edges(metal_polys, simplify_edge_min_angle_deg)

        #If there are any MultiPolygons, convert them into normal polygons...
        new_polys = ShapelyEx.shapely_to_list(metal_polys)

        substrate = ShapelyEx.rectangle((self._chip_centre[0]-self._chip_size_x/2)/unit_conv,
                                        (self._chip_centre[1]-self._chip_size_y/2)/unit_conv,
                                        (self._chip_centre[0]+self._chip_size_x/2)/unit_conv,
                                        (self._chip_centre[1]+self._chip_size_y/2)/unit_conv)
        dielectric_gaps = shapely.difference(substrate, metal_polys)
        dielectric_gaps = ShapelyEx.shapely_to_list(dielectric_gaps)

        return new_polys, dielectric_gaps
