# Copyright 2025 Prasanna Pakkiam
# SPDX-License-Identifier: Apache-2.0

from qiskit_metal.toolbox_metal.parsing import parse_value
from shapely.geometry import Polygon
import numpy as np
import shapely
import geopandas as gpd

import matplotlib.pyplot as plt

from SQDMetal.Utilities.QiskitShapelyRenderer import QiskitShapelyRenderer
from SQDMetal.Utilities.ShapelyEx import ShapelyEx

class PVD_Shadows:
    def __init__(self, design, chip_name='main'):
        self.design = design
        self._chip_name = chip_name
        self.update_pvd_profiles()

    def update_pvd_profiles(self):
        unit_conv, unit_conv_name = self.get_units()

        #Calculate chip-bounds
        cur_dict = parse_value(self.design.chips[self._chip_name].size, self.design.chips[self._chip_name].size.keys())
        self._chip_bounds = [((cur_dict['center_x']-cur_dict['size_x']/2)*unit_conv, (cur_dict['center_y']-cur_dict['size_y']/2)*unit_conv),
                             ((cur_dict['center_x']+cur_dict['size_x']/2)*unit_conv, (cur_dict['center_y']-cur_dict['size_y']/2)*unit_conv),
                             ((cur_dict['center_x']+cur_dict['size_x']/2)*unit_conv, (cur_dict['center_y']+cur_dict['size_y']/2)*unit_conv),
                             ((cur_dict['center_x']-cur_dict['size_x']/2)*unit_conv, (cur_dict['center_y']+cur_dict['size_y']/2)*unit_conv)]
        self._chip_bounds = shapely.Polygon(self._chip_bounds)

        #Setup dictionary indexed as layer index and then a list of the evaporation steps of which each step is a dictionary of required parameters...
        self.evap_profiles = {}
        if 'evaporations' in self.design.chips[self._chip_name]:
            for cur_layer in self.design.chips[self._chip_name].evaporations.keys():
                cur_dict = parse_value(self.design.chips[self._chip_name].evaporations[cur_layer], self.design.chips[self._chip_name].evaporations[cur_layer].keys())
                cur_layer_ind = int(cur_layer.lower().split("layer")[-1])
                pvds = sorted([x for x in cur_dict if x.startswith('pvd')])
                pvds = [cur_dict[x] for x in pvds]
                self.evap_profiles[cur_layer_ind] = self._get_cur_pvd_profile(cur_dict['bottom_layer'] * unit_conv,
                                                        cur_dict['top_layer'] * unit_conv,
                                                        cur_dict['undercut'] * unit_conv,
                                                        pvds, unit_conv)

    def _get_cur_pvd_profile(self, thickness_L, thickness_U, undercut, pvd_infos, unit_conv):
        ret_list = []
        for cur_pvd_info in pvd_infos:
            angle_phi = cur_pvd_info['angle_phi']/180*np.pi
            angle_theta = cur_pvd_info['angle_theta']/180*np.pi
            metal_thickness = cur_pvd_info.get('metal_thickness',0)*unit_conv
            #
            cur_profile = {}
            cur_profile['leading_dist'] = (thickness_L+thickness_U)*np.tan(angle_theta)
            cur_profile['tailing_dist'] = (thickness_L)*np.tan(angle_theta)
            cur_profile['tan_theta'] = np.tan(angle_theta)
            cur_profile['cur_evap_dir'] = -np.array([np.cos(angle_phi), np.sin(angle_phi)])
            cur_profile['metal_thickness'] = metal_thickness
            cur_profile['undercut'] = undercut
            ret_list.append(cur_profile)
        return ret_list

    def get_units(self):
        unit_conv = self.design.get_units()
        unit_conv_name = unit_conv
        if unit_conv == 'mm':
            unit_conv = 1e-3
        elif unit_conv == 'um':
            unit_conv = 1e-6
        elif unit_conv == 'nm':
            unit_conv = 1e-9
        else:
            assert False, f"Unrecognised units: {unit_conv}"
        return unit_conv, unit_conv_name

    def plot_all_layers(self, mode='separate', plot_overlap = False ,**kwargs):
        qmpl = QiskitShapelyRenderer(None, self.design, None)
        gsdf = qmpl.get_net_coordinates(kwargs.get('resolution',4))
        unit_conv, unit_conv_name = self.get_units()
        
        plot_mask = kwargs.get('plot_mask', True)

        layer_ids = np.unique(gsdf['layer'])
        plot_polys = []
        names = []
        new_metal_layers = [(0, shapely.Polygon(self._chip_bounds))]
        for layer_id in layer_ids:
            filt = gsdf.loc[(gsdf['layer'] == layer_id) & (~gsdf['subtract'])]
            if filt.shape[0] == 0:
                continue    #Shouldn't trigger...
            metal_polys = shapely.unary_union(filt['geometry'])
            metal_polys = shapely.affinity.scale(metal_polys, xfact=unit_conv, yfact=unit_conv, origin=(0,0))
            
            if layer_id in self.evap_profiles:
                if plot_mask:
                    plot_polys += [metal_polys]
                    names += [f"Layer {layer_id} Mask"]
                leShadows = self.get_all_shadows(metal_polys, layer_id, mode, list_evap_metals=new_metal_layers, **kwargs)
                plot_polys += leShadows
                names += [f"Layer {layer_id}, Step {x}" for x in range(len(leShadows))]
            else:
                plot_polys += [metal_polys]
                names += [f"Layer {layer_id}"]
        #split layers or steps and create geodataframe
        layers = [name.split(",")[0] for name in names]
        steps = [name.split(",")[1].strip() if "," in name else "Mask" for name in names]

        gdf = gpd.GeoDataFrame({'names':names, 'layers': layers, 'steps': steps}, geometry=plot_polys)
        gdf_filt = gdf[(gdf['steps'] != 'Mask') & (~gdf['layers'].str.contains('Mask'))] #ignore mask region to find the right overlap regions

        
        fig, ax = plt.subplots(1)
        unique_layers = gdf_filt['layers'].unique()
        unique_steps = gdf_filt['steps'].unique()

        assert len(unique_layers) == 2 or (len(unique_layers) == 1 and len(unique_steps) == 2), \
        "Error: The input data must contain either two unique layers or one layer with two unique steps."

        if plot_overlap:
            #separate the GeoDataFrame by layer or step
            if len(unique_layers) == 2: 
                layer1 = gdf[gdf['layers'] == unique_layers[0]]
                layer2 = gdf[gdf['layers'] == unique_layers[1]]
                overlap = gpd.overlay(layer1, layer2, how='intersection') #find the overlapping regions
            else:
                step1 = gdf[gdf['steps'] == unique_steps[0]]
                step2 = gdf[gdf['steps'] == unique_steps[1]]
                overlap = gpd.overlay(step1, step2, how='intersection')

            gdf.set_crs("EPSG:3857", inplace=True)
            total_overlap_area = overlap.geometry.area.sum()*1e12 #calculate the sum of area of the overlapping regions
            print(f"Total area of overlapping regions: {total_overlap_area:.6f} µm²")
            overlap.plot(ax=ax, facecolor='none', edgecolor='black', hatch='//', alpha=1, linestyle='--',linewidth=1.2, label='Overlap')

        gdf.boundary.plot(ax=ax, color='black', linestyle='--', linewidth=1, alpha=0.25)
        gdf.plot(ax = ax, column='layers', cmap='jet', alpha=0.35, categorical=True, legend=True)
        ax.set_xlabel('Position (m)')
        ax.set_ylabel('Position (m)')

    def plot_layer(self, layer_id, mode='separate', **kwargs):
        qmpl = QiskitShapelyRenderer(None, self.design, None)
        gsdf = qmpl.get_net_coordinates(kwargs.get('resolution',4))

        filt = gsdf.loc[(gsdf['layer'] == layer_id) & (~gsdf['subtract'])]
        if filt.shape[0] == 0:
            return
        metal_polys = shapely.unary_union(filt['geometry'])
        unit_conv, unit_conv_name = self.get_units()
        metal_polys = shapely.affinity.scale(metal_polys, xfact=unit_conv, yfact=unit_conv, origin=(0,0))

        plot_mask = kwargs.get('plot_mask', True)

        if layer_id in self.evap_profiles:
            plot_polys = []
            names = []
            if plot_mask:
                plot_polys += [metal_polys]
                names += [f"Layer {layer_id} Mask"]
            shad_polys = self.get_all_shadows(metal_polys, layer_id, mode, **kwargs)
            plot_polys += shad_polys
            names += [f"Layer {layer_id}, Step {x}" for x in range(len(shad_polys))]
        else:
            plot_polys = [metal_polys]
            names = [f"Layer {layer_id}"]
        
        gdf = gpd.GeoDataFrame({'names':names}, geometry=plot_polys)
        fig, ax = plt.subplots(1)
        gdf.plot(ax = ax, column='names', cmap='jet', alpha=0.2, categorical=True, legend=True)
        ax.set_xlabel(f'Position ({unit_conv_name})')
        ax.set_ylabel(f'Position ({unit_conv_name})')

    def get_shadows_for_component(self, qObj_name, mode='separate', **kwargs):
        comp_id = self.design.components[qObj_name].id

        qmpl = QiskitShapelyRenderer(None, self.design, None)
        gsdf = qmpl.get_net_coordinates(kwargs.get('resolution',4))

        filt = gsdf.loc[(gsdf['component'] == comp_id) & (~gsdf['subtract'])]
        if filt.shape[0] == 0:
            return
        metal_polys = shapely.unary_union(filt['geometry'])
        unit_conv, unit_conv_name = self.get_units()
        metal_polys = shapely.affinity.scale(metal_polys, xfact=unit_conv, yfact=unit_conv, origin=(0,0))

        plot_mask = kwargs.get('plot_mask', True)

        layer_id = int(self.design.components[qObj_name].options.layer)

        if layer_id in self.evap_profiles:
            plot_polys = []
            names = []
            if plot_mask:
                plot_polys += [metal_polys]
                names += [f"Layer {layer_id} Mask"]
            shad_polys = self.get_all_shadows(metal_polys, layer_id, mode, **kwargs)
            plot_polys += shad_polys
            names += [f"Layer {layer_id}, Step {x}" for x in range(len(shad_polys))]
        else:
            plot_polys = [metal_polys]
            names = [f"Layer {layer_id}"]
        
        return plot_polys

    def get_shadow_largest_interior_for_component(self, qObj_name, mode='separate', **kwargs):
        polys = self.get_shadows_for_component(qObj_name, mode, **kwargs)
        polys = ShapelyEx.fuse_polygons_threshold(polys, kwargs.get('threshold', 1e-9))
        if isinstance(polys, shapely.geometry.multipolygon.MultiPolygon):
            polys = [x for x in polys.geoms]
        else:
            polys = [polys]
        #
        cur_area = 0
        cur_max_int = None
        for poly in polys:
            for cur_int_line in poly.interiors:
                cur_int = shapely.Polygon(cur_int_line)
                if cur_int.area > cur_area:
                    cur_area = cur_int.area
                    cur_max_int = cur_int
        return np.array(cur_max_int.exterior.coords[:]), cur_max_int

    def get_all_shadows(self, cur_poly, layer_id, mode='merge', **kwargs):
        list_evap_metals = kwargs.get('list_evap_metals', [(0, shapely.Polygon(self._chip_bounds))])

        if layer_id not in self.evap_profiles:
            return cur_poly
        profiles = []
        for dict_evap_params in self.evap_profiles[layer_id]:
            profiles += [self._get_poly_shadow(cur_poly, dict_evap_params, list_evap_metals)]
            list_evap_metals = self._add_metallic_layer(list_evap_metals, profiles[-1], dict_evap_params['metal_thickness'])
        if mode == 'merge':
            return [shapely.unary_union(profiles)]
        elif mode == 'separate':
            return profiles
        elif mode == 'separate_delete_below':
            assert 'layer_trim_length' in kwargs, "Must supply layer_trim_length if using mode: separate_delete_below."
            trim_len = kwargs.get('layer_trim_length')
            for m in range(0,len(profiles)-1):
                profiles[m] = profiles[m].difference(profiles[m+1].buffer(trim_len, join_style=2, cap_style=3))
            return profiles
        elif mode == 'separate_delete_intersections':
            assert 'layer_trim_length' in kwargs, "Must supply layer_trim_length if using mode: separate_delete_below."
            trim_len = kwargs.get('layer_trim_length')
            for m in range(0,len(profiles)):
                for n in range(m+1,len(profiles)):
                    intersec = shapely.intersection(profiles[m], profiles[n]).buffer(trim_len, join_style=2, cap_style=3)
                    profiles[m] = profiles[m].difference(intersec)
                    profiles[n] = profiles[n].difference(intersec)
            return profiles
        elif mode == 'separate_plot_intersections':
            assert 'layer_trim_length' in kwargs, "Must supply layer_trim_length if using mode: separate_delete_below."
            trim_len = kwargs.get('layer_trim_length')
            for m in range(0,len(profiles)):
                for n in range(m+1,len(profiles)):
                    intersec = shapely.intersection(profiles[m], profiles[n]).buffer(trim_len, join_style=2, cap_style=3)

            return [intersec]

    def _get_poly_shadow(self, cur_poly, dict_evap_params, list_evap_metals):
        leading_dist = dict_evap_params['leading_dist']
        tailing_dist = dict_evap_params['tailing_dist']
        if leading_dist == 0:   #Basically a vertical evaporation...
            return cur_poly

        extSegs, intSegs = self.get_line_segments(cur_poly)

        undercut = dict_evap_params['undercut']
        allowed_region = self.get_allowed_region(cur_poly, undercut) #Mitre joint for sharp corners...
        
        cur_evap_dir = dict_evap_params['cur_evap_dir']
        bounds = list(allowed_region.bounds)
        #Pad the bounds a bit...
        pad_x = 0.1*(bounds[2] - bounds[0])
        pad_y = 0.1*(bounds[3] - bounds[1])
        bounds[0] -= pad_x
        bounds[1] -= pad_y
        bounds[2] += pad_x
        bounds[3] += pad_y
        allEdges = extSegs + intSegs
        light_profiles = []
        for cur_edge in allEdges:
            line_dir = np.diff(cur_edge, axis=0)[0]
            leading = (line_dir[0]*cur_evap_dir[1] - line_dir[1]*cur_evap_dir[0] > 0)
            if leading:
                shadow = leading_dist*cur_evap_dir
            else:
                shadow = tailing_dist*cur_evap_dir
            p0, p1 = np.array(cur_edge[0])+shadow, np.array(cur_edge[1])+shadow
            poly_light = self.create_light_profile(p0, p1, bounds, cur_evap_dir)
            if isinstance(poly_light, Polygon):
                if leading:
                    light_profiles += [(poly_light,1)]
                else:
                    light_profiles += [(poly_light,-1)]

        light_outlines = gpd.geoseries.GeoSeries([x[0] for x in light_profiles]).exterior.unary_union

        # gdf = gpd.GeoDataFrame({'names':[x[1] for x in light_profiles]}, geometry=gpd.geoseries.GeoSeries([x[0] for x in light_profiles]))
        # fig, ax = plt.subplots(1)
        # gdf.plot(ax = ax, column='names', cmap='jet', alpha=0.2, categorical=True, legend=True)
        # ax.set_xlabel(f'Position ({unit_conv_name})')
        # ax.set_ylabel(f'Position ({unit_conv_name})')

        new_polys = list(shapely.ops.polygonize(light_outlines))
        #Choose a fraction the minimum area as a threshold as the intersection command can glitch a bit with floating-point precision...
        min_area_thresh = np.min([y.area for y in new_polys])*0.75
        region_sums = [(new_polys[y], sum([x[1] for x in light_profiles if shapely.intersection(x[0],new_polys[y]).area > min_area_thresh])) for y in range(len(new_polys))]
        #Select regions that are lit
        total_light = [x[0] for x in region_sums if x[1] > 0]

        final_calc_metals = self.process_geometry(shapely.intersection(shapely.geometry.multipolygon.MultiPolygon(total_light), allowed_region))

        #Account for thickness of previously evaporated layers...
        if len(list_evap_metals) == 0:
            return final_calc_metals    #Surface is uniform in height...       
        allowed_polys = [self.process_geometry(x[1]) for x in list_evap_metals]

        # gdf = gpd.GeoDataFrame({'names':[x[0] for x in list_evap_metals]}, geometry=gpd.geoseries.GeoSeries([x[1] for x in list_evap_metals]))
        # fig, ax = plt.subplots(1)
        # gdf.plot(ax = ax, column='names', cmap='jet', alpha=0.2, categorical=True, legend=True)
        # ax.set_xlabel(f'Position ({unit_conv_name})')
        # ax.set_ylabel(f'Position ({unit_conv_name})')

        for m in range(len(list_evap_metals)-1, 0,-1):
            extSegs, intSegs = self.get_line_segments(list_evap_metals[m][1])
            lineSegs = extSegs + intSegs
            for n in range(m-1,-1,-1):
                height_diff = list_evap_metals[m][0] - list_evap_metals[n][0]
                dist_offset = height_diff * dict_evap_params['tan_theta']
                for cur_line_seg in lineSegs:
                    allowed_polys[n] = self._subtract_shadows(allowed_polys[n], cur_line_seg, cur_evap_dir, dist_offset)

        # gdf = gpd.GeoDataFrame({'names':[x[0] for x in list_evap_metals]}, geometry=gpd.geoseries.GeoSeries([x for x in allowed_polys]))
        # fig, ax = plt.subplots(1)
        # gdf.plot(ax = ax, column='names', cmap='jet', alpha=0.2, categorical=True, legend=True)
        
        allowed_polys = shapely.unary_union(allowed_polys)
        self.allowed_polys = final_calc_metals
        return self.process_geometry(shapely.intersection(allowed_polys, final_calc_metals))

    def _subtract_shadows(self, poly, line_seg, vec_light, dist_offset):
        vx, vy = vec_light
        p0, p1 = line_seg

        #Check light is not parallel with edge...
        vec_edge = (p1[0]-p0[0], p1[1]-p0[1])
        if np.abs(vx*vec_edge[1]-vy*vec_edge[0]) < 1e-9:
            return poly
        
        vx *= dist_offset
        vy *= dist_offset
        
        cut_poly = shapely.Polygon([(p0[0],p0[1]),(p0[0]+vx,p0[1]+vy),(p1[0]+vx,p1[1]+vy),(p1[0],p1[1])])
        culled_poly = shapely.difference(poly, cut_poly)

        # gdf = gpd.GeoDataFrame({'names':['culledPoly']}, geometry=gpd.geoseries.GeoSeries([culled_poly]))
        # fig, ax = plt.subplots(1)
        # gdf.plot(ax = ax, column='names', cmap='jet', alpha=0.2, categorical=True, legend=True)

        return culled_poly

    def process_geometry(self, geom):
        if geom.is_empty:
            return Polygon([])
        #Sometimes intersections return lines/points etc...
        if isinstance(geom, shapely.geometry.collection.GeometryCollection):
            ret_poly = shapely.unary_union([y for y in geom.geoms if (isinstance(y, Polygon) or isinstance(y, shapely.geometry.multipolygon.MultiPolygon) and y.area > 1e-25)])
        elif isinstance(geom, shapely.geometry.multipolygon.MultiPolygon):
            ret_poly = shapely.geometry.multipolygon.MultiPolygon([y for y in geom.geoms if y.area > 1e-25])
        elif isinstance(geom, Polygon) and geom.area > 1e-25:
            ret_poly = geom
        else:
            return Polygon([])
        #Sometimes intersection connects two polygons with an infinitely thin line...
        #https://gis.stackexchange.com/questions/280488/removing-thin-rectangles-from-a-shapely-polygon
        #Basically this trick attempts to remove such "features"...
        #Probably due to floating-point errors?        
        d = 1e-13 # distance
        return ret_poly.buffer(-d).buffer(d).intersection(ret_poly).simplify(d)

    def get_allowed_region(self, cur_poly, undercut):
        return cur_poly.buffer(undercut, join_style=2, cap_style=3)

    def get_line_segments(self, poly):
        #Fun one to try on a unit-test...
        # multiPoly = shapely.wkt.loads('''
        #             MULTIPOLYGON
        #             (((40 40, 20 45, 45 30, 40 40)),
        #             ((20 35, 10 30, 10 10, 30 5, 45 20, 20 35), (30 20, 20 15, 20 25, 30 20)))
        #             ''')
        extSegs = []
        intSegs = []
        if isinstance(poly, shapely.geometry.multipolygon.MultiPolygon):
            for x in poly.geoms:
                cur_poly = shapely.geometry.polygon.orient(x, 1)
                cur_ext_coords = cur_poly.exterior.coords[:]
                extSegs += [[cur_ext_coords[m-1], cur_ext_coords[m]] for m in range(1,len(cur_ext_coords))]
                cur_interiors = cur_poly.interiors
                if len(cur_interiors) == 0:
                    continue
                for y in cur_interiors:
                    cur_int_coords = y.coords[:]
                    intSegs += [[cur_int_coords[m-1], cur_int_coords[m]] for m in range(1,len(cur_int_coords))]
        else:
            poly = shapely.geometry.polygon.orient(poly, 1)
            cur_ext_coords = poly.exterior.coords[:]
            extSegs += [[cur_ext_coords[m-1], cur_ext_coords[m]] for m in range(1,len(cur_ext_coords))]
            cur_interiors = poly.interiors
            if len(cur_interiors) > 0:
                for y in cur_interiors:
                    cur_int_coords = y.coords[:]
                    intSegs += [[cur_int_coords[m-1], cur_int_coords[m]] for m in range(1,len(cur_int_coords))]
        return extSegs, intSegs

    def create_light_profile(self, p0, p1, bounds, vec_light):
        #bounds given as minx, miny, maxx, maxy
        #First figure out where p0 and p1 intersect onto the bounding box via ray vec_light
        xmin, ymin, xmax, ymax = bounds
        vx, vy = vec_light
        
        #Check light is not parallel with edge...
        vec_edge = (p1[0]-p0[0], p1[1]-p0[1])
        if np.abs(vx*vec_edge[1]-vy*vec_edge[0]) < 1e-9:
            return None
        
        #See if it's a trivial case:
        if np.abs(vx) < 1e-12:
            if vy > 0:
                return Polygon([p0, p1, [p1[0], ymax], [p0[0], ymax]])
            else:
                return Polygon([p0, p1, [p1[0], ymin], [p0[0], ymin]])
        if np.abs(vy) < 1e-12:
            if vx > 0:
                return Polygon([p0, p1, [xmax, p1[1]], [xmax, p0[1]]])
            else:
                return Polygon([p0, p1, [xmin, p1[1]], [xmin, p0[1]]])
        #Otherwise construct the polygon...
        intersecs = []
        intersec_inds = [] #0, 1, 2, 3 being left, bottom, right and top
        for p in [p0, p1]:
            #Check left side
            t = (xmin-p[0])/vx
            if t >= 0:
                y = p[1]+t*vy
                if y >= ymin and y <= ymax:
                    intersecs += [(xmin, y)]
                    intersec_inds += [0]
                    continue
            #Check bottom side
            t = (ymin-p[1])/vy
            if t >= 0:
                x = p[0]+t*vx
                if x >= xmin and x <= xmax:
                    intersecs += [(x, ymin)]
                    intersec_inds += [1]
                    continue
            #Check right side
            t = (xmax-p[0])/vx
            if t >= 0:
                y = p[1]+t*vy
                if y >= ymin and y <= ymax:
                    intersecs += [(xmax, y)]
                    intersec_inds += [2]
                    continue
            #It's the top side then...
            t = (ymax-p[1])/vy
            if t >= 0:
                x = p[0]+t*vx
                intersecs += [(x, ymax)]
                intersec_inds += [3]
            #Stuff it - floating-point issues; it's probably in one of the corners...
            corners = [[xmin,ymin],[xmax,ymin],[xmax,ymax],[xmin,ymax]]
            for cur_corn in corners:
                cand_vec = [cur_corn[0]-p[0],cur_corn[1]-p[1]]
                if np.abs(vx*cand_vec[1]-vy*cand_vec[0]) < 1e-9:
                    intersecs += [(cur_corn[0], cur_corn[1])]
                    intersec_inds += [-1]
                    break            
        if len(intersecs) <= 1:
            # print(p0, p1)
            # print(intersecs)
            return None
        intersecs = [intersecs[0], p0, p1, intersecs[1]]
        if intersec_inds[0] != intersec_inds[1]:
            intersec_inds = sorted(intersec_inds)
            if intersec_inds[0] == 0 and intersec_inds[1] == 1:
                intersecs.append((xmin, ymin))
            elif intersec_inds[0] == 1 and intersec_inds[1] == 2:
                intersecs.append((xmax, ymin))
            elif intersec_inds[0] == 2 and intersec_inds[1] == 3:
                intersecs.append((xmax, ymax))
            elif intersec_inds[0] == 0 and intersec_inds[1] == 3:
                intersecs.append((xmin, ymax))
        return Polygon(np.array(intersecs))

    def get_junction_area(self,qObj_name=None, layer_id=None, **kwargs):
        '''prints and returns the expected junction area in um^2

        The area is equal to the intersection area between the evaporations, i.e. the junction overlap.
        Must provide either qObj_name (find out by design.components) or layer_id

        IMPORTANT: make sure there are no other areas in the same object or layer that give a shadow
        (e.g. when using tapered electrodes), otherwise these will be included in the area.

        Does not include extra area from overlap between sides of first and second layer
        (relevant for small junctions)

        To do: add layer thickness contributions
        Added by: SZ, 11/9/2024
        '''
        if qObj_name is not None:

            comp_id = self.design.components[qObj_name].id

            qmpl = QiskitShapelyRenderer(None, self.design, None)
            gsdf = qmpl.get_net_coordinates(kwargs.get('resolution',4))

            filt = gsdf.loc[(gsdf['component'] == comp_id) & (~gsdf['subtract'])]
            if filt.shape[0] == 0:
                return
            metal_polys = shapely.unary_union(filt['geometry'])
            unit_conv, unit_conv_name = self.get_units()
            metal_polys = shapely.affinity.scale(metal_polys, xfact=unit_conv, yfact=unit_conv, origin=(0,0))

            layer_id = int(self.design.components[qObj_name].options.layer)
            if layer_id in self.evap_profiles:

                intersec = self.get_all_shadows(metal_polys, layer_id, 'separate_plot_intersections', **kwargs)
            else:
                print('Component does not have evap profile')
            print(f'Expected junction area for {qObj_name}: {intersec[0].area*1e12} um^2')
            return intersec[0].area*1e12

        elif layer_id is not None:

            qmpl = QiskitShapelyRenderer(None, self.design, None)
            gsdf = qmpl.get_net_coordinates(kwargs.get('resolution',4))

            filt = gsdf.loc[(gsdf['layer'] == layer_id) & (~gsdf['subtract'])]
            if filt.shape[0] == 0:
                return
            metal_polys = shapely.unary_union(filt['geometry'])
            unit_conv, unit_conv_name = self.get_units()
            metal_polys = shapely.affinity.scale(metal_polys, xfact=unit_conv, yfact=unit_conv, origin=(0,0))

            if layer_id in self.evap_profiles:
                intersec = self.get_all_shadows(metal_polys, layer_id, 'separate_plot_intersections', **kwargs)
            else:
                print('Component does not have evap profile')
            print(f'Expected junction area for layer {layer_id}: {intersec[0].area*1e12} um^2')
            return intersec[0].area*1e12

        else:
            print("Must supply layer_id or component_id!")

    def _add_metallic_layer(self, cur_metal_layers, new_metal_layer, new_metal_layer_height):
        #cur_metal_layers is a list of tuples where the first element is the height and the second element is the polygon...
        new_metal_layers = []
        for cur_height,cur_metal_poly in cur_metal_layers:
            polyInt = shapely.intersection(cur_metal_poly, new_metal_layer)
            non_int_orig_poly = shapely.difference(cur_metal_poly, polyInt)
            if non_int_orig_poly.area > 1e-25:
                new_metal_layers.append( (cur_height, non_int_orig_poly) )
            if polyInt.area > 1e-15:
                new_metal_layers.append((cur_height + new_metal_layer_height, polyInt))
        
        new_metal_layers = sorted(new_metal_layers, key=lambda x: x[0])

        return new_metal_layers
