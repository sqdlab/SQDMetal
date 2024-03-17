# from pint import UnitRegistry
import shapely
import numpy as np
from qiskit_metal.toolbox_python.utility_functions import bad_fillet_idxs
from SQDMetal.Utilities.PVD_Shadows import PVD_Shadows
from SQDMetal.Utilities.QiskitShapelyRenderer import QiskitShapelyRenderer
from qiskit_metal.qlibrary.terminations.launchpad_wb import LaunchpadWirebond
from SQDMetal.Utilities.ShapelyEx import ShapelyEx

class QUtilities:
    @staticmethod
    def get_units(design):
        unit_conv = design.get_units()
        if unit_conv == 'mm':
            unit_conv = 1e-3
        elif unit_conv == 'um':
            unit_conv = 1e-6
        elif unit_conv == 'nm':
            unit_conv = 1e-9
        else:
            assert False, f"Unrecognised units: {unit_conv}"
        return unit_conv

    @staticmethod
    def parse_value_length(strVal):
        #This is far too slow!?
        # ureg = UnitRegistry().Quantity
        # return ureg(strVal).to('m').magnitude

        #So do this instead...
        if isinstance(strVal, int) or isinstance(strVal, float):
            return strVal
        strVal = strVal.strip()
        assert len(strVal) > 1, f"Length \'{strVal}\' is invalid (no units?)."
        if strVal[-2:] == "mm":
            return float(strVal[:-2]+'e-3')
        elif strVal[-2:] == "um":
            return float(strVal[:-2]+'e-6')
        elif strVal[-2:] == "nm":
            return float(strVal[:-2]+'e-9')
        elif strVal[-2:] == "pm":
            return float(strVal[:-2]+'e-12')
        elif strVal[-2:] == "fm":
            return float(strVal[:-2]+'e-15')
        elif strVal[-1:] == "m":
            return float(strVal[:-2])
        else:
            assert len(strVal) > 1, f"Length \'{strVal}\' is invalid."

    @staticmethod
    def get_comp_bounds(design, objs, units_metres = False):
        if not (isinstance(objs, list) or isinstance(objs, np.ndarray)):
            objs = [objs]
        x_vals = []
        y_vals = []
        unit_conv = QUtilities.get_units(design)
        for cur_obj in objs:
            paths = design.components[cur_obj].qgeometry_table('path')
            for _, row in paths.iterrows():
                cur_minX, cur_minY, cur_maxX, cur_maxY = row['geometry'].buffer(row['width'] / 2, cap_style=shapely.geometry.CAP_STYLE.flat).bounds
                x_vals += [cur_minX, cur_maxX]
                y_vals += [cur_minY, cur_maxY]
            for cur_poly in design.components[cur_obj].qgeometry_list('poly'):
                cur_minX, cur_minY, cur_maxX, cur_maxY = cur_poly.bounds
                x_vals += [cur_minX, cur_maxX]
                y_vals += [cur_minY, cur_maxY]
            #In case the geometry is empty, just take the pos_x and pos_y identifiers (e.g. the Joint object)...
            if len(design.components[cur_obj].qgeometry_table('path')) == 0 and len(design.components[cur_obj].qgeometry_table('poly')) == 0:
                x_vals += [QUtilities.parse_value_length(design.components[cur_obj].options.pos_x)/unit_conv]
                y_vals += [QUtilities.parse_value_length(design.components[cur_obj].options.pos_y)/unit_conv]
        if units_metres:
            return unit_conv*min(x_vals), unit_conv*min(y_vals), unit_conv*max(x_vals), unit_conv*max(y_vals)
        else:
            return min(x_vals), min(y_vals), max(x_vals), max(y_vals)

    @staticmethod
    def calc_points_on_path(dists, design, component_name, trace_name='', trace_name_gap='', dists_are_fractional=False):
        '''
        Returns the points along a path taking into account the curved edges (i.e. fillets).

        Inputs:
            * dists - List of raw distances (or fractional distances from 0 to 1 if dists_are_fractional is made True) along the path. If dists
                      is given as a single value, then dists_are_fractional is ignored, and the returned path will be a bunch of points spaced
                      by the value given by dists. Units are in Qiskit-Metal design units (when not fractional).
            * component_name - Name of the QComponent in the design
            * trace_name - (Optional) If provided, then the path given by the trace_name is selected. Otherwise, the default path (typically
                           called 'trace' in normal wiring/routing objects) is selected.
            * trace_name_gap - (Optional) If provided, then the returned gap width is given by the trace_name_gap entity. Otherwise, the
                               default path (typically called 'cut' in normal wiring/routing objects) is taken for the gap width.
            * dists_are_fractional - (Defaults to True) If True, then the argument dists is taken as a list of  fractional distances from 0 to 1.
        
        Returns tuple (final_pts, normals, width, gap, total_dist) where:
            * final_pts - Numpy array of coordinates along the path given as (x,y) on every row.
            * normals   - Numpy array of unit-vectors for the right-hand normal vector for every point in final_pts given as (nx, ny) on every row.
            * width - CPW trace width.
            * gap   - CPW trace gap.
            * total_dist - Total length of path.
        Note that the units in the returned values are in Qiskit-Metal design units...
        '''

        df = design.qgeometry.tables['path']
        df = df[df['component'] == design.components[component_name].id]#['geometry'][0]
        if trace_name != '':
            dfTrace = df[df['name']==trace_name]
            if trace_name_gap != '':
                dfGap = df[df['name']==trace_name_gap]
            else:
                dfGap = df[df['name']==trace_name]
        else:
            dfTrace = df[df['subtract']==False]
            dfGap = df[df['subtract']==True]
        width = dfTrace['width'].iloc[0]
        if len(dfGap) > 0:
            gap = (dfGap['width'].iloc[0] - width)*0.5
        else:
            gap = 0
        assert len(dfTrace) > 0, f"The component \'{component_name}\' has no valid path."

        rFillet = dfTrace['fillet'].iloc[0]
        line_segs = QUtilities.calc_lines_and_fillets_on_path(np.array(dfTrace['geometry'].iloc[0].coords[:]), rFillet, design.template_options.PRECISION)

        total_dist = sum([x['dist'] for x in line_segs])
        
        if isinstance(dists, (list, tuple, np.ndarray)):
            dists = np.sort(dists)
            if dists_are_fractional:
                assert np.min(dists) >= 0 and np.max(dists) <= 1, "Fractional distances must be in the interval: [0,1]."
                dists = dists * total_dist
        else:
            dists = np.arange(0, total_dist, dists)
            if dists[-1] != total_dist:
                dists = np.concatenate([dists, [total_dist]])

        final_pts = []
        normals = []
        line_seg_ind = 0
        cur_seg_dist = 0
        for cur_dist in dists:
            while cur_dist > cur_seg_dist + line_segs[line_seg_ind]['dist']:
                cur_seg_dist += line_segs[line_seg_ind]['dist']
                line_seg_ind += 1
            if 'centre' in line_segs[line_seg_ind]:
                angleReq = (cur_dist - cur_seg_dist) / rFillet
                angleReq = line_segs[line_seg_ind]['angleStart'] + angleReq * np.sign(line_segs[line_seg_ind]['angleDelta'])
                norm_dir = np.array([np.cos(angleReq), np.sin(angleReq)])
                new_pt = line_segs[line_seg_ind]['centre'] + rFillet * norm_dir
                normals += [norm_dir * np.sign(line_segs[line_seg_ind]['angleDelta'])]
            else:
                new_pt = line_segs[line_seg_ind]['start'] + (cur_dist - cur_seg_dist) * line_segs[line_seg_ind]['dir']
                normals += [[line_segs[line_seg_ind]['dir'][1], -line_segs[line_seg_ind]['dir'][0]]]
            final_pts += [new_pt]
        return np.array(final_pts), np.array(normals), width, gap, total_dist

    @staticmethod
    def calc_lines_and_fillets_on_path(points, rFillet, precision):
        '''
        Returns the line segments and curved fillets on a given path of points.

        Inputs:
            * points - Numpy array of points given as (x,y) on every row.
            * rFillet - Radius of fillets on every corner.
            * precision - The numeric precision in which to calculate bad fillets taken usually as design.template_options.PRECISION

        The function returns list of segments each represented as dictionaries. Straight segments have the keys:
            * 'start' - Starting point of line segment
            * 'end'   - Ending point of line segment
            * 'dir'   - Unit-vector of the direction of the line segment
            * 'dist'  - Length of the line segment
        Curved segments (i.e. ones with fillets) have the keys:
            * 'start'  - Starting point of line segment
            * 'end'    - Ending point of line segment
            * 'centre' - Centre of the circle drawing the fillet arc
            * 'angleStart' - Starting angle (on Cartesian plane) of the fillet arc (in radians)
            * 'angleDelta' - Angle traversed by the fillet arc (in radians) using the typical polar angle direction convention in R2...
            * 'dist'   - Arclength of the fillet
        '''

        # Get list of vertices that can't be filleted
        no_fillet = bad_fillet_idxs(points, rFillet, precision)

        line_segs = []
        curPt = points[0]
        for (m, (corner, end)) in enumerate(zip(points[1:], points[2:])):
            a = corner - curPt
            b = end - corner
            aDist = np.linalg.norm(a)
            bDist = np.linalg.norm(b)
            aHat = a / aDist
            bHat = b / bDist
            theta = np.arccos(np.dot(-aHat,bHat))
            if m + 1 in no_fillet or np.abs(theta-np.pi) < 1e-9:
                line_segs += [{'start':curPt, 'end':corner, 'dir':aHat, 'dist':np.linalg.norm(corner-curPt)}]
                curPt = corner
            else:
                distTangentEdge = rFillet/np.tan(theta/2)
                ptFilletStart = corner - aHat*distTangentEdge
                ptFilletEnd = corner + bHat*distTangentEdge

                vCentre = bHat-np.dot(aHat,bHat)*aHat
                vCentre *= rFillet / np.linalg.norm(vCentre)
                ptCentre = ptFilletStart + vCentre

                angleStart = ptFilletStart - ptCentre
                angleStart = np.arctan2(angleStart[1], angleStart[0])
                #
                if aHat[0]*bHat[1] > aHat[1]*bHat[0]:
                    angleDelta = theta
                else:
                    angleDelta = -theta

                line_segs += [{'start':curPt, 'end':ptFilletStart, 'dir':aHat, 'dist':np.linalg.norm(ptFilletStart-curPt)}]
                line_segs += [{'start':ptFilletStart, 'end':ptFilletEnd, 'centre':ptCentre, 'angleStart':angleStart, 'angleDelta':angleDelta, 'dist':(np.pi-theta)*rFillet}]
                curPt = ptFilletEnd
        distLast = np.linalg.norm(points[-1]-curPt)
        line_segs += [{'start':curPt, 'end':points[-1], 'dir':(points[-1]-curPt)/distLast, 'dist':distLast}]

        return line_segs

    @staticmethod
    def get_metals_in_layer(design, layer_id, **kwargs):
        '''
        Partitions unique conductors from a layer in a Qiskit-Metal design object. If the particular layer has fancy PVD evaporation steps, the added
        metallic layer will account for said steps and merge the final result. In addition, all metallic elements that are contiguous are merged into
        single blobs.

        Inputs:
            - design - Qiskit-Metal deisgn object
            - layer_id - The index of the layer from which to take the metallic polygons
            - restrict_rect - (Optional) List highlighting the clipping rectangle [xmin,ymin,xmax,ymax]
            - resolution - (Optional) Defaults to 4. This is the number of points along a curved section of a path object in Qiskit-Metal.
            - threshold - (Optional) Defaults to -1. This is the threshold in metres, below which consecutive vertices along a given polygon are
                          combined into a single vertex. This simplification helps with meshing as COMSOL will not overdo the meshing. If this
                          argument is negative, the argument is ignored.
            - fuse_threshold - (Optional) Defaults to 1e-12. This is the minimum distance between metallic elements, below which they are considered
                               to be a single polygon and thus, the polygons are merged with the gap filled. This accounts for floating-point errors
                               that make adjacent elements fail to merge as a single element, due to infinitesimal gaps between them.
            - evap_mode - (Optional) Defaults to 'separate_delete_below'. These are the methods upon which to separate or merge overlapping elements
                          across multiple evaporation steps. See documentation on PVD_Shadows for more details on the available options.
            - group_by_evaporations - (Optional) Defaults to False. If set to True, if elements on a particular evaporation step are separated due
                                      to the given evap_mode, they will still be selected as a part of the same conductor (useful for example, in
                                      capacitance matrix simulations).
            - evap_trim - (Optional) Defaults to 20e-9. This is the trimming distance used in certain evap_mode profiles. See documentation on
                          PVD_Shadows for more details on its definition.
            - unit_conv - (Optional) Unit conversion factor to convert the Qiskit-Metal design units. Defaults to converting the units into metres.
        
        Outputs a tuple containing:
            - metal_polys_all - All separate metallic islands
            - metal_sel_ids   - Selection indices for each of the metallic islands (mostly relevant when using group_by_evaporations)
        '''
        #Fresh update on PVD profiles...
        pvd_shadows = PVD_Shadows(design)

        thresh = kwargs.get('threshold',  -1)
        resolution = kwargs.get('resolution', 4)

        qmpl = QiskitShapelyRenderer(None, design, None)
        gsdf = qmpl.get_net_coordinates(resolution)

        filt = gsdf.loc[(gsdf['layer'] == layer_id) & (gsdf['subtract'] == False)]
        if filt.shape[0] == 0:
            return

        unit_conv = kwargs.get('unit_conv', QUtilities.get_units(design))
        
        fuse_threshold = kwargs.get('fuse_threshold', 1e-12)

        #Merge the metallic elements
        metal_polys = shapely.unary_union(filt['geometry'])
        metal_polys = shapely.affinity.scale(metal_polys, xfact=unit_conv, yfact=unit_conv, origin=(0,0))
        metal_polys = ShapelyEx.fuse_polygons_threshold(metal_polys, fuse_threshold)
        restrict_rect = kwargs.get('restrict_rect', None)
        if isinstance(restrict_rect, list):
            metal_polys = shapely.clip_by_rect(metal_polys, *restrict_rect)
        #Calculate the individual evaporated elements if required
        evap_mode = kwargs.get('evap_mode', 'separate_delete_below')
        group_by_evaporations = kwargs.get('group_by_evaporations', False)
        if group_by_evaporations and evap_mode != 'merge':
            metal_evap_polys_separate = pvd_shadows.get_all_shadows(metal_polys, layer_id, 'separate')
            #Convert all MultiPolygons into individual polygons...
            if not isinstance(metal_evap_polys_separate, list):
                metal_evap_polys_separate = [metal_evap_polys_separate]
        #Calculate evaporated shadows
        evap_trim = kwargs.get('evap_trim', 20e-9)
        metal_evap_polys = pvd_shadows.get_all_shadows(metal_polys, layer_id, evap_mode, layer_trim_length=evap_trim)
        #Convert all MultiPolygons into individual polygons...
        if not isinstance(metal_evap_polys, list):
            metal_evap_polys = [metal_evap_polys]

        metal_polys_all = []
        metal_sel_ids = []
        for m, cur_poly in enumerate(metal_evap_polys):
            if isinstance(cur_poly, shapely.geometry.multipolygon.MultiPolygon):
                temp_cur_metals = [x for x in cur_poly.geoms]
                metal_polys_all += temp_cur_metals
                num_polys = len(temp_cur_metals)
            else:
                temp_cur_metals = [cur_poly] #i.e. it's just a lonely Polygon object...
                metal_polys_all += temp_cur_metals
                num_polys = 1

            if group_by_evaporations and evap_mode != 'merge':
                #Collect the separate polygons that live in the current evaporation layer
                cur_polys_separate = metal_evap_polys_separate[m]
                if isinstance(cur_polys_separate, shapely.geometry.multipolygon.MultiPolygon):
                    cur_polys_separate = [x for x in cur_polys_separate.geoms]
                #Find the separate polygon in which the given polygon fits... 
                cur_sel_inds = []
                for cur_metal in temp_cur_metals:
                    for sep_piece_ind, cur_metal_piece in enumerate(cur_polys_separate):
                        if shapely.intersection(cur_metal, cur_metal_piece).area >= 0.99*cur_metal.area:
                            cur_sel_inds += [len(metal_sel_ids) + sep_piece_ind]
                            break
                metal_sel_ids += cur_sel_inds
            else:
                #Just enumerate to all separate metallic pieces...
                metal_sel_ids += [x for x in range(len(metal_sel_ids), len(metal_sel_ids) + num_polys)]

        return metal_polys_all, metal_sel_ids

    @staticmethod
    def get_perimetric_polygons(design, comp_names, **kwargs):
        thresh = kwargs.get('threshold',  -1)
        resolution = kwargs.get('resolution', 4)

        qmpl = QiskitShapelyRenderer(None, design, None)
        gsdf = qmpl.get_net_coordinates(resolution)

        ids = [design.components[x].id for x in comp_names]
        filt = gsdf[gsdf['component'].isin(ids)]

        if filt.shape[0] == 0:
            return
    
        unit_conv = kwargs.get('unit_conv', QUtilities.get_units(design))
        fuse_threshold = kwargs.get('fuse_threshold', 1e-12)

        metal_polys = shapely.unary_union(filt['geometry'])
        metal_polys = shapely.affinity.scale(metal_polys, xfact=unit_conv, yfact=unit_conv, origin=(0,0))
        metal_polys = ShapelyEx.fuse_polygons_threshold(metal_polys, fuse_threshold)
        restrict_rect = kwargs.get('restrict_rect', None)
        if isinstance(restrict_rect, list):
            metal_polys = shapely.clip_by_rect(metal_polys, *restrict_rect)

        if isinstance(metal_polys, shapely.geometry.multipolygon.MultiPolygon):
            lePolys = list(metal_polys.geoms)
        else:
            lePolys = [metal_polys]
        
        return lePolys

    @staticmethod
    def _get_LauncherWB_params(design, launcher_name, unit_conv_extra=1):        
        launcher_len = QUtilities.parse_value_length(design.components[launcher_name].options['pad_height']) + QUtilities.parse_value_length(design.components[launcher_name].options['taper_height']) + QUtilities.parse_value_length(design.components[launcher_name].options['lead_length'])
        unit_conv = QUtilities.get_units(design)
        startPt = design.components[launcher_name].pins['tie']['middle']*unit_conv - design.components[launcher_name].pins['tie']['normal']*launcher_len
        padDir = design.components[launcher_name].pins['tie']['normal']*1.0
        padWid = QUtilities.parse_value_length(design.components[launcher_name].options['pad_width'])
        padGap = QUtilities.parse_value_length(design.components[launcher_name].options['pad_gap'])

        return startPt * unit_conv_extra, padDir, padWid * unit_conv_extra, padGap * unit_conv_extra

    @staticmethod
    def _get_Route_params(design, route_name, pin_name, unit_conv_extra=1):
        unit_conv = QUtilities.get_units(design)

        startPt = design.components[route_name].pins[pin_name]['middle']*unit_conv
        padDir = -1.0*design.components[route_name].pins[pin_name]['normal']
        padWid = QUtilities.parse_value_length(design.components[route_name].options.trace_width)
        padGap = QUtilities.parse_value_length(design.components[route_name].options.trace_gap)

        return startPt * unit_conv_extra, padDir, padWid * unit_conv_extra, padGap * unit_conv_extra

    @staticmethod
    def get_RFport_CPW_coords_Launcher(design, qObjName, len_launch = 20e-6, unit_conv_extra = 1):
        qObj = design.components[qObjName]
        if isinstance(qObj, LaunchpadWirebond):
            vec_ori, vec_launch, cpw_wid, cpw_gap = QUtilities._get_LauncherWB_params(design, qObjName, unit_conv_extra)
        else:
            assert False, f"\'{qObjName}\' is an unsupported object type."

        vec_perp = np.array([vec_launch[1],-vec_launch[0]])
        vec_launch *= len_launch
        
        launchesA = [vec_ori + vec_perp * cpw_wid*0.5, vec_ori + vec_perp * (cpw_wid*0.5+cpw_gap),
                    vec_ori + vec_launch + vec_perp * (cpw_wid*0.5+cpw_gap), vec_ori + vec_launch + vec_perp * cpw_wid*0.5]
        launchesA = [[p[0],p[1]] for p in launchesA]

        launchesB = [vec_ori - vec_perp * cpw_wid*0.5, vec_ori - vec_perp * (cpw_wid*0.5+cpw_gap),
                    vec_ori + vec_launch - vec_perp * (cpw_wid*0.5+cpw_gap), vec_ori + vec_launch - vec_perp * cpw_wid*0.5]
        launchesB = [[p[0],p[1]] for p in launchesB]

        return launchesA, launchesB, vec_perp

    @staticmethod
    def get_RFport_CPW_coords_Route(design, qObjName, pin_name, len_launch = 20e-6, unit_conv_extra = 1):
        qObj = design.components[qObjName]

        #TODO: Add in type-checking somehow here?
        vec_ori, vec_launch, cpw_wid, cpw_gap = QUtilities._get_Route_params(design, qObjName, pin_name, unit_conv_extra)

        vec_perp = np.array([vec_launch[1],-vec_launch[0]])
        vec_launch *= len_launch

        launchesA = [vec_ori + vec_perp * cpw_wid*0.5, vec_ori + vec_perp * (cpw_wid*0.5+cpw_gap),
                    vec_ori + vec_launch + vec_perp * (cpw_wid*0.5+cpw_gap), vec_ori + vec_launch + vec_perp * cpw_wid*0.5]
        launchesA = [[p[0],p[1]] for p in launchesA]

        launchesB = [vec_ori - vec_perp * cpw_wid*0.5, vec_ori - vec_perp * (cpw_wid*0.5+cpw_gap),
                    vec_ori + vec_launch - vec_perp * (cpw_wid*0.5+cpw_gap), vec_ori + vec_launch - vec_perp * cpw_wid*0.5]
        launchesB = [[p[0],p[1]] for p in launchesB]

        return launchesA, launchesB, vec_perp

