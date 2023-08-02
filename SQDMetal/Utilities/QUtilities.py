# from pint import UnitRegistry
import shapely
import numpy as np
from qiskit_metal.renderers.renderer_mpl.mpl_renderer import QMplRenderer
from qiskit_metal.toolbox_python.utility_functions import bad_fillet_idxs

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
        strVal = strVal.strip()
        assert len(strVal) > 1, f"Length \'{strVal}\' is invalid."
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
    def chk_within(chk_geom, main_geom, thresh=0.99):
        return shapely.intersection(chk_geom, main_geom).area / chk_geom.area > thresh

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
            * dists - List of raw distances (or fractional distances from 0 to 1 if dists_are_fractional is made True) along the path.
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
        '''
        dists = np.sort(dists)

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
        if dists_are_fractional:
            assert np.min(dists) >= 0 and np.max(dists) <= 1, "Fractional distances must be in the interval: [0,1]."
            dists = dists * total_dist

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
    def fuse_polygons_threshold(polys, threshold=1e-12):
        if not isinstance(polys, list):
            polys = [polys]
        lePolys = []
        for cur_poly in polys:
            if isinstance(cur_poly, shapely.geometry.multipolygon.MultiPolygon):
                lePolys += list(cur_poly.geoms)
            else:
                lePolys += [cur_poly]
        lePolys = [p.buffer(threshold, join_style=2, cap_style=3) for p in lePolys]
        return shapely.unary_union(lePolys).buffer(-threshold, join_style=2, cap_style=3)
