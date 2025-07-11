import shapely
import numpy as np

class ShapelyEx:
    @staticmethod
    def chk_within(chk_geom, main_geom, thresh=0.99):
        return shapely.intersection(chk_geom, main_geom).area / chk_geom.area > thresh

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
    
    @staticmethod
    def rectangle(x1,y1,x2,y2, make_poly=True):
        x_1 = min(x1,x2)
        x_2 = max(x1,x2)
        y_1 = min(y1,y2)
        y_2 = max(y1,y2)
        coords = [(x_1, y_2), (x_2, y_2), (x_2, y_1), (x_1, y_1)]
        if make_poly:
            return shapely.Polygon(coords)
        else:
            return coords

    @staticmethod
    def rectangle_from_line(pos1,pos2, rect_width, make_poly=True):
        v_parl = pos2-pos1
        v_perp = np.array([-v_parl[1], v_parl[0]])
        v_perp /= np.linalg.norm(v_perp)
        v_perp *= rect_width*0.5
        coords = np.array([pos1+v_perp, pos1-v_perp, pos1-v_perp+v_parl, pos1+v_perp+v_parl])
        if make_poly:
            return shapely.Polygon(coords)
        else:
            return coords

    @staticmethod
    def shapely_to_list(shapely_obj):
        if isinstance(shapely_obj, shapely.geometry.multipolygon.MultiPolygon) or isinstance(shapely_obj, shapely.geometry.multilinestring.MultiLineString):
            return [x for x in shapely_obj.geoms]
        else:
            return [shapely_obj] #i.e. it's just a lonely Polygon or LineString object...

    @staticmethod
    def get_points_uniform_in_polygon(shapely_poly, spacing_x, spacing_y):
        min_x, min_y, max_x, max_y = shapely_poly.bounds
        smpl_x = np.arange(min_x,max_x+spacing_x, spacing_x)
        smpl_y = np.arange(min_y,max_y+spacing_y, spacing_y)
        ret_array = []
        for x in smpl_x:
            for y in smpl_y:
                if shapely_poly.contains(shapely.geometry.Point(x,y)):
                    ret_array.append([x, y])
        return np.array(ret_array)

