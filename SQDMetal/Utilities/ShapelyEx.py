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

    @staticmethod
    def simplify_line_edges(coords, min_angle_deg):
        '''
        Give coordinates as a list of (x,y) pairs. It must be a closed shell (i.e. coords[0]==coords[-1]). The
        min_angle_deg specifies the angle upon which to discard a corner if below this angle...
        '''
        min_angle = min_angle_deg/180*np.pi
        cur_coords = coords[:-1]    #Remove the repeated closure point...
        deleted = True
        while deleted:
            deleted = False
            numPts = len(cur_coords)
            chosen_coords = np.zeros(numPts, dtype=int)
            for m in range(numPts):
                if chosen_coords[m] == 1:
                    continue
                vec_a = np.array([cur_coords[m][0]-cur_coords[m-1][0], cur_coords[m][1]-cur_coords[m-1][1]])
                vec_b = np.array([cur_coords[(m+1)%numPts][0]-cur_coords[m][0], cur_coords[(m+1)%numPts][1]-cur_coords[m][1]])
                vec_a_norm = np.linalg.norm(vec_a)
                vec_b_norm = np.linalg.norm(vec_b)
                len_fac = vec_a_norm/vec_b_norm if vec_a_norm>vec_b_norm else vec_b_norm/vec_a_norm
                phi = np.arccos(np.dot(vec_a, vec_b)/(vec_a_norm*vec_b_norm))
                if phi* len_fac <= min_angle:
                    chosen_coords[m-1] = 1
                    chosen_coords[(m+1)%numPts] = 1
                    # chosen_coords[m] = 0
                    deleted = True
                else:
                    chosen_coords[m] = 1
            cur_coords = [cur_coords[m] for m in range(numPts) if chosen_coords[m] == 1]
        cur_coords.append(cur_coords[0])    #Close the loop again...
        return cur_coords

    @staticmethod
    def simplify_poly_edges(shapely_poly, min_angle_deg):
        list_shapely_polys = ShapelyEx.shapely_to_list(shapely_poly)

        filtered_polys = []
        for cur_poly in list_shapely_polys:
            exteriors = ShapelyEx.simplify_line_edges(cur_poly.exterior.coords[:], min_angle_deg)
            interiors = [ShapelyEx.simplify_line_edges(cur_int.coords[:], min_angle_deg) for cur_int in cur_poly.interiors]
            filtered_polys.append(shapely.Polygon(exteriors, interiors))

        filtered_polys = shapely.unary_union(filtered_polys)

        return filtered_polys
