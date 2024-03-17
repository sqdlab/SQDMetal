import shapely

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
