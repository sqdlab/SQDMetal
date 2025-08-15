# Copyright 2025 Prasanna Pakkiam
# SPDX-License-Identifier: Apache-2.0

import numpy as np
import shapely
from SQDMetal.Utilities.ShapelyEx import ShapelyEx

class GeomBase:
    @property
    def ParserType(self):
        raise NotImplementedError()

    @property
    def chip_size_x(self):
        raise NotImplementedError()

    @property
    def chip_size_y(self):
        raise NotImplementedError()

    @property
    def chip_size_z(self):
        raise NotImplementedError()

    @property
    def chip_centre(self):
        raise NotImplementedError()

    def process_layers(self, metallic_layers, ground_plane, **kwargs):
        raise NotImplementedError()

    def create_CPW_feed_via_point(self, pt_near_centre_end, len_launch, metallic_layers, ground_plane, **kwargs):
        metals, gaps = self.process_layers(metallic_layers, ground_plane, **kwargs)
        poly_all = shapely.unary_union(metals)

        thresh_90_degrees = kwargs.get('thresh_90_degrees', 0.1)

        thresh_90 = np.cos((90-thresh_90_degrees)/180*np.pi)
        ptPort = shapely.Point(pt_near_centre_end[0], pt_near_centre_end[1])
        is_inside = poly_all.contains(ptPort)

        end_point = shapely.shortest_line(poly_all.boundary, ptPort).coords[0]

        poly_boundary_lines = ShapelyEx.shapely_to_list(poly_all.boundary)

        corner_points = []
        for cur_line in poly_boundary_lines:
            for m in range(len(cur_line.coords)-1):  #Ignore last point...
                if m > 0:
                    coord_1 = cur_line.coords[m-1]
                else:
                    coord_1 = cur_line.coords[-2]
                coord_2 = cur_line.coords[m+1]
                #
                vec1 = np.array([coord_1[0]-cur_line.coords[m][0], coord_1[1]-cur_line.coords[m][1]])
                vec2 = np.array([coord_2[0]-cur_line.coords[m][0], coord_2[1]-cur_line.coords[m][1]])
                #
                norm1 = np.linalg.norm(vec1)
                if norm1 < 1e-16:
                    continue
                norm2 = np.linalg.norm(vec2)
                if norm2 < 1e-16:
                    continue
                vec1 /= norm1
                vec2 /= norm2
                #
                if np.abs(np.dot(vec1, vec2)) < thresh_90:
                    corner_points.append(cur_line.coords[m])
        corner_points = np.array(corner_points)

        #Filter corner points that are collinear to the edge segment...
        vec_perp = np.array([end_point[0]-pt_near_centre_end[0], end_point[1]-pt_near_centre_end[1]])
        if is_inside:
            vec_perp = -vec_perp
        vec_perp /= np.linalg.norm(vec_perp)
        vec_parallel = np.array([vec_perp[1], -vec_perp[0]])
        collinear_corners_1 = []
        collinear_corners_2 = []
        for cur_corner in corner_points:
            vec_cur_corner = np.array([cur_corner[0]-end_point[0], cur_corner[1]-end_point[1]])
            vec_cur_corner_norm = np.linalg.norm(vec_cur_corner)
            vec_cur_corner /= vec_cur_corner_norm
            if np.abs(np.dot(vec_perp, vec_cur_corner)) < thresh_90:
                if np.dot(vec_cur_corner, vec_parallel) > 0:
                    collinear_corners_1.append((vec_cur_corner_norm, cur_corner))
                else:
                    collinear_corners_2.append((vec_cur_corner_norm, cur_corner))

        #Find 4 closest collinear corners to the edge point...
        pad_pts_1, pad_pts_2 = sorted(collinear_corners_1)[:2], sorted(collinear_corners_2)[:2]

        assert len(pad_pts_1) == 2, f"Seems that a CPW couldn't be formed. Check that the line-segment closest to the point {pt_near_centre_end} is indeed the edge of the central conductor of the CPW feed."
        assert len(pad_pts_2) == 2, f"Seems that a CPW couldn't be formed. Check that the line-segment closest to the point {pt_near_centre_end} is indeed the edge of the central conductor of the CPW feed."

        #Draw the pads...
        pad_1 = [
            pad_pts_1[0][1],
            pad_pts_1[1][1],
            pad_pts_1[1][1]+vec_perp*len_launch,
            pad_pts_1[0][1]+vec_perp*len_launch
        ]
        pad_1 = [[p[0], p[1]] for p in pad_1]
        #
        pad_2 = [
            pad_pts_2[0][1],
            pad_pts_2[1][1],
            pad_pts_2[1][1]+vec_perp*len_launch,
            pad_pts_2[0][1]+vec_perp*len_launch
        ]
        pad_2 = [[p[0], p[1]] for p in pad_2]

        #Use for debugging...
        # p = gpd.GeoSeries(poly_all)
        # p.plot()
        # plt.plot([pt_near_centre_end[0]],[pt_near_centre_end[1]],'bx')
        # plt.plot([end_point[0]],[end_point[1]],'bx')
        # plt.plot(corner_points[:,0], corner_points[:,1], 'ro')
        # plt.plot([x[1][0] for x in collinear_corners_1], [x[1][1] for x in collinear_corners_1], 'gx')
        # plt.plot([x[1][0] for x in collinear_corners_2], [x[1][1] for x in collinear_corners_2], 'cx')
        # plt.plot([x[0] for x in pad_1], [x[1] for x in pad_1], 'm')
        # plt.plot([x[0] for x in pad_2], [x[1] for x in pad_2], 'm')
        # plt.show()

        return pad_1, pad_2, vec_perp.tolist()