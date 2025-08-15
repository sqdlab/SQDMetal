# Copyright 2025 Prasanna Pakkiam
# SPDX-License-Identifier: Apache-2.0

import shapely
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt

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
        new_poly = shapely.unary_union(lePolys).buffer(-threshold, join_style=2, cap_style=3)
        #Now that gaps are joined, handle edge-case where the interior hole infinitesimally touches exterior outlines. Sometimes
        #it doesn't cut the polygon into a MultiPolygon...
        return new_poly.buffer(-threshold, join_style=2, cap_style=3).buffer(threshold, join_style=2, cap_style=3)
    
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

    @staticmethod
    def get_smallest_feature_size(shapely_poly):
        lePolys = ShapelyEx.shapely_to_list(shapely_poly)
        minx, miny, maxx, maxy = shapely_poly.bounds

        smallest_dist2 = max(maxx-minx,maxy-miny)**2
        for cur_poly in lePolys:
            for m in range(1,len(cur_poly.exterior.coords)):
                cur_dist2 = (cur_poly.exterior.coords[m][0]-cur_poly.exterior.coords[m-1][0])**2 + (cur_poly.exterior.coords[m][1]-cur_poly.exterior.coords[m-1][1])**2
                smallest_dist2 = min(smallest_dist2, cur_dist2)
            for cur_int in cur_poly.interiors:
                for m in range(1,len(cur_int.coords)):
                    cur_dist2 = (cur_int.coords[m][0]-cur_int.coords[m-1][0])**2 + (cur_int.coords[m][1]-cur_int.coords[m-1][1])**2
                    smallest_dist2 = min(smallest_dist2, cur_dist2)
        
        return np.sqrt(smallest_dist2)

    @staticmethod
    def rasterise_shapely_poly(shapely_poly, size_x, size_y):
        px = 1/plt.rcParams['figure.dpi']  # pixel in inches
        
        fig, ax = plt.subplots(1, figsize=(size_x*px, size_y*px))
        p = gpd.GeoSeries(shapely_poly)
        p.plot(ax=ax)

        fig.tight_layout()
        ax.autoscale()
        ax.set_axis_off()

        left = ax.figure.subplotpars.left
        right = ax.figure.subplotpars.right
        top = ax.figure.subplotpars.top
        bottom = ax.figure.subplotpars.bottom
        figw = float(size_x*px)/(right-left)
        figh = float(size_y*px)/(top-bottom)
        ax.figure.set_size_inches(figw, figh)

        fig.canvas.draw()
        width, height = fig.canvas.get_width_height()
        bbox = ax.get_tightbbox(fig.canvas.get_renderer())
        image_data = fig.canvas.copy_from_bbox(bbox)
        #
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        #
        plt.close(fig)

        img_array = np.asarray(image_data)
        img_array = np.clip(255-img_array[:,:,0], 0,1)
        return img_array[::-1,:].T, [xlim[0], ylim[0], xlim[1], ylim[1]]

    @staticmethod
    def get_shapely_polygon_skeleton(shapely_polyon : shapely.Polygon, **kwargs):
        def get_bisecting_rays(coords_ring: list, orient_CCW=True):
            rays = []
            for m in range(0, len(coords_ring)-1):    #Assuming repeated initial coordinate...
                p0 = coords_ring[m-1] if m > 0 else coords_ring[-2]
                p1 = coords_ring[m]
                p2 = coords_ring[m+1]
                p0 = np.array(p0) 
                p1 = np.array(p1) 
                p2 = np.array(p2)
                #
                v2 = p2-p1
                v2 /= np.linalg.norm(v2)
                v1 = p1-p0
                v1 /= np.linalg.norm(v1)
                #
                bisec_ray = v2-v1
                bisec_ray_norm = np.linalg.norm(bisec_ray)
                #Check if the points are collinear...
                if bisec_ray_norm < 1e-16:
                    bisec_ray = np.array([-v1[1], v1[0]])
                else:
                    bisec_ray /= bisec_ray_norm
                if v1[0]*v2[1]-v1[1]*v2[0] < 0:
                    if orient_CCW:
                        bisec_ray = -bisec_ray
                elif not orient_CCW:
                    bisec_ray = -bisec_ray
                rays.append({
                    'dir': bisec_ray,
                    'pt' : p1,
                    'prev_edge' : -v1,
                    'next_edge' : v2,
                    'child': False
                })
            return rays

        def find_intersec_collapse_ray(shapely_poly, rays: list, overall_rays: list, overall_rays_stilts: list):
            minx,miny,maxx,maxy = shapely_poly.bounds
            biggest_distance =  maxx-minx + maxy-miny

            intersect_pts = []  # noqa: F841 # abhishekchak52: unused variable intersect_pts
            lowest_dist = 1e15
            lowest_dist_ind = -1
            dist_intersection_m = np.zeros(len(rays))+biggest_distance
            dist_intersection_m_1 = np.zeros(len(rays))+biggest_distance
            intersections = np.zeros((len(rays),2))
            for m in range(0,len(rays)):
                a, a0 = rays[m-1]['dir'], rays[m-1]['pt']
                b, b0 = rays[m]['dir'], rays[m]['pt']
                #Writing it as l1 = a0 + a*t and l2 = b0 + b*u
                soln_test = a[1]*b[0] - a[0]*b[1]
                #Disregard the case where the two lines intersect each other directly as that isn't possible for neighbouring rays...
                if np.abs(soln_test) > 1e-9:
                    t = 1/soln_test * (-b[1]*(b0[0]-a0[0]) + b[0]*(b0[1]-a0[1]))
                    if t >= 0:
                        new_pt = a0+a*t
                        # if shapely_poly.contains(shapely.LineString([new_pt,a0])) and shapely_poly.contains(shapely.LineString([new_pt,b0])):
                        if shapely.Point(*new_pt).within(shapely_poly):
                            dist_m_1 = np.linalg.norm(new_pt-a0)
                            dist_m = np.linalg.norm(new_pt-b0)
                            dist_intersection_m[m] = dist_m
                            dist_intersection_m_1[m] = dist_m_1
                            intersections[m] = new_pt
            for m in range(0,len(rays)):
                dist_m = dist_intersection_m[m]
                dist_m_1 = dist_intersection_m_1[m]
                new_pt = intersections[m]
                if dist_intersection_m_1[(m+1)%len(rays)] < dist_intersection_m[m]:
                    continue
                if dist_m < lowest_dist:
                    lowest_dist = dist_m
                    lowest_dist_pt = new_pt
                    lowest_dist_ind = m
                if dist_m_1 < lowest_dist:
                    lowest_dist = dist_m_1
                    lowest_dist_pt = new_pt
                    lowest_dist_ind = m

            if lowest_dist_ind >= 0:
                if rays[lowest_dist_ind-1]['child']:
                    overall_rays.append((rays[lowest_dist_ind-1]['pt'], lowest_dist_pt))
                else:
                    overall_rays_stilts.append((rays[lowest_dist_ind-1]['pt'], lowest_dist_pt))
                if rays[lowest_dist_ind]['child']:
                    overall_rays.append((rays[lowest_dist_ind]['pt'], lowest_dist_pt))
                else:
                    overall_rays_stilts.append((rays[lowest_dist_ind]['pt'], lowest_dist_pt))
                if lowest_dist_ind == 0:
                    ray1 = rays.pop(-1)
                    ray2 = rays.pop(0)
                    ins_ind = 0
                else:
                    ray1 = rays.pop(lowest_dist_ind-1)
                    ray2 = rays.pop(lowest_dist_ind-1)
                    ins_ind = lowest_dist_ind-1
                #ray2 is ahead of ray1
                child_vec = (ray2['next_edge'] + ray1['prev_edge'])/2
                child_vec /= np.linalg.norm(child_vec)
                #Check if the points are collinear...
                # if bisec_ray_norm < 1e-16:
                #     bisec_ray = np.array([-v1[1], v1[0]])
                # else:
                rays.insert(ins_ind, {
                    'dir': child_vec,
                    'pt' : lowest_dist_pt,
                    'prev_edge' : ray1['prev_edge'],
                    'next_edge' : ray2['next_edge'],
                    'child' : True
                })
                return True
            else:
                #Edge-case, so handle it differently...
                min_dist = None
                for m in range(0,len(rays)):
                    target_ind = m
                    if rays[m]['child'] and not rays[m-1]['child']:
                        target_pt = rays[m]['pt']
                    elif not rays[m]['child'] and rays[m-1]['child']:
                        target_pt = rays[m-1]['pt']
                    else:
                        continue
                    cur_dist = np.linalg.norm(rays[m]['pt'] - rays[m-1]['pt'])
                    if min_dist is None or cur_dist < min_dist:
                        min_dist = cur_dist
                        min_targ_pt = target_pt
                        min_targ_ind = target_ind
                if min_dist is None:
                    return False

                child_vec = (rays[min_targ_ind]['next_edge'] + rays[min_targ_ind-1]['prev_edge'])/2
                assert np.linalg.norm(child_vec) > 0

                if min_targ_ind == 0:
                    ray1 = rays.pop(-1)
                    ray2 = rays.pop(0)
                    ins_ind = 0
                else:
                    ray1 = rays.pop(min_targ_ind-1)
                    ray2 = rays.pop(min_targ_ind-1)
                    ins_ind = min_targ_ind-1

                child_vec = (ray2['next_edge'] + ray1['prev_edge'])/2
                child_vec /= np.linalg.norm(child_vec)
                #Check if the points are collinear...
                # if bisec_ray_norm < 1e-16:
                #     bisec_ray = np.array([-v1[1], v1[0]])
                # else:
                rays.insert(min_targ_ind-1, {
                    'dir': child_vec,
                    'pt' : min_targ_pt,
                    'prev_edge' : ray1['prev_edge'],
                    'next_edge' : ray2['next_edge'],
                    'child' : True
                })
                return True


        shapely_poly = shapely.geometry.polygon.orient(shapely_polyon, sign=1.0)

        overall_rays = []
        overall_rays_stilts = []
        rays = get_bisecting_rays(shapely_poly.exterior.coords)

        debug_plots = kwargs.get('debug_plots', False)
        plot_result = kwargs.get('plot_result', False)
        if plot_result or debug_plots:
            vecMag = ShapelyEx.get_smallest_feature_size(shapely_poly)/2
        limit_iters = kwargs.get('limit_iterations', 1e9)

        iters = 0
        while(len(rays) > 2):
            iters += 1
            if iters > limit_iters:
                break

            res = find_intersec_collapse_ray(shapely_poly, rays, overall_rays, overall_rays_stilts)

            if debug_plots:
                fig, ax = plt.subplots(1)
                gdf = gpd.GeoDataFrame(geometry=[shapely_poly])
                gdf.plot(ax=ax)
                for cur_ray in rays:
                    vec, pt = cur_ray['dir'], cur_ray['pt']
                    plt.plot([pt[0], pt[0]+vec[0]*vecMag], [pt[1], pt[1]+vec[1]*vecMag], 'r')
                for cur_ray in overall_rays:
                    plt.plot([cur_ray[0][0], cur_ray[1][0]], [cur_ray[0][1], cur_ray[1][1]], 'b-')
                for cur_ray in overall_rays_stilts:
                    plt.plot([cur_ray[0][0], cur_ray[1][0]], [cur_ray[0][1], cur_ray[1][1]], 'g-')

            if not res:
                break
        
        if plot_result:
            fig, ax = plt.subplots(1)
            gdf = gpd.GeoDataFrame(geometry=[shapely_poly])
            gdf.plot(ax=ax)
            for cur_ray in rays:
                vec, pt = cur_ray['dir'], cur_ray['pt']
                plt.plot([pt[0], pt[0]+vec[0]*vecMag], [pt[1], pt[1]+vec[1]*vecMag], 'r')
            for cur_ray in overall_rays:
                plt.plot([cur_ray[0][0], cur_ray[1][0]], [cur_ray[0][1], cur_ray[1][1]], 'b-')
            for cur_ray in overall_rays_stilts:
                plt.plot([cur_ray[0][0], cur_ray[1][0]], [cur_ray[0][1], cur_ray[1][1]], 'g-')

        return {
            'skeleton': overall_rays,
            'stilts': overall_rays_stilts
        }

