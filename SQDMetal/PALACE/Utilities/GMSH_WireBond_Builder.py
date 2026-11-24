# Copyright 2025 Prasanna Pakkiam
# SPDX-License-Identifier: Apache-2.0
#
# Construction of 3D wire-bond geometries in GMSH for Palace simulations.
#
# The builder works in the same (mm) units as GMSH_Geometry_Builder. It has three jobs:
#   1. Generate the 2D centreline (along-span coordinate s vs height z) for a given bridge profile.
#   2. Classify each bond by what its two feet land on (ground plane, CPW trace/other conductor,
#      bare dielectric or off-chip) and what it crosses.
#   3. Build the OCC solid (round wire / rectangular ribbon) or zero-thickness sheet in GMSH.

import gmsh
import numpy as np
import shapely
import warnings


class WireBondProfiles:
    """
    Centreline generators. Every profile returns (s, z) where s runs from 0 to span (along the bond,
    horizontal) and z is the height above the straight baseline joining the two feet. All values in
    the same length units as the inputs.
    """

    SUPPORTED = ['semicircle', 'arc', 'elliptical', 'parabolic', 'trapezoid', 'ball_wedge', 'custom']

    @staticmethod
    def get(profile, span, height, num_pts, params=None):
        params = params or {}
        if callable(profile) or isinstance(profile, (list, tuple, np.ndarray)):
            return WireBondProfiles.custom(profile, span, height, num_pts)
        profile = profile.lower()
        if profile == 'semicircle':
            return WireBondProfiles.arc(span, span/2, num_pts)
        elif profile in ('arc', 'shallow', 'shallow_arc', 'circular'):
            return WireBondProfiles.arc(span, height, num_pts)
        elif profile in ('elliptical', 'ellipse'):
            return WireBondProfiles.elliptical(span, height, num_pts)
        elif profile in ('parabolic', 'parabola'):
            return WireBondProfiles.parabolic(span, height, num_pts)
        elif profile in ('trapezoid', 'square', 'box'):
            rise = params.get('rise_fraction', 0.15)
            r_c = params.get('corner_radius', None)
            pts = [(0, 0), (rise*span, height), ((1-rise)*span, height), (span, 0)]
            return WireBondProfiles.filleted_polyline(pts, r_c if r_c is not None else 0.3*min(height, rise*span*2), num_pts)
        elif profile in ('ball_wedge', 'standard', 'ball'):
            #Vertical neck above the ball (first bond), kink at the loop top, then a near-straight descent to the wedge.
            pk = params.get('peak_position', 0.2)
            neck = params.get('neck_fraction', 0.75)
            r_c = params.get('corner_radius', None)
            pts = [(0, 0), (0, neck*height), (pk*span, height), (span, 0)]
            return WireBondProfiles.filleted_polyline(pts, r_c if r_c is not None else 0.25*height, num_pts)
        assert False, f"Wire-bond profile '{profile}' is unsupported. Use one of {WireBondProfiles.SUPPORTED} or a callable."

    @staticmethod
    def arc(span, height, num_pts):
        #Circular arc through (0,0), (span/2,height), (span,0). height > span/2 gives a 'major' (omega-like) arc.
        assert height > 0, "Wire-bond height must be positive."
        R = (height**2 + (span/2)**2) / (2*height)
        zc = height - R
        phiA = np.arctan2(-zc, -span/2)
        if phiA < 0:
            phiA += 2*np.pi
        phis = np.linspace(phiA, np.pi - phiA, num_pts)
        s = span/2 + R*np.cos(phis)
        z = zc + R*np.sin(phis)
        s[0], z[0], s[-1], z[-1] = 0.0, 0.0, span, 0.0
        return s, z

    @staticmethod
    def elliptical(span, height, num_pts):
        th = np.linspace(np.pi, 0, num_pts)
        return span/2*(1 + np.cos(th)), height*np.sin(th)

    @staticmethod
    def parabolic(span, height, num_pts):
        t = 0.5*(1 - np.cos(np.linspace(0, np.pi, num_pts)))    #Cluster points near the feet
        s = t*span
        return s, 4*height*s*(span - s)/span**2

    @staticmethod
    def custom(profile, span, height, num_pts):
        if callable(profile):
            #f(t) with t in [0,1] returning the height (same units as span), f(0)=f(1)=0
            t = 0.5*(1 - np.cos(np.linspace(0, np.pi, num_pts)))
            z = np.array([profile(x) for x in t], dtype=float)
            return t*span, z
        pts = np.array(profile, dtype=float)
        assert pts.ndim == 2 and pts.shape[1] == 2, "Custom profiles must be given as Nx2 array of (s_fraction, z)."
        return pts[:, 0]*span, pts[:, 1]

    @staticmethod
    def filleted_polyline(pts, radius, num_pts):
        #Polyline with circular fillets on interior vertices, densely sampled so that a spline through it is faithful.
        pts = np.array(pts, dtype=float)
        segs = [pts[0][None, :]]
        n_arc = max(6, num_pts//6)
        for m in range(1, len(pts)-1):
            a, b, c = pts[m-1], pts[m], pts[m+1]
            u1 = (a - b)/np.linalg.norm(a - b)
            u2 = (c - b)/np.linalg.norm(c - b)
            ang = np.arccos(np.clip(np.dot(u1, u2), -1, 1))
            if ang > np.pi - 1e-6 or radius <= 0:
                segs.append(b[None, :])
                continue
            d = radius/np.tan(ang/2)
            dmax = 0.45*min(np.linalg.norm(a - b), np.linalg.norm(c - b))
            r = radius if d <= dmax else dmax*np.tan(ang/2)
            d = min(d, dmax)
            t1, t2 = b + u1*d, b + u2*d
            bis = (u1 + u2)/np.linalg.norm(u1 + u2)
            ctr = b + bis*r/np.sin(ang/2)
            a1 = np.arctan2(*(t1 - ctr)[::-1])
            a2 = np.arctan2(*(t2 - ctr)[::-1])
            da = (a2 - a1 + np.pi) % (2*np.pi) - np.pi
            angs = a1 + da*np.linspace(0, 1, n_arc)
            segs.append(np.c_[ctr[0] + r*np.cos(angs), ctr[1] + r*np.sin(angs)])
        segs.append(pts[-1][None, :])
        #Densify straight sections
        out = [segs[0][0]]
        for seg in segs[1:]:
            for p in seg:
                prev = out[-1]
                L = np.linalg.norm(p - prev)
                n_sub = int(np.ceil(L / (np.ptp(pts[:, 0])/num_pts + 1e-30)))
                for k in range(1, max(n_sub, 1)):
                    out.append(prev + (p - prev)*k/n_sub)
                out.append(p)
        out = np.array(out)
        keep = np.r_[True, np.linalg.norm(np.diff(out, axis=0), axis=1) > 1e-9*(np.ptp(pts[:, 0]) + 1)]
        out = out[keep]
        return out[:, 0], out[:, 1]


class GMSH_WireBond_Builder:
    #Materials treated as superconducting (i.e. PEC by default) at cryogenic temperatures
    SUPERCONDUCTORS = ['al', 'aluminium', 'aluminum', 'alsi', 'al1%si', 'nb', 'niobium', 'in', 'indium', 'nbti', 'ta', 'tantalum', 'pec']

    def __init__(self, metals, chip_bounds, z_surface, metal_thickness, extents, ground_plane_included=True):
        """
        Parameters (all lengths in GMSH units - i.e. mm):
            metals      - list of shapely Polygons of the (fused) metallic regions, indexed as in GMSH_Geometry_Builder
            chip_bounds - (xmin, ymin, xmax, ymax) of the chip
            z_surface   - z-coordinate of the top of the substrate
            metal_thickness - thickness of the planar metals (0 for quasi-3D sheets)
            extents     - simulation air-box extents [(xmin,xmax),(ymin,ymax),(zmin,zmax)]
        """
        self.metals = metals
        self.chip_poly = shapely.box(*chip_bounds)
        self.z_surf = z_surface
        self.t_metal = metal_thickness
        self.extents = np.array(extents)
        self.ground_ids = self._identify_ground_metals(ground_plane_included)
        self._footprints = []

    def _identify_ground_metals(self, ground_plane_included):
        #Ground = largest conductor touching the chip perimeter, plus any other large conductors touching the perimeter
        #(e.g. ground islands separated by a feedline). Cropped CPW traces also touch the perimeter but are small.
        if len(self.metals) == 0:
            return []
        perim = self.chip_poly.exterior.buffer(1e-6)
        areas = np.array([x.area for x in self.metals])
        touching = [m for m, x in enumerate(self.metals) if x.intersects(perim)]
        if len(touching) == 0:
            m = int(np.argmax(areas))
            return [m] if areas[m] > 0.3*self.chip_poly.area else []
        a_max = max(areas[touching])
        return [m for m in touching if areas[m] >= 0.05*a_max]

    @staticmethod
    def _rect(centre, u, length, width):
        n = np.array([-u[1], u[0]])
        c = np.array(centre)
        return shapely.Polygon([c - u*length/2 - n*width/2, c + u*length/2 - n*width/2,
                                c + u*length/2 + n*width/2, c - u*length/2 + n*width/2])

    def _locate_foot(self, pt, u, foot_len, foot_wid):
        if not self.chip_poly.contains(shapely.Point(pt)):
            return {'location': 'off_chip', 'metal_index': None, 'coverage': 0.0}
        fp = self._rect(pt, u, max(foot_len, 1e-9), max(foot_wid, 1e-9))
        overlaps = np.array([fp.intersection(x).area for x in self.metals]) if len(self.metals) else np.array([0.0])
        #The conductor directly under the bond point wins; otherwise fall back to the largest overlap with the foot
        under = [m for m, x in enumerate(self.metals) if x.covers(shapely.Point(pt))]
        m = under[0] if len(under) > 0 else int(np.argmax(overlaps))
        if overlaps[m] <= 0:
            return {'location': 'dielectric', 'metal_index': None, 'coverage': 0.0}
        return {'location': 'ground' if m in self.ground_ids else 'trace', 'metal_index': m,
                'coverage': float(overlaps[m]/fp.area), 'centre_on_metal': len(under) > 0}

    def classify(self, wb):
        """
        Classify a wire-bond (dict in mm units, already normalised by the Model) based on its landing sites.

        Returns a dict with keys:
            category  - 'ground' (both ends on ground / ground-to-off-chip) or 'cpw_trace' (any end on a non-ground
                        conductor) or 'floating'.
            bond_type - 'ground_stitch', 'ground_island_bridge', 'ground_offchip', 'trace_jumper', 'trace_stitch',
                        'trace_offchip', 'trace_to_ground' or 'floating_foot'
            feet      - per-foot location dictionaries
            crosses   - metal indices crossed by the bond's span (excluding the landing conductors)
        """
        pA, pB = np.array(wb['pos_start']), np.array(wb['pos_end'])
        span = np.linalg.norm(pB - pA)
        assert span > 0, f"Wire-bond {wb['name']} has coincident start and end positions."
        u = (pB - pA)/span
        fl, fw = wb['foot_length'], wb['foot_width']
        fA = self._locate_foot(pA, u, fl, fw)
        fB = self._locate_foot(pB, u, fl, fw)
        locs = (fA['location'], fB['location'])
        assert locs != ('off_chip', 'off_chip'), f"Wire-bond {wb['name']}: both ends are off-chip."

        issues = []
        if 'dielectric' in locs:
            bond_type, category = 'floating_foot', 'floating'
            issues.append("a foot lands on bare dielectric (no conductor under it)")
        elif locs == ('ground', 'ground'):
            category = 'ground'
            bond_type = 'ground_stitch' if fA['metal_index'] == fB['metal_index'] else 'ground_island_bridge'
        elif set(locs) == {'ground', 'off_chip'}:
            bond_type, category = 'ground_offchip', 'ground'
        elif locs == ('trace', 'trace'):
            category = 'cpw_trace'
            bond_type = 'trace_stitch' if fA['metal_index'] == fB['metal_index'] else 'trace_jumper'
        elif set(locs) == {'trace', 'off_chip'}:
            bond_type, category = 'trace_offchip', 'cpw_trace'
        else:
            bond_type, category = 'trace_to_ground', 'cpw_trace'
            issues.append("it shorts a non-ground conductor (e.g. CPW centre trace) to ground")

        for f, lbl in [(fA, 'start'), (fB, 'end')]:
            if f['metal_index'] is not None and not f['centre_on_metal']:
                issues.append(f"the {lbl} foot is centred on a dielectric gap (assigned to metal {f['metal_index']} by largest overlap)")
            if f['metal_index'] is not None and f['coverage'] < 0.999:
                issues.append(f"the {lbl} foot overhangs its conductor ({100*f['coverage']:.0f}% on metal)")

        #Conductors crossed by the span (e.g. the CPW centre trace under a ground stitch)
        line = shapely.LineString([pA, pB])
        landing = {fA['metal_index'], fB['metal_index']}
        crosses = [m for m, x in enumerate(self.metals) if m not in landing and line.intersects(x)]

        user_type = wb.get('bond_type', 'auto')
        if user_type not in ('auto', None):
            if user_type in ('ground', 'cpw_trace'):
                if user_type != category:
                    issues.append(f"declared as '{user_type}' but detected as '{category}' ({bond_type})")
            elif user_type != bond_type:
                issues.append(f"declared as '{user_type}' but detected as '{bond_type}'")

        for msg in issues:
            warnings.warn(f"Wire-bond {wb['name']}: {msg}.")

        return {'category': category, 'bond_type': bond_type, 'feet': [fA, fB], 'crosses': crosses,
                'span': span, 'issues': issues}
    
    def _foot_z(self, foot, wb):
        if foot['location'] == 'off_chip':
            if wb['offchip_termination'] == 'floor':
                return self.extents[2][0]
            return self.z_surf + wb['offchip_z']
        if foot['location'] == 'dielectric':
            return self.z_surf
        return self.z_surf + self.t_metal

    def centreline(self, wb, cls):
        pA, pB = np.array(wb['pos_start'], dtype=float), np.array(wb['pos_end'], dtype=float)
        zA, zB = self._foot_z(cls['feet'][0], wb), self._foot_z(cls['feet'][1], wb)
        span = cls['span']
        u = (pB - pA)/span
        height = wb['height'] if wb['height'] is not None else 0.3*span
        s, z = WireBondProfiles.get(wb['profile'], span, height, wb['num_path_points'], wb['profile_params'])
        base = self._baseline(s, span, pA, pB, zA, zB, cls)
        pts = np.c_[pA[0] + u[0]*s, pA[1] + u[1]*s, base + z]
        return pts, (zA, zB), u

    def _baseline(self, s, span, pA, pB, zA, zB, cls):
        '''
        Height of the straight 'baseline' under the bridge profile. When one end is off-chip at a different height
        (e.g. the enclosure floor or a recessed PCB), the baseline stays level with the on-chip foot until the bond has
        cleared the chip edge and only then descends/ascends to the off-chip landing height - so that the wire never
        cuts through the substrate.
        '''
        if abs(zA - zB) < 1e-12:
            return np.full_like(s, zA)
        locs = [f['location'] for f in cls['feet']]
        if 'off_chip' not in locs:
            return zA + (zB - zA)*s/span
        crossing = shapely.LineString([pA, pB]).intersection(self.chip_poly.boundary)
        pts = [np.array(g.coords[0]) for g in getattr(crossing, 'geoms', [crossing]) if not g.is_empty]
        if len(pts) == 0:
            return zA + (zB - zA)*s/span
        s_edges = [np.linalg.norm(x - pA) for x in pts]
        if locs[1] == 'off_chip':
            s0 = min(s_edges)
            return np.where(s <= s0, zA, zA + (zB - zA)*(s - s0)/max(span - s0, 1e-30))
        s0 = max(s_edges)
        return np.where(s >= s0, zB, zA + (zB - zA)*s/max(s0, 1e-30))

    def _check_centreline(self, wb, pts, r):
        ext = self.extents
        tol = 1e-9
        top = pts[:, 2].max() + r
        assert top < ext[2][1] - tol, f"Wire-bond {wb['name']} (top at z={top:.4g}) pokes out of the air-box (z_max={ext[2][1]:.4g}). Increase the z-boundary."
        assert pts[:, 0].min() - r > ext[0][0] and pts[:, 0].max() + r < ext[0][1] and \
               pts[:, 1].min() - r > ext[1][0] and pts[:, 1].max() + r < ext[1][1], \
               f"Wire-bond {wb['name']} leaks outside the x/y boundaries of the air-box."
        #Interior points over (or within a wire radius of) the chip must clear the chip surface
        chip_zone = self.chip_poly.buffer(r) if r > 0 else self.chip_poly
        for p in pts[1:-1]:
            if chip_zone.contains(shapely.Point(p[:2])) and p[2] < self.z_surf - tol:
                assert False, f"Wire-bond {wb['name']} intersects the substrate at {p}. Check offchip_z/height or move the off-chip end further from the chip edge."

    def check_footprint_overlap(self, wb, cls):
        """Returns the names of previously registered bonds whose plan-view footprint overlaps this one (candidates for a
        3D intersection check - crossing at different heights is fine)."""
        line = shapely.LineString([wb['pos_start'], wb['pos_end']]).buffer(max(wb['width'], wb['foot_width'])/2 + 1e-9)
        hits = [other_name for other_name, other_fp in self._footprints if line.intersects(other_fp)]
        self._footprints.append((wb['name'], line))
        return hits

    @staticmethod
    def check_3D_intersection(name_a, dimtags_a, name_b, dimtags_b):
        """Asserts that two built wire-bonds do not physically intersect (OCC boolean intersection of copies)."""
        for da in dimtags_a:
            for db in dimtags_b:
                if da[0] == 2 and db[0] == 2:
                    continue    #Two sheets can only meet along curves, which fragmentation handles
                ca = gmsh.model.occ.copy([da])
                cb = gmsh.model.occ.copy([db])
                res, _ = gmsh.model.occ.intersect(ca, cb, removeObject=True, removeTool=True)
                gmsh.model.occ.synchronize()
                size = sum(gmsh.model.occ.getMass(*x) for x in res if x[0] == min(da[0], db[0]))
                if len(res) > 0:
                    gmsh.model.occ.remove(res, recursive=True)
                gmsh.model.occ.synchronize()
                assert size <= 1e-12, f"Wire-bonds {name_a} and {name_b} physically intersect. Move them apart, or change their heights/profiles so that one passes over the other."

    def build(self, wb, cls):
        """
        Builds the wire-bond in GMSH (OCC kernel). Returns list of (dim, tag): 3D volumes for 'round'/'ribbon'
        cross-sections, 2D surfaces for 'sheet' (plus any 2D off-chip landing pads).
        """
        xs = wb['cross_section']
        w = wb['width']
        t = wb['thickness'] if xs == 'ribbon' else w
        r = 0.5*max(w, t) if xs != 'sheet' else 0.0
        pts, (zA, zB), u = self.centreline(wb, cls)
        self._check_centreline(wb, pts, r)
        n_h = np.array([-u[1], u[0], 0.0])  #Horizontal normal to the bond plane

        #Extend the ends into the clip region along the end tangents so that the clipped ends are full cross-sections
        ext = 1.5*r + wb['foot_thickness'] + 1e-4*cls['span']
        t0 = pts[1] - pts[0]; t0 /= np.linalg.norm(t0)
        t1 = pts[-1] - pts[-2]; t1 /= np.linalg.norm(t1)
        if xs != 'sheet':
            pts = np.vstack([pts[0] - t0*ext, pts, pts[-1] + t1*ext])
        p_tags = [gmsh.model.occ.addPoint(*p) for p in pts]
        spline = gmsh.model.occ.addSpline(p_tags)
        wire = gmsh.model.occ.addWire([spline])

        tan0 = pts[1] - pts[0]; tan0 /= np.linalg.norm(tan0)
        p0 = pts[0]
        if xs == 'round':
            prof = [(2, gmsh.model.occ.addDisk(p0[0], p0[1], p0[2], r, r, zAxis=tan0.tolist(), xAxis=n_h.tolist()))]
        elif xs == 'ribbon':
            b = np.cross(tan0, n_h); b /= np.linalg.norm(b)
            corners = [p0 - n_h*w/2 - b*t/2, p0 + n_h*w/2 - b*t/2, p0 + n_h*w/2 + b*t/2, p0 - n_h*w/2 + b*t/2]
            prof = [(2, self._plane_from_pts(corners))]
        elif xs == 'sheet':
            la = gmsh.model.occ.addPoint(*(p0 - n_h*w/2))
            lb = gmsh.model.occ.addPoint(*(p0 + n_h*w/2))
            prof = [(1, gmsh.model.occ.addLine(la, lb))]
        else:
            assert False, f"Wire-bond cross_section '{xs}' unsupported (use 'round', 'ribbon' or 'sheet')."

        swept = gmsh.model.occ.addPipe(prof, wire, wb['trihedron'])
        gmsh.model.occ.remove(prof, recursive=True)
        gmsh.model.occ.remove([(1, wire)], recursive=True)
        gmsh.model.occ.synchronize()

        if xs == 'sheet':
            return swept + self._offchip_pads(wb, cls, (zA, zB), u)

        #Feet (flattened wedge/ball deformation) fused onto the wire
        feet = []
        if wb['include_feet'] and wb['foot_length'] > 0 and wb['foot_thickness'] > 0:
            for p, z in [(wb['pos_start'], zA), (wb['pos_end'], zB)]:
                feet.append((3, self._oriented_box(p, z, u, wb['foot_length'], wb['foot_width'], wb['foot_thickness'])))
        vol = swept
        if feet:
            vol, _ = gmsh.model.occ.fuse(swept, feet, removeObject=True, removeTool=True)

        #Clip everything below each landing plane
        clip = []
        half = 0.5*max(wb['foot_length'], wb['foot_width']) + 3*r + ext
        for p, z in [(wb['pos_start'], zA), (wb['pos_end'], zB)]:
            depth = 4*ext + 4*r
            clip.append((3, gmsh.model.occ.addBox(p[0] - half, p[1] - half, z - depth, 2*half, 2*half, depth)))
        vol, _ = gmsh.model.occ.cut(vol, clip, removeObject=True, removeTool=True)
        gmsh.model.occ.synchronize()
        vol = [x for x in vol if x[0] == 3]
        assert len(vol) > 0, f"Wire-bond {wb['name']}: GMSH failed to build the solid (check width vs. profile curvature)."
        return vol + self._offchip_pads(wb, cls, (zA, zB), u)

    AIRBRIDGE_PROFILES = ['arch', 'flat', 'trapezoid', 'elliptical', 'parabolic']
    AIRBRIDGE_ALIASES = {'arc': 'arch', 'circular': 'arch', 'reflow': 'arch', 'box': 'flat', 'rectangular': 'flat',
                         'pedestal': 'flat', 'ellipse': 'elliptical', 'parabola': 'parabolic'}

    @staticmethod
    def airbridge_profile_name(profile):
        if callable(profile) or isinstance(profile, (list, tuple, np.ndarray)):
            return 'custom'
        name = str(profile).lower()
        return GMSH_WireBond_Builder.AIRBRIDGE_ALIASES.get(name, name)

    def build_airbridge(self, wb, cls):
        """
        Builds a 3D air-bridge (OCC kernel). The bridge lies in the vertical plane containing its two feet: the side profile
        (feet + deck) is drawn in that plane and extruded across the deck width. Feet (the I-bar pad sections) wider than the
        deck are fused on as boxes. Returns a list of (dim, tag): 3D volumes for 'solid' and 2D surfaces for 'sheet'.

        Geometry (s along the bridge, z above the landing surface; feet of length foot_length centred on the end points):
            * 'arch' (reflowed resist), 'trapezoid', 'elliptical', 'parabolic', custom: a conformal film of (vertical)
              thickness ``thickness`` whose underside follows the profile with apex ``height`` over the inner span.
            * 'flat' (two-layer resist pedestal): vertical posts filling the feet up to a flat deck whose underside is at
              ``height`` and whose thickness is ``thickness``.
        """
        locs = [f['location'] for f in cls['feet']]
        assert 'off_chip' not in locs, f"Air-bridge {wb['name']}: both feet must lie on the chip."
        xs = wb['cross_section']
        w, t, h, Lf = wb['width'], wb['thickness'], wb['height'], wb['foot_length']
        span = cls['span']
        pA, pB = np.array(wb['pos_start'], dtype=float), np.array(wb['pos_end'], dtype=float)
        u = (pB - pA)/span
        n_h = np.array([-u[1], u[0], 0.0])
        zA, zB = self._foot_z(cls['feet'][0], wb), self._foot_z(cls['feet'][1], wb)
        s1, s2 = Lf/2, span - Lf/2
        L = s2 - s1
        assert L > 0, f"Air-bridge {wb['name']}: the feet overlap - the distance between the foot centres ({span*1e3:.3g}um) must exceed foot_length ({Lf*1e3:.3g}um)."
        sA, sB = -Lf/2, span + Lf/2

        prof = self.airbridge_profile_name(wb['profile'])
        flat = prof == 'flat'
        if not flat:
            name = {'arch': 'arc'}.get(prof, prof)
            s, z = WireBondProfiles.get(wb['profile'] if name == 'custom' else name, L, h, wb['num_path_points'], wb['profile_params'])
            s = np.array(s, dtype=float) + s1
            lower = np.c_[s, zA + (zB - zA)*(s - s1)/L + np.array(z, dtype=float)]
            lower[0], lower[-1] = (s1, zA), (s2, zB)

        p0 = pA - n_h[:2]*w/2
        to3 = lambda s, z: (p0[0] + u[0]*s, p0[1] + u[1]*s, z)
        top = (max(zA, zB) + h + (t if xs == 'solid' else 0)) if flat else lower[:, 1].max() + (t if xs == 'solid' else 0)
        ext = self.extents
        assert top < ext[2][1] - 1e-9, f"Air-bridge {wb['name']} (top at z={top:.4g}) pokes out of the air-box (z_max={ext[2][1]:.4g})."
        r = 0.5*max(w, wb['foot_width'])
        for p in [pA, pB]:
            assert p[0] - r - Lf > ext[0][0] and p[0] + r + Lf < ext[0][1] and p[1] - r - Lf > ext[1][0] and p[1] + r + Lf < ext[1][1], \
                f"Air-bridge {wb['name']} leaks outside the x/y boundaries of the air-box."

        if xs == 'sheet':
            #Zero-thickness deck (no feet); its bottom edges sit on the landing conductors
            if flat:
                segs = [[(s1, zA), (s1, zA + h)], [(s1, zA + h), (s2, zB + h)], [(s2, zB + h), (s2, zB)]]
            else:
                segs = [lower]
            curves = self._add_chain(segs, to3, closed=False)
            out = gmsh.model.occ.extrude([(1, c) for c in curves], *(n_h*w))
            gmsh.model.occ.synchronize()
            surfs = [x for x in out if x[0] == 2]
            if len(surfs) > 1:
                surfs, _ = gmsh.model.occ.fuse(surfs[:1], surfs[1:], removeObject=True, removeTool=True)
                gmsh.model.occ.synchronize()
                surfs = [x for x in surfs if x[0] == 2]
            return surfs

        assert xs == 'solid', f"Air-bridge cross_section '{xs}' unsupported (use 'solid' or 'sheet')."
        if flat:
            hA, hB = zA + h + t, zB + h + t
            segs = [[(sA, zA), (s1, zA)], [(s1, zA), (s1, zA + h)], [(s1, zA + h), (s2, zB + h)], [(s2, zB + h), (s2, zB)],
                    [(s2, zB), (sB, zB)], [(sB, zB), (sB, hB)], [(sB, hB), (sA, hA)], [(sA, hA), (sA, zA)]]
            foot_h = h + t
        else:
            upper = lower[::-1] + np.array([0, t])
            segs = [[(sA, zA), (s1, zA)], lower, [(s2, zB), (sB, zB)], [(sB, zB), (sB, zB + t)], [(sB, zB + t), (s2, zB + t)],
                    upper, [(s1, zA + t), (sA, zA + t)], [(sA, zA + t), (sA, zA)]]
            foot_h = t
        curves = self._add_chain(segs, to3, closed=True)
        surf = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop(curves)])
        out = gmsh.model.occ.extrude([(2, surf)], *(n_h*w))
        gmsh.model.occ.synchronize()
        vol = [x for x in out if x[0] == 3]
        if wb['foot_width'] > w*(1 + 1e-9):
            feet = [(3, self._oriented_box(p, z, u, Lf, wb['foot_width'], foot_h)) for p, z in [(pA, zA), (pB, zB)]]
            vol, _ = gmsh.model.occ.fuse(vol, feet, removeObject=True, removeTool=True)
            gmsh.model.occ.synchronize()
            vol = [x for x in vol if x[0] == 3]
        assert len(vol) > 0, f"Air-bridge {wb['name']}: GMSH failed to build the solid."
        return vol

    @staticmethod
    def _add_chain(segs, to3, closed):
        '''
        Adds a chain of curves through the (s,z) segments (2 points -> line, more -> spline). Consecutive segments share
        their end points; for a closed chain, the last segment must end on the first point.
        '''
        first, prev, curves = None, None, []
        for k, seg in enumerate(segs):
            seg = np.asarray(seg, dtype=float)
            start = prev if prev is not None else gmsh.model.occ.addPoint(*to3(*seg[0]))
            if first is None:
                first = start
            end = first if (closed and k == len(segs) - 1) else gmsh.model.occ.addPoint(*to3(*seg[-1]))
            mids = [gmsh.model.occ.addPoint(*to3(*p)) for p in seg[1:-1]]
            curves.append(gmsh.model.occ.addLine(start, end) if len(mids) == 0 else gmsh.model.occ.addSpline([start] + mids + [end]))
            prev = end
        return curves

    def _offchip_pads(self, wb, cls, zs, u):
        pads = []
        if wb['offchip_termination'] != 'pad':
            return pads
        for f, p, z in zip(cls['feet'], [wb['pos_start'], wb['pos_end']], zs):
            if f['location'] == 'off_chip':
                L, W = wb['offchip_pad_size']
                rect = self._rect(p, u, L, W)
                bx = rect.bounds
                assert bx[0] > self.extents[0][0] and bx[2] < self.extents[0][1] and bx[1] > self.extents[1][0] and bx[3] < self.extents[1][1], \
                    f"Wire-bond {wb['name']}: the off-chip landing pad leaks outside the air-box. Increase the x/y boundaries or reduce offchip_pad_size."
                assert not rect.intersects(self.chip_poly), f"Wire-bond {wb['name']}: the off-chip landing pad overlaps the chip; move the end further off-chip or reduce offchip_pad_size."
                pads.append((2, self._plane_from_pts([[c[0], c[1], z] for c in rect.exterior.coords[:-1]])))
        return pads

    @staticmethod
    def _plane_from_pts(corners):
        pt = [gmsh.model.occ.addPoint(*c) for c in corners]
        ln = [gmsh.model.occ.addLine(pt[m], pt[(m + 1) % len(pt)]) for m in range(len(pt))]
        return gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop(ln)])

    @staticmethod
    def _oriented_box(p, z, u, length, width, height):
        box = gmsh.model.occ.addBox(p[0] - length/2, p[1] - width/2, z, length, width, height)
        ang = np.arctan2(u[1], u[0])
        if abs(ang) > 1e-12:
            gmsh.model.occ.rotate([(3, box)], p[0], p[1], z, 0, 0, 1, ang)
        return box

    @staticmethod
    def resolve_material(material, boundary='auto'):
        """
        Returns (boundary_type, MaterialConductor or None). Superconductors default to PEC; normal metals
        (e.g. Au, Cu) default to a finite-conductivity (surface impedance) boundary in Palace.
        """
        from SQDMetal.Utilities.Materials import MaterialConductor
        if isinstance(material, MaterialConductor):
            return ('pec' if boundary == 'pec' else 'conductivity'), material
        name = str(material).lower().replace(' ', '').replace('_', '').replace('-', '')
        if boundary == 'pec' or (boundary == 'auto' and name in GMSH_WireBond_Builder.SUPERCONDUCTORS):
            return 'pec', None
        alias = {'cu': 'ofhc', 'copper': 'ofhc'}
        try:
            return 'conductivity', MaterialConductor(alias.get(name, material))
        except AssertionError:
            assert False, f"No conductivity is known for '{material}' (a finite-conductivity boundary needs one). Pass a " \
                          f"MaterialConductor object, or use boundary='pec' for superconductors."
