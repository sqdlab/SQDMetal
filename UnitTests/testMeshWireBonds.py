"""
Tests for the 3D wire-bonds and air-bridges in the Palace simulations. Nothing here needs Palace itself: the fast tests
check the profiles, landing-site classification, materials and the simulation API; the integration tests build the GMSH
mesh and Palace config file and check what ends up in them.
"""
from qiskit_metal import designs, Dict
from qiskit_metal.qlibrary.terminations.launchpad_wb import LaunchpadWirebond
from qiskit_metal.qlibrary.tlines.straight_path import RouteStraight
from qiskit_metal.qlibrary.sample_shapes.rectangle import Rectangle

from SQDMetal.Comps.WireBonds import WireBond, WireBondsCPW
from SQDMetal.Comps.AirBridges import AirBridgeCPW
from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation
from SQDMetal.PALACE.Capacitance_Simulation import PALACE_Capacitance_Simulation
from SQDMetal.PALACE.Utilities.GMSH_WireBond_Builder import GMSH_WireBond_Builder, WireBondProfiles
from SQDMetal.PALACE.Utilities.GMSH_Geometry_Builder import GMSH_Geometry_Builder
from SQDMetal.Utilities.Materials import MaterialInterface

import contextlib
import io
import json
import os
import shutil
import tempfile
import unittest
import warnings

import gmsh
import numpy as np
import shapely

#The mesh/config files of the integration tests are kept here for inspection (e.g. open the .msh files in Gmsh). Each test
#class's folder is cleared when it starts, so its contents always come from the latest run.
OUTPUT_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "wirebond_test_outputs")


def output_dir(*parts):
    path = os.path.join(OUTPUT_ROOT, *parts)
    shutil.rmtree(path, ignore_errors=True)
    os.makedirs(path)
    return path


def save_summary(sim, path):
    #Human-readable classification of every bond/bridge, next to the mesh
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        sim.get_wirebond_summary()
    with open(path, 'w') as f:
        f.write(buf.getvalue())


def build_feedline_design(island=False):
    """2.5mm x 1.2mm chip with a CPW feedline between two launch pads (and optionally an isolated metal island)."""
    design = designs.DesignPlanar({}, True)
    design.chips.main.size['size_x'] = '2.5mm'
    design.chips.main.size['size_y'] = '1.2mm'
    LaunchpadWirebond(design, 'lpL', options=dict(pos_x='-0.9mm', pos_y='0mm', orientation='0'))
    LaunchpadWirebond(design, 'lpR', options=dict(pos_x='0.9mm', pos_y='0mm', orientation='180'))
    RouteStraight(design, 'feed', Dict(trace_width='10um', trace_gap='6um',
                  pin_inputs=Dict(start_pin=Dict(component='lpL', pin='tie'), end_pin=Dict(component='lpR', pin='tie'))))
    if island:
        Rectangle(design, 'island_cut', options=dict(pos_x='0mm', pos_y='0.35mm', width='200um', height='150um', subtract=True))
        Rectangle(design, 'island', options=dict(pos_x='0mm', pos_y='0.35mm', width='100um', height='60um'))
    design.rebuild()
    return design


def launch_pad_centre(design, name):
    pad = design.qgeometry.get_component(name)['poly']
    c = pad[pad['name'] == 'launch_pad']['geometry'].iloc[0].centroid
    return c.x, c.y


def quiet(func, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return func(*args, **kwargs)


class _TempDirTest(unittest.TestCase):
    def setUp(self):
        self._dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self._dir, ignore_errors=True)


# ------------------------------------------------------------------ Profiles

class TestWireBondProfiles(unittest.TestCase):
    SPAN, HEIGHT, N = 1.0, 0.3, 41

    def check_ends(self, s, z, span):
        self.assertAlmostEqual(s[0], 0, places=9)
        self.assertAlmostEqual(s[-1], span, places=9)
        self.assertAlmostEqual(z[0], 0, places=9)
        self.assertAlmostEqual(z[-1], 0, places=9)
        self.assertTrue(np.all(z >= -1e-12))

    def test_apex_heights(self):
        for name in ['arc', 'elliptical', 'parabolic']:
            s, z = WireBondProfiles.get(name, self.SPAN, self.HEIGHT, self.N)
            self.check_ends(s, z, self.SPAN)
            self.assertAlmostEqual(z.max(), self.HEIGHT, delta=1e-3, msg=name)
        s, z = WireBondProfiles.get('semicircle', self.SPAN, self.HEIGHT, self.N)
        self.check_ends(s, z, self.SPAN)
        self.assertAlmostEqual(z.max(), self.SPAN/2, delta=1e-3)

    def test_filleted_profiles(self):
        #Fillets round off the corners, so the apex is slightly below the requested height
        for name in ['trapezoid', 'ball_wedge']:
            s, z = WireBondProfiles.get(name, self.SPAN, self.HEIGHT, self.N)
            self.check_ends(s, z, self.SPAN)
            self.assertLessEqual(z.max(), self.HEIGHT + 1e-9, msg=name)
            self.assertGreater(z.max(), 0.8*self.HEIGHT, msg=name)

    def test_major_arc(self):
        #Taller than half the span gives an omega-like loop that bulges past the feet
        s, z = WireBondProfiles.get('arc', self.SPAN, 0.8, self.N)
        self.check_ends(s, z, self.SPAN)
        self.assertLess(s.min(), 0)
        self.assertAlmostEqual(z.max(), 0.8, delta=1e-3)

    def test_custom(self):
        s, z = WireBondProfiles.get(lambda t: 0.2*np.sin(np.pi*t), 2.0, None, self.N)
        self.check_ends(s, z, 2.0)
        self.assertAlmostEqual(z.max(), 0.2, delta=1e-3)
        s, z = WireBondProfiles.get(np.array([[0, 0], [0.5, 0.1], [1, 0]]), 2.0, None, self.N)
        np.testing.assert_allclose(s, [0, 1, 2])
        np.testing.assert_allclose(z, [0, 0.1, 0])

    def test_unsupported(self):
        with self.assertRaises(AssertionError):
            WireBondProfiles.get('zigzag', 1.0, 0.3, 11)


# ------------------------------------------------------------------ Classification

class TestWireBondClassification(unittest.TestCase):
    """Units are mm (as inside the GMSH builders). Metals: 0 = ground, 1 = CPW trace, 2 = isolated island."""

    @classmethod
    def setUpClass(cls):
        chip = shapely.box(-1, -0.5, 1, 0.5)
        trace = shapely.box(-0.8, -0.005, 0.8, 0.005)
        island = shapely.box(0.3, 0.2, 0.35, 0.25)
        ground = chip.difference(shapely.box(-0.9, -0.011, 0.9, 0.011)).difference(shapely.box(0.28, 0.18, 0.37, 0.27))
        extents = [(-1.5, 1.5), (-1, 1), (-0.5, 0.5)]
        cls.builder = GMSH_WireBond_Builder([ground, trace, island], (-1, -0.5, 1, 0.5), 0.0, 0.0, extents, True)

    def classify(self, pA, pB, foot=(0.01, 0.005), bond_type='auto'):
        wb = dict(name='test', pos_start=np.array(pA), pos_end=np.array(pB), foot_length=foot[0], foot_width=foot[1], bond_type=bond_type)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            cls = self.builder.classify(wb)
        return cls, [str(x.message) for x in w]

    def test_ground_identification(self):
        self.assertEqual(self.builder.ground_ids, [0])

    def test_ground_stitch(self):
        cls, warns = self.classify((0, -0.05), (0, 0.05))
        self.assertEqual((cls['category'], cls['bond_type']), ('ground', 'ground_stitch'))
        self.assertEqual(cls['crosses'], [1])
        self.assertAlmostEqual(cls['span'], 0.1)
        self.assertEqual(warns, [])

    def test_trace_jumper(self):
        cls, warns = self.classify((0.5, 0), (0.325, 0.225))
        self.assertEqual((cls['category'], cls['bond_type']), ('cpw_trace', 'trace_jumper'))
        self.assertEqual([f['metal_index'] for f in cls['feet']], [1, 2])
        self.assertIn(0, cls['crosses'])

    def test_trace_to_ground_warns(self):
        cls, warns = self.classify((0, 0), (0, 0.05))
        self.assertEqual(cls['bond_type'], 'trace_to_ground')
        self.assertTrue(any('shorts' in x for x in warns))

    def test_offchip(self):
        cls, _ = self.classify((0, 0.3), (0, 0.7))
        self.assertEqual((cls['category'], cls['bond_type']), ('ground', 'ground_offchip'))
        self.assertEqual(cls['feet'][1]['location'], 'off_chip')
        cls, _ = self.classify((-0.7, 0), (-1.2, 0))
        self.assertEqual((cls['category'], cls['bond_type']), ('cpw_trace', 'trace_offchip'))
        with self.assertRaises(AssertionError):
            self.classify((0, 0.7), (0, 0.8))

    def test_dielectric_foot(self):
        cls, warns = self.classify((0, 0.008), (0, 0.05), foot=(0.002, 0.001))
        self.assertEqual((cls['category'], cls['bond_type']), ('floating', 'floating_foot'))
        self.assertTrue(any('dielectric' in x for x in warns))

    def test_foot_overhang_warns(self):
        #A wide foot on the 10um trace spills over into the gaps
        cls, warns = self.classify((0.5, 0), (0.325, 0.225), foot=(0.03, 0.03))
        self.assertLess(cls['feet'][0]['coverage'], 0.999)
        self.assertTrue(any('overhangs' in x for x in warns))

    def test_declared_type_mismatch_warns(self):
        _, warns = self.classify((0.5, 0), (0.325, 0.225), bond_type='ground')
        self.assertTrue(any("declared as 'ground'" in x for x in warns))
        _, warns = self.classify((0, -0.05), (0, 0.05), bond_type='ground')
        self.assertEqual(warns, [])

    def test_coincident_ends(self):
        with self.assertRaises(AssertionError):
            self.classify((0, 0.1), (0, 0.1))


# ------------------------------------------------------------------ Materials and conductor merging

class TestWireBondMaterials(unittest.TestCase):
    def test_superconductors_are_pec(self):
        for m in ['Al', 'aluminium', 'Nb', 'In', 'NbTi']:
            self.assertEqual(GMSH_WireBond_Builder.resolve_material(m), ('pec', None), msg=m)

    def test_normal_metals_use_conductivity(self):
        for m in ['Au', 'Cu']:
            bc, mat = GMSH_WireBond_Builder.resolve_material(m)
            self.assertEqual(bc, 'conductivity', msg=m)
            self.assertGreater(mat.conductivity, 0)

    def test_boundary_override(self):
        self.assertEqual(GMSH_WireBond_Builder.resolve_material('Au', 'pec'), ('pec', None))
        #Superconductors have no tabulated conductivity, so forcing a conductivity boundary must say so clearly
        with self.assertRaisesRegex(AssertionError, 'No conductivity is known'):
            GMSH_WireBond_Builder.resolve_material('Al', 'conductivity')
        from SQDMetal.Utilities.Materials import MaterialConductor
        self.assertEqual(GMSH_WireBond_Builder.resolve_material(MaterialConductor('gold'))[0], 'conductivity')


class TestConductorMerging(unittest.TestCase):
    @staticmethod
    def meta(name, pg, feet_idx):
        return {'name': name, 'physical_group': pg, 'pos_start': np.array([0, 0]), 'pos_end': np.array([1, 0]), 'width': 0.1,
                'classification': {'feet': [{'metal_index': x} for x in feet_idx]}}

    def test_union_find(self):
        metals = [shapely.box(0, 0, 1, 1), shapely.box(2, 0, 3, 1), shapely.box(4, 0, 5, 1), shapely.box(6, 0, 7, 1)]
        mapping = [[10], [11], [12], [13]]
        wbs = [self.meta('a', 20, [2, 1]), self.meta('b', 21, [3, 2]), self.meta('c', 22, [None, None])]
        new_map, polys, groups = GMSH_Geometry_Builder._merge_wirebonded_conductors(None, mapping, metals, wbs)
        self.assertEqual(groups, [[0], [1, 2, 3], []])
        self.assertEqual(new_map, [[10], [11, 12, 13, 20, 21], [22]])
        self.assertEqual(len(polys), 3)
        self.assertAlmostEqual(polys[1].area, 3.0)
        self.assertEqual([x['conductor_index'] for x in wbs], [1, 1, 2])


# ------------------------------------------------------------------ Simulation API (no meshing)

class TestWireBondAPI(_TempDirTest):
    @classmethod
    def setUpClass(cls):
        cls.design = build_feedline_design()
        AirBridgeCPW(cls.design, 'abs', options=dict(pathObj='feed', dist_air_bridges='300um', dist_air_bridges_ends='200um'))
        WireBondsCPW(cls.design, 'stitches', options=dict(pathObj='feed', dist_wirebonds='400um', dist_wirebonds_ends='150um'))
        quiet(cls.design.rebuild)

    def make_sim(self, cls=PALACE_Eigenmode_Simulation):
        return quiet(cls, name='api_test', metal_design=self.design, sim_parent_directory=self._dir + '/', mode='simPC',
                     meshing='GMSH', create_files=False, user_options={'mesh_max': 100e-6})

    def test_option_validation(self):
        sim = self.make_sim()
        with self.assertRaises(AssertionError):
            sim.add_wirebond([0, 0], [1e-4, 0], diameter=25e-6)
        with self.assertRaises(AssertionError):
            sim.add_airbridge([0, 0], [1e-4, 0], cross_section='round', foot_thickness=1e-6)
        with self.assertRaises(AssertionError):
            sim.set_wirebond_defaults(colour='gold')
        with self.assertRaises(AssertionError):
            sim.set_airbridge_defaults(offchip_termination='floor')

    def test_names(self):
        sim = self.make_sim()
        sim.add_wirebond([0, 0], [1e-4, 0])
        sim.add_airbridge([0, 0], [1e-4, 0])
        sim.add_airbridge([0, 0], [1e-4, 0])
        self.assertEqual([x['name'] for x in sim._wirebonds], ['wb_0', 'ab_0', 'ab_1'])
        self.assertEqual([x['kind'] for x in sim._wirebonds], ['wirebond', 'airbridge', 'airbridge'])
        with self.assertRaises(AssertionError):
            sim.add_wirebond([0, 0], [1e-4, 0], name='ab_0')
        for bad in ['metal_bond', 'x_neg', 'air_box', 'rf_port_1']:
            with self.assertRaises(AssertionError, msg=bad):
                sim.add_airbridge([0, 0], [1e-4, 0], name=bad)

    def test_resolve_defaults(self):
        sim = self.make_sim()
        sim.set_wirebond_defaults(width=20e-6)
        sim.set_airbridge_defaults(height=4e-6)
        sim.add_wirebond([0, 0], [1e-4, 0], name='round')
        sim.add_wirebond([0, 0], [1e-4, 0], name='ribbon', cross_section='ribbon', width=30e-6, thickness=4e-6)
        sim.add_wirebond([0, 0], [1e-4, 0], name='au', material='Au')
        sim.add_airbridge([0, 0], [1e-4, 0], name='ab', width=6e-6)
        sim.add_airbridge([0, 0], [1e-4, 0], name='ab_flat', profile='pedestal', foot_width=15e-6)
        res = {x['name']: x for x in sim._resolve_wirebonds()}
        self.assertEqual(res['round']['width'], 20e-6)
        np.testing.assert_allclose([res['round'][k] for k in ['foot_length', 'foot_width', 'foot_thickness']], [50e-6, 30e-6, 8e-6])
        np.testing.assert_allclose([res['ribbon'][k] for k in ['foot_length', 'foot_width', 'foot_thickness']], [60e-6, 30e-6, 4e-6])
        self.assertEqual(res['round']['mesh_max'], 100e-6)
        self.assertEqual(res['round']['boundary'], 'pec')
        self.assertEqual(res['au']['boundary'], 'conductivity')
        self.assertEqual(res['ab']['height'], 4e-6)
        self.assertEqual((res['ab']['foot_length'], res['ab']['foot_width']), (6e-6, 6e-6))
        self.assertEqual(res['ab']['mesh_min'], 3e-6)
        self.assertEqual(res['ab_flat']['foot_width'], 15e-6)
        self.assertEqual(res['ab']['kind'], 'airbridge')

    def test_resolve_rejects_bad_options(self):
        cases = [('wb', dict(cross_section='square')), ('wb', dict(width=0)), ('wb', dict(profile='zigzag')),
                 ('wb', dict(offchip_termination='pcb')), ('wb', dict(cross_section='ribbon', thickness=0)),
                 ('ab', dict(cross_section='round')), ('ab', dict(height=0)), ('ab', dict(thickness=0)),
                 ('ab', dict(profile='ball_wedge'))]
        for kind, opts in cases:
            sim = self.make_sim()
            (sim.add_wirebond if kind == 'wb' else sim.add_airbridge)([0, 0], [1e-4, 0], **opts)
            with self.assertRaises(AssertionError, msg=f"{kind} {opts}"):
                sim._resolve_wirebonds()
        sim = self.make_sim()
        sim.add_airbridge([0, 0], [1e-4, 0], cross_section='sheet', thickness=0)
        sim._resolve_wirebonds()

    def test_airbridge_specs_from_component(self):
        specs = self.design.components['abs'].get_airbridge_specs()
        self.assertEqual(len(specs), 5)
        #Feet centred on the I-bar pad sections: trace/2 + gap + pad_gap + bridge_pad_gap + bridge_pad_width/2 = 23.5um
        for s in specs:
            np.testing.assert_allclose(np.abs(s['pos_start'][1]), 23.5e-6, atol=1e-12)
            np.testing.assert_allclose(s['pos_end'][1], -s['pos_start'][1], atol=1e-12)
            np.testing.assert_allclose(s['pos_end'][0], s['pos_start'][0], atol=1e-12)
            self.assertAlmostEqual(s['width'], 8e-6)
            self.assertAlmostEqual(s['foot_length'], 5e-6)
            self.assertAlmostEqual(s['foot_width'], 15e-6)
            self.assertEqual(s['bond_type'], 'ground')
        np.testing.assert_allclose(np.diff([s['pos_start'][0] for s in specs]), 300e-6)

    def test_from_design(self):
        sim = self.make_sim()
        self.assertEqual(sim.add_airbridges_from_design(height=5e-6), 5)
        self.assertEqual(sim.add_wirebonds_from_design(), len(self.design.components['stitches'].wirebond_feet))
        self.assertTrue(all(x['height'] == 5e-6 for x in sim._wirebonds if x['kind'] == 'airbridge'))
        with self.assertRaises(AssertionError):
            sim.add_airbridges_from_design(component_names=['feed'])
        with self.assertRaises(AssertionError):
            sim.add_wirebonds_from_design(component_names=['abs'])


# ------------------------------------------------------------------ Mesh + config (GMSH only)

def physical_group_extents(msh_path, prefixes):
    """Returns {name: (xyz_min, xyz_max)} in metres for the 2D physical groups whose names start with the given prefixes."""
    gmsh.initialize()
    gmsh.option.setNumber('General.Verbosity', 0)
    try:
        gmsh.open(msh_path)
        tags, coords, _ = gmsh.model.mesh.getNodes()
        coords = np.array(coords).reshape(-1, 3)
        idx = {t: i for i, t in enumerate(tags)}
        ret = {}
        for dim, pg in gmsh.model.getPhysicalGroups(2):
            name = gmsh.model.getPhysicalName(dim, pg)
            if not name.startswith(tuple(prefixes)):
                continue
            nodes = set()
            for e in gmsh.model.getEntitiesForPhysicalGroup(dim, pg):
                _, _, node_tags = gmsh.model.mesh.getElements(dim, e)
                for x in node_tags:
                    nodes.update(x)
            pts = coords[[idx[t] for t in nodes]]*1e-3
            ret[name] = (pts.min(axis=0), pts.max(axis=0))
        return ret
    finally:
        gmsh.finalize()


class TestWireBondEigenmodeMesh(unittest.TestCase):
    """One eigenmode simulation with every kind of bond and bridge; the mesh/config are built once and then inspected."""

    @classmethod
    def setUpClass(cls):
        cls._dir = output_dir(cls.__name__)
        design = build_feedline_design()
        xL, yL = launch_pad_centre(design, 'lpL')
        WireBond(design, 'feed_bond', options=dict(start_x=f'{xL}mm', start_y=f'{yL}mm', end_x='-1.45mm', end_y='0mm',
                                                   profile='ball_wedge', height='150um'))
        WireBondsCPW(design, 'stitches', options=dict(pathObj='feed', dist_wirebonds='400um', dist_wirebonds_ends='150um'))
        WireBond(design, 'gp', options=dict(start_x='0.6mm', start_y='0.5mm', end_x='0.6mm', end_y='0.75mm', material='Au',
                                            profile='elliptical', height='120um'))
        AirBridgeCPW(design, 'abs', options=dict(pathObj='feed', dist_air_bridges='500um', dist_air_bridges_ends='300um'))
        quiet(design.rebuild)

        sim = quiet(PALACE_Eigenmode_Simulation, name='wb_eig', metal_design=design, sim_parent_directory=cls._dir + '/',
                    mode='simPC', meshing='GMSH', create_files=True,
                    user_options={'mesh_max': 200e-6, 'mesh_min': 20e-6, 'number_of_freqs': 2, 'starting_freq': 4e9})
        sim.add_metallic(1)
        sim.add_ground_plane()
        sim.set_xBoundary_as_proportion(0.2)
        sim.set_yBoundary_as_proportion(0.4)
        sim.set_wirebond_defaults(width=25e-6, offchip_termination='pad')
        sim.add_wirebonds_from_design(component_names=['feed_bond', 'stitches'])
        sim.add_wirebonds_from_design(component_names=['gp'], offchip_termination='floor')
        sim.add_wirebond([0.3e-3, -45e-6], [0.3e-3, 45e-6], name='ribbon', cross_section='ribbon', width=30e-6,
                         thickness=5e-6, profile='arc', height=40e-6)
        sim.set_airbridge_defaults(height=3e-6, thickness=0.5e-6)
        sim.add_airbridges_from_design()
        sim.add_airbridge([0.2e-3, -30e-6], [0.2e-3, 30e-6], name='ab_flat', profile='flat', foot_length=10e-6, foot_width=15e-6)
        sim.add_airbridge([-0.25e-3, -30e-6], [-0.25e-3, 30e-6], name='ab_sheet', cross_section='sheet')
        sim.setup_EPR_interfaces(metal_air=MaterialInterface('Aluminium-Vacuum'), substrate_air=MaterialInterface('Silicon-Vacuum'),
                                 substrate_metal=MaterialInterface('Silicon-Aluminium'))
        quiet(sim.prepare_simulation)
        cls.summary = {x['name']: x for x in sim.get_wirebond_summary(print_summary=False)}
        save_summary(sim, os.path.join(cls._dir, 'wb_eig', 'wirebond_summary.txt'))
        cls.num_cpw_bridges = len(design.components['abs'].get_airbridge_specs())
        with open(os.path.join(cls._dir, 'wb_eig', 'wb_eig.json')) as f:
            cls.config = json.load(f)
        cls.extents = physical_group_extents(os.path.join(cls._dir, 'wb_eig', 'wb_eig.msh'), ['wb_', 'ab_', 'ribbon'])

    def test_classification(self):
        s = self.summary
        self.assertEqual(s['wb_feed_bond']['bond_type'], 'trace_offchip')
        self.assertEqual(s['wb_gp']['bond_type'], 'ground_offchip')
        self.assertTrue(s['wb_gp']['grounds_to_enclosure'])
        self.assertEqual(s['ribbon']['bond_type'], 'ground_stitch')
        stitches = [x for x in s if x.startswith('wb_stitches_')]
        bridges = [x for x in s if x.startswith('ab_abs_')]
        self.assertGreater(len(stitches), 0)
        self.assertGreater(len(bridges), 0)
        self.assertEqual(len(bridges), self.num_cpw_bridges)
        for name in stitches + bridges + ['ab_flat', 'ab_sheet']:
            self.assertEqual(s[name]['bond_type'], 'ground_stitch', msg=name)
            self.assertEqual(s[name]['issues'], [], msg=name)
        self.assertEqual({s[x]['kind'] for x in bridges + ['ab_flat', 'ab_sheet']}, {'airbridge'})
        self.assertEqual({s[x]['kind'] for x in stitches + ['wb_feed_bond', 'wb_gp', 'ribbon']}, {'wirebond'})
        #All ground bonds join the ground conductor; the feed bond joins the trace
        self.assertEqual(len({s[x]['conductor_index'] for x in stitches + bridges + ['wb_gp', 'ribbon']}), 1)
        self.assertNotEqual(s['wb_feed_bond']['conductor_index'], s['ribbon']['conductor_index'])

    def test_config_boundaries(self):
        bnds = self.config['Boundaries']
        pec = bnds['PEC']['Attributes']
        cond = {x['Attributes'][0]: x for x in bnds.get('Conductivity', [])}
        self.assertEqual(len(pec), len(set(pec)))
        for name, x in self.summary.items():
            pg = x['physical_group']
            if name == 'wb_gp':
                self.assertNotIn(pg, pec)
                self.assertIn(pg, cond)
                self.assertGreater(cond[pg]['Conductivity'], 1e7)
            else:
                self.assertIn(pg, pec, msg=name)
                self.assertNotIn(pg, cond, msg=name)

    def test_epr_interfaces(self):
        bond_pgs = {x['physical_group'] for x in self.summary.values()}
        diel = self.config['Boundaries']['Postprocessing']['Dielectric']
        ma = [x for x in diel if x.get('Type', '') == 'MA']
        self.assertGreater(len(ma), 0)
        self.assertTrue(bond_pgs.issubset(set(sum([x['Attributes'] for x in ma], []))))
        for x in diel:
            if x.get('Type', '') != 'MA':
                self.assertTrue(bond_pgs.isdisjoint(x['Attributes']), msg=x.get('Type', ''))

    def test_mesh_has_every_bond(self):
        self.assertEqual(set(self.extents), set(self.summary))

    def test_airbridge_geometry(self):
        #Landing surface is z=0 (planar metals); solid bridges top out at height + thickness, the sheet at height
        for name, (lo, hi) in self.extents.items():
            if not name.startswith('ab_'):
                continue
            self.assertAlmostEqual(lo[2], 0, delta=1e-9, msg=name)
            self.assertAlmostEqual(hi[2], 3e-6 if name == 'ab_sheet' else 3.5e-6, delta=1e-9, msg=name)
        lo, hi = self.extents['ab_abs_0']
        self.assertAlmostEqual(hi[1] - lo[1], 52e-6, delta=1e-9)   #Feet centred at +/-23.5um, 5um long
        self.assertAlmostEqual(hi[0] - lo[0], 15e-6, delta=1e-9)   #I-bar pad sections are 15um wide
        lo, hi = self.extents['ab_flat']
        self.assertAlmostEqual(hi[1] - lo[1], 70e-6, delta=1e-9)
        self.assertAlmostEqual(hi[0] - lo[0], 15e-6, delta=1e-9)
        lo, hi = self.extents['ab_sheet']
        self.assertAlmostEqual(hi[1] - lo[1], 52e-6, delta=1e-9)   #Sheet runs between the inner edges of 8um feet
        self.assertAlmostEqual(hi[0] - lo[0], 8e-6, delta=1e-9)

    def test_wirebond_geometry(self):
        lo, hi = self.extents['ribbon']
        self.assertAlmostEqual(lo[2], 0, delta=1e-9)
        self.assertGreater(hi[2], 40e-6)
        self.assertLess(hi[2], 50e-6)
        #The floor-terminated bond reaches down to the air-box floor below the substrate
        lo, hi = self.extents['wb_gp']
        self.assertLess(lo[2], -100e-6)


class TestWireBondCapacitance(unittest.TestCase):
    """Capacitance terminals: a jumper merges two conductors and a bond to the enclosure floor grounds a conductor."""

    def run_sim(self, with_bonds):
        design = build_feedline_design(island=True)
        name = 'wb_cap_bonds' if with_bonds else 'wb_cap_no_bonds'
        sim = quiet(PALACE_Capacitance_Simulation, name=name, metal_design=design, sim_parent_directory=self._dir + '/',
                    mode='simPC', meshing='GMSH', create_files=True, user_options={'mesh_max': 200e-6, 'mesh_min': 20e-6})
        sim.add_metallic(1)
        sim.add_ground_plane()
        sim.set_xBoundary_as_proportion(0.2)
        sim.set_yBoundary_as_proportion(0.4)
        if with_bonds:
            xL, yL = launch_pad_centre(design, 'lpL')
            sim.add_wirebond([xL*1e-3, yL*1e-3], [0, 0.35e-3], name='jumper', height=150e-6)
            sim.add_wirebond([0.6e-3, 0.5e-3], [0.6e-3, 0.75e-3], name='to_floor', offchip_termination='floor')
            sim.add_airbridge([0.5e-3, -30e-6], [0.5e-3, 30e-6], name='bridge')
        quiet(sim.prepare_simulation)
        if with_bonds:
            save_summary(sim, os.path.join(self._dir, name, 'wirebond_summary.txt'))
        with open(os.path.join(self._dir, name, name + '.json')) as f:
            bnds = json.load(f)['Boundaries']
        return sim, bnds

    def setUp(self):
        self._dir = output_dir(type(self).__name__, self._testMethodName)

    def test_terminals(self):
        _, bnds = self.run_sim(with_bonds=False)
        self.assertEqual(len(bnds['Terminal']), 3)      #Ground plane, feedline and island

        sim, bnds = self.run_sim(with_bonds=True)
        s = {x['name']: x for x in sim.get_wirebond_summary(print_summary=False)}
        self.assertEqual(s['jumper']['bond_type'], 'trace_jumper')
        self.assertEqual(s['to_floor']['bond_type'], 'ground_offchip')
        self.assertEqual(s['bridge']['bond_type'], 'ground_stitch')
        self.assertEqual(s['bridge']['conductor_index'], s['to_floor']['conductor_index'])
        #Feedline + island merged into one terminal; the ground plane (with its bonds) is moved into Ground
        self.assertEqual(len(bnds['Terminal']), 1)
        self.assertEqual(bnds['Terminal'][0]['Index'], 1)
        term = bnds['Terminal'][0]['Attributes']
        ground = bnds['Ground']['Attributes']
        self.assertIn(s['jumper']['physical_group'], term)
        for name in ['to_floor', 'bridge']:
            self.assertIn(s[name]['physical_group'], ground)
            self.assertNotIn(s[name]['physical_group'], term)
        self.assertEqual(len(sim._cur_cap_terminals), 1)


class TestTrenchedAirBridges(unittest.TestCase):
    """Substrate trenching with lumped ports on the launch pads and 3D air-bridges. Regression test: the trench walls used to
    end exactly on the port edges, leaving unmerged vertices at the port corners ('1D mesh not forming a closed loop')."""

    def test_trenched_mesh(self):
        out = output_dir(type(self).__name__, self._testMethodName)
        design = build_feedline_design()
        AirBridgeCPW(design, 'abs', options=dict(pathObj='feed', layer=4, layer_bridge=5, dist_air_bridges='500um', dist_air_bridges_ends='300um'))
        quiet(design.rebuild)
        sim = quiet(PALACE_Eigenmode_Simulation, name='wb_trench', metal_design=design, sim_parent_directory=out + '/',
                    mode='simPC', meshing='GMSH', create_files=True,
                    user_options={'mesh_max': 200e-6, 'mesh_min': 20e-6, 'number_of_freqs': 1, 'starting_freq': 4e9})
        sim.add_metallic(1)
        sim.add_ground_plane()
        sim.create_port_CPW_on_Launcher('lpL', 20e-6)
        sim.create_port_CPW_on_Launcher('lpR', 20e-6)
        sim.enforce_full_3D_simulation(0, substrate_trenching=10e-6)
        #Fine enough at the launch pads to put mesh nodes inside the 20um x 60um ports (checked below)
        sim.fine_mesh_components(['lpL', 'lpR'], min_size=5e-6, max_size=200e-6, taper_dist_min=10e-6, metals_only=False)
        sim.set_airbridge_defaults(cross_section='sheet', height=3e-6)
        num_bridges = sim.add_airbridges_from_design()
        self.assertGreater(num_bridges, 0)
        sim.setup_EPR_interfaces(metal_air=MaterialInterface('Aluminium-Vacuum'), substrate_air=MaterialInterface('Silicon-Vacuum'),
                                 substrate_metal=MaterialInterface('Silicon-Aluminium'))
        quiet(sim.prepare_simulation)
        save_summary(sim, os.path.join(out, 'wb_trench', 'wirebond_summary.txt'))

        with open(os.path.join(out, 'wb_trench', 'wb_trench.json')) as f:
            bnds = json.load(f)['Boundaries']
        self.assertEqual(len(bnds['LumpedPort']), 2)
        summary = sim.get_wirebond_summary(print_summary=False)
        self.assertEqual(len(summary), num_bridges)
        for x in summary:
            self.assertIn(x['physical_group'], bnds['PEC']['Attributes'])

        msh = os.path.join(out, 'wb_trench', 'wb_trench.msh')
        ext = physical_group_extents(msh, ['ab_', 'rf_port_'])
        self.assertEqual(len([x for x in ext if x.startswith('ab_')]), num_bridges)
        self.assertEqual(len([x for x in ext if x.startswith('rf_port_')]), 4)
        for name, (lo, hi) in ext.items():
            self.assertAlmostEqual(lo[2], 0, delta=1e-9, msg=name)
        #The substrate-air EPR interface is the trench floors and walls - it must not reference empty physical groups
        #(Palace: 'Unknown boundary postprocessing attributes')
        sa = [x for x in bnds['Postprocessing']['Dielectric'] if x['Type'] == 'SA'][0]['Attributes']
        gaps = physical_group_extents(msh, ['dielectric_gaps'])
        self.assertEqual(len(gaps), 1)
        lo, hi = gaps['dielectric_gaps']
        self.assertAlmostEqual(lo[2], -10e-6, delta=1e-9)
        self.assertAlmostEqual(hi[2], 0, delta=1e-9)
        gmsh.initialize()
        gmsh.option.setNumber('General.Verbosity', 0)
        try:
            gmsh.open(msh)
            for attr in sa:
                ents = gmsh.model.getEntitiesForPhysicalGroup(2, attr)
                self.assertGreater(sum(sum(len(t) for t in gmsh.model.mesh.getElements(2, e)[1]) for e in ents), 0, msg=f"SA attribute {attr}")
        finally:
            gmsh.finalize()
        #The substrate is trenched (10um) under the CPW gaps, including under the lumped ports
        gmsh.initialize()
        gmsh.option.setNumber('General.Verbosity', 0)
        try:
            gmsh.open(msh)
            sub = [pg for d, pg in gmsh.model.getPhysicalGroups(3) if gmsh.model.getPhysicalName(d, pg) == 'dielectric_substrate'][0]
            zmax_under_port = []
            for e in gmsh.model.getEntitiesForPhysicalGroup(3, sub):
                tags, coords, _ = gmsh.model.mesh.getNodes(3, e, includeBoundary=True)
                c = np.array(coords).reshape(-1, 3)
                lo, hi = ext['rf_port_0a']
                sel = (c[:, 0] > lo[0]*1e3 + 1e-4) & (c[:, 0] < hi[0]*1e3 - 1e-4) & (c[:, 1] > lo[1]*1e3 + 1e-4) & (c[:, 1] < hi[1]*1e3 - 1e-4)
                zmax_under_port += c[sel, 2].tolist()
        finally:
            gmsh.finalize()
        #No substrate above the trench floor under the port (it would sit at z=0 without trenching)
        self.assertGreater(len(zmax_under_port), 0)
        self.assertLessEqual(max(zmax_under_port)*1e-3, -10e-6 + 1e-9)


if __name__ == '__main__':
    unittest.main()
