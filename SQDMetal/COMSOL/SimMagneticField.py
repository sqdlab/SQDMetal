from SQDMetal.COMSOL.Model import COMSOL_Simulation_Base
from SQDMetal.Utilities.QUtilities import QUtilities

import mph
import jpype.types as jtypes
import geopandas as gpd
import shapely
import matplotlib.pyplot as plt
import numpy as np

class COMSOL_Simulation_MagneticField(COMSOL_Simulation_Base):
    
    #terminal names
    terminal_names = []
    
    def __init__(self, model):
        self.model = model
        self.jc = model._get_java_comp()
        self.dset_name = ""
        self.phys_ecis = ""
        self.phys_mf = ""
        self._study = ""
        self._sub_study = "capMatStationary"
        self._soln = ""
        self.current_feed_UBend = None
        self._bfield_int = None

    def _prepare_simulation(self):
        self.jc.param().set("PortName", jtypes.JInt(1))    #Just use default port name...
        self.jc.study().create(self._study)
        self.jc.study(self._study).create(self._sub_study,"Stationary")
        self._soln = self.model._add_solution("solBField")
        self.jc.sol().create(self._soln)
        self.jc.sol(self._soln).createAutoSequence(self._study)
        self.dset_name = self.model._get_dset_name()

    def _run_premesh(self):
        select_3D_name = "condPlusUBend"
        self.model._model.java.component("comp1").geom("geom1").create(select_3D_name, "UnionSelection")
        self.model._model.java.component("comp1").geom("geom1").feature(select_3D_name).label(select_3D_name)
        self.model._model.java.component("comp1").geom("geom1").feature(select_3D_name).set("entitydim", jtypes.JInt(2))
        conds = ["condAll"]
        if self.current_feed_UBend:
            conds += ["condAll", self.current_feed_UBend]
        self.model._model.java.component("comp1").geom("geom1").feature(select_3D_name).set("input", jtypes.JArray(jtypes.JString)(conds))
        self.model._model.java.component("comp1").geom("geom1").run("condAll")

        self._study = self.model._add_study("stdCapMat")

        self.phys_ecis = self.model._add_physics("ecis")
        self.jc.component("comp1").physics().create(self.phys_ecis, "ElectricCurrentsShell", "geom1")
        self.jc.component("comp1").physics(self.phys_ecis).selection().named("geom1_condPlusUBend")
        self.jc.component("comp1").physics(self.phys_ecis).feature("csh1").create("bterm1", "BoundaryTerminal", 1)
        self.jc.component("comp1").physics(self.phys_ecis).feature("csh1").feature("bterm1").selection().named(self._sel_term_src)
        self.jc.component("comp1").physics(self.phys_ecis).feature("csh1").feature("bterm1").set("I0", jtypes.JFloat(1.0))
        self.jc.component("comp1").physics(self.phys_ecis).feature("csh1").create("bgnd1", "BoundaryGround", 1)
        self.jc.component("comp1").physics(self.phys_ecis).feature("csh1").feature("bgnd1").selection().named(self._sel_term_gnd)

        self.phys_mf = self.model._add_physics("mf")
        self.jc.component("comp1").physics().create(self.phys_mf, "InductionCurrents", "geom1")
        self.jc.component("comp1").physics(self.phys_mf).create("scu1", "SurfaceCurrent", 2)
        self.jc.component("comp1").physics(self.phys_mf).feature("scu1").selection().named("geom1_condPlusUBend")
        self.jc.component("comp1").physics(self.phys_mf).feature("scu1").set("Js0", jtypes.JArray(jtypes.JString)(["ecis.JsX", "ecis.JsY", "ecis.JsZ"]))
        self.jc.component("comp1").physics(self.phys_mf).create("mi2", "MagneticInsulation", 2)
        self.jc.component("comp1").physics(self.phys_mf).feature("mi2").selection().named(f"geom1_{self.current_feed_Pad}")

        self.jc.component("comp1").material("Metal").selection().named("geom1_condPlusUBend")

    def _get_required_physics(self):
        return [self.phys_ecis, self.phys_mf]

    def set_current_feed_on_CPW_on_Route(self, qObjName, pin_name='end', width_U = 20e-6, src_gnd_gap=10e-6):
        assert self.current_feed_UBend == None, "Cannot have more than one current feed..."

        qObj = self.model.design.components[qObjName]
        # if isinstance(qObj, LaunchpadWirebond):
        vec_ori, vec_launch, cpw_wid, cpw_gap = self._get_Route_params(qObjName, pin_name)
        # else:
        #     assert False, f"\'{qObjName}\' is an unsupported object type."

        vec_perp = np.array([-vec_launch[1],vec_launch[0]])
        vec_launch = -vec_launch

        pol_name = "current_feed_"
        
        uBend = [vec_ori + vec_perp*(cpw_wid*0.5+cpw_gap),
                 vec_ori + vec_perp*(cpw_wid*0.5+cpw_gap) - vec_launch*src_gnd_gap,
                 vec_ori + vec_perp*(cpw_wid*0.5+cpw_gap) - vec_launch*src_gnd_gap - vec_perp*(cpw_wid+2*cpw_gap),
                 vec_ori - vec_perp*(cpw_wid*0.5+cpw_gap),
                 vec_ori - vec_perp*(cpw_wid*0.5+cpw_gap+width_U),
                 vec_ori - vec_perp*(cpw_wid*0.5+cpw_gap+width_U) - vec_launch*(src_gnd_gap+width_U),
                 vec_ori + vec_perp*(cpw_wid*0.5+cpw_gap+width_U) - vec_launch*(src_gnd_gap+width_U),
                 vec_ori + vec_perp*(cpw_wid*0.5+cpw_gap+width_U)]
        uBend = [[p[0],p[1]] for p in uBend]
        self.model._create_poly(pol_name + "uBend", uBend)
        self.current_feed_UBend = self.model._setup_selection_boundaries("FeedLnU", pol_name + "uBend")
        self._sel_term_src = self.model._create_boundary_selection_sphere(min(cpw_wid, width_U)*0.1, vec_ori[0], vec_ori[1], is_edge = True)

        feedPad = [vec_ori + vec_perp*(cpw_wid*0.5),
                   vec_ori - vec_perp*(cpw_wid*0.5),
                   vec_ori - vec_perp*(cpw_wid*0.5) - vec_launch*src_gnd_gap,
                   vec_ori + vec_perp*(cpw_wid*0.5) - vec_launch*src_gnd_gap]
        feedPad = [[p[0],p[1]] for p in feedPad]
        self.model._create_poly(pol_name + "feedPad", feedPad)
        self.current_feed_Pad = self.model._setup_selection_boundaries("FeedLnP", pol_name + "feedPad")
        vec_gnd = vec_ori - vec_launch*src_gnd_gap
        self._sel_term_gnd = self.model._create_boundary_selection_sphere(min(cpw_wid, width_U)*0.1, vec_gnd[0], vec_gnd[1], is_edge = True)

    def set_current_feed_point_point_verticalU(self, x1, y1, x2, y2, uWidth):
        wpName = "wpCurrentFeed"
        self.jc.component("comp1").geom("geom1").feature().create(wpName, "WorkPlane")
        self.jc.component("comp1").geom("geom1").feature(wpName).set("unite", jtypes.JBoolean(True))
        self.jc.component("comp1").geom("geom1").feature(wpName).set("planetype", "normalvector")
        self.jc.component("comp1").geom("geom1").feature(wpName).set("normalvector", jtypes.JArray(jtypes.JDouble)([y1-y2, x2-x1, 0.0]))
        self.jc.component("comp1").geom("geom1").feature(wpName).set("normalcoord", jtypes.JArray(jtypes.JDouble)([(x1+x2)*0.5, (y1+y2)*0.5, 0.0]))
        self.jc.component("comp1").geom("geom1").feature(wpName).set("rot", jtypes.JDouble(90))

        dist = np.sqrt((x2-x1)**2+(y2-y1)**2)
        poly_coords = [[dist/2+uWidth/2, 0],
                       [dist/2+uWidth/2, uWidth*2],
                       [-dist/2-uWidth/2, uWidth*2],
                       [-dist/2-uWidth/2, 0],
                       [-dist/2+uWidth/2, 0],
                       [-dist/2+uWidth/2, uWidth],
                       [dist/2-uWidth/2, uWidth],
                       [dist/2-uWidth/2, 0]
                       ]

        polyUname = "currentFeedU"
        selPolyNameWP = "cSelCurrentFeedU"
        self.jc.component("comp1").geom("geom1").feature(wpName).geom().create(polyUname, "Polygon")
        self.jc.component("comp1").geom("geom1").feature(wpName).geom().feature(polyUname).set('selresult', 'on') #To generate automatic selections...
        self.jc.component("comp1").geom("geom1").feature(wpName).geom().feature(polyUname).set('selresultshow', 'bnd')
        self.jc.component("comp1").geom("geom1").feature(wpName).geom().feature(polyUname).set("source", "table")
        self.jc.component("comp1").geom("geom1").feature(wpName).geom().feature(polyUname).set('table', jtypes.JArray(jtypes.JDouble,2)(poly_coords))    #This is much faster...
        #
        self.jc.component("comp1").geom("geom1").feature(wpName).geom().selection().create(selPolyNameWP, "CumulativeSelection")
        self.jc.component("comp1").geom("geom1").feature(wpName).geom().feature(polyUname).set("contributeto", selPolyNameWP)
        #
        select_3D_name = f"selCurrentFeedVerticalU"
        self.jc.component("comp1").geom("geom1").create(select_3D_name, "UnionSelection")
        self.jc.component("comp1").geom("geom1").feature(select_3D_name).label(select_3D_name)
        self.jc.component("comp1").geom("geom1").feature(select_3D_name).set("entitydim", jtypes.JInt(2))
        self.jc.component("comp1").geom("geom1").feature(select_3D_name).set("input", f"{wpName}_{selPolyNameWP}")

        self._sel_term_src = self.model._create_boundary_selection_sphere(uWidth*0.1, x1, y1, is_edge=True)
        self._sel_term_gnd = self.model._create_boundary_selection_sphere(uWidth*0.1, x2, y2, is_edge=True)
        self.current_feed_UBend = None
        self.current_feed_Pad = select_3D_name

    def set_Bfield_integration_area(self, coords):
        uBend = [[p[0],p[1]] for p in coords]
        self.model._create_poly("bfieldIntArea", uBend)
        self._bfield_int = self.model._setup_selection_boundaries("BFieldIntegration", "bfieldIntArea")

    def display_conductor_indices(self):
        '''
        Plots a coloured visualisation of the metallic conductors and their corresponding row/column indices of the capacitance matrix.
        '''
        minX = (self.model.chip_centre[0] - self.model.chip_len*0.5)
        maxX = (self.model.chip_centre[0] + self.model.chip_len*0.5)
        minY = (self.model.chip_centre[1] - self.model.chip_wid*0.5)
        maxY = (self.model.chip_centre[1] + self.model.chip_wid*0.5)
        chip_bounding_poly = shapely.LineString([(minX,minY),(maxX,minY),(maxX,maxY),(minX,maxY),(minX,minY)])
        
        leGeoms = [chip_bounding_poly] + [x[0] for x in self.model._cond_polys]
        leNames = ["Chip"] + [f"Cond{x[1]}" for x in self.model._cond_polys]
        gdf = gpd.GeoDataFrame({'names':leNames}, geometry=leGeoms)
        fig, ax = plt.subplots(1)
        gdf.plot(ax = ax, column='names', cmap='jet', alpha=0.5, categorical=True, legend=True)
        ax.set_xlabel(f'Position (m)')
        ax.set_ylabel(f'Position (m)')

    def run(self):
        assert self._bfield_int != None, "Must define an integration surface"

        self.jc.result().numerical().create("int1", "IntSurface")
        self.jc.result().numerical("int1").set("intvolume", jtypes.JBoolean(True))
        self.jc.result().numerical("int1").selection().named("geom1_" + self._bfield_int)
        self.jc.result().numerical("int1").set("expr", jtypes.JArray(jtypes.JString)(["mf.Bz"]))
        self.jc.result().numerical("int1").set("descr", jtypes.JArray(jtypes.JString)(["Magnetic flux density, z-component"]))
        self.jc.result().numerical("int1").set("unit", jtypes.JArray(jtypes.JString)(["Wb"]))
        #
        self.jc.sol(self._soln).runAll()
        #
        flux = np.array(self.jc.result().numerical("int1").computeResult())[0]
        self.jc.result().numerical().remove("int1")
    
        return flux[0]

    def _get_Route_params(self, route_name, pin_name):
        design = self.model.design
        unit_conv = QUtilities.get_units(design)

        startPt = design.components[route_name].pins[pin_name]['middle']*unit_conv
        padDir = -1.0*design.components[route_name].pins[pin_name]['normal']
        padWid = QUtilities.parse_value_length(design.components[route_name].options.trace_width)
        padGap = QUtilities.parse_value_length(design.components[route_name].options.trace_gap)

        return startPt, padDir, padWid, padGap
