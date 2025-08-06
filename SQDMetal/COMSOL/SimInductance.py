from SQDMetal.COMSOL.Model import COMSOL_Simulation_Base
from SQDMetal.Utilities.QUtilities import QUtilities
from SQDMetal.Utilities.QiskitShapelyRenderer import QiskitShapelyRenderer

import jpype.types as jtypes
import shapely
import numpy as np

class COMSOL_Simulation_Inductance(COMSOL_Simulation_Base):

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
        self.jc.component("comp1").physics(self.phys_ecis).feature("csh1").feature("bterm1").set("I0", jtypes.JFloat(1.0))
        self.jc.component("comp1").physics(self.phys_ecis).feature("csh1").feature("bterm1").selection().named(self._loop_closure_start_end[0])
        self.jc.component("comp1").physics(self.phys_ecis).feature("csh1").create("bgnd1", "BoundaryGround", 1)
        self.jc.component("comp1").physics(self.phys_ecis).feature("csh1").feature("bgnd1").selection().named(self._loop_closure_start_end[1])

        self.phys_mf = self.model._add_physics("mf")
        self.jc.component("comp1").physics().create(self.phys_mf, "InductionCurrents", "geom1")
        self.jc.component("comp1").physics(self.phys_mf).create("scu1", "SurfaceCurrent", 2)
        self.jc.component("comp1").physics(self.phys_mf).feature("scu1").selection().named("geom1_condAll")
        self.jc.component("comp1").physics(self.phys_mf).feature("scu1").set("Js0", jtypes.JArray(jtypes.JString)(["ecis.JsX", "ecis.JsY", "ecis.JsZ"]))
        self.jc.component("comp1").physics(self.phys_mf).create("mi2", "MagneticInsulation", 2)
        self.jc.component("comp1").physics(self.phys_mf).feature("mi2").selection().named(f"geom1_{self._sel_loop_close}")

        self.jc.component("comp1").material("Metal").selection().named("geom1_condPlusUBend")

    def _get_required_physics(self):
        return [self.phys_ecis, self.phys_mf]

    def set_return_loop_on_route(self, qObjName, **kwargs):
        thresh = kwargs.get('threshold', -1)
        sel_feed_rad = kwargs.get('feed_edge_min_size', 1e-9)

        #Presumes that it is simply-closed polygon
        qmpl = QiskitShapelyRenderer(None, self.model.design, None)
        gsdf = qmpl.get_net_coordinates()
        comp_id = self.model.design.components[qObjName].id
        metal_polys_all = gsdf.loc[(gsdf["component"] == comp_id) & (~gsdf["subtract"])]['geometry'].iloc[0]
        #
        unit_conv = kwargs.get("unit_conv", QUtilities.get_units(self.model.design))
        metal_polys_all = shapely.affinity.scale( metal_polys_all, xfact=unit_conv, yfact=unit_conv, origin=(0, 0) )

        pol_name = "loopClose"+str(self.model._num_polys)
        #Convert coordinates into metres...
        cur_poly = np.array(metal_polys_all.exterior.coords[:])
        cur_poly = self.model._simplify_geom(cur_poly, thresh)
        sel_x, sel_y, sel_r = self.model._create_poly(pol_name, cur_poly[:-1])
        poly_interiors = metal_polys_all.interiors
        if len(poly_interiors) > 0:
            pol_name_ints = []
            for ind,cur_int in enumerate(poly_interiors):
                pol_name_int = "loopClose"
                pol_name_ints.append(pol_name_int)
                #Convert coordinates into metres...
                cur_poly_int = np.array(cur_int.coords[:])
                cur_poly_int = self.model._simplify_geom(cur_poly_int, thresh)
                sel_x, sel_y, sel_r = self.model._create_poly(pol_name_int, cur_poly_int[:-1])
            #Subtract interiors from the polygon
            # abhishekchak52: undefined variable m
            diff_name = f"difInt{self.model._num_polys+m}"
            self.model._model.java.component("comp1").geom("geom1").feature("wp1").geom().create(diff_name, "Difference")
            self.model._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(diff_name).selection("input").set(pol_name)
            self.model._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(diff_name).selection("input2").set(*pol_name_ints)
            select_obj_name = diff_name
        else:
            select_obj_name = pol_name

            self.loop_closure = select_obj_name
        
        start_pt = QUtilities._get_Route_params(self.model.design, qObjName, 'start')[0]
        end_pt = QUtilities._get_Route_params(self.model.design, qObjName, 'end')[0]
        self._loop_closure_start_end = ( self.model._create_boundary_selection_sphere(sel_feed_rad, *start_pt, is_edge=True),
                                         self.model._create_boundary_selection_sphere(sel_feed_rad, *end_pt, is_edge=True) )

        cur_sel_index = len(self.model._conds)
        self._sel_loop_close = self.model._setup_selection_boundaries(cur_sel_index, select_obj_name, sel_3D_prefix="selLoopClose")

    def run(self):
        self.jc.result().numerical().create("gev1", "EvalGlobal")
        self.jc.result().numerical("gev1").setIndex("expr", "2*mf.intWm/ecis.I0_1^2", 0)
        self.jc.result().numerical("gev1").setIndex("expr", "ecis.V0_1/ecis.I0_1", 1)

        self.jc.sol(self._soln).runAll()
        #
        res = np.array(self.jc.result().numerical("gev1").computeResult())
        res = res[0][0]

        inductance, resistance = res

        return inductance, resistance

    def _get_Route_params(self, route_name, pin_name):
        design = self.model.design
        unit_conv = QUtilities.get_units(design)

        startPt = design.components[route_name].pins[pin_name]['middle']*unit_conv
        padDir = -1.0*design.components[route_name].pins[pin_name]['normal']
        padWid = QUtilities.parse_value_length(design.components[route_name].options.trace_width)
        padGap = QUtilities.parse_value_length(design.components[route_name].options.trace_gap)

        return startPt, padDir, padWid, padGap
