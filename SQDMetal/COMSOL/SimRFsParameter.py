from SQDMetal.COMSOL.Model import COMSOL_Simulation_Base
from SQDMetal.Utilities.QUtilities import QUtilities
from SQDMetal.Utilities.ShapelyEx import ShapelyEx

import jpype.types as jtypes
import shapely
import numpy as np
import scipy.optimize

import os


class COMSOL_Simulation_RFsParameters(COMSOL_Simulation_Base):
    def __init__(self, model, adaptive='None', modal_min_freq_num_eigs = (1e9,7), relative_tolerance=0.01, include_loss_tangents=False):
        self.model = model
        self.jc = model._get_java_comp()
        self.dset_name = ""
        self.phys_emw = ""
        self._study = ""
        if adaptive == 'Multiple':
            self._sub_study = ["eig", "frmod"]
        else:
            self._sub_study = "freq"
        self._soln = ""
        #
        self._ports = []
        self.lumped_elems = []
        assert adaptive == 'None' or adaptive == 'Single' or adaptive == 'Multiple', "Adaptive Simulation must be None, Single or Multiple"
        self.adaptive = adaptive
        self.modal_min_freq_num_eigs = modal_min_freq_num_eigs
        self.relative_tolerance = relative_tolerance
        self._include_loss_tangents = include_loss_tangents
        self._current_sources = []
        self._electric_point_dipoles = []
        self._discard_PEC_BC = False

    def _prepare_simulation(self):
        self._study = self.model._add_study("stdRFsparams")
        self.jc.study().create(self._study)
        if self.adaptive == 'None':
            self.jc.study(self._study).create(self._sub_study, "Frequency")
            self.jc.study(self._study).feature(self._sub_study).set("solnum", "auto")
            self.jc.study(self._study).feature(self._sub_study).set("notsolnum", "auto")
            self.jc.study(self._study).feature(self._sub_study).set("savesolsref", jtypes.JBoolean(False))
            self.jc.study(self._study).feature(self._sub_study).set("ngen", "5")
        elif self.adaptive == 'Single':
            self.jc.study(self._study).create(self._sub_study, "FrequencyAdaptive")
            self.jc.study(self._study).feature(self._sub_study).set("rtol", self.relative_tolerance)
        elif self.adaptive == 'Multiple':
            self.jc.study(self._study).create(self._sub_study[0], "Eigenfrequency")
            self.jc.study(self._study).feature(self._sub_study[0]).set("eigwhich", "lr")
            self.jc.study(self._study).feature(self._sub_study[0]).set("conrad", "1")
            self.jc.study(self._study).feature(self._sub_study[0]).activate("emw", jtypes.JBoolean(True))
            self.jc.study(self._study).create(self._sub_study[1], "Frequencymodal")
            self.jc.study(self._study).feature(self._sub_study[1]).set("solnum", "auto")
            self.jc.study(self._study).feature(self._sub_study[1]).set("notsolnum", "auto")
            self.jc.study(self._study).feature(self._sub_study[1]).activate("emw", jtypes.JBoolean(True))
            self.jc.study(self._study).feature(self._sub_study[0]).set("shift", f"{self.modal_min_freq_num_eigs[0]*1e-9}[GHz]")
            self.jc.study(self._study).feature(self._sub_study[0]).set("neigsactive", jtypes.JBoolean(True))
            self.jc.study(self._study).feature(self._sub_study[0]).set("neigs", f"{self.modal_min_freq_num_eigs[1]}")
            #Only store the port solutions for Multi-Modal...
            lst_ports = ["geom1_"+cur_port['polys'][0] for cur_port in self._ports]
            self.jc.study(self._study).feature(self._sub_study[1]).set("usestoresel", "selection")
            self.jc.study(self._study).feature(self._sub_study[1]).set("storesel", jtypes.JArray(jtypes.JString)(lst_ports))
        self._soln = self.model._add_solution("solRFsparams")
        self.jc.sol().create(self._soln)
        if True:#self.adaptive != 'None':
            self.jc.sol(self._soln).createAutoSequence(self._study)
        else:
            self.jc.sol(self._soln).study(self._study)
            self._temporary_freq_domain_AutoSequence()
        # if self.adaptive == 'None':
        #     self.jc.sol(self._soln).feature("s1").set("stol", self.relative_tolerance)
        self.dset_name = self.model._get_dset_name()
        
        #Setup surface integrals (e.g. for EPR calculations)
        int_name = "intMetals"
        self.jc.result().numerical().create(int_name, "IntSurface")
        self.jc.result().numerical(int_name).set("intvolume", jtypes.JBoolean(True))
        self.jc.result().numerical(int_name).selection().named("geom1_condAll")
        self.jc.result().numerical(int_name).set("expr", jtypes.JArray(jtypes.JString)(["emw.Ez"]))
        self.jc.result().numerical(int_name).set("descr", jtypes.JArray(jtypes.JString)(["Electric field, z-component"]))
        #
        int_name = "intDielectric"
        self.jc.result().numerical().create(int_name, "IntSurface")
        self.jc.result().numerical(int_name).set("intvolume", jtypes.JBoolean(True))
        self.jc.result().numerical(int_name).selection().named("geom1_dielectricSurface")
        self.jc.result().numerical(int_name).set("expr", jtypes.JArray(jtypes.JString)(["emw.Ez"]))
        self.jc.result().numerical(int_name).set("descr", jtypes.JArray(jtypes.JString)(["Electric field, z-component"]))


    def _temporary_freq_domain_AutoSequence(self):
        self.jc.sol(self._soln).create("st1", "StudyStep")
        self.jc.sol(self._soln).feature("st1").set("study", self._study)
        self.jc.sol(self._soln).feature("st1").set("studystep", "freq")
        self.jc.sol(self._soln).create("v1", "Variables")
        self.jc.sol(self._soln).feature("v1").set("control", "freq")
        self.jc.sol(self._soln).create("s1", "Stationary")
        self.jc.sol(self._soln).feature("s1").set("stol", self.relative_tolerance)
        self.jc.sol(self._soln).feature("s1").create("p1", "Parametric")
        self.jc.sol(self._soln).feature("s1").feature().remove("pDef")
        self.jc.sol(self._soln).feature("s1").feature("p1").set("pname", jtypes.JArray(jtypes.JString)(["freq"]))
        self.jc.sol(self._soln).feature("s1").feature("p1").set("plistarr", jtypes.JArray(jtypes.JString)(["1.0, 2.0, 3.0"]))
        self.jc.sol(self._soln).feature("s1").feature("p1").set("punit", jtypes.JArray(jtypes.JString)(["GHz"]))
        self.jc.sol(self._soln).feature("s1").feature("p1").set("pcontinuationmode", "no")
        self.jc.sol(self._soln).feature("s1").feature("p1").set("preusesol", "auto")
        self.jc.sol(self._soln).feature("s1").feature("p1").set("pdistrib", "off")
        self.jc.sol(self._soln).feature("s1").feature("p1").set("plot", "on")
        self.jc.sol(self._soln).feature("s1").feature("p1").set("plotgroup", "Default")
        self.jc.sol(self._soln).feature("s1").feature("p1").set("probesel", "all")
        self.jc.sol(self._soln).feature("s1").feature("p1").set("probes", jtypes.JArray(jtypes.JString)([]))
        self.jc.sol(self._soln).feature("s1").feature("p1").set("control", "freq")
        self.jc.sol(self._soln).feature("s1").set("linpmethod", "sol")
        self.jc.sol(self._soln).feature("s1").set("linpsol", "zero")
        self.jc.sol(self._soln).feature("s1").set("control", "freq")
        self.jc.sol(self._soln).feature("s1").feature("aDef").set("complexfun", jtypes.JBoolean(True))
        self.jc.sol(self._soln).feature("s1").create("fc1", "FullyCoupled")
        self.jc.sol(self._soln).feature("s1").create("i1", "Iterative")
        self.jc.sol(self._soln).feature("s1").feature("i1").set("linsolver", "gmres")
        self.jc.sol(self._soln).feature("s1").feature("i1").set("prefuntype", "right")
        self.jc.sol(self._soln).feature("s1").feature("i1").set("itrestart", "300")
        self.jc.sol(self._soln).feature("s1").feature("i1").label("Suggested Iterative Solver (emw)")
        self.jc.sol(self._soln).feature("s1").feature("i1").create("mg1", "Multigrid")
        self.jc.sol(self._soln).feature("s1").feature("i1").feature("mg1").set("iter", "1")
        self.jc.sol(self._soln).feature("s1").feature("i1").feature("mg1").feature("pr").create("sv1", "SORVector")
        self.jc.sol(self._soln).feature("s1").feature("i1").feature("mg1").feature("pr").feature("sv1").set("prefun", "sorvec")
        self.jc.sol(self._soln).feature("s1").feature("i1").feature("mg1").feature("pr").feature("sv1").set("iter", jtypes.JInt(2))
        self.jc.sol(self._soln).feature("s1").feature("i1").feature("mg1").feature("pr").feature("sv1").set("relax", jtypes.JInt(1))
        self.jc.sol(self._soln).feature("s1").feature("i1").feature("mg1").feature("pr").feature("sv1").set("sorvecdof", jtypes.JArray(jtypes.JString)(["comp1_E"]))
        self.jc.sol(self._soln).feature("s1").feature("i1").feature("mg1").feature("po").create("sv1", "SORVector")
        self.jc.sol(self._soln).feature("s1").feature("i1").feature("mg1").feature("po").feature("sv1").set("prefun", "soruvec")
        self.jc.sol(self._soln).feature("s1").feature("i1").feature("mg1").feature("po").feature("sv1").set("iter", jtypes.JInt(2))
        self.jc.sol(self._soln).feature("s1").feature("i1").feature("mg1").feature("po").feature("sv1").set("relax", jtypes.JInt(1))
        self.jc.sol(self._soln).feature("s1").feature("i1").feature("mg1").feature("po").feature("sv1").set("sorvecdof", jtypes.JArray(jtypes.JString)(["comp1_E"]))
        self.jc.sol(self._soln).feature("s1").feature("i1").feature("mg1").feature("cs").create("d1", "Direct")
        self.jc.sol(self._soln).feature("s1").feature("i1").feature("mg1").feature("cs").feature("d1").set("linsolver", "pardiso")
        self.jc.sol(self._soln).feature("s1").feature("fc1").set("linsolver", "i1")
        self.jc.sol(self._soln).feature("s1").feature().remove("fcDef")
        self.jc.sol(self._soln).attach("stdRFsparams")

    def _run_premesh(self):
        #Create physics and study for RF analysis (e.g. s-parameters)
        self.phys_emw = self.model._add_physics("emw")
        self.jc.component("comp1").physics().create(self.phys_emw, "ElectromagneticWaves", "geom1")
        #Create any lumped element ports
        for cur_elem_id, cur_lelem in enumerate(self.lumped_elems):
            elem_name = "lelem" + str(cur_elem_id)
            self.jc.component("comp1").physics(self.phys_emw).create(elem_name, "LumpedElement", 2)
            elem_bnds = self.model._get_selection_boundaries(cur_lelem[0])
            self.jc.component("comp1").physics(self.phys_emw).feature(elem_name).selection().set(jtypes.JArray(jtypes.JInt)(elem_bnds))
            self.jc.component("comp1").physics(self.phys_emw).feature(elem_name).set("LumpedElementType", cur_lelem[2])
            if cur_lelem[2] == 'Inductor':
                self.jc.component("comp1").physics(self.phys_emw).feature(elem_name).set('Lelement', jtypes.JDouble(cur_lelem[1]))
            elif cur_lelem[2] == 'LCparallel':
                self.jc.component("comp1").physics(self.phys_emw).feature(elem_name).set('Lelement', jtypes.JDouble(cur_lelem[1][0]))
                self.jc.component("comp1").physics(self.phys_emw).feature(elem_name).set('Celement', jtypes.JDouble(cur_lelem[1][1]))
            elif cur_lelem[2] == 'RLCparallel':
                self.jc.component("comp1").physics(self.phys_emw).feature(elem_name).set('Lelement', jtypes.JDouble(cur_lelem[1][0]))
                self.jc.component("comp1").physics(self.phys_emw).feature(elem_name).set('Celement', jtypes.JDouble(cur_lelem[1][1]))
                self.jc.component("comp1").physics(self.phys_emw).feature(elem_name).set('Relement', jtypes.JDouble(cur_lelem[1][2]))


        if self.model._include_loss_tangents:
            self.jc.component("comp1").physics("emw").feature("wee1").set("DisplacementFieldModel", "LossTangent")

        #Note that PEC1 is the default exterior boundary condition
        if not self._discard_PEC_BC:
            self.jc.component("comp1").physics(self.phys_emw).create("pec2", "PerfectElectricConductor", 2)
            self.jc.component("comp1").physics(self.phys_emw).feature("pec2").selection().named("geom1_condAll")
        #Create the excitation ports
        self.port_names = []
        for cur_port_id,cur_port in enumerate(self._ports):
            port_name = "lport" + str(cur_port_id)
            self.jc.component("comp1").physics(self.phys_emw).create(port_name, "LumpedPort", 2)
            if cur_port['type'] == 'CPW':
                self.jc.component("comp1").physics(self.phys_emw).feature(port_name).set('PortType', 'MultiElementUniform')
                self.jc.component("comp1").physics(self.phys_emw).feature(port_name).selection().named('geom1_'+cur_port['polys'][0])
                self.jc.component("comp1").physics(self.phys_emw).feature(port_name).feature('ue1').selection().named('geom1_'+cur_port['polys'][1])
                self.jc.component("comp1").physics(self.phys_emw).feature(port_name).feature('ue2').selection().named('geom1_'+cur_port['polys'][2])
                self.jc.component("comp1").physics(self.phys_emw).feature(port_name).feature('ue1').set('ahUniformElement', jtypes.JArray(jtypes.JDouble)([cur_port['vec_perp'][0],cur_port['vec_perp'][1],0.0]))
                self.jc.component("comp1").physics(self.phys_emw).feature(port_name).feature('ue2').set('ahUniformElement', jtypes.JArray(jtypes.JDouble)([-cur_port['vec_perp'][0],-cur_port['vec_perp'][1],0.0]))
            else:
                self.jc.component("comp1").physics(self.phys_emw).feature(port_name).set('PortType', 'Uniform')
                self.jc.component("comp1").physics(self.phys_emw).feature(port_name).selection().named('geom1_'+cur_port['polys'][0])
            #
            self.port_names.append(port_name)
        #Create any current sources
        for m,cur_current_src in enumerate(self._current_sources):
            assert False, "Do not add surface current sources for now"
            #TODO: Finish this feature later...
            region, extrusion_depth, J_z = cur_current_src[1], cur_current_src[2], cur_current_src[3]  # noqa: F841
            
            if cur_current_src[0] == 'region':
                cur_extrusion = f"ext{m}"
                self.jc.component("comp1").geom("geom1").feature().create(cur_extrusion, "Extrude")
                self.jc.component("comp1").geom("geom1").feature(cur_extrusion).set("extrudefrom", "faces")
                self.jc.component("comp1").geom("geom1").feature(cur_extrusion).selection("inputface").named("condAll")
                self.jc.component("comp1").geom("geom1").feature(cur_extrusion).setIndex("distance", str(extrusion_depth), 0)

                #It'll have to be:
                # - Create workplane + extrusion separately for each polygon
                # - Contribute entire thing to a selection
                # - Create selection spheres on top and bottom sections to slice out the top/bottom
                # - Create difference-selection to isolate side-walls
                # - Create S.C. B.C. on side-walls and PEC on top/bottom sections...

                # self.jc.component("comp1").physics(self.phys_emw).create(f"scu{m}", "SurfaceCurrent", 2)
                # 
                # if region == 'metals':
                #     self.jc.component("comp1").physics(self.phys_emw).feature(f"scu{m}").selection().named("geom1_condAll")
                # else:
                #     self.jc.component("comp1").physics(self.phys_emw).feature(f"scu{m}").selection().named("geom1_dielectricSurface")
                # self.jc.component("comp1").physics(self.phys_emw).feature(f"scu{m}").set("Js0", jtypes.JArray(jtypes.JDouble)([0, 0, J_z]))
        #Create electric point dipoles
        for m,cur_dipole in enumerate(self._electric_point_dipoles):
            sel_name, dipole_moment_p, dipole_moment_vecdir = cur_dipole
            edp_name = f"edp{m}"
            self.jc.component("comp1").physics(self.phys_emw).create(edp_name, "ElectricPointDipole", 0)
            self.jc.component("comp1").physics(self.phys_emw).feature(edp_name).selection().named(f"geom1_{sel_name}_pnt")  #Why does it require this _pnt?
            self.jc.component("comp1").physics(self.phys_emw).feature(edp_name).set("enpI", jtypes.JArray(jtypes.JDouble)(dipole_moment_vecdir.tolist()))
            self.jc.component("comp1").physics(self.phys_emw).feature(edp_name).set("normpI", jtypes.JFloat(dipole_moment_p))

    def create_port_2_conds(self, qObjName1, pin1, qObjName2, pin2, rect_width=20e-6):
        unit_conv = QUtilities.get_units(self.model.design)
        pos1 = self.model.design.components[qObjName1].pins[pin1]['middle'] * unit_conv
        pos2 = self.model.design.components[qObjName2].pins[pin2]['middle'] * unit_conv
        self.create_port_2_conds_by_position(pos1, pos2, rect_width)

    def create_port_2_conds_by_position(self, pos1, pos2, rect_width=20e-6):
        coords = ShapelyEx.rectangle_from_line(pos1, pos2, rect_width, False)

        sel_x, sel_y, sel_r = self.model._create_poly(f"port{len(self._ports)}", coords)
        select_3D_name = self.model._setup_selection_boundaries(len(self._ports), f"port{len(self._ports)}", 'port')
        cur_port = {'type':'Single', 'polys': [select_3D_name]}
        self._ports += [cur_port]

    def create_port_CPW_on_Launcher(self, qObjName, len_launch = 20e-6):
        '''
        Creates an RF port on a CPW inlet. The elements form fins to ground from the central CPW stripline. Note that the first port will automatically be
        the 50Ohm excitation port in the E-field plots while the second port is a 50Ohm ground. The s-parameters will calculate S11 and S21 if 2 such ports
        are defined.
        Inputs:
            - CPW_obj - A CPW object that has the attributes: start, end, width, gap
            - is_start - If True, then the port is attached to the start of the CPW, while False attaches the port to the end of the CPW
            - len_launch - (Default: 20e-6) Length of the inlet port fins along the the CPW. It is a good idea to keep it thin w.r.t. CPW gap 
        '''
        pol_name = "port_launch" + str(len(self._ports))
        
        launchesA, launchesB, vec_perp = QUtilities.get_RFport_CPW_coords_Launcher(self.model.design, qObjName, len_launch)

        cur_port = {'type':'CPW'}

        sel_x, sel_y, sel_r = self.model._create_poly(pol_name + "a", launchesA)
        select_3D_nameA = self.model._setup_selection_boundaries(len(self._ports), pol_name + "a", 'porta')

        sel_x, sel_y, sel_r = self.model._create_poly(pol_name + "b", launchesB)
        select_3D_nameB = self.model._setup_selection_boundaries(len(self._ports), pol_name + "b", 'portb')

        sel_all_name = f'uniselPort{len(self._ports)}'
        self.jc.component("comp1").geom("geom1").create(sel_all_name, "UnionSelection")
        self.jc.component("comp1").geom("geom1").feature(sel_all_name).label(sel_all_name)
        self.jc.component("comp1").geom("geom1").feature(sel_all_name).set("entitydim", jtypes.JInt(2))
        self.jc.component("comp1").geom("geom1").feature(sel_all_name).set("input", jtypes.JArray(jtypes.JString)([select_3D_nameA, select_3D_nameB]))

        cur_port['polys'] = [sel_all_name, select_3D_nameA, select_3D_nameB]

        #Each port is defined as: [portA-selection-name, portB-selection-name, vec_CPW2GND_1] where vec_CPW2GND_1 is a db.DVector pointing in the direction
        #of ground from the CPW for portA.
        cur_port['vec_perp'] = vec_perp
        self._ports += [cur_port]
    
    def create_port_CPW_on_Route(self, qObjName, pin_name='end', len_launch = 20e-6):
        '''
        Creates an RF port on a CPW inlet. The elements form fins to ground from the central CPW stripline. Note that the first port will automatically be
        the 50Ohm excitation port in the E-field plots while the second port is a 50Ohm ground. The s-parameters will calculate S11 and S21 if 2 such ports
        are defined.
        Inputs:
            - CPW_obj - A CPW object that has the attributes: start, end, width, gap
            - is_start - If True, then the port is attached to the start of the CPW, while False attaches the port to the end of the CPW
            - len_launch - (Default: 20e-6) Length of the inlet port fins along the the CPW. It is a good idea to keep it thin w.r.t. CPW gap 
        '''

        pol_name = "port_launch" + str(len(self._ports))

        launchesA, launchesB, vec_perp = QUtilities.get_RFport_CPW_coords_Route(self.model.design, qObjName, pin_name, len_launch)

        cur_port = {'type':'CPW', 'polys': []}

        sel_x, sel_y, sel_r = self.model._create_poly(pol_name + "a", launchesA)
        select_3D_nameA = self.model._setup_selection_boundaries(len(self._ports), pol_name + "a", 'porta')

        sel_x, sel_y, sel_r = self.model._create_poly(pol_name + "b", launchesB)
        select_3D_nameB = self.model._setup_selection_boundaries(len(self._ports), pol_name + "b", 'portb')

        sel_all_name = f'uniselPort{len(self._ports)}'
        self.jc.component("comp1").geom("geom1").create(sel_all_name, "UnionSelection")
        self.jc.component("comp1").geom("geom1").feature(sel_all_name).label(sel_all_name)
        self.jc.component("comp1").geom("geom1").feature(sel_all_name).set("entitydim", jtypes.JInt(2))
        self.jc.component("comp1").geom("geom1").feature(sel_all_name).set("input", jtypes.JArray(jtypes.JString)([select_3D_nameA, select_3D_nameB]))

        cur_port['polys'] = [sel_all_name, select_3D_nameA, select_3D_nameB]

        #Each port is defined as: [portA-selection-name, portB-selection-name, vec_CPW2GND_1] where vec_CPW2GND_1 is a db.DVector pointing in the direction
        #of ground from the CPW for portA.
        cur_port['vec_perp'] = vec_perp
        self._ports += [cur_port]

    def create_RFport_CPW_groundU_Launcher_inplane(self, qObjName, thickness_side=20e-6, thickness_back=20e-6, separation_gap=0e-6):
        Uclip = QUtilities.get_RFport_CPW_groundU_Launcher_inplane(self.model.design, qObjName, thickness_side, thickness_back, separation_gap)
        self.model._add_cond(Uclip)

    def create_RFport_CPW_groundU_Route_inplane(self, route_name, pin_name, thickness_side=20e-6, thickness_back=20e-6, separation_gap=0e-6):
        Uclip = QUtilities.get_RFport_CPW_groundU_Route_inplane(self.model.design, route_name, pin_name, thickness_side, thickness_back, separation_gap)
        self.model._add_cond(Uclip)

    def create_lumped_inductance(self, qObj1, pin1, qObj2, pin2, width, value):
        unit_conv = QUtilities.get_units(self.model.design)
        pos1 = self.model.design.components[qObj1].pins[pin1]['middle'] * unit_conv
        pos2 = self.model.design.components[qObj2].pins[pin2]['middle'] * unit_conv
        v_parl = pos2-pos1
        v_perp = np.array([-v_parl[1], v_parl[0]])
        v_perp /= np.linalg.norm(v_perp)
        v_perp *= width*0.5

        sel_x, sel_y, sel_r = self.model._create_poly(f"lumpedPort{len(self.lumped_elems)}", np.array([pos1+v_perp, pos1-v_perp, pos1-v_perp+v_parl, pos1+v_perp+v_parl]))

        self.lumped_elems.append([self.model._create_boundary_selection_sphere(sel_r, sel_x, sel_y), value, 'Inductor'])

    def create_port_JosephsonJunction(self, qObjName, **kwargs):
        junction_index = kwargs.get('junction_index', 0)

        comp_id = self.model.design.components[qObjName].id
        gsdf = self.model.design.qgeometry.tables['junction']
        gsdf = gsdf.loc[gsdf["component"] == comp_id]
        ls = gsdf['geometry'].iloc[junction_index]
        assert isinstance(ls, shapely.geometry.linestring.LineString), "The junction must be defined as a LineString object. Check the code for the part to verify this..."

        coords = [x for x in ls.coords]
        assert len(coords) == 2, "Currently only junctions with two points (i.e. a rectangle) are supported."
        rect_width = gsdf['width'].iloc[junction_index]

        if 'E_J_Hertz' in kwargs:
            elem_charge = 1.60218e-19
            hbar = 1.05457e-34
            h = hbar*2*np.pi
            phi0 = hbar/(2*elem_charge)
            E_J = kwargs.pop('E_J_Hertz') * h
            L_ind = phi0**2 / E_J
        elif 'L_J' in kwargs:
            L_ind = kwargs.pop('L_J')
        else:
            assert False, "Must supply either E_J_Hertz or L_J to define a Josephson Junction port."

        C_J = kwargs.get('C_J', 0)
        R_J = kwargs.get('R_J', 0)

        port_name = "rf_port_" + str(len(self._ports))  # noqa: F841

        unit_conv = QUtilities.get_units(self.model.design)
        pos1 = np.array(coords[0]) * unit_conv
        pos2 = np.array(coords[1]) * unit_conv
        v_parl = pos2-pos1
        v_parl /= np.linalg.norm(v_parl)

        portCoords = ShapelyEx.rectangle_from_line(pos1, pos2, rect_width * unit_conv, False)
        portCoords = [x for x in portCoords]

        sel_x, sel_y, sel_r = self.model._create_poly(f"lumpedPort{len(self.lumped_elems)}", np.array(portCoords))
        if C_J > 0:
            if R_J > 0:
                self.lumped_elems.append([self.model._create_boundary_selection_sphere(sel_r, sel_x, sel_y), (L_ind, C_J, R_J), 'RLCparallel'])
            else:
                self.lumped_elems.append([self.model._create_boundary_selection_sphere(sel_r, sel_x, sel_y), (L_ind, C_J), 'LCparallel'])
        else:
            self.lumped_elems.append([self.model._create_boundary_selection_sphere(sel_r, sel_x, sel_y), L_ind, 'Inductor'])

    def set_freq_range(self, freq_start, freq_end, num_points, use_previous_solns = False):
        '''
        Set the frequency range to sweep in the s-parameter and E-field simulation. THE SIMULATION MUST BE BUILT BY THIS POINT (i.e. running build_geom_mater_elec_mesh).
        Inputs:
            - freq_start, freq_end - Start and end (inclusive) frequencies in units of Hertz
            - num_points - Number of points to use in the sweep (must be above 1)
            - use_previous_solns - (Default: False) If True, it will use the previous frequency value as a starting point for the current point; it can have the tendency
                                   to falsely find local minima and cause the solution to have 'discrete steps' - so it's better to keep it False if accuracy is important.
                                   Note that this is ignored if using Frequency Modal simulations (i.e. setting adaptive to 'Multiple')
        '''
        if num_points > 1:
            str_freqs = "range({0}[GHz],({1}[GHz]-({0}[GHz]))/{2},{1}[GHz])".format(freq_start*1e-9, freq_end*1e-9, num_points-1)
        else:
            str_freqs = "{0}[GHz]*1^range(1,1)".format(freq_start*1e-9, )
        if self.adaptive == 'Multiple':
            self.jc.study(self._study).feature(self._sub_study[1]).set("plist", str_freqs)
        else:
            self.jc.study(self._study).feature(self._sub_study).set("errestandadap", "none")
            self.jc.study(self._study).feature(self._sub_study).set("plist", str_freqs)
            # if use_previous_solns:
            #     self.jc.study(self._study).feature(self._sub_study).set("preusesol", "yes")
            # else:
            #     self.jc.study(self._study).feature(self._sub_study).set("preusesol", "no")

    def set_freq_values(self, list_freqs, use_previous_solns = False):
        '''
        Set the frequency range to sweep in the s-parameter and E-field simulation. THE SIMULATION MUST BE BUILT BY THIS POINT (i.e. running build_geom_mater_elec_mesh).
        Inputs:
            - freq_start, freq_end - Start and end (inclusive) frequencies in units of Hertz
            - num_points - Number of points to use in the sweep (must be above 1)
            - use_previous_solns - (Default: False) If True, it will use the previous frequency value as a starting point for the current point; it can have the tendency
                                   to falsely find local minima and cause the solution to have 'discrete steps' - so it's better to keep it False if accuracy is important.
                                   Note that this is ignored if using Frequency Modal simulations (i.e. setting adaptive to 'Multiple')
        '''
        str_freqs = ", ".join(str(x*1e-9) for x in list_freqs)
        if self.adaptive == 'Multiple':
            self.jc.study(self._study).feature(self._sub_study[1]).set("plist", str_freqs)
        else:
            self.jc.study(self._study).feature(self._sub_study).set("errestandadap", "none")
            self.jc.study(self._study).feature(self._sub_study).set("plist", str_freqs)

    def add_surface_current_source_region(self, region="", extrusion_depth=1e-6, J_z=1.0):
        assert self.adaptive == 'None', "Cannot use eigenmode-based simulations when setting current sources."

        assert region in ['metals', 'dielectric'], "Must set region to either 'metals' or 'dielectric'"
        self._current_sources.append(('region', region, extrusion_depth, J_z))
        if region == 'metals':
            self._discard_PEC_BC = True

    def add_electric_point_dipole(self, vec_pos, dipole_moment_p, dipole_moment_vecdir):
        pt_name = f"pt{len(self._electric_point_dipoles)}"
        sel_name = "sel_" + pt_name
        self.jc.component("comp1").geom("geom1").create(pt_name, "Point")
        self.jc.component("comp1").geom("geom1").feature(pt_name).setIndex("p", str(vec_pos[0]), 0)
        self.jc.component("comp1").geom("geom1").feature(pt_name).setIndex("p", str(vec_pos[1]), 1)
        self.jc.component("comp1").geom("geom1").feature(pt_name).setIndex("p", str(vec_pos[2]), 2)
        self.jc.component("comp1").geom("geom1").run(pt_name)
        self.jc.component("comp1").geom("geom1").selection().create(sel_name, "CumulativeSelection")
        self.jc.component("comp1").geom("geom1").selection(sel_name).label(pt_name)
        self.jc.component("comp1").geom("geom1").feature(pt_name).set("contributeto", sel_name)
        #
        dipole_moment_vecdir = np.array(dipole_moment_vecdir)
        self._electric_point_dipoles.append((sel_name, dipole_moment_p, dipole_moment_vecdir/np.linalg.norm(dipole_moment_vecdir)))

    def _get_required_physics(self):
        return [self.phys_emw]

    def run(self, recompute=True, is_dB=False):
        '''
        Run simulation to get s-parameters after running the simulation. Returns a 3-row array in which the rows are: frequency values, S11s, S21s.
        MUST RUN AFTER COMPILING SIMULATION (i.e. after running build_geom_mater_elec_mesh...)
        
        Inputs:
            - recompute - (Default True) If true, the solution result is recomputed
        '''
        if (recompute):
            self.jc.sol(self._soln).runAll()
        
        self.jc.result().numerical().create("ev1", "Eval")
        self.jc.result().numerical("ev1").set("data", self.dset_name)
        self.jc.result().numerical("ev1").set("expr", self._sub_study)
        freqs = self.jc.result().numerical("ev1").getData()
        freqs = np.array(freqs[0])[:,1]
        s1Remains = []
        for cur_ind in range(0,len(self._ports)):
            if is_dB:
                self.jc.result().numerical("ev1").set("expr", f"emw.S{cur_ind+1}1dB")
                cur_sVals = self.jc.result().numerical("ev1").getData()
                #For some reason the returned arrays are rows of the same data apparently repeated across the columns...
                cur_sVals = np.array(cur_sVals[0])[:,1]
            else:
                #For some reason the returned arrays are rows of the same data apparently repeated across the columns...
                self.jc.result().numerical("ev1").set("expr", f"real(emw.S{cur_ind+1}1)")
                cur_sValsR = self.jc.result().numerical("ev1").getData()
                self.jc.result().numerical("ev1").set("expr", f"imag(emw.S{cur_ind+1}1)")
                cur_sValsI = self.jc.result().numerical("ev1").getData()
                #
                cur_sVals = np.array(cur_sValsR[0])[:,0] + 1j*np.array(cur_sValsI[0])[:,0]
            s1Remains += [cur_sVals]
        self.jc.result().numerical().remove("ev1")

        return np.vstack([freqs] + s1Remains)

    def eval_fields_over_mesh(self):
        x,y,z,Ex,Ey,Ez,Bx,By,Bz,Dx,Dy,Dz,Hx,Hy,Hz = self.model._model.evaluate(['x', 'y', 'z', 'emw.Ex', 'emw.Ey', 'emw.Ez', 'emw.Bx', 'emw.By', 'emw.Bz', 'emw.Dx', 'emw.Dy', 'emw.Dz', 'emw.Hx', 'emw.Hy', 'emw.Hz'])
        return {'coords': np.vstack([x,y,z]).T,
                'E': np.vstack([Ex,Ey,Ez]).T,
                'B': np.vstack([Bx,By,Bz]).T,
                'D': np.vstack([Dx,Dy,Dz]).T,
                'H': np.vstack([Hx,Hy,Hz]).T}

    def eval_field_at_pts(self, expr, coords: np.ndarray, freq_index=1):
        allowed_exprs = ['Ex', 'Ey', 'Ez', 'E',
                         'Bx', 'By', 'Bz', 'B',
                         'Dx', 'Dy', 'Dz', 'D',
                         'Hx', 'Hy', 'Hz', 'H']
        assert expr in allowed_exprs, f"Must only ask for: {str(allowed_exprs)}"

        if self.model._dset_exists("cpt1"):
            self.jc.result().dataset().remove("cpt1")
        if self.model._numeval_exists("pev1"):
            self.jc.result().numerical().remove("pev1")
        
        self.jc.result().dataset().create("cpt1", "CutPoint3D")

        np.savetxt("temp_data.csv", np.real(coords), delimiter=',')

        #Loading via Java is slow...
        # self.jc.result().dataset("cpt1").set("pointx", jtypes.JArray(jtypes.JDouble)(np.real(coords[:,0]).tolist()))
        # self.jc.result().dataset("cpt1").set("pointy", jtypes.JArray(jtypes.JDouble)(np.real(coords[:,1]).tolist()))
        # self.jc.result().dataset("cpt1").set("pointz", jtypes.JArray(jtypes.JDouble)(np.real(coords[:,2]).tolist()))
        self.jc.result().dataset("cpt1").set("method", "file")
        self.jc.result().dataset("cpt1").set("filename", "temp_data.csv")

        self.jc.result().dataset("cpt1").set("data", self.dset_name)
        
        self.jc.result().numerical().create("pev1", "EvalPoint")
        self.jc.result().numerical("pev1").set("data", "cpt1")
        self.jc.result().numerical("pev1").setIndex("looplevelinput", "manualindices", 0)
        self.jc.result().numerical("pev1").setIndex("looplevelindices", jtypes.JInt(freq_index), 0)

        if expr in ['Ex', 'Ey', 'Ez', 'Bx', 'By', 'Bz', 'Dx', 'Dy', 'Dz', 'Hx', 'Hy', 'Hz']:
            self.jc.result().numerical("pev1").setIndex("expr", f"emw.{expr}", 0)
            cur_vals_R = self.jc.result().numerical("pev1").getReal()
            cur_vals_I = self.jc.result().numerical("pev1").getImag()
            ret_vals = (np.array(cur_vals_R) + 1j*np.array(cur_vals_I))[:,0]
        elif expr in ['E', 'B', 'D', 'H']:
            self.jc.result().numerical("pev1").setIndex("expr", f"emw.{expr}x", 0)
            cur_vals_R = self.jc.result().numerical("pev1").getReal()
            cur_vals_I = self.jc.result().numerical("pev1").getImag()
            x_vals = (np.array(cur_vals_R) + 1j*np.array(cur_vals_I))[:,0]
            #
            self.jc.result().numerical("pev1").setIndex("expr", f"emw.{expr}y", 0)
            cur_vals_R = self.jc.result().numerical("pev1").getReal()
            cur_vals_I = self.jc.result().numerical("pev1").getImag()
            y_vals = (np.array(cur_vals_R) + 1j*np.array(cur_vals_I))[:,0]
            #
            self.jc.result().numerical("pev1").setIndex("expr", f"emw.{expr}z", 0)
            cur_vals_R = self.jc.result().numerical("pev1").getReal()
            cur_vals_I = self.jc.result().numerical("pev1").getImag()
            z_vals = (np.array(cur_vals_R) + 1j*np.array(cur_vals_I))[:,0]
            #
            ret_vals = np.vstack([x_vals, y_vals, z_vals]).T

        self.jc.result().numerical().remove("pev1")
        self.jc.result().dataset().remove("cpt1")
        if os.path.exists("temp_data.csv"):
            os.remove("temp_data.csv")

        return ret_vals

    def run_only_eigenfrequencies(self, return_dofs=False):
        #TODO: Check it is in eigenmode type and find a better solution than this... This basically evaluates a table and then extracts the first column
        #which happens to hold the simulated eigenfrequencies...
        self.jc.sol(self._soln).runFromTo("st1", "su1")
        self.jc.result().numerical().create("gev1", "EvalGlobal")
        self.jc.result().numerical("gev1").set("data", self.dset_name)
        self.jc.result().numerical("gev1").set("expr", jtypes.JArray(jtypes.JString)(["numberofdofs"]))

        self.jc.result().table().create("tbl1", "Table")
        self.jc.result().numerical("gev1").set("table", "tbl1")
        self.jc.result().numerical("gev1").computeResult()
        self.jc.result().numerical("gev1").setResult()

        freqs = [1e9*np.complex128(str(y).replace('i','j')) for y in [x[0] for x in self.jc.result().table("tbl1").getTableData(True)]]
        dofs = int(np.array(self.jc.result().numerical("gev1").computeResult())[0,0,0])

        self.jc.result().table().remove("tbl1")
        self.jc.result().numerical().remove("gev1")

        if return_dofs:
            return freqs, dofs
        else:
            return freqs

    def _cost_func(self, freq, param_index=1, find_max=False):
        self.set_freq_range(freq[0], freq[0],1)
        data = self.run()
        
        self.optimised_data += [data]
        if find_max:
            return -data[param_index]
        else:
            return data[param_index]

    def find_peak(self, initial_guess, search_bnd_min_freq, search_bnd_max_freq, s_param_index, find_max=False, freq_accuracy=1e6):
        self.optimised_data = []
        sol = scipy.optimize.minimize(lambda x: self._cost_func(x, s_param_index, find_max), [initial_guess], method='Nelder-Mead', bounds=(scipy.optimize.Bounds(search_bnd_min_freq, search_bnd_max_freq)), options={'xatol': freq_accuracy} )
        self.optimised_data = np.hstack(self.optimised_data).T
        if find_max:
            self.optimal_value = (sol.x, -sol.fun)
        else:
            self.optimal_value = (sol.x, sol.fun)
        return self.optimal_value
