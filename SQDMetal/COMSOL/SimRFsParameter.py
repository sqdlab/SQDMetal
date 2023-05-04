from SQDMetal.COMSOL.Model import COMSOL_Simulation_Base
from SQDMetal.Utilities.QUtilities import QUtilities

import mph
import jpype.types as jtypes
import geopandas as gpd
import shapely
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize

from qiskit_metal.qlibrary.terminations.launchpad_wb import LaunchpadWirebond

class COMSOL_Simulation_RFsParameters(COMSOL_Simulation_Base):
    class Lumped_Element:
        def __init__(self, selection, elem_type, value):
            self._sel = selection
            self._elem_type = elem_type
            self.value = value
            self._setter = None
        
        def set_value(self, val):
            self.value = val
            if not self._setter is None:
                self._setter.set('Lelement', jtypes.JDouble(val))

    def __init__(self, model, adaptive='Single', modal_min_freq_num_eigs = (1e9,7), relative_tolerance=0.01):
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
            lst_ports = [item for sublist in self._ports for item in sublist[:2]]
            self.jc.study(self._study).feature(self._sub_study[1]).set("usestoresel", "selection");
            self.jc.study(self._study).feature(self._sub_study[1]).set("storesel", jtypes.JArray(jtypes.JString)(lst_ports))
        self._soln = self.model._add_solution("solRFsparams")
        self.jc.sol().create(self._soln)
        self.jc.sol(self._soln).createAutoSequence(self._study)
        # if self.adaptive == 'None':
        #     self.jc.sol(self._soln).feature("s1").set("stol", self.relative_tolerance)
        self.dset_name = self.model._get_dset_name()

    def _run_premesh(self):
        #Create physics and study for RF analysis (e.g. s-parameters)
        self.phys_emw = self.model._add_physics("emw")
        self.jc.component("comp1").physics().create(self.phys_emw, "ElectromagneticWaves", "geom1")
        #Create any lumped element ports
        for cur_elem_id, cur_lelem in enumerate(self.lumped_elems):
            elem_name = "lelem" + str(cur_elem_id)
            self.jc.component("comp1").physics(self.phys_emw).create(elem_name, "LumpedElement", 2)
            elem_bnds = self.model._get_selection_boundaries(cur_lelem._sel[0])
            self.jc.component("comp1").physics(self.phys_emw).feature(elem_name).selection().set(jtypes.JArray(jtypes.JInt)(elem_bnds))
            self.jc.component("comp1").physics(self.phys_emw).feature(elem_name).set("LumpedElementType", cur_lelem._elem_type)
            cur_lelem._setter = self.jc.component("comp1").physics(self.phys_emw).feature(elem_name)
            cur_lelem.set_value(cur_lelem.value)

        #Note that PEC1 is the default exterior boundary condition
        self.jc.component("comp1").physics(self.phys_emw).create("pec2", "PerfectElectricConductor", 2)
        self.jc.component("comp1").physics(self.phys_emw).feature("pec2").selection().named("geom1_condAll")
        #Create the excitation ports
        for cur_port_id,cur_port in enumerate(self._ports):
            port_name = "lport" + str(cur_port_id)
            self.jc.component("comp1").physics(self.phys_emw).create(port_name, "LumpedPort", 2)
            self.jc.component("comp1").physics(self.phys_emw).feature(port_name).set('PortType', 'MultiElementUniform')
            port_bndsA = self.model._get_selection_boundaries(cur_port[0])
            port_bndsB = self.model._get_selection_boundaries(cur_port[1])
            self.jc.component("comp1").physics(self.phys_emw).feature(port_name).selection().set(jtypes.JArray(jtypes.JInt)(port_bndsA+port_bndsB))
            self.jc.component("comp1").physics(self.phys_emw).feature(port_name).feature('ue1').selection().set(jtypes.JArray(jtypes.JInt)(port_bndsA))
            self.jc.component("comp1").physics(self.phys_emw).feature(port_name).feature('ue2').selection().set(jtypes.JArray(jtypes.JInt)(port_bndsB))
            self.jc.component("comp1").physics(self.phys_emw).feature(port_name).feature('ue1').set('ahUniformElement', jtypes.JArray(jtypes.JDouble)([cur_port[2][0],cur_port[2][1],0.0]))
            self.jc.component("comp1").physics(self.phys_emw).feature(port_name).feature('ue2').set('ahUniformElement', jtypes.JArray(jtypes.JDouble)([-cur_port[2][0],-cur_port[2][1],0.0]))

    def create_port_on_CPW(self, qObjName, is_start=True, len_launch = 20e-6):
        '''
        Creates an RF port on a CPW inlet. The elements form fins to ground from the central CPW stripline. Note that the first port will automatically be
        the 50Ohm excitation port in the E-field plots while the second port is a 50Ohm ground. The s-parameters will calculate S11 and S21 if 2 such ports
        are defined.
        Inputs:
            - CPW_obj - A CPW object that has the attributes: start, end, width, gap
            - is_start - If True, then the port is attached to the start of the CPW, while False attaches the port to the end of the CPW
            - len_launch - (Default: 20e-6) Length of the inlet port fins along the the CPW. It is a good idea to keep it thin w.r.t. CPW gap 
        '''
        
        qObj = self.model.design.components[qObjName]
        if isinstance(qObj, LaunchpadWirebond):
            vec_ori, vec_launch, cpw_wid, cpw_gap = self._get_LauncherWB_params(qObjName)
        else:
            assert False, f"\'{qObjName}\' is an unsupported object type."

        vec_perp = np.array([vec_launch[1],-vec_launch[0]])
        vec_launch *= len_launch

        pol_name = "port_launch" + str(len(self._ports))
        
        launches = [vec_ori + vec_perp * cpw_wid*0.5, vec_ori + vec_perp * (cpw_wid*0.5+cpw_gap),
                    vec_ori + vec_launch + vec_perp * (cpw_wid*0.5+cpw_gap), vec_ori + vec_launch + vec_perp * cpw_wid*0.5]
        launches = [[p[0],p[1]] for p in launches]
        sel_x, sel_y, sel_r = self.model._create_poly(pol_name + "a", launches)
        cur_launch = [self.model._create_boundary_selection_sphere(sel_r, sel_x, sel_y)]

        launches = [vec_ori - vec_perp * cpw_wid*0.5, vec_ori - vec_perp * (cpw_wid*0.5+cpw_gap),
                    vec_ori + vec_launch - vec_perp * (cpw_wid*0.5+cpw_gap), vec_ori + vec_launch - vec_perp * cpw_wid*0.5]
        launches = [[p[0],p[1]] for p in launches]
        sel_x, sel_y, sel_r = self.model._create_poly(pol_name + "b", launches)
        cur_launch += [self.model._create_boundary_selection_sphere(sel_r, sel_x, sel_y)]

        #Each port is defined as: [portA-selection-name, portB-selection-name, vec_CPW2GND_1] where vec_CPW2GND_1 is a db.DVector pointing in the direction
        #of ground from the CPW for portA.
        cur_launch += [vec_perp]
        self._ports += [cur_launch]

    def create_lumped_element(self, vec_start, vec_end, width, value, elem_type='Inductor'):
        '''
        Creates a lumped element. Note that it must be vertical or horizontal between conductors...
        Inputs:
            - vec_start - A DPoint for the start of the rectangle (in Klayout units)
            - vec_end   - A DPoint for the end of the rectangle (in Klayout units)
            - width     - Width of the rectangle in metres
            - value     - Value of the component raw units (e.g. Henries for Inductors)
        '''
        vec_start *= self.model.kLy2metre
        vec_end *= self.model.kLy2metre
        
        vec_elem = vec_end - vec_start
        vec_perp = db.DVector(vec_elem.y,-vec_elem.x) / vec_elem.length()

        pol_name = "lumped_elem" + str(len(self.lumped_elems))
        
        launches = [vec_start - vec_perp * width*0.5, vec_start + vec_perp * width*0.5,
                    vec_start + vec_perp * width*0.5 + vec_elem, vec_start - vec_perp * width*0.5 + vec_elem]
        launches = [[p.x,p.y] for p in launches]
        sel_x, sel_y, sel_r = self.model._create_poly(pol_name, launches)
        cur_lumped_elem = [self.model._create_boundary_selection_sphere(sel_r, sel_x, sel_y)]

        new_elem = COMSOL_Simulation_RFsParameters.Lumped_Element(cur_lumped_elem, elem_type, value)
        self.lumped_elems += [new_elem]
        return new_elem

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
            str_freqs = "{0}[GHz]*1^range(1,1)".format(freq_start*1e-9, freq_end*1e-9, num_points-1)
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


    def _get_required_physics(self):
        return [self.phys_emw]

    def run(self, recompute=True):
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
            self.jc.result().numerical("ev1").set("expr", f"emw.S{cur_ind+1}1dB")
            cur_sVals = self.jc.result().numerical("ev1").getData()
            #For some reason the returned arrays are rows of the same data apparently repeated across the columns...
            cur_sVals = np.array(cur_sVals[0])[:,1]
            s1Remains += [cur_sVals]
        self.jc.result().numerical().remove("ev1")

        return np.vstack([freqs] + s1Remains)

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
    
    def _get_LauncherWB_params(self, launcher_name):
        design = self.model.design
        
        launcher_len = QUtilities.parse_value_length(design.components[launcher_name].options['pad_height']) + QUtilities.parse_value_length(design.components[launcher_name].options['taper_height']) + QUtilities.parse_value_length(design.components[launcher_name].options['lead_length'])
        unit_conv = QUtilities.get_units(design)
        startPt = design.components[launcher_name].pins['tie']['middle']*unit_conv - design.components[launcher_name].pins['tie']['normal']*launcher_len
        padDir = design.components[launcher_name].pins['tie']['normal']*1.0
        padWid = QUtilities.parse_value_length(design.components[launcher_name].options['pad_width'])
        padGap = QUtilities.parse_value_length(design.components[launcher_name].options['pad_gap'])

        return startPt, padDir, padWid, padGap
