from SQDMetal.COMSOL.Model import COMSOL_Simulation_Base

import jpype.types as jtypes
import geopandas as gpd
import shapely
import matplotlib.pyplot as plt
import numpy as np

class COMSOL_Simulation_CapMats(COMSOL_Simulation_Base):
    
    #terminal names
    terminal_names = []
    
    def __init__(self, model):
        self.model = model
        self.jc = model._get_java_comp()
        self.dset_name = ""
        self.phys_es = ""
        self._study = ""
        self._sub_study = "capMatStationary"
        self._soln = ""

    def _prepare_simulation(self):
        self.jc.param().set("PortName", jtypes.JInt(1))    #Just use default port name...
        self.jc.study().create(self._study)
        self.jc.study(self._study).create(self._sub_study,"Stationary")
        self._soln = self.model._add_solution("solCapMat")
        self.jc.sol().create(self._soln)
        self.jc.sol(self._soln).createAutoSequence(self._study)
        self.dset_name = self.model._get_dset_name()

    def _run_premesh(self):
        #Create physics and study for electrostatics (e.g. generating capacitance matrices)
        self._study = self.model._add_study("stdCapMat")
        self.phys_es = self.model._add_physics("es")
        self.jc.component("comp1").physics().create(self.phys_es, "Electrostatics", "geom1")
        self.jc.component("comp1").physics(self.phys_es).prop("PortSweepSettings").set("useSweep", jtypes.JBoolean(True))
        self.jc.component("comp1").physics(self.phys_es).prop("PortSweepSettings").set("PortParamName", "PortName")

        #Note that PEC1 is the default exterior boundary condition
        #Create terminals for capacitance matrix simulations
        for cur_term in range(len(self.model._conds)):
            term_name = "term"+str(cur_term)
            self.terminal_names.append(term_name)
            self.jc.component("comp1").physics(self.phys_es).create(term_name, "Terminal", 2)           
            self.jc.component("comp1").physics(self.phys_es).feature(term_name).selection().named("".join(self.model._conds[cur_term]))
            self.jc.component("comp1").physics(self.phys_es).feature(term_name).set("TerminalType", "Voltage")
            self.jc.component("comp1").physics(self.phys_es).feature(term_name).set("TerminalName", jtypes.JInt(cur_term+1))
        #Set the exterior boundaries to ground - otherwise capacitances will no longer have self-capacitances...
        self.jc.component("comp1").physics(self.phys_es).create("gnd1", "Ground", 2)
        self.jc.component("comp1").physics(self.phys_es).feature("gnd1").selection().named(self.model._sel_ext_boundaries)
        
    def _get_required_physics(self):
        return [self.phys_es]

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
        ax.set_xlabel('Position (m)')
        ax.set_ylabel('Position (m)')

    def run(self):
        '''
        Runs the simulation and returns a capacitance matrix. MUST RUN AFTER COMPILING SIMULATION (i.e. after running build_geom_mater_elec_mesh...)
        '''
        num_ports = len(self.model._conds)
        capMatFull = np.zeros([num_ports,num_ports])
        
        if num_ports == 1:
            self.jc.result().numerical().create("gev1", "EvalGlobal")
            self.jc.result().numerical("gev1").set("data", self.dset_name)
            self.jc.result().numerical("gev1").set("expr", "es.C11")

            self.jc.param().set("PortName", jtypes.JInt(1))
            #Evaluate the single self-capacitance
            self.jc.sol(self._soln).runAll()
            capMatFull = np.array(self.jc.result().numerical("gev1").computeResult())[0]
            self.jc.result().numerical().remove("gev1")
        else:
            #Setup temporary results dataset 
            self.jc.result().numerical().create("gmev1", "EvalGlobalMatrix")
            self.jc.result().numerical("gmev1").set("data", self.dset_name)
            self.jc.result().numerical("gmev1").set("expr", "es.C")
            
            for cur_port in range(num_ports):
                self.jc.param().set("PortName", jtypes.JInt(cur_port+1))
                #Evaluate column of capacitance matrix
                self.jc.sol(self._soln).runAll()
                #Extract column from result
                capCol = self.jc.result().numerical("gmev1").computeResult()
                capCol = np.array(capCol[0])
                capMatFull[:,cur_port] = capCol[:,cur_port]
            
            self.jc.result().numerical().remove("gmev1")
        return capMatFull