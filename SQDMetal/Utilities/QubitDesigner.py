import numpy as np
import scipy.optimize
from SQDMetal.Utilities.QUtilities import QUtilities

class ResonatorBase:
    def get_res_capacitance(self):
        raise NotImplementedError()
    def get_res_inductance(self):
        raise NotImplementedError()
    def get_res_impedance(self):
        raise NotImplementedError()
    def get_res_frequency(self):
        raise NotImplementedError()
    def print(self):
        print("Resonator:")
        print(f"\tFrequency: {QUtilities.add_units(self.get_res_frequency())}Hz")
        print(f"\tInductance: {QUtilities.add_units(self.get_res_inductance())}H")
        print(f"\tCapacitance: {QUtilities.add_units(self.get_res_capacitance())}F")

class ResonatorHalfWave(ResonatorBase):
    def __init__(self, f0, shorted=False, impedance=50):
        self.f0 = f0
        self.Zr = impedance
        self.shorted = shorted
    
    def get_res_capacitance(self):
        if self.shorted:
            #Equivalent to a SERIES RLC
            return 2/(np.pi * self.Zr * 2*np.pi*self.f0)
        else:
            #Equivalent to a PARALLEL RLC
            return np.pi/(2 * self.Zr * 2*np.pi*self.f0)
    
    def get_res_inductance(self):
        if self.shorted:
            #Equivalent to a SERIES RLC
            return (self.Zr*np.pi)/(2 * 2*np.pi*self.f0)
        else:
            #Equivalent to a PARALLEL RLC
            return (2*self.Zr)/(np.pi * 2*np.pi*self.f0)
    
    def get_res_impedance(self):
        return self.Zr
    
    def get_res_frequency(self):
        return self.f0

class ResonatorQuarterWave(ResonatorBase):
    def __init__(self, f0, shorted=False, impedance=50):
        self.f0 = f0
        self.Zr = impedance
        self.shorted = shorted
    
    def get_res_capacitance(self):
        if self.shorted:
            #Equivalent to a PARALLEL RLC
            return np.pi/(4 * self.Zr * 2*np.pi*self.f0)
        else:
            #Equivalent to a SERIES RLC
            return 4/(np.pi * self.Zr * 2*np.pi*self.f0)
    
    def get_res_inductance(self):
        if self.shorted:
            #Equivalent to a PARALLEL RLC
            return (4*self.Zr)/(np.pi * 2*np.pi*self.f0)
        else:
            #Equivalent to a SERIES RLC
            return (self.Zr*np.pi)/(4 * 2*np.pi*self.f0)
    
    def get_res_impedance(self):
        return self.Zr
    
    def get_res_frequency(self):
        return self.f0

class ResonatorLC(ResonatorBase):
    def __init__(self, f0=None, L=None, C=None):
        if f0 == None:
            assert L != None and C!= None, "Must supply L and C if omitting f0."
            self.L = L
            self.C = C
            self.f0 = 1/(2*np.pi * np.sqrt(L*C))
        elif L == None:
            assert f0 != None and C!= None, "Must supply f0 and C if omitting L."
            self.f0 = f0
            self.C = C
            self.L = 1/((2*np.pi*f0)**2 * C)
        elif C == None:
            assert f0 != None and L!= None, "Must supply f0 and L if omitting C."
            self.f0 = f0
            self.L = L
            self.C = 1/((2*np.pi*f0)**2 * L)
        else:
            assert np.abs(2*np.pi*f0 - 1/np.sqrt(L*C)) < 1e-9, "System overconstrained. Just supply 2 of L, C and F0..."
    
    def get_res_capacitance(self):
        if self.shorted:
            #Equivalent to a PARALLEL RLC
            return np.pi/(4 * self.Zr * 2*np.pi*self.f0)
        else:
            #Equivalent to a SERIES RLC
            return 4/(np.pi * self.Zr * 2*np.pi*self.f0)
    
    def get_res_inductance(self):
        if self.shorted:
            #Equivalent to a PARALLEL RLC
            return (4*self.Zr)/(np.pi * 2*np.pi*self.f0)
        else:
            #Equivalent to a SERIES RLC
            return (self.Zr*np.pi)/(4 * 2*np.pi*self.f0)
    
    def get_res_impedance(self):
        return np.sqrt(self.L/self.C)
    
    def get_res_frequency(self):
        return self.f0


class XmonDesigner:
    def __init__(self, resonator):
        self.resonator = resonator
    
    def get_free_params(self):
        return ['fQubit', 'C_g', 'C_J', 'chi']
    
    @staticmethod
    def _eff_CJ(C_J, C_g, C_r):
        return 0.5 * (2*C_J + 2*C_g/C_r*(C_J+C_r)) / (1+C_g/C_r)
    @staticmethod
    def _eff_Cr(C_J, C_g, C_r):
        return 0.5 * (2*C_r + 2*C_g/C_J*(C_J+C_r)) / (1+C_g/C_J)
    @staticmethod
    def _EJ_from_QubitFreqHertz(QubiFreqHertz, C_J, C_g, C_r):
        elem = 1.60217663e-19
        h = 6.62607015e-34
        CJeff = XmonDesigner._eff_CJ(C_J, C_g, C_r)
        Ec = (2*elem)**2 / (2*CJeff)
        return (QubiFreqHertz * h + Ec/4)**2 / (2*Ec)
    @staticmethod
    def _g_hertz(f_res, C_g, C_J, QubiFreqHertz, C_r, L_r):
        elem = 1.60217663e-19
        h = 6.62607015e-34
        Z_0 = 376.730313668
        alpha = 0.0072973525693
        CJeff = XmonDesigner._eff_CJ(C_J, C_g, C_r)
        Ec = (2*elem)**2 / (2*CJeff) / h
        E_J_hertz = XmonDesigner._EJ_from_QubitFreqHertz(QubiFreqHertz, C_J, C_g, C_r) / h

        Creff = XmonDesigner._eff_Cr(C_J, C_g, C_r)
        Z_r = np.sqrt(L_r/Creff)

        return (2*np.pi*f_res * C_g)/(C_J+C_g/C_r*(C_J+C_r)) * (2*E_J_hertz/Ec)**(1/4) * np.sqrt(Z_r/Z_0) * np.sqrt(2*np.pi*alpha) / (2*np.pi)

    @staticmethod
    def _chi_hertz(f_res, C_g, C_J, QubiFreqHertz, C_r, L_r):
        elem = 1.60217663e-19
        h = 6.62607015e-34
        g_hertz = XmonDesigner._g_hertz(f_res, C_g, C_J, QubiFreqHertz, C_r, L_r)
        CJeff = XmonDesigner._eff_CJ(C_J, C_g, C_r)
        Ec = (2*elem)**2 / (2*CJeff)
        delta = QubiFreqHertz - f_res

        return -(Ec/(24*h) * g_hertz**2) / (delta * (delta - Ec/(24*h)))

    @staticmethod
    def _anharmonicity_hertz(C_g, C_J, C_r):
        elem = 1.60217663e-19
        h = 6.62607015e-34
        CJeff = XmonDesigner._eff_CJ(C_J, C_g, C_r)
        Ec = (2*elem)**2 / (2*CJeff)
        return Ec/(24*h)

    def optimise(self, param_constraints):
        func = lambda x: abs(XmonDesigner._chi_hertz(self.resonator.get_res_frequency(), x[1], x[2], x[0], self.resonator.get_res_capacitance(), self.resonator.get_res_inductance()) - x[3])

        params = self.get_free_params()
        num_params = len(params)
        constrs = []
        x0 = []
        for m, cur_param in enumerate(params):
            assert cur_param in params, f"Must provide {cur_param} parameter constraint."
            if isinstance(param_constraints[cur_param], tuple) or isinstance(param_constraints[cur_param], list) or isinstance(param_constraints[cur_param], np.ndarray):
                cons = np.array(param_constraints[cur_param])
                assert cons.size == 2, "Must provide the constraint as either a single value (equality) or a tuple (bounded region)."
                constrs.append((cons[0], cons[1]))
                x0.append(0.5*(cons[0]+cons[1]))
            else:
                constrs.append((param_constraints[cur_param], param_constraints[cur_param]))
                x0.append(param_constraints[cur_param])
        sol = scipy.optimize.minimize(func, x0, bounds=constrs, method='Nelder-Mead')

        g = XmonDesigner._g_hertz(self.resonator.get_res_frequency(), sol.x[1], sol.x[2], sol.x[0], self.resonator.get_res_capacitance(), self.resonator.get_res_inductance())
        chi = XmonDesigner._chi_hertz(self.resonator.get_res_frequency(), sol.x[1], sol.x[2], sol.x[0], self.resonator.get_res_capacitance(), self.resonator.get_res_inductance())
        anh = XmonDesigner._anharmonicity_hertz(sol.x[1], sol.x[2], self.resonator.get_res_capacitance())

        self.resonator.print()
        print("Qubit:")
        print(f"\tFrequency: {QUtilities.add_units(sol.x[0])}Hz")
        print(f"\tAnharmonicity: {QUtilities.add_units(anh)}Hz")
        print(f"\tg: {QUtilities.add_units(g)}Hz")
        print(f"\tDelta: {QUtilities.add_units(sol.x[0]-self.resonator.get_res_frequency())}Hz")
        print(f"\tchi: {QUtilities.add_units(chi)}Hz")
        print(f"\tCg: {QUtilities.add_units(sol.x[1])}F")
        print(f"\tCJ: {QUtilities.add_units(sol.x[2])}F")
        

if __name__ == '__main__':
    XmonDesigner(ResonatorHalfWave(10.5e9)).optimise({'fQubit':9.5e9, 'C_g':(0.5e-15,10e-15), 'C_J':(10e-15,100e-15), 'chi':(-3e6,-1.7e6)})
