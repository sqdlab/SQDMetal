# Copyright 2025 Prasanna Pakkiam
# SPDX-License-Identifier: Apache-2.0

import numpy as np
import scipy.optimize
from SQDMetal.Utilities.GenUtilities import GenUtilities


##################################################################################
####################################RESONATORS####################################

class ResonatorBase:
    def get_res_capacitance(self):
        raise NotImplementedError()
    def get_res_inductance(self):
        raise NotImplementedError()
    def get_res_impedance(self):
        raise NotImplementedError()
    def get_res_frequency(self):
        raise NotImplementedError()
    def is_res_parallelLC(self):
        raise NotImplementedError()
    def print(self):
        print("Resonator:")
        print(f"\tFrequency: {GenUtilities.add_units(self.get_res_frequency())}Hz")
        print(f"\tInductance: {GenUtilities.add_units(self.get_res_inductance())}H")
        print(f"\tCapacitance: {GenUtilities.add_units(self.get_res_capacitance())}F")

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
    
    def is_res_parallelLC(self):
        return not self.shorted

class ResonatorQuarterWave(ResonatorBase):
    def __init__(self, f0, shorted=True, impedance=50):
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
    
    def is_res_parallelLC(self):
        return self.shorted

class ResonatorLC(ResonatorBase):
    def __init__(self, f0=None, L=None, C=None):
        if f0 is None:
            assert L is not None and C is not None, "Must supply L and C if omitting f0."
            self.L = L
            self.C = C
            self.f0 = 1/(2*np.pi * np.sqrt(L*C))
        elif L is None:
            assert f0 is not None and C is not None, "Must supply f0 and C if omitting L."
            self.f0 = f0
            self.C = C
            self.L = 1/((2*np.pi*f0)**2 * C)
        elif C is None:
            assert f0 is not None and L is not None, "Must supply f0 and L if omitting C."
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
    
    def is_res_parallelLC(self):
        return True


##################################################################################
##################################TRANSMON QUBIT##################################

class TransmonBase:
    def _eff_CJ(self, C_J, C_g, C_r):
        return 0.5 * (2*C_J + 2*C_g/C_r*(C_J+C_r)) / (1+C_g/C_r)
    
    def _eff_Cr(self, C_J, C_g, C_r):
        return 0.5 * (2*C_r + 2*C_g/C_J*(C_J+C_r)) / (1+C_g/C_J)
    
    def _EJ_from_QubitFreqHertz(self, QubiFreqHertz, C_J, C_g, C_r):
        elem = 1.60217663e-19
        h = 6.62607015e-34
        CJeff = self._eff_CJ(C_J, C_g, C_r)
        Ec = elem**2 / (2*CJeff)
        return (QubiFreqHertz * h + Ec)**2 / (8*Ec)
    
    def _g_hertz(self, C_g, C_J, QubiFreqHertz, C_r, L_r):
        elem = 1.60217663e-19
        h = 6.62607015e-34
        Z_0 = 376.730313668
        alpha = 0.0072973525693
        CJeff = self._eff_CJ(C_J, C_g, C_r)
        Ec = elem**2 / (2*CJeff) / h
        E_J_hertz = self._EJ_from_QubitFreqHertz(QubiFreqHertz, C_J, C_g, C_r) / h

        Creff = self._eff_Cr(C_J, C_g, C_r)
        Z_r = np.sqrt(L_r/Creff)
        w_res = 1/np.sqrt(L_r*Creff)

        return (w_res * C_g)/(C_J+C_g/C_r*(C_J+C_r)) * (0.5*E_J_hertz/Ec)**(1/4) * np.sqrt(Z_r/Z_0) * np.sqrt(2*np.pi*alpha) / (2*np.pi)
    
    def _chi_hertz(self, C_g, C_J, QubiFreqHertz, f_res, C_r, L_r):
        elem = 1.60217663e-19
        h = 6.62607015e-34
        g_hertz = self._g_hertz(C_g, C_J, QubiFreqHertz, C_r, L_r)
        CJeff = self._eff_CJ(C_J, C_g, C_r)
        Ec = (2*elem)**2 / (2*CJeff)
        delta = QubiFreqHertz - f_res

        return -(Ec/(4*h) * g_hertz**2) / (delta * (delta - Ec/(4*h)))

    def EJonEC(self, CJeff, fQubitHertz):
        elem = 1.60217663e-19
        h = 6.62607015e-34
        Ec = (elem)**2 / (2*CJeff)
        return fQubitHertz / (Ec/h)

    def _anharmonicity_hertz(self, C_g, C_J, C_r):
        elem = 1.60217663e-19
        h = 6.62607015e-34
        CJeff = self._eff_CJ(C_J, C_g, C_r)
        Ec = (elem)**2 / (2*CJeff)
        return Ec/h

    def get_free_params(self):
        raise NotImplementedError()

    def parse_params(self, param_constraints):
        default_params = self.get_free_params()
        x0 = []
        constrs = []
        mapped_params = {}
        for m, cur_param in enumerate(default_params):
            cur_param_constraints = param_constraints if cur_param in param_constraints else default_params
            #
            if isinstance(cur_param_constraints[cur_param], tuple) or isinstance(cur_param_constraints[cur_param], list) or isinstance(cur_param_constraints[cur_param], np.ndarray):
                cons = np.array(cur_param_constraints[cur_param])
                assert cons.size == 2, "Must provide the constraint as either a single value (equality) or a tuple (bounded region)."
                constrs.append((cons[0], cons[1]))
                x0.append(0.5*(cons[0]+cons[1]))
            else:
                constrs.append((cur_param_constraints[cur_param], cur_param_constraints[cur_param]))
                x0.append(cur_param_constraints[cur_param])
            mapped_params[cur_param] = m
        return x0, constrs, mapped_params
    
    def chk_constr_resp(self, constr_name, x, constrs, mps):
        constrs = constrs[mps[constr_name]]
        if isinstance(constrs, (list, tuple, np.ndarray)):
            leVal = x[mps[constr_name]]
            respected = leVal / constrs[0] - 1 >= -1e-6  and leVal / constrs[1] - 1 <= 1e-6
        else:
            respected = np.abs(constrs / x[mps[constr_name]] - 1) < 1e-6
        if respected:
            return ""
        else:
            return f"[FAIL - should be {constrs}]"

class XmonDesigner(TransmonBase):
    def __init__(self, resonator):
        assert resonator.is_res_parallelLC(), "Resonator must be in a parallel-coupled arrangment."
        self.resonator = resonator
    
    def get_free_params(self):
        return {'fQubit':(1,20e9), 'C_g':(1e-18,1e-9), 'C_J':(1e-18,1e-9), 'chi':(-10e9,-1), 'Ej/Ec':(0.01,1000)}

    def optimise(self, param_constraints):
        x0, constrs, mps = self.parse_params(param_constraints)

        fres, Cres, Lres = self.resonator.get_res_frequency(), self.resonator.get_res_capacitance(), self.resonator.get_res_inductance()

        def func(x):
            return abs(self._chi_hertz(x[mps['C_g']], x[mps['C_J']], x[mps['fQubit']], fres, Cres, Lres) / x[mps['chi']] - 1) + abs(self.EJonEC(self._eff_CJ(x[mps['C_J']],x[mps['C_g']], Cres), x[mps['fQubit']]) / x[mps['Ej/Ec']] - 1)

        sol = scipy.optimize.minimize(func, x0, bounds=constrs, method='TNC', tol=1e-5)
        x = sol.x


        g = self._g_hertz(x[mps['C_g']], x[mps['C_J']], x[mps['fQubit']], Cres, Lres)
        chi = self._chi_hertz(x[mps['C_g']], x[mps['C_J']], x[mps['fQubit']], fres, Cres, Lres)
        anh = self._anharmonicity_hertz(sol.x[1], sol.x[2], Cres)

        x[mps['chi']] = chi
        x[mps['Ej/Ec']] = self.EJonEC(self._eff_CJ(x[mps['C_J']],x[mps['C_g']], Cres), x[mps['fQubit']])

        sigfigs = 5

        print(f"Cost Function Error: {func(x)/2*100}%")
        self.resonator.print()
        print("Qubit:")
        print(f"\tFrequency: {GenUtilities.add_units(sol.x[0],sigfigs)}Hz")
        print(f"\tAnharmonicity: {GenUtilities.add_units(anh,sigfigs)}Hz")
        print(f"\tg: {GenUtilities.add_units(g,sigfigs)}Hz")
        print(f"\tDelta: {GenUtilities.add_units(x[mps['fQubit']]-fres,sigfigs)}Hz")
        print(f"\tchi: {GenUtilities.add_units(chi,sigfigs)}Hz")
        print(f"\tCg: {GenUtilities.add_units(x[mps['C_g']],sigfigs)}F")
        print(f"\tCJ: {GenUtilities.add_units(x[mps['C_J']],sigfigs)}F")
        print(f"\tEj/Ec: {GenUtilities.add_units(x[mps['Ej/Ec']],sigfigs)}")

class FloatingTransmonDesigner(TransmonBase):
    def __init__(self, resonator):
        assert resonator.is_res_parallelLC(), "Resonator must be in a parallel-coupled arrangment."
        self.resonator = resonator
    """
    Symbol conventions
    C_g1 : Pad 1 to resonator
    C_g2 : Pad 2 to resonator
    C_q1 : Pad 1 to ground
    C_q2 : Pad 2 to ground
    C_J  : Pad mutual inductance + JJ capacitance 
    """

    def get_free_params(self):
        return {'fQubit':(1,20e9),
                'C_q1':(1e-18,1e-9), 'C_q2':(1e-18,1e-9),
                'C_g1':(1e-18,1e-9), 'C_g2':(1e-18,1e-9),
                'C_J':(1e-18,1e-9),
                'chi':(-10e9,-1), 'Ej/Ec':(0.01,1000),
                'C_sigma': (1e-18,1e-6),
                'beta': (0,1)}
    
    def _eff_Cg(self, C_q1, C_q2, C_g1, C_g2):
        return (C_g1*C_q2 - C_g2*C_q1) / (C_q1 + C_q2 + C_g1 + C_g2)
    def _eff_Cq(self, C_q1, C_q2, C_g1, C_g2):
        return (C_g1*C_g2 + 2*C_g2*C_q1 + C_q1*C_q2) / (C_q1 + C_q2 + C_g1 + C_g2)
    def _Csigma(self, C_q1, C_q2, C_g1, C_g2, C_J):
        return (C_q1+C_g1)*(C_q2+C_g2) / (C_q1 + C_q2 + C_g1 + C_g2) + C_J
    def _beta(self, C_q1, C_q2, C_g1, C_g2, C_J):
        return (C_g1*C_q2 - C_g2*C_q1) / ((C_q1+C_g1)*(C_q2+C_g2) + (C_q1 + C_q2 + C_g1 + C_g2)*C_J)

    def optimise(self, param_constraints, print_results=True):
        x0, constrs, mps = self.parse_params(param_constraints)

        fres, Cres, Lres = self.resonator.get_res_frequency(), self.resonator.get_res_capacitance(), self.resonator.get_res_inductance()

        def func(x):
            return abs(self._chi_hertz(self._eff_Cg(*x[1:5]), self._eff_Cq(*x[1:5])+x[5], x[0], fres, Cres, Lres) / x[mps['chi']] - 1) \
                      + abs(self.EJonEC(self._eff_CJ(self._eff_Cq(*x[1:5])+x[5], self._eff_Cg(*x[1:5]), Cres), x[mps['fQubit']]) / x[mps['Ej/Ec']] - 1) \
                      + abs(self._Csigma(*x[1:6]) / x[mps['C_sigma']] - 1) \
                      + abs(self._beta(*x[1:6]) / x[mps['beta']] - 1)

        x0[mps['C_sigma']] = self._Csigma(*x0[1:6])
        x0[mps['beta']] = self._beta(*x0[1:6])
        sol = scipy.optimize.minimize(func, x0, bounds=constrs, method='Nelder-Mead')
        x = sol.x

        g = self._g_hertz(self._eff_Cg(*sol.x[1:5]), self._eff_Cq(*sol.x[1:5])+sol.x[5], sol.x[0], Cres, Lres)
        chi = self._chi_hertz(self._eff_Cg(*sol.x[1:5]), self._eff_Cq(*sol.x[1:5])+sol.x[5], sol.x[0], fres, Cres, Lres)
        anh = self._anharmonicity_hertz(self._eff_Cg(*sol.x[1:5]), self._eff_Cq(*sol.x[1:5])+sol.x[5], Cres)

        x[mps['beta']] = self._beta(*x[1:6])
        x[mps['C_sigma']] = self._Csigma(*x[1:6])
        x[mps['Ej/Ec']] = self.EJonEC(self._eff_CJ(self._eff_Cq(*x[1:5])+x[5], self._eff_Cg(*x[1:5]), Cres), x[mps['fQubit']])

        sigfigs = 5
        cost_err = func(sol.x)/2*100
        if print_results:
            print(f"Cost Function Error: {cost_err}%")
            self.resonator.print()
            print("Qubit:")
            print(f"\tFrequency: {GenUtilities.add_units(sol.x[0],sigfigs)}Hz" + self.chk_constr_resp('fQubit', x,constrs,mps))
            print(f"\tAnharmonicity: {GenUtilities.add_units(anh,sigfigs)}Hz")
            print(f"\tg: {GenUtilities.add_units(g,sigfigs)}Hz")
            print(f"\tDelta: {GenUtilities.add_units(x[mps['fQubit']]-fres,sigfigs)}Hz")
            print(f"\tchi: {GenUtilities.add_units(chi,sigfigs)}Hz")
            print(f"\tCq1: {GenUtilities.add_units(x[mps['C_q1']],sigfigs)}F" + self.chk_constr_resp('C_q1', x,constrs,mps))
            print(f"\tCq2: {GenUtilities.add_units(x[mps['C_q2']],sigfigs)}F" + self.chk_constr_resp('C_q2', x,constrs,mps))
            print(f"\tCg1: {GenUtilities.add_units(x[mps['C_g1']],sigfigs)}F" + self.chk_constr_resp('C_g1', x,constrs,mps))
            print(f"\tCg2: {GenUtilities.add_units(x[mps['C_g2']],sigfigs)}F" + self.chk_constr_resp('C_g2', x,constrs,mps))
            print(f"\tCJ: {GenUtilities.add_units(x[mps['C_J']],sigfigs)}F" + self.chk_constr_resp('C_J', x,constrs,mps))
            print(f"\tEj/Ec: {GenUtilities.add_units(x[mps['Ej/Ec']],sigfigs, True)}" + self.chk_constr_resp('Ej/Ec', x,constrs,mps))
            print(f"\tCΣ: {GenUtilities.add_units(x[mps['C_sigma']],sigfigs)}F" + self.chk_constr_resp('C_sigma', x,constrs,mps))
            print(f"\tbeta: {GenUtilities.add_units(x[mps['beta']],sigfigs, True)}" + self.chk_constr_resp('beta', x,constrs,mps))

        results = {
            "cost_error_%": {cost_err},
            "fQubit_Hz": x[mps["fQubit"]],
            "anh_Hz": anh,
            "g_Hz": g,
            "Delta_Hz": x[mps["fQubit"]] - fres,
            "chi_Hz": chi,
            "C_q1_F": x[mps["C_q1"]],
            "C_q2_F": x[mps["C_q2"]],
            "C_g1_F": x[mps["C_g1"]],
            "C_g2_F": x[mps["C_g2"]],
            "C_J_F": x[mps["C_J"]],
            "Ej_Ec_ratio": x[mps["Ej/Ec"]],
            "C_sigma_F": x[mps["C_sigma"]],
            "beta": x[mps["beta"]],
        }
        # return results as dictionary
        return results

if __name__ == '__main__':
    XmonDesigner(ResonatorHalfWave(10.5e9)).optimise({
    'fQubit':5.5e9,
    'C_g':8.1e-15,#(0.1e-15,8.1e-15),
    'C_J':22e-15,
    'chi':(-10e6,-0.0e6),
    'Ej/Ec':(1,100)})
    FloatingTransmonDesigner(ResonatorHalfWave(7.5e9)).optimise({'fQubit':(1e9, 10e9), 
                                                                  'C_q1':(30.862e-15),
                                                                  'C_q2':(31.974e-15),
                                                                  'C_g1':(10.910e-15),
                                                                  'C_g2':(0.916e-15),
                                                                  'C_J':(42.461e-15, 55e-15),
                                                                  'chi':(-0.9e6,-0.2e6),
                                                                  'Ej/Ec':20})
