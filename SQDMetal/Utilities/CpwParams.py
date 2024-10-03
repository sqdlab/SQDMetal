import numpy as np
import scipy.special
import scipy.optimize
from SQDMetal.Utilities.QUtilities import QUtilities
from SQDMetal.Utilities.Materials import Material

class CpwParams:
    def __init__(self, rel_permittivity, dielectric_thickness):
        self.rel_permittivity = rel_permittivity
        self.dielectric_thickness = dielectric_thickness
    
    @classmethod
    def fromQDesign(cls, design, chip_name='main'):
        if isinstance(design.chips[chip_name].material, str):
            matr = design.chips[chip_name].material
        else:
            matr = design.chips[chip_name].material[0]
        er = Material(matr).permittivity
        h = QUtilities.parse_value_length(design.chips[chip_name].size['size_z'])
        return CpwParams(er, h)
    
    def fromMaterial(cls, material_obj, h):
        #That is, it can be called as: CpwParams.fromMaterial("sapphire", 500e-6) or CpwParams.fromMaterial(Material(...), 500e-6) etc...
        if isinstance(material_obj, str):
            Er = Material(material_obj).permittivity
        elif isinstance(material_obj, Material):
            Er = material_obj.permittivity
        else:
            assert False, "The first argument \"material_obj\" must be a Material object type or a string for the material name."
        return CpwParams(Er, h)

    @staticmethod
    def calc_impedance(tr_wid, tr_gap, er, h):
        a = tr_wid
        b = tr_wid + 2.0*tr_gap
        k = a/b
        k1 = np.tanh(np.pi*a/(4.0*h)) / np.tanh(np.pi*b/(4.0*h))
        kp = np.sqrt(1-k**2)
        k1p = np.sqrt(1-k1**2)

        kp_k1_on_k_k1p = scipy.special.ellipk(kp**2)*scipy.special.ellipk(k1**2) / (scipy.special.ellipk(k**2)*scipy.special.ellipk(k1p**2))
        e_eff = (1.0 + er*kp_k1_on_k_k1p)/(1.0 + kp_k1_on_k_k1p)

        z0 = 60.0*np.pi/np.sqrt(e_eff) * 1.0 / (scipy.special.ellipk(k**2)/scipy.special.ellipk(kp**2) + scipy.special.ellipk(k1**2)/scipy.special.ellipk(k1p**2))

        return z0

    def get_gap_from_width(self, trace_width, target_impedance=50):
        er = self.rel_permittivity
        h = self.dielectric_thickness

        x0 = [trace_width]

        res = scipy.optimize.minimize(lambda x: np.abs(CpwParams.calc_impedance(trace_width, x[0], er, h)-target_impedance), x0, method='Nelder-Mead', tol=1e-6)

        return res.x[0]

    def get_width_from_gap(self, trace_gap, target_impedance=50):
        er = self.rel_permittivity
        h = self.dielectric_thickness

        x0 = [trace_gap]

        res = scipy.optimize.minimize(lambda x: np.abs(CpwParams.calc_impedance(x[0], trace_gap, er, h)-target_impedance), x0, method='Nelder-Mead', tol=1e-6)

        return res.x[0]

    def get_width_from_space(self, gnd_plane_gap, target_impedance=50):
        #Space is the distance between the CPW ground planes

        er = self.rel_permittivity
        h = self.dielectric_thickness

        x0 = [gnd_plane_gap*0.5]

        res = scipy.optimize.minimize(lambda x: np.abs(CpwParams.calc_impedance(x[0], 0.5*(gnd_plane_gap-x[0]), er, h)-target_impedance), x0, method='Nelder-Mead', tol=1e-6)

        return res.x[0]
