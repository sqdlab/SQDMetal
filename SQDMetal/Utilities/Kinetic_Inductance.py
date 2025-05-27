import scipy
import scipy.integrate as integrate
import numpy as np

class Kinetic_Inductance:
    
    def __init__(self):
        pass

    def calculate_kinetic_inductance(self, coherence_len, delta_0, normal_cond, film_thickness, mean_free_path, t_crit, temp, omega):
        
        #constants
        k_b = 1.380649e-23               #boltzman constant
        hbar = 1.05457e-34               #reduced plank's constant, h/2*pi
        q = 1.60218e-19                  #fundamental charge
        mu_0 = 1.256637e-6               #vacuum permeability

        #calculate quantities needed for numerical integation
        fermi_velocity = np.pi * coherence_len * (delta_0*q) / hbar #fermi velovity, m/s
        tau = mean_free_path / fermi_velocity #scattering time
        delta = (delta_0 * q) * np.sqrt( np.cos((np.pi/2) * (temp/t_crit)**2)) #temperature dependent gap energy - related to the gap energy at zero temp (delta_0)

        #Calculate complex conductivity:
        x = (hbar * omega) / (2 * delta)
        y = hbar / (2*tau*(delta))   # Impurity scattering
        tt = temp/t_crit  # T / Tc
        sigma = self.calc_complex_cond(x, y, tt, delta_0, t_crit)
        
        #calculate effective penetration depth
        lambda_eff = np.sqrt(1/(mu_0 * omega * sigma.imag * normal_cond))
        
        #calculate surface impedance and surface inductance
        Z_s = mu_0 * omega * lambda_eff * (1 / np.tanh(film_thickness/lambda_eff))
        L_s = Z_s / omega

        #print values
        print(f"Complex conductivity: σ = {sigma} S")
        print(f"Surface Impedance: Z = {Z_s} Ω/□")
        print(f"Surface Inductance: L = {L_s} H/□" )
        print(f"Effective Penetrattion Depth: λ_eff = {lambda_eff/1e-9} nm")

        return L_s, Z_s, sigma

    def calc_complex_cond(self, x, y, tt, delta_0, t_crit):
        """
        Calculation of complex conductivity sigma(q=0, omega) for superconductors
        using BCS-theory.
        
        Input:
            x  = omega / (2 * Delta)
            y  = l / (2 * Delta * tau)
            tt = temperature / Tc
            delta_0 = gap energy at 0 K
            t_crit = critical temperature
        
        Output:
            s = sigma = sigma1 - j * sigma2
        """
        M = 500
        k_b = 1.380649e-23
        q = 1.60218e-19
        dl = 1.0 / M
        dx = 1.0 / int(M * max(1.0, np.sqrt(x)))
        doverk = (delta_0*q) / (k_b * t_crit) 
        t = tt / (doverk * 2 * np.sqrt(1 - tt) * (0.9963 + 0.7733 * tt))
        
        sl = 0j
        s2 = 0j
        s3 = 0j
        
        # First integral loop
        u_values = np.arange(dx * 0.5, 1.0 , dx)
        for u in u_values:
            s2 += self.GK(0.5 + (u/(1 - u) )** 2, x, y, t, 2) * u / ((1 - u) ** 3)
        s = s2 * dx * 2.0
        
        if x < 1:
            for u in u_values:
                sl += self.GK(0.5 + x * u ** 2 * (3 - u - u), x, y, t, 1) * u * (1 - u)
            s += sl * dx * 6.0 * x
        else:
            for u in u_values:
                s3 += self.GK(0.5 + (x - 1) * u ** 2 * (3 - u - u), x, y, t, 3) * u * (1 - u)
            #sl = 0j  # Reset sl
            for u in np.arange(dl * 0.5, 1.0 , dl):
                sl += self.GK(x - 0.5 + u ** 2 * (3 - u - u), x, y, t, 1) * u * (1 - u)
            s += (s3 * dx * (x - 1) + sl * dl) * 6.0
        
        s *= complex(0.0, y) * 0.5 / x
        return s
    
    def GK(self, e, x, y, t, k):
        """
        Computes the integral function GK for the given parameters.
        """
        cy = complex(0.0, y)
        
        if k == 1:
            p4 = complex(0.0, np.sqrt(0.25 - (e - x) ** 2))
            p2 = np.sqrt(e ** 2 - 0.25)
            c42 = (0.25 + e * (e - x)) / (p4 * p2 + 1e-20)
            th = np.tanh(e / (t + t + 0.001))
            return th * ((1 - c42) / (p4 + p2 + cy) - (1 + c42) / (p4 - p2 + cy))
        
        elif k == 2:
            p1 = np.sqrt((e + x) ** 2 - 0.25)
            p2 = np.sqrt(e ** 2 - 0.25)
            c12 = (0.25 + e * (e + x)) / (p1 * p2 + 1e-20)
            th = np.tanh(e / (t + t + 0.001))
            return (np.tanh((e + x) / (t + t + 0.001)) * ((1 + c12) / (p1 - p2 + cy) - (1 - c12) / (-p1 - p2 + cy))
                    + th * ((1 - c12) / (p1 + p2 + cy) - (1 + c12) / (p1 - p2 + cy)))
        
        elif k == 3:
            p3 = np.sqrt((e - x) ** 2 - 0.25)
            p2 = np.sqrt(e ** 2 - 0.25)
            c32 = (0.25 + e * (e - x)) / (p3 * p2 + 1e-20)
            th = np.tanh(e / (t + t + 0.001))
            return th * ((1 - c32) / (p3 + p2 + cy) - (1 + c32) / (p3 - p2 + cy))
        
        return 0j
