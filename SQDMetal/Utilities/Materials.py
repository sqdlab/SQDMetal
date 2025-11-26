# Copyright 2025 Prasanna Pakkiam
# SPDX-License-Identifier: Apache-2.0

# Author: Prasanna Pakkiam
# Creation Date: 04/07/2023
# Description: A simple structure class to house the parameters/properties of different materials/substrates.
#              Currently, it just holds the permittivities and permeabilities.

class Material:
    def __init__(self, name="", **kwargs):
        if name == "":
            self.permittivity = kwargs.pop('permittivity', 1)   #Relative Permittivity
            self.permeability = kwargs.pop('permeability', 1)   #Relative Permeability
            self.loss_tangent = kwargs.pop('loss_tangent', 0)   #Dielectric Loss Tangent
            for m in kwargs:
                assert False, f"Option {m} is not recognised"
        else:
            #The names can have capitals and can be delimited via spaces, underscores or camelcase
            orig_name = name
            name = name.lower().replace(" ", "").replace("_", "")
            ################################
            #
            #SILICON
            #
            #References: https://ieeexplore.ieee.org/document/1717770
            if name == "silicon" or name == "siliconcold" or name == "siliconcryo":
                self.permittivity = 11.45
                self.permeability = 1
                self.loss_tangent = 2.7e-6 #PHYS. REV. APPLIED 18, 034013 (2022) Measurement of the Low-Temperature Loss Tangent of High-Resistivity Silicon Using a High-Q Superconducting Resonator
            elif name == "siliconwarm":
                self.permittivity = 11.7
                self.permeability = 1
                self.loss_tangent = 2.7e-6  #TODO: Change this to actual value!
            ################################
            #
            #SAPPHIRE
            #
            #References: https://indico.fnal.gov/event/4162/attachments/54289/64791/SapphireAndOtherMaterialsAt_100K.pdf
            #TODO: Look into what the default orientation of sapphire should be considered.
            elif name == "sapphire" or name == "sapphireperp" or name == "sapphireperpendicular":
                self.permittivity = 9.4
                self.permeability = 1
                self.loss_tangent = 1.2e-6  #TODO: Change this to actual value!
            elif name == "sapphireparallel":
                self.permittivity = 11.35
                self.permeability = 1
                self.loss_tangent = 1.2e-6  #TODO: Change this to actual value!
            ################################
            else:
                assert False, f"Material \"{orig_name}\" unrecognised."

class MaterialInterface:
    def __init__(self, name="", **kwargs):
        if name == "":
            self.permittivity = kwargs.pop('permittivity', 1)   #Relative Permittivity
            self.permeability = kwargs.pop('permeability', 1)   #Relative Permeability
            self.loss_tangent = kwargs.pop('loss_tangent', 0)   #Dielectric Loss Tangent
            for m in kwargs:
                assert False, f"Option {m} is not recognised"
        else:
            #The names can have capitals and can be delimited via spaces, underscores or camelcase
            orig_name = name
            name = name.lower().replace(" ", "").replace("_", "").replace("-", "")
            ################################
            #
            #SILICON-AIR
            #
            #References
            #[1] https://www.epjap.org/articles/epjap/pdf/2004/10/ap04128.pdf
            #[2] https://ieeexplore.ieee.org/abstract/document/6419
            #[3] https://link.springer.com/article/10.1134/S106378341402022X
            if name == "siliconair" or name == "siliconvacuum" or name == "airsilicon" or name == "vacuumsilicon" or name == "siliconoxide":
                self.permittivity = 3.85 #10.0
                self.permeability = 1
                self.loss_tangent = 0.006 #0.9e-3
            ################################
            #
            #ALUMINIUM-OXIDE
            #
            elif name == "aluminiumair" or name == "aluminiumvacuum" or name == "airaluminium" or name == "vacuumaluminium" or name == "aluminiumoxide":
                self.permittivity = 10.0
                self.permeability = 1
                self.loss_tangent = 0.5e-3  #TODO: Change this to actual value!
            ################################
            #
            #SILICON-ALUMINIUM
            #
            elif name == "siliconaluminium" or name == "aluminiumsilicon":
                self.permittivity = 10.0
                self.permeability = 1
                self.loss_tangent = 0.3e-3
            ################################
            #
            #SILICON-TANTALUM
            #
            #References
            #[1] https://www.epjap.org/articles/epjap/pdf/2004/10/ap04128.pdf
            #[2] https://ieeexplore.ieee.org/abstract/document/6419781
            #[3] https://link.springer.com/article/10.1134/S106378341402022X
            elif name == "silicontantalum" or name == "tantalumsilicon":
                self.permittivity = 3.85 # assuming this interface is silicon oxide
                self.permeability = 1
                self.loss_tangent = 0.006
            ################################
            #
            #TANTALUM-OXIDE
            #
            #References
            #[1] https://pubs.aip.org/aip/jap/article/54/11/6502/13154/Electrical-properties-of-amorphous-tantalum
            #[2] https://ieeexplore.ieee.org/document/1717770
            elif name == "tantalumair" or name == "tantalumvacuum" or name == "airtantalum" or name == "vacuumtantalum" or name == "tantalumoxide":
                self.permittivity = 20.0
                self.permeability = 1
                self.loss_tangent = 0.005
            else:
                assert False, f"Material \"{orig_name}\" unrecognised."

class MaterialConductor:
    def __init__(self, name="", **kwargs):
        if name == "":
            self.conductivity = kwargs.pop('conductivity', 1)   #In S/m
            self.permeability = kwargs.pop('permeability', 1)   #Relative Permeability
            for m in kwargs:
                assert False, f"Option {m} is not recognised"
        else:
            #The names can have capitals and can be delimited via spaces, underscores or camelcase
            orig_name = name
            name = name.lower().replace(" ", "").replace("_", "").replace("-", "")
            ################################
            #
            #OFHC (DEOXYGENATED COPPER)
            #
            #References
            #[1] 
            if name == "deoxygenatedcopper" or name == "ofhc" or name == "oxygenfreecopper":
                self.conductivity = 1/2e-9
                self.permeability = 1
            ################################
            else:
                assert False, f"Material \"{orig_name}\" unrecognised."

