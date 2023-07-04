# Author: Prasanna Pakkiam
# Creation Date: 04/07/2023
# Description: A simple structure class to house the parameters/properties of different materials/substrates.
#              Currently, it just holds the permittivities and permeabilities.

class Material:
    def __init__(self, name="", **kwargs):
        if name == "":
            self.permittivity = kwargs.pop('permittivity', 1)   #Relative Permittivity
            self.permeability = kwargs.pop('permeability', 1)   #Relative Permeability
            for m in kwargs:
                assert False, f"Option {m} is not recognised"
        else:
            #The names can have capitals and can be delimited via spaces, underscores or camelcase
            orig_name = name
            name = name.tolower().replace(" ", "").replace("_", "")
            ################################
            #
            #SILICON
            #
            #References: https://ieeexplore.ieee.org/document/1717770
            if name == "silicon" or name == "siliconcold" or name == "siliconcryo":
                self.permittivity = 11.45
                self.permeability = 1
            elif name == "siliconwarm":
                self.permittivity = 11.7
                self.permeability = 1
            ################################
            #
            #SAPPHIRE
            #
            #References: https://indico.fnal.gov/event/4162/attachments/54289/64791/SapphireAndOtherMaterialsAt_100K.pdf
            #TODO: Look into what the default orientation of sapphire should be considered.
            elif name == "sapphire" or name == "sapphireperp" or name == "sapphireperpendicular":
                self.permittivity = 9.4
                self.permeability = 1
            elif name == "sapphireparallel":
                self.permittivity = 11.35
                self.permeability = 1
            ################################
            else:
                assert False, f"Material \"{orig_name}\" unrecognised."

