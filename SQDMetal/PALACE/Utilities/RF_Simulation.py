from Utilities.QUtilities import QUtilities
import gmsh
import numpy as np

class RF_Simulation:

    config_ports = []
    ports_index = 1

    def __init__(self, name, ports_dict, user_options, jj_dict, hpc_options = {}):
        self.name = name
        self.ports_dict = ports_dict
        self.user_options = user_options
        self.jj_dict = jj_dict
        self.hpc_options = hpc_options


    def prepare_simulation(self):

        #dictionary for physical groups/bundary conditions
        physical_groups = {}

        physical_group_ids = gmsh.model.get_physical_groups()
        for _,group in enumerate(physical_group_ids):
            group_name = gmsh.model.get_physical_name(group[0], group[1])
            if group_name == 'metals':
                physical_groups['metals'] = group[1]
            elif group_name == 'dielectric_gaps':
                physical_groups['dielectric_gaps'] = group[1]
            elif group_name == 'ground_plane':
                physical_groups['ground_plane'] = group[1]
            elif group_name == 'dielectric_substrate':
                physical_groups['dielectric_substrate'] = group[1]
            elif group_name == 'air_box':
                physical_groups['air_box'] = group[1]
            elif group_name == 'far_field':
                physical_groups['far_field']  = group[1]

        return physical_groups


    def _process_ports_for_config_file(self):
        #Assumes that ports is a dictionary that contains the port names (with separate keys with suffixes a and b for multi-element ports)
        #where each value is a list of element IDs corresponding to the particular port...

        #If there are ports for the launch pads in the dictionary process them
        if self.ports_dict:
            for _, port in enumerate(self.ports_dict):
                elements = list(self.ports_dict[port].keys())
                
                ports_json = {}

                ports_json['Index'] = self.ports_index
                self.ports_index += 1

                ports_json['R'] = 50

                ports_json['Elements'] = [
                        {
                        "Attributes": [self.ports_dict[port][elements[0]][0]],
                        "Direction": self.ports_dict[port][elements[0]][1] #vec_field + [0]    #TODO: Change this after SPACK updates Palace beyond August 2023
                        },
                        {
                        "Attributes": [self.ports_dict[port][elements[1]][0]],
                        "Direction": self.ports_dict[port][elements[1]][1] #[-x for x in vec_field] + [0]
                        }
                    ]
                
                if self.ports_index == 1: #Excitation is always on first port in the ports dictionary
                    ports_json['Excitation'] = True

                self.config_ports.append(ports_json)
        
        #If there are junctions in the dictionary process them
        if self.jj_dict:
            for jj in self.jj_dict.keys():
            
                jj_physical_group = self.jj_dict[jj][0]
                jj_inductance = self.jj_dict[jj][1][0]

                jj_json = {}
                jj_json['Index'] = self.ports_index
                self.ports_index += 1
                jj_json['L'] = jj_inductance           #jj inductance in nH
                jj_json['Attributes'] = [jj_physical_group]
                jj_json['Direction'] = "+Y"
                self.config_ports.append(jj_json)

        return self.config_ports
    
        
        
        

    
    
    
        

                
