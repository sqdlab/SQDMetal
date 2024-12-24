import gmsh

class GMSH_Navigator:
    def __init__(self, mesh_file):
        gmsh.initialize()
        gmsh.model.remove()
        gmsh.finalize()
        gmsh.initialize()
        gmsh.open(mesh_file)
        #
        lePhysGroups = gmsh.model.getPhysicalGroups()
        self.physical_group_names_entities = [(gmsh.model.get_physical_name(*x), gmsh.model.getEntitiesForPhysicalGroup(*x), x) for x in lePhysGroups]

    def open_GUI(self, show_metals=True, hide_box=True, hide_substrate=True, hide_substrate_gaps=True):
        for cur_entity in self.physical_group_names_entities:
            name, tags, physical_group = cur_entity
            leEntityTags = [(physical_group[0], x) for x in tags]
            #
            if name.startswith('metal'):
                visibility = 1 if show_metals else 0
            elif name == 'far_field' or name == 'air_box':
                visibility = 0 if hide_box else 1
            elif name == 'dielectric_substrate':
                visibility = 0 if hide_substrate else 1
            elif name == 'dielectric_gaps':
                visibility = 0 if hide_substrate_gaps else 1
            else:
                visibility = 0
            #
            gmsh.model.setVisibility(leEntityTags, visibility, True)
        gmsh.fltk.run()

