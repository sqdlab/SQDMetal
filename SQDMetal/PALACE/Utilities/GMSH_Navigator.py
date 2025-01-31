import gmsh
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import collections  as mc
import numpy as np

class GMSH_Navigator:
    def __init__(self, mesh_file, units=0.001, gmsh_verbosity=1):
        gmsh.initialize()
        gmsh.option.setNumber('General.Verbosity', gmsh_verbosity)
        gmsh.option.setNumber("General.Terminal", gmsh_verbosity)
        gmsh.model.remove()
        gmsh.finalize()
        gmsh.initialize()
        gmsh.option.setNumber('General.Verbosity', gmsh_verbosity)
        gmsh.option.setNumber("General.Terminal", gmsh_verbosity)
        gmsh.open(mesh_file)

        #
        self._file_path = mesh_file
        self._units = units
        #
        lePhysGroups = gmsh.model.getPhysicalGroups()
        self.physical_group_names_entities = [(gmsh.model.get_physical_name(*x), gmsh.model.getEntitiesForPhysicalGroup(*x), x) for x in lePhysGroups]
        #
        self.leCols = plt.rcParams['axes.prop_cycle'].by_key()['color']
        self.leColsRGB = [tuple(int(round(v*255)) for v in matplotlib.colors.to_rgb(c)) for c in self.leCols]

    def open_GUI(self, show_metals=True, hide_box=True, hide_substrate=True, hide_substrate_gaps=True):
        for m, cur_entity in enumerate(self.physical_group_names_entities):
            name, tags, physical_group = cur_entity
            leEntityTags = [(physical_group[0], x) for x in tags]
            #
            if name.startswith('metal'):
                visibility = 1 if show_metals else 0
            elif name.startswith('rf_port'):
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
            gmsh.model.setColor(leEntityTags, *self.leColsRGB[m%len(self.leColsRGB)])
        gmsh.option.setNumber('Mesh.SurfaceFaces', 1)
        gmsh.fltk.run()

    def export_to_png(self, filt_names = ['metal', 'rf_port', 'dielectric_gaps'], file_path=""):
        #Get mesh-coordinates using code from:
        #   - https://stackoverflow.com/questions/73882211/how-to-get-all-edges-and-faces-triangles-in-a-mesh-and-their-nodes-3d-co-ordi
        #   - https://bthierry.pages.math.cnrs.fr/tutorial/gmsh/api/detail/
        #
        nodeTags, nodeCoords, _ = gmsh.model.mesh.getNodes()
        elementType = gmsh.model.mesh.getElementType("triangle", 1)
        # faceNodes = gmsh.model.mesh.getElementFaceNodes(elementType, 3)
        edgeNodes = gmsh.model.mesh.getElementEdgeNodes(elementType)
        #
        nodes = np.reshape(nodeCoords, (int(len(nodeCoords)/3), 3))
        # faces = np.reshape(faceNodes, (int(len(faceNodes)/3), 3))
        edges = np.reshape(edgeNodes, (int(len(edgeNodes)/2), 2))

        leNodes = {nodeTags[x]:nodes[x] for x in range(nodeTags.size)}

        fig, ax = plt.subplots(1)
        for m, cur_phys_group in enumerate(self.physical_group_names_entities):
            found_match = False
            for filt_name in filt_names:
                if filt_name in cur_phys_group[0]:
                    found_match = True
                    break                
            if not found_match:
                continue

            filtNodeTags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(*cur_phys_group[2])
            filtEdges = [[leNodes[x[0]][:2]*self._units, leNodes[x[1]][:2]*self._units] for x in edges if x[0] in filtNodeTags and x[1] in filtNodeTags]
            #Trick adapted from: https://stackoverflow.com/questions/21352580/plotting-numerous-disconnected-line-segments-with-different-colors
            lc = mc.LineCollection(filtEdges, linewidths=1, color=self.leCols[m % len(self.leCols)])
            ax.add_collection(lc)
            ax.autoscale()
        
        if file_path == "":
            file_path = '.'.join(self._file_path.split('.')[:-1]) + '.png'
        fig.savefig(file_path)
        plt.close(fig)
