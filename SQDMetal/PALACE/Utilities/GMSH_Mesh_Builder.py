# Copyright 2025 Prasanna Pakkiam
# SPDX-License-Identifier: Apache-2.0

import gmsh

class GMSH_Mesh_Builder:

    def __init__(self, fine_meshes, user_options):
        
        self._fine_meshes = fine_meshes
        self.user_options = user_options


    def build_mesh(self):

        #turn off these parametes as we will determine the mesh element size from a mesh field
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

        #Idea is to:
        #   - Use the distance field to calculate distance from surface
        #   - Combine that with a threshold field to have minimum size when inside surface and slowly grow when outside
        #   - Take Minimum of all threshold fields to calculate mesh element size in current location!
        #
        #Note:
        #   - The algorithm needs to calculate how far it is from a surface - thus, it approximates it via a discretisation
        #   - The Sampling parameter only affects this discretisation - NOT mesh density! Lowering it may decrease mesh generation time
        #
        #c.f. https://gitlab.onelab.info/gmsh/gmsh/blob/gmsh_4_10_0/tutorials/python/t10.py
        #     https://gmsh.info/doc/texinfo/gmsh.html#Gmsh-mesh-size-fields

        fields_to_min = []
        for m, cur_fine_mesh in enumerate(self._fine_meshes):
            #create distance field and threshold field - see Gmsh documentation - python tutorial 10 for more detail
            gmsh.model.mesh.field.add("Distance", 2*m+1)
            #
            if 'region' in cur_fine_mesh:
                gmsh.model.mesh.field.setNumbers(2*m+1, "SurfacesList", [x[1] for x in cur_fine_mesh['region']])
            elif 'path' in cur_fine_mesh:
                gmsh.model.mesh.field.setNumbers(2*m+1, "CurvesList", cur_fine_mesh['path'])
            else:
                assert False, "The fine-mesh parameter does not have a supported object to reference/mesh - e.g. a region or path key..."
            #
            gmsh.model.mesh.field.setNumber(2*m+1, "Sampling", self.user_options.get('gmsh_dist_func_discretisation', 150))

            #A threshold field built in conjunction with the distance field
            thresh_field_id = 2*(m+1)
            gmsh.model.mesh.field.add("Threshold", thresh_field_id)
            gmsh.model.mesh.field.setNumber(thresh_field_id, "InField", 2*m+1)
            gmsh.model.mesh.field.setNumber(thresh_field_id, "SizeMin", cur_fine_mesh['mesh_min'])              #minimum size of mesh element
            gmsh.model.mesh.field.setNumber(thresh_field_id, "SizeMax", cur_fine_mesh['mesh_max'])              #maximum size of mesh element
            gmsh.model.mesh.field.setNumber(thresh_field_id, "DistMin", cur_fine_mesh['taper_dist_min'])        #distance from the gmsh surface to keep minimum mesh element size 
            gmsh.model.mesh.field.setNumber(thresh_field_id, "DistMax", cur_fine_mesh['taper_dist_max'])        #how far from the element until the max element size can be implemented 
            
            fields_to_min.append(thresh_field_id)
        min_field_id = thresh_field_id+1
        gmsh.model.mesh.field.add("Min", min_field_id)
        gmsh.model.mesh.field.setNumbers(min_field_id, "FieldsList", fields_to_min)


        #set mesh algorithm Mesh.Algorithm3D to HXT (option - 10) or Delaunay 3D (option - 1)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10) #HXT seems to be producing less slivers near curvature in the design
        gmsh.model.mesh.field.setAsBackgroundMesh(min_field_id)
        gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
        gmsh.model.mesh.generate(3)

       