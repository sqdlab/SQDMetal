import gmsh

class GMSH_Mesh_Builder:

    def __init__(self, surfaces, user_options):
        
        self.surfaces = surfaces
        self.user_options = user_options


    def build_mesh(self):

        #turn off these parametes as we will determine the mesh element size from a mesh field
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

        #create list of surfaces to mesh finely
        surfaces_to_mesh = [x[2] for x in self.surfaces]

        print('Finely meshing surfaces:', surfaces_to_mesh)

        #create distance field and threshold field - see Gmsh documentation - python tutorial 10 for more detail
        gmsh.model.mesh.field.add("Distance", 1)
        gmsh.model.mesh.field.setNumbers(1, "SurfacesList", surfaces_to_mesh)
        gmsh.model.mesh.field.setNumber(1, "Sampling", self.user_options['mesh_sampling'])

        #A threshold field built in conjunction with the distance field
        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", self.user_options['mesh_min'])        #minimum size of mesh element
        gmsh.model.mesh.field.setNumber(2, "SizeMax", self.user_options['mesh_max'])        #maximum size of mesh element
        gmsh.model.mesh.field.setNumber(2, "DistMin", 30e-3)        #distance from the gmsh surface to keep minimum mesh element size 
        gmsh.model.mesh.field.setNumber(2, "DistMax", 200e-3)       #how far from the element until the max element size can be implemented 
        
        #set mesh algorithm Mesh.Algorithm3D to HXT (option - 10) or Delaunay 3D (option - 1)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10) #HXT seems to be producing less slivers near curvature in the design
        gmsh.model.mesh.field.setAsBackgroundMesh(2)
        gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
        gmsh.model.mesh.generate(3)
        
        print('Mesh successfully built in Gmsh.')

       