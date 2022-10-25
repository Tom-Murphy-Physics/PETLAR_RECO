"""
Written by Tom Murphy
An stl file is a type of file used in CAD
It makes triangles in 3D space to create a "Mesh"
Meshes can be viewed like any other object in CAD

Each triangle needs:
    3 Verticies (in 3D cartesian coordinates)
    Vector normal to the surface
    
This file takes the verticies of a cube and prdocues an stl file 
Note that it will only take 7 arguments 
    The file it will write to
    6 verticies 
    
The vertices will be on oppisate corners of the square
"""

import numpy as np
import struct

#   Writes the first line of the stl file 
#   Decalares the object as a solid
def intro(F):
    F.write("solid ASCII\n")
    
#   Takes a cartesian direction and returns the stl code for that direction
def get_normal(n):
    if n == "x":
        return "  facet normal 1.000000e+00 0.000000e+00 0.000000e+00\n"
    if n == "y":
        return "  facet normal 0.000000e+00 1.000000e+00 0.000000e+00\n"
    if n == "z":
        return "  facet normal 0.000000e+00 0.000000e+00 1.000000e+00\n"
    if n == "-x":
        return "  facet normal -1.000000e+00 0.000000e+00 0.000000e+00\n"
    if n == "-y":
        return "  facet normal 0.000000e+00 -1.000000e+00 0.000000e+00\n"
    if n == "-z":
        return "  facet normal 0.000000e+00 0.000000e+00 -1.000000e+00\n"

#   Takes the File and the vertex and writes it to the stl file    
def make_vertex(F,a,b,c):
    string = "      vertex   " + str(a) + " " + str(b) + " " + str(c) + "\n"
    F.write(string)    # Vertex 1
    
#   Makes the cube in stl format    
def make_cube(F,a,b,c,d,e,f):
    a = np.single(a/20)
    b = np.single(b/20)
    c = np.single(c/20)
    d = np.single(d/20)
    e = np.single(e/20)
    f = np.single(f/20)
    
    #   Front Face
    
    F.write(get_normal("-y"))       # Triangle 1
    F.write("    outer loop\n")  
    make_vertex(F,a,b,c)
    make_vertex(F,d,b,c)
    make_vertex(F,d,b,f)
    F.write("    endloop\n")
    F.write("  endfacet\n")
    
    F.write(get_normal("-y"))        # Triangle 2
    F.write("    outer loop\n")  
    make_vertex(F,a,b,c)
    make_vertex(F,a,b,f)
    make_vertex(F,d,b,f)
    F.write("    endloop\n")
    F.write("  endfacet\n")
    
    #   Right Face
    
    F.write(get_normal("x"))        # Triangle 3
    F.write("    outer loop\n")  
    make_vertex(F,d,b,c)
    make_vertex(F,d,e,c)
    make_vertex(F,d,e,f)
    F.write("    endloop\n")
    F.write("  endfacet\n")
    
    F.write(get_normal("x"))        # Triangle 4
    F.write("    outer loop\n")  
    make_vertex(F,d,b,c)
    make_vertex(F,d,b,f)
    make_vertex(F,d,e,f)
    F.write("    endloop\n")
    F.write("  endfacet\n")
    
    #   Back Face
    
    F.write(get_normal("y"))        # Triangle 5
    F.write("    outer loop\n")  
    make_vertex(F,d,e,c)
    make_vertex(F,d,e,f)
    make_vertex(F,a,e,f)
    F.write("    endloop\n")
    F.write("  endfacet\n")
    
    F.write(get_normal("y"))        # Triangle 6
    F.write("    outer loop\n")  
    make_vertex(F,d,e,c)
    make_vertex(F,a,e,c)
    make_vertex(F,a,e,f)
    F.write("    endloop\n")
    F.write("  endfacet\n")
    
    #   Left Face 
    
    F.write(get_normal("-x"))        # Triangle 7
    F.write("    outer loop\n")  
    make_vertex(F,a,e,c)
    make_vertex(F,a,b,c)
    make_vertex(F,a,b,f)
    F.write("    endloop\n")
    F.write("  endfacet\n")
    
    F.write(get_normal("-x"))        # Triangle 8
    F.write("    outer loop\n")  
    make_vertex(F,a,e,c)
    make_vertex(F,a,e,f)
    make_vertex(F,a,b,f)
    F.write("    endloop\n")
    F.write("  endfacet\n")
    
    #   Top Face
    
    F.write(get_normal("z"))        # Triangle 9
    F.write("    outer loop\n")  
    make_vertex(F,a,b,f)
    make_vertex(F,d,b,f)
    make_vertex(F,d,e,f)
    F.write("    endloop\n")
    F.write("  endfacet\n")
    
    F.write(get_normal("z"))        # Triangle 10
    F.write("    outer loop\n")  
    make_vertex(F,a,b,f)
    make_vertex(F,a,e,f)
    make_vertex(F,d,e,f)
    F.write("    endloop\n")
    F.write("  endfacet\n")
    
    #   Bottom Face
    
    F.write(get_normal("-z"))        # Triangle 11
    F.write("    outer loop\n")  
    make_vertex(F,a,b,c)
    make_vertex(F,d,b,c)
    make_vertex(F,d,e,c)
    F.write("    endloop\n")
    F.write("  endfacet\n")
    
    F.write(get_normal("-z"))        # Triangle 12
    F.write("    outer loop\n")  
    make_vertex(F,a,b,c)
    make_vertex(F,a,e,c)
    make_vertex(F,d,e,c)
    F.write("    endloop\n")
    F.write("  endfacet\n")

#   Ends the file
def end(F):
    F.write("endsolid")