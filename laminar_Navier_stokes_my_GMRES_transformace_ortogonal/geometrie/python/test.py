

#The import of Kernel module is necessary in order to be able tu use
#the basic geometric objects

from CGAL.Kernel import *

# we want to use a Delaunay_triangulation_2

from CGAL.Triangulations_2 import *
from random import *

dt = Delaunay_triangulation_2()   # we construct the triangulation here
vertices = {}
for i in range(20):
    vertices[i] = dt.insert(Point_2(random(),random()))  # we generate a random points and insert them in the dt triangulation

cir_faces = dt.incident_faces(vertices[0]) # use of the incident_faces function
                                           # cir_faces is circulator over the incident faces of the Vertex vertices[0]
# compute the finite faces incident to vertices[0]
finites_faces = []
f1 = cir_faces.next()
if dt.is_infinite(f1) == False:
    finites_faces.append(f1)
    
for f in cir_faces:
    if f == f1:
        break   
    finites_faces.append(f) 

    

# the finite edges of the triangulation
for e in dt.edges:
    print e.vertex()
    #draw_edge(e)

# the finite vertices
for v in dt.vertices:
    print v.point().x(),v.point().y()

#the finite points:
for p in dt.points:
    print p.x(),p.y()

     

# use of line_walk function
# first, we have to define a line that will cross the triangulation
p1 = Point_2(0,0)
p2 = Point_2(1,1)

faces = dt.line_walk(p1,p2)   # this function return a list of faces
for f in faces:
    if dt.is_infinite(f) == False:
        #draw_face(f)
        print "this face is finite"