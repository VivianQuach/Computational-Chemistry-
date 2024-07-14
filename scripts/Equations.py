import math 
#Pythagoream Theorem 
# c^2 = a^2 + b^2 (in a xy plane)
#This function will be in a xyz plane 
#This is used to help find the bond length 
#The parameter represent the atom 
def pythagoream(a1,a2):
    return math.sqrt((a1[0] - a2[0])**2 + (a1[1] - a2[1])**2 + (a1[2] - a2[2])**2)

#Vector 
# V(a,b) = (x(b) - x(a))x + (y(b) - y(a))t + (z(b) - z(a))z
#The outer x,y,z have " ^ " on top 
#This will give you the vector between two atoms 
def vector(a1,a2):
    x = a2[1] - a1[1]
    y = a2[2] - a1[2]
    z = a2[3] - a1[3]
    return [x,y,z]

#Magnitude of Vector 
# ||V(a,b)|| = sqrt((x(b) - x(a))^2 + (y(b) - y(a))^2 + (z(b) - z(a))^2)
#This is the distance between 2 3d cartesian coordinates 
#(aka distance between two vectors)
def mag_vector(a1,a2):
    vec = vector(a1,a2)
    return math.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)

#Unit Vector 
# uv(a,b) = (V(a,b))/ ||V(a,b)|| 
#This is between two atoms 
def unit_vector(a1,a2):
    vec = vector(a1,a2)
    mag_vec = mag_vector(a1,a2)
    unit_vec = []
    for i in range(3):
        unit_vec.append(vec[i]/mag_vec)
    return unit_vec

#Dot Product 
# a.b = a[x]*b[x] + a[y]*b[y] + a[z]*b[z]
# a and b are both vector 
#This will be between 3 atoms 
def dot_product(a,b):
    x = a[0]*b[0]
    y = a[1]*b[1]
    z = a[2]*b[2]
    return (x + y + z)

#Cross Product 
# r(ji) x r(jk) = (||r(ji)||)*(||r(jk||)sin(θijk) = sin(θijk)
#This cross product is from two unit vectors 
# a x b = (a[1]*b[2] - a[2]*b[1]) - (a[0]*b[2] - a[2]*b[0]) + (a[0]*b[1] - a[1]*b[0])
def cross_product(a,b):
    i = a[1]*b[2] - a[2]*b[1]
    j = a[0]*b[2] - a[2]*b[0]
    k = a[0]*b[1] - a[1]*b[0]
    return [i,j,k]

