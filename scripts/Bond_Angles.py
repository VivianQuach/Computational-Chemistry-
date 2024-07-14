import os,sys,math,numpy
from Equations import * 

#All the molecules covalence radius 
cov_rads = {  'H' : 0.37, 'C' : 0.77, 'O' : 0.73, 'N' : 0.75, 'F' : 0.71,
'P' : 1.10, 'S' : 1.03, 'Cl': 0.99, 'Br': 1.14, 'I' : 1.33, 'He': 0.30,
'Ne': 0.84, 'Ar': 1.00, 'Li': 1.02, 'Be': 0.27, 'B' : 0.88, 'Na': 1.02,
'Mg': 0.72, 'Al': 1.30, 'Si': 1.18, 'K' : 1.38, 'Ca': 1.00, 'Sc': 0.75,
'Ti': 0.86, 'V' : 0.79, 'Cr': 0.73, 'Mn': 0.67, 'Fe': 0.61, 'Co': 0.64,
'Ni': 0.55, 'Cu': 0.46, 'Zn': 0.60, 'Ga': 1.22, 'Ge': 1.22, 'As': 1.22,
'Se': 1.17, 'Kr': 1.03, 'X' : 0.00}

#What we are looking for: What atom is bonded to which atom. 
#1st: Grab the xyz file 
#2nd: Caculate the bond length between atoms to see if they are less than the threshold of the covalent radii added together
#3rd: See what the bond angle are between the molecule Note: If their is less than 3, it should not have an angle

#See if the file exists or not 
def is_xyz_file(xyz_file):
    try: 
        file = open(xyz_file, "r")
    except IOError: 
        print("File does not exist")
        sys.exit() # terminates the program as the file doesn't exist RIP 
    lines = file.readlines() 
    file.close()
    return lines; 

#####
#Bond Length 
######

#Checking if they are covalently bonded together 
def is_covalently_bond(a,b):
    #from 1 to len(vector) -> Since at position 0, its the symbol and that can't be caculated
    rij_eq = pythagoream(a[1::],b[1:len(b)])
    k = 1.2 # Threshold factor 
    #looking for the radius based off their symbol 
    a_rad = cov_rads[a[0]]
    b_rad = cov_rads[b[0]]
    if(rij_eq <= k*(a_rad + b_rad)):
       return True
    else: 
        return False

#Bond Connection 
#Which atom is connected together
def bond_connection(xyz): 
    #Note that the xyz, the first two lines are invalid here 
    xyz = xyz[2::]
    #Stores all the bonds that are connected together 
    bond_len = []
    bonds = []
    for i in range(len(xyz)):
        is_connect = []
        for j in range(len(xyz)):
            check_bond = is_covalently_bond(xyz[i], xyz[j])
            if (check_bond and i != j):
                rij_eq = pythagoream(xyz[i][1::], xyz[j][1::])
                bond_len.append([i, j, xyz[i], xyz[j], rij_eq])
                is_connect.append(j)
        bonds.append([xyz[i][0], i, is_connect])
    return (bonds, bond_len); 

######
#Bond Angles 
######

#Gets the angle between the three atoms
def angle(c1,c2,c3):
    #gets the unit vector inverse cosine(r(ji) . r(jk))
    uv1 = unit_vector(c2,c1)
    uv2 = unit_vector(c2,c3)
    dot = dot_product(uv1,uv2)
    return (180.0/math.pi) * math.acos(dot)

#Gets the angle of all the connects of the molecules 
def get_angle(xyz, bonds):
    xyz = xyz[2::]
    angles = [] 
    for j in range (len(xyz)):
        for a in range (len(bonds)): 
            #this loops make sure that the atom is connected to the j on the left side (i)
            if (a < len(bonds[j][2])): 
                i = bonds[j][2][a]
                for b in range (a+1, len(bonds)):
                    #this loops make sure that the atom is connected to the j on the right side (k)
                    if(b < len(bonds[j][2])):
                        k = bonds[j][2][b]
                        connection = i , '-' , (j) , '-', k
                        angles.append([connection,angle(xyz[i],xyz[j],xyz[k])])
    return angles; 

######
#Dihedral/Torison Angles 
######

#Calculating the angle between 4 atoms 
def dihedral(c1,c2,c3,c4):
    #1st. Get unit vector between c2 and c1, c2 and c3, c3 and c2 and c3 and c4
    ji = unit_vector(c2,c1)
    jk = unit_vector(c2,c3)
    kj = unit_vector(c3,c2)
    kl = unit_vector(c3,c4)
    #2nd. Get the cross product between ji and jk, kj and kl 
    ji_jk = cross_product(ji,jk)
    kj_kl = cross_product(kj,kl)
    #3rd. Get the dot product between the two cross product 
    dot = dot_product(ji_jk, kj_kl)
    #4th. Get the sign of the torison angle. Equation (normal vector of ijk) dot product to the vector of kl
    sign = 2 * float(dot_product(ji_jk, kl) < 0) - 1 
    #Using math.acos(dot) = cos^(-1)(cos(theta)) = theta 
    return (180.0 / math.pi) * sign * math.acos(dot) 

def get_dihedral(xyz,bonds):
    xyz = xyz[2::]
    angles = []
    for i in range(len(xyz)):
        if (len(bonds[i][2]) <= 1):
            continue
        i_bonds = bonds[i][2]
        for j in i_bonds:
            if (len(bonds[j][2]) <= 1):
                    continue
            j_bonds = bonds[j][2]
            for k in i_bonds:
                #we can do this way because the loop is going through the same thing 
                if (j == k): 
                    continue
                for x in j_bonds: 
                    if (i == x or x == k):
                        continue
                    angles.append([bonds[k][1], bonds[i][1],bonds[j][1],bonds[x][1],dihedral(xyz[k], xyz[i], xyz[j], xyz[x])])
    return angles
