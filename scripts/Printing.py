from Bond_Angles import *


#xyz file 
def print_xyz_file(xyz_file):
    xyz = [] 
    for x in xyz_file: 
        xyz.append(x.split()) # this is so that the coordinates are split as well
    #Heading 
    print("XYZ coordinates of the molecule")
    print()
    for x in xyz: 
        print(" " .join(x))
    print() 

    #Making all the xyz coordinates be floats instead of strings
    for x in range(2, len(xyz)):
        xyz[x][1] = float(xyz[x][1])
        xyz[x][2] = float(xyz[x][2])
        xyz[x][3] = float(xyz[x][3])

    #Since the file has how many molecule there is and energy of it, i've decide to return it, for other purposes
    return xyz; 


#Bond length
def print_bond_length(bonds): 
    #Don't want any repeats
    non_repeats = []
    for x in range(len(bonds)):
        for y in range(x+1, len(bonds)): 
            if(x != y and bonds[x][0] == bonds[y][1] and bonds[x][1] == bonds[y][0]):
                non_repeats.append(bonds[x])
                break; 
    print("Bond length between two atoms in Angstrom")
    for x in non_repeats:
        print(x[0], "-",x[1], "(", x[2][0], "-", x[3][0], ")", round(x[4],3))
    print()

#Angle between 3 atoms 
def print_angles(angle):
    print("Bond Angles")
    for i in range (len(angle)):
        print(''.join(str(i) for i in angle[i][0]), ":", angle[i][1])
    print()

#Dihedral Angle
def print_dihedral(dihedral):
    if(len(dihedral) == 0):
        print("There is no dihedral angles")
        print("")
        return
    no_repeats = [] 
    for i in range(len(dihedral)):
        add = False
        for j in range(len(no_repeats)):
            if (no_repeats[j][0] != dihedral[i][3] and no_repeats[j][1] != dihedral[i][2]): 
                add = True
        if(add or len(no_repeats) == 0):
            no_repeats.append(dihedral[i])
        
    print("Dihedral Angle")
    for i in range(len(no_repeats)): 
        display = "-".join(str(no_repeats[i][x]) for x in range(4))
        print(display, ":", no_repeats[i][4])
    print()

#Out of Plane Angle
def print_outOfPlane(outOfPlane):
    if(len(outOfPlane) == 0):
        print("There is no out of plane angles")
        print("")
        return
    print("Out of Plane Angle")
    for i in range(len(outOfPlane)):
        display = "-".join(str(outOfPlane[i][x]) for x in range(4))
        print(display, ":", round(outOfPlane[i][4],3))
    print()

#Center of Mass 
def print_com(com):
    print("Center of Mass")
    print("x:",round(com[0],3), " y:", round(com[1],3), " z:", round(com[2],3))
    print()

def print_recenter(xyz):
    print("Recentering the xyz file based off of Center of Mass")
    print()
    for x in xyz: 
        print(x[0] + ": " +  " , ".join(str(round(x[y],3)) for y in range(1,4)))
    print() 

#Moment of Inertia 
def print_moi(moi):
    print("Moment of Inertia Tensor")
    print("\t x \t y \t z ")
    print("x: \t" + str(round(moi[0][0],3)) + "\t" +  str(round(moi[0][1],3)) + "\t" +  str(round(moi[0][2],3)))
    print("y: \t" + str(round(moi[1][0],3)) + "\t" +  str(round(moi[1][1],3)) + "\t" +  str(round(moi[1][2],3)))
    print("z: \t" + str(round(moi[2][0],3)) + "\t" +  str(round(moi[2][1],3)) + "\t" +  str(round(moi[2][2],3)))
    print()