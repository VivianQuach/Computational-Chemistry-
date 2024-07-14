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