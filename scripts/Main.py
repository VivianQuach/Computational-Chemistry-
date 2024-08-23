
from Bond_Angles import *
from Printing import * 


file = "../Molecules/xyz/" + input() + ".xyz"
xyz_file = is_xyz_file(file)
#Return with the xyz coordinates as float instead of string 
xyz = print_xyz_file(xyz_file)

bonds, bond_len = bond_connection(xyz) 
print_bond_length(bond_len)

angle = get_angle(xyz,bonds)
print_angles(angle)

dihedral = get_dihedral(xyz,bonds)
print_dihedral(dihedral)


outOfPlane = get_out_of_plane(xyz,bonds)
print_outOfPlane(outOfPlane)

com = get_com(xyz)
print_com(com)

coords = recenter(xyz,com)
print_recenter(coords)

moi = get_moi(xyz)
print_moi(moi)