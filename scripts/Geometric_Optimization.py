# Lennard Jones Potential Equation 
#Appromixates the interaction between two molecules
# Eplison: How strongly does the two atoms attract each other 
# Sigma: How close two nonbonded particals can get (Van der Waals Radius)
# r: the distance between the two molecules center point 
def lennard_jones_potential(r, epsilon = 1.0, sigma = 1.0):
    return 4 * epsilon( (sigma/r)**12 - (sigma/r)**6)

