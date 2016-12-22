"""
Created on 12/01/16 by Mihir Metkar
This script is adapted from https://sourceforge.net/p/pymol/mailman/message/25744364/
This script doesnt work with,
resn (str): the residue name
resi (str): the residue identifier (residue number) as a string, including optional insertion code
resv (int): the residue identifier (residue number) as an integer, excluding insertion code
It works with,
ID (int): PDB atom id (not guaranteed to be unique)
To test if the distances are correct, use the following 2 commands in Pymol.

dist /4ug0_5S///2/P, /4ug0_5S///120/P
dist id 80206, id 82724

Where, 80206 is id for resi 2, name P of 4ug0_5S and 82724 is id for resi 121, name P of 4ug0_5S

To get distance between 2 different chains/ molecules, make a 2nd list b instead of a2.

"""

from pymol import stored

stored.atom1 = []
stored.atom2 = []
cmd.iterate_state(1, "4ug0_5.8S and name P", "stored.atom1.append((resv, ID))")
cmd.iterate_state(1, "4ug0_28S and name P", "stored.atom2.append((resv, ID))")


outFile = open("4ug0_5-8S_28S_distances.txt", 'w')

outFile.write("resi1(4ug0_5.8S)\tresi2(4ug0_28S)\tpdb_id_1\tpdb_id_2\tdistance\n")

for i,a1 in enumerate(stored.atom1):
    a1_resv = a1[0]
    a1_ID = a1[1]

    for j,a2 in enumerate(stored.atom2):
        a2_resv = a2[0]
        a2_ID = a2[1]

        outFile.write( "%d\t%d\t%s\t%s\t%s\n" % (a1_resv, a2_resv, a1_ID, a2_ID, cmd.get_distance( "id %s" % a1_ID, "id %s" % a2_ID)))
print (a1_resv, a2_resv)
print ("Done")

outFile.close()
