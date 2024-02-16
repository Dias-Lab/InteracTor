#!/usr/bin/env python

"""
Routines to calculate the Accessible Surface Area of a set of atoms.
The algorithm is adapted from the Rose lab's chasa.py, which uses
the dot density technique found in:

Shrake, A., and J. A. Rupley. "Environment and Exposure to Solvent
of Protein Atoms. Lysozyme and Insulin." JMB (1973) 79:351-371.
"""


import math
from vector3d import pos_distance, Vector3d, pos_distance_sq


def generate_sphere_points(n):
    """
    Returns list of 3d coordinates of points on a sphere using the
    Golden Section Spiral algorithm.
    """
    points = []
    inc = math.pi * (3 - math.sqrt(5))
    offset = 2 / float(n)
    for k in range(int(n)):
        y = k * offset - 1 + (offset / 2)
        r = math.sqrt(1 - y*y)
        phi = k * inc
        points.append([math.cos(phi)*r, y, math.sin(phi)*r])
    return points


def find_neighbor_indices(atoms, probe, k):
    """
    Returns list of indices of atoms within probe distance to atom k. 
    """
    neighbor_indices = []
    atom_k = atoms[k]
    radius = atom_k.radius + probe + probe
    indices = range(k)
    indices = range(k+1, len(atoms))
    for i in indices:
        atom_i = atoms[i]
        dist = pos_distance(atom_k.pos, atom_i.pos)
        if dist < radius + atom_i.radius:
            neighbor_indices.append(i)
    return neighbor_indices


def calculate_asa(atoms, probe, n_sphere_point=960):
    """
    Returns list of accessible surface areas of the atoms, using the probe
    and atom radius to define the surface.
    """
    sphere_points = generate_sphere_points(n_sphere_point)

    const = 4.0 * math.pi / len(sphere_points)
    test_point = Vector3d()
    areas = []
    for i, atom_i in enumerate(atoms):

        neighbor_indices = find_neighbor_indices(atoms, probe, i)
        n_neighbor = len(neighbor_indices)
        j_closest_neighbor = 0
        radius = probe + atom_i.radius

        n_accessible_point = 0
        for point in sphere_points:
            is_accessible = True

            test_point.x = point[0]*radius + atom_i.pos.x
            test_point.y = point[1]*radius + atom_i.pos.y
            test_point.z = point[2]*radius + atom_i.pos.z

            cycled_indices = range(j_closest_neighbor, n_neighbor)
            cycled_indices = (range(j_closest_neighbor))

            for j in cycled_indices:
                atom_j = atoms[neighbor_indices[j]]
                r = atom_j.radius + probe
                diff_sq = pos_distance_sq(atom_j.pos, test_point)
                if diff_sq < r*r:
                    j_closest_neighbor = j
                    is_accessible = False
                    break
            if is_accessible:
                n_accessible_point += 1

        area = const*n_accessible_point*radius*radius 
        areas.append(area)
    return areas


def main():
  import sys
  import getopt
  import molecule
  

  usage = \
  """

  Copyright (c) 2007 Bosco Ho
  
  Calculates the total Accessible Surface Area (ASA) of atoms in a 
  PDB file. 

  Usage: asa.py -s n_sphere in_pdb [out_pdb]
  
  - out_pdb    PDB file in which the atomic ASA values are written 
               to the b-factor column.
               
  -s n_sphere  number of points used in generating the spherical
               dot-density for the calculation (default=960). The 
               more points, the more accurate (but slower) the 
               calculation. 

  """

  opts, args = getopt.getopt(sys.argv[1:], "n:")
  if len(args) < 1:
    print (usage)
    return
    
  mol = molecule.Molecule(args[0])
  atoms = mol.atoms()
  molecule.add_radii(atoms)

  n_sphere = 960
  for o, a in opts:
    if '-n' in o:
      n_sphere = int(a)
      print ("Points on sphere: ", n_sphere)
  asas = calculate_asa(atoms, 1.4, n_sphere)
  print ("%.1f angstrom squared." % sum(asas))

  if len(args) > 1:
    for asa, atom in zip(asas, atoms):
      atom.bfactor = asa
    mol.write_pdb(args[1])
  
  
if __name__ == "__main__":
  main()


    
