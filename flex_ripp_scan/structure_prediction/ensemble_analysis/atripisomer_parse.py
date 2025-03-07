import argparse
import gzip
import numpy as np

def get_dihedral(atom1, atom2, atom3, atom4):
    #https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
    """Praxeolitic formula
    1 sqrt, 1 cross product"""
    p0 = np.array(atom1)
    p1 = np.array(atom2)
    p2 = np.array(atom3)
    p3 = np.array(atom4)

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    angle = np.degrees(np.arctan2(y, x))
    return angle

def compute_plane(point1, point2, point3):
    """
    Computes the formula for a plane given three non-collinear points.

    Parameters:
    point1, point2, point3 (tuple): Three tuples each representing the (x, y, z) coordinates of a point.

    Returns:
    tuple: The coefficients (A, B, C, D) of the plane equation Ax + By + Cz + D = 0.
    """
    # Create vectors from point1 to point2 and from point1 to point3
    vector1 = np.array(point2) - np.array(point1)
    vector2 = np.array(point3) - np.array(point1)
    
    # Compute the cross product of vector1 and vector2 to get the normal vector
    normal_vector = np.cross(vector1, vector2)
    
    # The normal vector components are the coefficients A, B, and C
    A, B, C = normal_vector
    
    # Compute D using the normal vector and point1
    D = -np.dot(normal_vector, np.array(point1))
    
    return (A, B, C, D)

def check_point_side(plane_params, test_point):
    """
    Determines which side of a plane a point is on.
    
    Parameters:
    plane_normal (tuple): The normal vector of the plane (A, B, C).
    plane_point (tuple): A known point on the plane (x0, y0, z0).
    test_point (tuple): The point to test (x, y, z).
    
    Returns:
    int: -1 if the point is on the side opposite to the normal,
         1 if on the same side as the normal,
         0 if the point is on the plane.
    """
    A, B, C, D = plane_params
    #x0, y0, z0 = plane_point
    x, y, z = test_point
    
    # Calculate D using the known point on the plane
    #D = - (A * x0 + B * y0 + C * z0)
    
    # Substitute the test point into the plane equation
    result = A * x + B * y + C * z + D
    
    #if result > 0:
    #    return 1
    #elif result < 0:
    #    return -1
    #else:
    #    return 0
    return result/abs(result)

def extract_xyz_coordinates(pdb_filename):
    """
    Extracts the XYZ coordinates from a PDB file.

    Parameters:
    pdb_filename (str): The path to the PDB file.

    Returns:
    dict: A dictionary with atom serial numbers as keys and (x, y, z) tuples as values.
    """
    coordinates = {}#TODO check for gzip
    if '.gz' in pdb_filename:
        with gzip.open(pdb_filename, 'rb') as file:
            for line in file:
                line = line.decode('utf-8')
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    # Extract fields from the line
                    resi = line[22:26].strip()
                    atom = line[12:16].strip()
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    # Store the coordinates in the dictionary
                    coordinates[f'{resi}_{atom}'] = (x, y, z)
    else:
        with open(pdb_filename, 'r') as file:
            for line in file:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    # Extract fields from the line
                    resi = line[22:26].strip()
                    atom = line[12:16].strip()
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    # Store the coordinates in the dictionary
                    coordinates[f'{resi}_{atom}'] = (x, y, z)
    return coordinates

def filter_obs(pdb_filenames, dihedrals):
    observations = {}
        
    for pdb_filename in pdb_filenames:
        coords = extract_xyz_coordinates(pdb_filename)
        dih_angs = []
        for dih in dihedrals:
            #parse the dihedral values
            dih_atoms = dih.split('-')
            dih_xyz = []
            for name in dih_atoms:
                dih_xyz.append(coords[name])
            plane_params = compute_plane(dih_xyz[0], dih_xyz[1], dih_xyz[2])
            dih_angs.append(check_point_side(plane_params, dih_xyz[3]))
        # Write the PDB name and the signs of the dihedrals
        #row = [str(angle) for angle in dih_angs] + [pdb_filename]
        
        #report stats
        obs = tuple(dih_angs)
        if obs in observations.keys():
            observations[obs] = observations[obs] + 1
        else:
            observations[obs] = 1

    #print(observations)
    max_value = max(observations.values())
    if list(observations.values()).count(max_value) > 1:
        print(f'error max value {max_value} is present multiple times')
        exit(1)
    else:
        max_key = max(observations, key=observations.get)
        return max_key

def main():
    parser = argparse.ArgumentParser(description='Calculate dihedral angles for PDB files and output the signs to a CSV file.')
    parser.add_argument('-p', '--pdb_list', type=str, help='File containing a list of PDB filenames')
    parser.add_argument('-d', '--dihedrals', nargs='+', help='List of dihedral angles to calculate') #Using plane instead of dihedral now
    parser.add_argument('-o', '--output', type=str, default='', help='Output CSV file name')
    parser.add_argument('-r', '--report', action='store_true', default=False, help='If specified, report occupency stats')
    parser.add_argument('-f', '--filt', default='', help='filter based on most common type in the ensemble')

    args = parser.parse_args()

    if args.output == '':
        output_csv = args.pdb_list.replace('.txt', '') + '.csv'
    else:
        output_csv = args.output
    
    with open(args.pdb_list, 'r') as file:
        pdb_filenames = file.read().splitlines()

    if 'score' in pdb_filenames[0]:
        pdb_filenames = pdb_filenames[1:]

    top_group = None
    if args.filt != '':
        filt_filenames = []
        with open(args.filt, 'r') as file:
            filt_filenames = file.read().splitlines()
 
        if 'score' in filt_filenames[0]:
            filt_filenames = filt_filenames[1:]
        top_group = filter_obs(filt_filenames, args.dihedrals)
    #print(top_group)

    #print(args.dihedrals)
    observations = {}

    with open(output_csv, 'w') as csvfile:
        # Write the header row
        dih_string = ''
        for dih in args.dihedrals:
            dih_string += dih + ','
        csvfile.write(dih_string  + 'description' + '\n')
        
        for pdb_filename in pdb_filenames:
            coords = extract_xyz_coordinates(pdb_filename)
            #dihedral_signs = get_dihedral(pdb_filename, dihedrals)
            dih_angs = []
            for dih in args.dihedrals:
                #parse the dihedral values
                dih_atoms = dih.split('-')
                dih_xyz = []
                for name in dih_atoms:
                    #resi = name.split('_')[0]
                    #atom = name.split('_')[1]
                    dih_xyz.append(coords[name])
                #temp_dih = get_dihedral(dih_xyz[0], dih_xyz[1], dih_xyz[2], dih_xyz[3])
                plane_params = compute_plane(dih_xyz[0], dih_xyz[1], dih_xyz[2])
                #print(dih, temp_dih, dih_xyz)
                #print(dih_xyz, plane_params)
                #dih_angs.append(temp_dih/abs(temp_dih))
                dih_angs.append(check_point_side(plane_params, dih_xyz[3]))
            #dihedral_signs = get_dihedral(pdb_filename, dihedrals)
            # Write the PDB name and the signs of the dihedrals
            obs = tuple(dih_angs)
            #print(obs)
            description = pdb_filename.split('.pdb')[0]
            if '/' in description:
                description = description.split('/')[-1]
            row = [str(angle) for angle in dih_angs] + [description]
            #print(','.join(row))
            if args.filt != '':
                if obs == top_group:
                    csvfile.write(','.join(row) + '\n')
            else:
                csvfile.write(','.join(row) + '\n')
            
            #report stats
            if obs in observations.keys():
                observations[obs] = observations[obs] + 1
            else:
                observations[obs] = 1

    if args.report:
        print(observations)

if __name__ == '__main__':
    main()

