from pymol import cmd

def customize_view():
    cmd.bg_color('white')
    cmd.set('ray_trace_mode', 1)
    cmd.set('ray_opaque_background', 1)
    cmd.set('ray_shadows', 'off')
    cmd.set('surface_mode', 1)
    cmd.set('antialias', 4)
    cmd.set('spec_reflect', 0)
    cmd.set('cartoon_oval_length', 0.5)

def split_and_copy_first_state():
    # Get the names of all objects in the session
    objects = cmd.get_object_list()

    if objects:
        first_object = objects[0]
        # Split states of the first object
        cmd.split_states(first_object)

        # Copy the first state to a new object named 'obj01'
        cmd.create('obj01', f'{first_object}_0001')

        # Parse the connection file
        conn_file_path = f'/home/tydingcw/Documents/EGFR_antibodies/RiPP_design/renamed_pdb/{first_object}_conn.txt'
        try:
            with open(conn_file_path, 'r') as conn_file:
                for line in conn_file:
                    parts = line.strip().split(',')
                    if len(parts) >= 4:
                        residue1 = parts[0]
                        residue2 = parts[3]

                        # Show sidechains for these residues
                        #cmd.show('sticks', f'(resi {residue1} or resi {residue2}&(sc.|(n. CA|n. N&r. PRO)))')
                        #cmd.show('sticks', f'(byres(resi {residue1} or resi {residue2})&(sc.|(n. CA|n. N&r. PRO)))')
                        cmd.show('sticks', f'(resi {residue1}&(sc.|(n. CA|n. N&r. PRO)))')
                        cmd.show('sticks', f'(resi {residue2}&(sc.|(n. CA|n. N&r. PRO)))')
            cmd.hide("(all and hydro)")
        except FileNotFoundError:
            print(f"Connection file not found: {conn_file_path}")
        
        # Define a list of colors to cycle through
        colors = [
            'blue', 'tv_blue', 'marine', 'slate', 'lightblue', 'skyblue', 
            'purpleblue', 'deepblue', 'density', 'cyan', 'palecyan', 
            'aquamarine', 'greencyan', 'teal', 'lightteal', 'deepteal'
        ]

        # Iterate over the objects starting from the second one and apply colors
        for i, obj in enumerate(objects[1:]):
            color = colors[i % len(colors)]  # Cycle through colors
            cmd.color(color, obj)
        
        util.cnc("all",_self=cmd) #color by element

customize_view()
split_and_copy_first_state()

# To ensure the script is executed when loaded
if __name__ == "__main__":
    customize_view()
    split_and_copy_first_state()
