from pymol import cmd, util

all_objects = cmd.get_names('objects')

four_letter = [name for name in all_objects if len(name) == 4]

blues = ['blue', 'purpleblue', 'deepblue', 'tv_blue', 'lightblue', 'skyblue', 'bluewhite', 'slate', 'marine', 'density', 'cyan', 'palecyan', 'greencyan', 'teal', 'deepteal', 'lightteal', 'aquamarine']
num_blues = len(blues)

for obj in four_letter:
    cmd.split_states(obj)
    cmd.create(None,f"{obj}_0001")
    base = obj.upper()
    count = 0
    for name in all_objects:
        if f'_{base}_' in name:
            cmd.color_deep(blues[count%num_blues], name, 0)
            count += 1

residue_names = 'CYS+DAL+DBB+DBR+DBS'
cmd.select('my_selection', f'resn {residue_names}')
cmd.show("sticks","((byres (my_selection))&(sc.|(n. CA|n. N&r. PRO)))")
cmd.hide("(my_selection and hydro)")
util.cnc("all",_self=cmd)

#for obj in all_objects:
#    if '6VHJ' in obj:
#        cmd.align("{} and name CA".format(obj), "obj01 and name CA")
#        cmd.align("{} and name CA".format(obj), "1aj1 and name CA")
#        cmd.align(f"{obj} and name CA", f"1aj1 and name CA")

#all_objects = cmd.get_names('objects')
#for obj in all_objects:
#    if '1AJ1' in obj:
#        cmd.align(f"{obj} and name CA", f"1aj1 and name CA")
