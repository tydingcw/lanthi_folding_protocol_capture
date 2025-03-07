from pymol import cmd

all_objects = cmd.get_names('objects')
for obj in all_objects:
    if '6VHJ' in obj:
        cmd.align("{} and name CA and resi 3-7".format(obj), "obj01 and name CA and resi 3-7")
#        cmd.align("{} and name CA".format(obj), "1aj1 and name CA")
#        cmd.align(f"{obj} and name CA", f"1aj1 and name CA")

#all_objects = cmd.get_names('objects')
#for obj in all_objects:
#    if '1AJ1' in obj:
#        cmd.align(f"{obj} and name CA", f"1aj1 and name CA")
