import argparse
import random

#need to identify the number of unique xml needed
#based on (num and type of lanthionine) and number of mutations
#make another document containing the params for a Rosetta run (variable names)

def make_ripp_topo(ripp_len_dict, renamed_path, write_cst=False):
    lan_name_dict = {}
    lan_name_dict['DAL'] = 'lanthionine'
    lan_name_dict['ALA'] = 'lanthionine'
    lan_name_dict['DBS'] = 'methyllanthionine'
    lan_name_dict['DBR'] = 'methyllanthionine'

    variations_set = set()
    ripp_topology_dict = {}

    #things that can vary are the number, type, and length of the rings and the linkers
    #lets exlude tail linkers for now, but could presumably use randomize by rama* or shear

    cst_line = 'AtomPair SG CYS CA ALA HARMONIC 2.8 0.1'
    cst_chainbrk = 'AtomPair C CYS N ALA HARMONIC 1.33 0.1'

    #This section creates the ripp_topology dict
    for ripp in ripp_len_dict.keys():
        with open(renamed_path+ripp+'_conn.txt', 'r') as conn_file:
            cst_out = ''
            rama_list = [x for x in range (1, ripp_len_dict[ripp]+1)]
            ring_list = []
            conn_list = []
            ringlen_list = []
            linker_list = []
            for raw_line in conn_file:
                split = raw_line.rstrip().split(',')
                ring_st = int(split[0])
                ring_end = int(split[3])
                ring_len = abs(ring_st - ring_end) + 1
                conn_list.append([split[0], split[3], ring_len, lan_name_dict[split[1]]])
                ringlen_list.append(ring_len)
                cst_out += cst_line.replace('ALA', str(ring_st)).replace('CYS', str(ring_end)) + '\n'
                for x in range(min(ring_st, ring_end), max(ring_st, ring_end)+1):
                    ring_list.append(x)
            #which residues are not involved in a ring?
            rama_list = [x for x in rama_list if x not in ring_list]
            #determine the linker grops
            while len(rama_list) > 0:
                temp_list = []
                temp_list.append(rama_list[0])
                for x in range(len(rama_list)):
                    if temp_list[-1] + 1 in rama_list:
                        temp_list.append(temp_list[-1] + 1)
                    else:
                        rama_list = [x for x in rama_list if x not in temp_list]
                        if 1 not in temp_list and ripp_len_dict[ripp] not in temp_list: #don't include flexible termini
                            linker_list.append((len(temp_list), min(temp_list), max(temp_list)))
                        break
            if len(linker_list) >= 1:
                for i in range(len(linker_list)):
                    ringlen_list.append(-1)
            ringlen_list.sort()
            variations_set.add(tuple(ringlen_list))
            ripp_topology_dict[ripp] = [conn_list, linker_list]
        if write_cst:
            with open(ripp+'.cst', 'w') as cst_file:
                if ripp in ['6PQG', '7JVF', '7JU9']:
                    for i in range(1, ripp_len_dict[ripp]):
                        cst_out += cst_chainbrk.replace('ALA', str(i)).replace('CYS', str(i+1)) + '\n'
                cst_file.write(cst_out)
    return ripp_topology_dict

#6vlj and 6vhj appear to be the only ones requiring a linker and each only have two chains

#header will include the score functions to use
#variable section with the selectors and linkers
#simple metrics and task opations
#mover section, start by defining the crosslinker
#static functions are the fast relax movers and constraint movers to collect nmr agreement info
#genkic will rely on parsed protocols to connect randomize by rama will be used for the linkers between sections
#for rings of 10 AA or less, do full prediction
#for larger rings, break up the ring and predict 4 AA on each side of the linker
#the protocols secion will vary based on the number of crosslinker movers

#how to decide when to randomize a linker?
#the first time that a ring next to the linker is closed

header = '''<ROSETTASCRIPTS>
    <SCOREFXNS>
    	<ScoreFunction name="r15" weights="ref2015"/>
    	<ScoreFunction name="r15_cst" weights="ref2015_cst">
           <Reweight scoretype="chainbreak" weight="40.0" />
        </ScoreFunction>
    	<ScoreFunction name="r15_cart" weights="ref2015_cart_cst">
            <Reweight scoretype="angle_constraint" weight="3.0" />
        </ScoreFunction>
    	<ScoreFunction name="bb_only" weights="empty.wts" >
            <Reweight scoretype="fa_rep" weight="0.1" />
            <Reweight scoretype="fa_atr" weight="0.2" />
            <Reweight scoretype="hbond_sr_bb" weight="2.0" />
            <Reweight scoretype="hbond_lr_bb" weight="2.0" />
            <Reweight scoretype="rama_prepro" weight="0.45" />
            <Reweight scoretype="omega" weight="0.4" />
            <Reweight scoretype="p_aa_pp" weight="0.6" />
        </ScoreFunction>
    </SCOREFXNS>
    <TASKOPERATIONS>
        <IncludeCurrent name="incl_curr" />
    </TASKOPERATIONS>
    <RESIDUE_SELECTORS>
        <Index name="peptide" resnums="%%start%%-%%end%%"/>
'''

phase1_1 = '''    </RESIDUE_SELECTORS>
    <SIMPLE_METRICS>
        <RMSDMetric name="calc_rmsd_full" residue_selector="peptide" rmsd_type="rmsd_all_heavy" use_native="true" super="true"/>
        <RMSDMetric name="calc_rmsd_bb" residue_selector="peptide" rmsd_type="rmsd_protein_bb_heavy_including_O" use_native="true" super="true"/>
        <TotalEnergyMetric name="cst_eng" scorefxn="r15_cst" scoretype="atom_pair_constraint" />
        <TotalEnergyMetric name="eng_no_cst" scorefxn="r15" />
        <SasaMetric name="SASA" />
'''

phase1_2 = '''    </SIMPLE_METRICS>
'''

chainbreak = '''    <FILTERS>
        <ChainBreak name="chainbrk" />
    </FILTERS>
'''

phase1_3 = '''    <MOVERS>
        <FastRelax name="init_frlx" repeats="1" scorefxn="r15_cst" min_type="dfpmin"
            task_operations="incl_curr" />
        <RandomizeBBByRamaPrePro name="ramarand_pep" residue_selector="peptide"/>
        <AddChainBreak name="findbrk" find_automatically="true" distance_cutoff="1.5" />
            '''
phase1_4 = '''        <FastRelax name="final_frlx" repeats="3" scorefxn="r15_cst" min_type="dfpmin"
            task_operations="incl_curr" />
'''

cart_relax = '''        <FastRelax name="cart_frlx" repeats="1" scorefxn="r15_cart" min_type="dfpmin"
            task_operations="incl_curr" cartesian="true" >
                <MoveMap name="cart_mm">
                    <Span begin="1" end="999" chi="1" bb="1" bondlength="1" bondangle="1" />
                </MoveMap>
            </FastRelax>
'''

phase_old = '''    </RESIDUE_SELECTORS>
    <MOVERS>
    <FastRelax name="init_frlx" repeats="1" scorefxn="r15_cst" min_type="dfpmin"
            task_operations="incl_curr" FastRelax />
    <ConstraintSetMover name="add_nmr_cst" add_constraints="true" cst_file="%%cst_file%%"/>
    <ClearConstraintsMover name="rm_cst" />


'''
crosslink = '''        <CrosslinkerMover name="cyclize1" linker_name="lanthionine"
            residue_selector="lan1" scorefxn="r15_cst"
            filter_by_constraints_energy="false"
            filter_by_sidechain_distance="false"
            pack_and_minimize_linker_and_sidechains="true"
        />
'''

crosslink_nomin = '''        <CrosslinkerMover name="cyclize1" linker_name="lanthionine"
            residue_selector="lan1" scorefxn="r15_cst"
            filter_by_constraints_energy="false"
            filter_by_sidechain_distance="false"
            pack_and_minimize_linker_and_sidechains="false"
        />
'''

genkic_head_child = '''        <GeneralizedKIC name="genkic_1" closure_attempts="100"
            stop_when_n_solutions_found="1" 
            selector="lowest_energy_selector" selector_scorefunction="r15_cst">
'''

genkic_head_parent = '''        <GeneralizedKIC name="genkic_fin" closure_attempts="10000"
            stop_when_n_solutions_found="50" pre_selection_mover="genkic_steps"
            selector="lowest_energy_selector" selector_scorefunction="bb_only">
'''

genkic_pivots = '''            <SetPivots res1="%%piv0_1%%" res2="%%piv0_2%%" res3="%%piv0_3%%" atom1="CA" atom2="CA" atom3="CA" /> #first pivot must be second atom in chain, last is 2nd to last
'''

genkic_closebond = '''            <CloseBond res1="%%lan_stX%%" atom1="CB" res2="%%lan_endX%%" atom2="SG" bondlength="1.807" angle2="115.230" angle1="103.176" />
'''
#lan_end is cys resi

parsed_head = '''        <ParsedProtocol name="genkic_steps">
            <Add mover="genkic_1" />
'''

parsed_head_individual = '''        <ParsedProtocol name="genkic_steps">
'''

parsed_tail = '''            <Add mover="init_frlx" />
        </ParsedProtocol>
'''

dihedral_cst = '''        <AddConstraints name="add_omega_cst" >
            <DihedralConstraintGenerator name="add_from_current" dihedral="omega" residue_selector="peptide" dihedral_angle="180" />
        </AddConstraints>
'''

phase2 = '''        <ConstraintSetMover name="add_cst" add_constraints="true" cst_file="ripp.cst"/>
        <ConstraintSetMover name="add_nmr_cst" add_constraints="true" cst_file="nmr.cst"/>
        <ClearConstraintsMover name="rm_cst" />
        <MinMover name="cart_min" cartesian="true" bondlength="0" bondangle="1" bb="false" omega="false" chi="true" scorefxn="r15_cart" />
    </MOVERS>
    <PROTOCOLS>
'''

footer = '''        <Add metrics="calc_rmsd_full" labels="rmsd_full"/>
        <Add metrics="calc_rmsd_bb" labels="rmsd_bb"/>
        <Add metrics="SASA" labels="SASA"/>
        <Add metrics="run_metrics" labels="run_metrics"/>
        <Add metrics="eng_no_cst" labels="eng_no_cst" />
    </PROTOCOLS>
    <OUTPUT scorefxn="r15_cst" />
</ROSETTASCRIPTS>'''

atom_tree = '        <AtomTree name="new_tree" simple_ft="1" /> \n'

declare_bond = '        <DeclareBond name="connect_termini" atom1="C" res1="cres" atom2="N" res2="nres" add_termini="false" /> \n'

#for getting the different orders that rings can be closed
def get_permutations(topo):
    #print(topo)
    if len(topo) == 0:
        return [[]]

    perm = []
    for i in range(len(topo)):
        current = topo[i]
        remaining = topo[:i] + topo[i+1:]
        sub_perm = get_permutations(remaining)
        for sub in sub_perm:
            perm.append([current] + sub)
    return perm

def get_ring_indicies(link1, link2):
    start = int(link1)
    end = int(link2)
    if end < start:
        temp = start
        start = end
        end = temp
    return start, end

#classify the ring type
def det_ring_type(index, topology):
    #print(index, topology)
    start, end = get_ring_indicies(topology[index-1][0], topology[index-1][1])

    contained = False
    overlap = False
    start_link = False
    end_link = False

    for i in range(len(topology)):
        if i != index:
            check_start, check_end = get_ring_indicies(topology[i][0], topology[i][1])
            if check_start < start < check_end:
                start_link = True
            if check_start < end < check_end:
                end_link = True
            if check_start < start < end < check_end: #if within a second ring, add link cst and relax prior to genkic
                #if start < check_start < check_end < end: #contains a second ring, outer ring is different
                contained = True
            elif check_start < start < check_end or check_start < end < check_end:
                overlap = True

    if contained:
        return 'contained' #ring inside another ring
    elif start_link and end_link:
        return 'linker' #ring that overlaps with at least two rings
    elif overlap:
        return 'overlap' #ring that overlaps with another
    else:
        return 'free'

def get_genkic_steps(ripp, ring_num, ring_len, ring_type, findbrk, tails, individual, randomize, link_rlx=False):
    i = ring_num
    genkic_xml = ''

    if individual:
        #make genkic mover
        if not (ring_type == 'linker' and link_rlx): #default is to always do genkic prediction
            genkic_xml += genkic_head_child.replace('genkic_1', f'genkic_{i}')
            if ring_len <= 10:
                genkic_xml += f'            <AddResidue res_index="%%piv{i}_1%%" />\n'
                for j in range(1, ring_len - 2):
                    genkic_xml += f'            <AddResidue res_index="%%int{i}_{j}%%" />\n'
                genkic_xml += f'            <AddResidue res_index="%%piv{i}_3%%" />\n'

                # print(base,tails,ring_len,rings, i)
                if ring_type != 'overlap' and ring_type != 'linker':  # for overlapping rings, adding tails will make a ton of chainbreaks.
                    for j in range(tails):
                        genkic_xml += f'            <AddTailResidue res_index="%%tail{i}_{j + 1}%%" />\n'
            else:
                print(f'long ring {ripp}')
                genkic_xml += f'            <AddResidue res_index="%%piv{i}_1%%" />\n'
                for j in range(1, 7):
                    genkic_xml += f'            <AddResidue res_index="%%int{i}_{j}%%" />\n'
                genkic_xml += f'            <AddResidue res_index="%%piv{i}_3%%" />\n'
                for j in range(tails):
                    genkic_xml += f'            <AddTailResidue res_index="%%tail{i}_{j + 1}%%" />\n'
            genkic_xml += genkic_pivots.replace('0', str(i))
            genkic_xml += genkic_closebond.replace('X', str(i))  # lan_st number set to i
            genkic_xml += '        </GeneralizedKIC>\n'

        # make parsed protocol
        genkic_xml += parsed_head_individual.replace('genkic_steps', f'genkic_steps{i}')
        if not randomize:
            genkic_xml += f'            <Add mover="cyclize{i}" />\n'
        if not (ring_type == 'linker' and link_rlx):
            genkic_xml += f'            <Add mover="genkic_{i}" />\n' #otherwise relax ring to close
        if findbrk:
            genkic_xml += '            <Add mover="findbrk" />\n'
        genkic_xml += '            <Add mover="init_frlx" />\n'
        genkic_xml += '        </ParsedProtocol>\n'

    else: #not individual
        if i == 1:
            genkic_xml += genkic_head_child
        else:
            # make parsed protocol
            genkic_xml += parsed_head.replace('genkic_steps', f'genkic_steps{i}').replace('genkic_1', f'genkic_{i - 1}')
            if findbrk:
                genkic_xml += '            <Add mover="findbrk" />\n'
            genkic_xml += '            <Add mover="init_frlx" />\n'
            if not randomize:
                genkic_xml += f'            <Add mover="cyclize{i}" />\n'
            genkic_xml += '        </ParsedProtocol>\n'
            # genkic parent
            # for individual, use genkic_head_child
            genkic_xml += genkic_head_parent.replace('genkic_steps', f'genkic_steps{i}').replace('genkic_fin', f'genkic_{i}')
        if ring_len <= 10:
            genkic_xml += f'            <AddResidue res_index="%%piv{i}_1%%" />\n'
            for j in range(1, ring_len - 2):
                genkic_xml += f'            <AddResidue res_index="%%int{i}_{j}%%" />\n'
            genkic_xml += f'            <AddResidue res_index="%%piv{i}_3%%" />\n'

            if ring_type != 'overlap' and ring_type != 'linker':  # for overlapping rings, adding tails will make a ton of chainbreaks.
                for j in range(tails):
                    genkic_xml += f'            <AddTailResidue res_index="%%tail{i}_{j + 1}%%" />\n'
        else:
            print(f'long ring {ripp}')
            genkic_xml += f'            <AddResidue res_index="%%piv{i}_1%%" />\n'
            for j in range(1, 7):
                genkic_xml += f'            <AddResidue res_index="%%int{i}_{j}%%" />\n'
            genkic_xml += f'            <AddResidue res_index="%%piv{i}_3%%" />\n'
            for j in range(tails):
                genkic_xml += f'            <AddTailResidue res_index="%%tail{i}_{j + 1}%%" />\n'
        genkic_xml += genkic_pivots.replace('0', str(i))
        genkic_xml += genkic_closebond.replace('X', str(i))  # lan_st number set to i
        genkic_xml += '        </GeneralizedKIC>\n'

    return genkic_xml


#now, decide how the stucture will be predicted
def make_ripp_xml(ripp_topology_dict, ripp_len_dict, randomize=True, individual=True, link_rlx=False, rand_link=False):
    if link_rlx and not individual:
        print('link_rlx requires individual')
        exit(1)
    if randomize and rand_link:
        print('randomize and rand_link are both true, setting rand_link to false')
        rand_link=False
    for ripp in ripp_topology_dict:
        #print('ripp', ripp)
        metric_names = ''
        metric_nmr = ''
        ripp_set = set()
        add_linker = False
        overlapping_ring = False
        if rand_link and len(ripp_topology_dict[ripp][1]) == 1:
            add_linker = True
        if rand_link and (len(ripp_topology_dict[ripp][1]) > 1 or (len(ripp_topology_dict[ripp][1]) == 1 and len(ripp_topology_dict[ripp][0]) > 3)):
            print('not set up to handle multiple linkers or linkers with more than 2 chains')
            exit(1)
        else:
            perms = get_permutations(list(range(len(ripp_topology_dict[ripp][0]))))
            for perm in perms:
                if add_linker:
                    ripp_set.add(tuple(perm+[ripp_len_dict[ripp]]+[0])) #0 was for linker initially
                ripp_set.add(tuple(perm))
            for topo in ripp_set:
                    base = ''
                    base += f'_{ripp}'
                    for x in topo:
                        base += '_' + str(ripp_topology_dict[ripp][0][x][0]) #str(x)
                    rings = len(ripp_topology_dict[ripp][0])
                    temp_xml = header
                    #now add the residue selectors (one for each conjugation, one for linker if needed, one for full peptide)
                    for i in range(1, rings+1):
                        temp_xml += f'        <Index name="lan{i}" resnums="%%lan_st{i}%%,%%lan_end{i}%%"/>\n'
                        temp_xml += f'        <Index name="ring{i}" resnums="%%ring_st{i}%%-%%ring_end{i}%%"/>\n'
                    if add_linker:
                        temp_xml += '        <Index name="linker" resnums="%%link_st%%-%%link_end%%"/>\n'
                    #simple metrics section and relax mover
                    temp_xml += phase1_1
                    metric_names = ''
                    metric_nmr = ''
                    with open(renamed_path+ripp+'_conn.txt', 'r') as conn_file:
                        i=1
                        for line in conn_file:
                            split = line.split(',')
                            res1 = int(split[0])
                            res2 = int(split[3])
                            start = min(res1, res2)
                            end = max(res1, res2)
                            temp_xml += f'        <DihedralDistanceMetric name="ring{i}_dih" use_native="true" residue_selector="ring{i}"/>\n'
                            temp_xml += f'        <TotalEnergyMetric name="ring{i}_eng" scorefxn="r15" scoretype="total_score" residue_selector="ring{i}"/>\n'
                            temp_xml += f'        <PeptideInternalHbondsMetric name="ring{i}_hbond" residue_selector="ring{i}"/>\n'
                            temp_xml += f'        <TotalEnergyMetric name="ring{i}_nmr" scorefxn="r15_cst" scoretype="atom_pair_constraint" residue_selector="ring{i}"/>\n'
                            temp_xml += f'        <RMSDMetric name="ring{i}_rmsd_full" residue_selector="ring{i}" rmsd_type="rmsd_all_heavy" use_native="true" super="true"/>\n'
                            temp_xml += f'        <RMSDMetric name="ring{i}_rmsd_bb" residue_selector="ring{i}" rmsd_type="rmsd_protein_bb_heavy_including_O" use_native="true" super="true"/>\n'
                            metric_names += f'ring{i}_dih,ring{i}_eng,ring{i}_hbond,ring{i}_rmsd_full,ring{i}_rmsd_bb,'
                            metric_nmr += f'ring{i}_nmr,'
                            i += 1
                    temp_xml += phase1_2
                    temp_xml += chainbreak
                    temp_xml += phase1_3
                    temp_xml += phase1_4 #Consider r15_cart term here
                    temp_xml += cart_relax #Consider r15_cart term here
                    temp_xml += dihedral_cst
                    #movers (one for each conjugation, linker if needed (rama rand), genkic)
                    for i in range(1, rings+1):
                        temp_xml += crosslink.replace('lan1', f'lan{i}').replace('cyclize1', f'cyclize{i}')
                        temp_xml += crosslink_nomin.replace('lan1', f'lan{i}').replace('cyclize1', f'cyclize_nomin{i}')
                    #randomize whole peptide at the start instead
                    if add_linker:
                        temp_xml += '        <RandomizeBBByRamaPrePro name="rand_link" residue_selector="linker"/>\n'

                    #iterating through all of the rings and writing out genkic steps for them
                    for i in range(1,rings+1):
                        ring_num = i
                        ring_len = ripp_topology_dict[ripp][0][topo[i-1]][2] #topo[i-1]
                        #what kind of connection is it? adjacent/ std, overlapping, encompased
                        ring_type = det_ring_type(i, ripp_topology_dict[ripp][0])
                        findbrk = ripp in ['1AJ1', '6PQG', '7JVF', '7JU9'] #overlapping ring or linker
                        tails = ripp_len_dict[ripp] - ring_len
                        if ring_type == 'overlap' or ring_type == 'linker':
                            overlapping_ring = True
                        temp_xml += get_genkic_steps(ripp, ring_num, ring_len, ring_type, findbrk, tails,
                                                     individual=individual, randomize=randomize, link_rlx=link_rlx)

                    if overlapping_ring:
                        temp_xml += atom_tree
                        for i in range(2, ripp_len_dict[ripp]+1):
                            temp_xml += declare_bond.replace('connect_termini', f'termini{i}').replace('cres', str(i - 1)).replace('nres', str(i))
                        temp_xml += '        <ParsedProtocol name="fix_peptide_bond">\n'
                        for i in range(2, ripp_len_dict[ripp]+1):
                            temp_xml += f'            <Add mover="termini{i}" />\n'
                        temp_xml += '        </ParsedProtocol>\n'

                    #nmr_constraint
                    temp_xml += phase2.replace('ripp', ripp).replace('nmr', ripp + '_nmr').replace('nmr.cst', 'noe.cst')
                    for i in range(1, rings+1):
                        temp_xml += f'        <Add mover="cyclize{i}" />\n'
                    if not randomize:
                        temp_xml += '        <Add mover="rm_cst" />\n'

                    if randomize:
                        temp_xml += '        <Add mover="ramarand_pep" />\n'
                    if add_linker:
                        temp_xml += '        <Add mover="rand_link" />\n'
                    temp_xml += '        <Add mover="add_omega_cst" />\n'
                    if individual:
                        for i in range(1, rings+1):
                            temp_xml += f'        <Add mover="genkic_steps{i}" />\n'  # others in parsed protocol
                    else:
                        if not randomize:
                            temp_xml += f'        <Add mover="cyclize1" />\n' #others in parsed protocol
                        temp_xml += f'        <Add mover="genkic_{rings}" />\n'

                    temp_xml += '        <Add mover="init_frlx" />\n'
                    temp_xml += '        <Add mover="findbrk" />\n'
                    temp_xml += '        <Add mover="rm_cst" />\n'
                    for i in range(1, rings+1):
                        temp_xml += f'        <Add mover="cyclize{i}" />\n'
                    #output section
                    if overlapping_ring:
                        temp_xml += '        <Add mover="fix_peptide_bond" />\n'
                        temp_xml += '        <Add mover="init_frlx" />\n'
                        temp_xml += '        <Add mover="new_tree" />\n'
                    if overlapping_ring:
                        temp_xml += '        <Add mover="fix_peptide_bond" />\n'
                    if overlapping_ring:
                        temp_xml += '        <Add mover="fix_peptide_bond" />\n'
                    temp_xml += '        <Add mover="rm_cst" />\n'
                    for i in range(1, rings+1):
                        temp_xml += f'        <Add mover="cyclize{i}" />\n'
                    temp_xml += '        <Add mover="findbrk" />\n'
                    temp_xml += '        <Add mover="final_frlx" />\n'
                    temp_xml += '        <Add filter="chainbrk" />\n'
                    temp_xml += '        <Add mover="rm_cst" />\n'
                    temp_xml += f'        <Add mover="add_{ripp}_nmr_cst" />\n'
                    temp_xml += '        <Add metrics="cst_eng" labels="nmr_cst"/>\n'
                    temp_xml += '        <Add metrics="metric_nmr" labels="metric_nmr"/>\n'.replace('metric_nmr', metric_nmr[:-1])
                    temp_xml += '        <Add mover="rm_cst" />\n'
                    for i in range(1, rings+1):
                        temp_xml += f'        <Add mover="cyclize{i}" />\n'
                    temp_xml += '        <Add mover="cart_min" />\n'
                    temp_xml += footer.replace('run_metrics', metric_names[:-1])
                    #save to file
                    if randomize:
                        base += '_rand'
                    if individual:
                        base += '_ind'
                    if link_rlx:
                        base += '_lr'
                    if add_linker:
                        base += '_rl'
                    with open('ripp_pred_xml' + base +'.xml', 'w') as out_file:
                        out_file.write(temp_xml)
                    #need to write out script variables for each ripp into a csv
                    #for each unique loop permutation will need a .sh file
                    #for each ringsize in the permutation
                    #for each ring that matches the size select all possible variable closure cominations
                    #and exclude the ring from being chosen again
                    #write out all combinations to a file

def make_ripp_options(ripp_topology_dict, renamed_path, randomize=True, individual=True, link_rlx=False, rand_link=False):
    options_head = '''-parser:protocol ${XML}
-in:file:s $PDB
-in:file:native $PDB
-nstruct 100
-out:file:scorefile 6VLJ_scan.sc
-use_input_sc true
-in:file:fullatom
-ignore_zero_occupancy false
-linmem_ig 10
-out:pdb_gz
-out:prefix $PREFIX
-in:detect_disulf
-mute all
-in:file:extra_res_path /home/tydingcw/Documents/EGFR_antibodies/RiPP_design/NCAA/uff_params/final_params/
-parser:script_vars'''
    if link_rlx and not individual:
        print('link_rlx requires individual')
        exit(1)
    if randomize and rand_link:
        print('randomize and rand_link are both true, setting rand_link to false')
        rand_link=False
    #generate the instructions for folding
    for ripp in ripp_topology_dict:
        add_linker = False
        instruction_list = []
        if rand_link and len(ripp_topology_dict[ripp][1]) == 1:
            add_linker = True
        if len(ripp_topology_dict[ripp][1]) > 1 or (len(ripp_topology_dict[ripp][1]) == 1 and len(ripp_topology_dict[ripp][0]) > 3):
            print('not set up to handle multiple linkers or linkers with more than 2 chains')
        else:
            perms = get_permutations([str(x[2])+'_'+x[0]+'_'+x[1]+'_'+x[3] for x in ripp_topology_dict[ripp][0]])
            perm_count = 1
            for perm in perms:
                # make the topology
                #now add the residue selectors (one for each conjugation, one for linker if needed, one for full peptide)
                rings = len(perm)
                options_str = options_head
                options_str += ' start=1 '
                options_str += f'end={ripp_len_dict[ripp]} '
                with open(renamed_path+ripp+'_conn.txt', 'r') as conn_file:
                    i=1
                    for line in conn_file:
                        split = line.split(',')
                        res1 = int(split[0])
                        res2 = int(split[3])
                        start = min(res1, res2)
                        end = max(res1, res2)
                        options_str += f'ring_st{i}='+str(start)+ ' ' #ripp_topology_dict[ripp][0][i-1][0]
                        options_str += f'ring_end{i}='+str(end)+ ' ' #ripp_topology_dict[ripp][0][i-1][1]
                        i += 1
                for i in range(1, len(perm)+1): # these are permutations of the rings
                    ring = perm[i-1]
                    res1 = int(ring.split('_')[1])
                    res2 = int(ring.split('_')[2])
                    options_str += f'lan_st{i}='+str(res1)+ ' ' #ripp_topology_dict[ripp][0][i-1][0]
                    options_str += f'lan_end{i}='+str(res2)+ ' ' #ripp_topology_dict[ripp][0][i-1][1]
                if add_linker:
                    options_str += 'link_st='+str(ripp_topology_dict[ripp][1][0][1]) + ' '
                    options_str += 'link_end='+str(ripp_topology_dict[ripp][1][0][2]) + ' '
                options_perm = options_str
                #for i in range(1,rings+1):
                #for each ring, generate possible ring closures, add to a list
                #sample all combinations of anchors, random middle pivot
                #unless ring is >10, then sample all middle pivots
                #for each ring, add params to list
                #for the next ring, take those params and copy them over the number of anchor times
                #then add the new params to that list
                #write those out to options files
                options_list = []
                for i in range(1, len(perm)+1):
                    ring = perm[i-1]
                    ring_len = ring_len = int(ring.split('_')[0])
                    if options_list != []:
                        #don't let anchor be lanthionine resi
                        anchors = ring_len-2
                        if anchors > 8: #if more than 10 residues, don't close the whole ring, only 6 resi can be middle pivot
                            anchors = 6
                        options_len = len(options_list)
                        temp_list = []
                        for k in range(options_len): #the options_list parameters are copied over and duplicated like 1,1,1,2,2,2
                            for j in range(anchors):
                                temp_list.append(options_list[k])
                        options_list = temp_list
                    max_ring = max(int(ring.split('_')[1]), int(ring.split('_')[2]))
                    min_ring = min(int(ring.split('_')[1]), int(ring.split('_')[2]))
                    full_ring = list(range(min_ring, max_ring+1))
                    if ring_len <= 10:
                        anchor_list = list(range(min_ring+1, max_ring)) #need to choose an anchor
                        #anchor
                        if options_list == []: #if list is empty, intialize to have an item for each anchor
                            for k in range(len(anchor_list)):
                                options_list.append('')
                        for k in range(len(anchor_list)):
                            to_add_list = list(range(k, len(options_list), len(anchor_list)))
                            anchor = anchor_list[k]
                            random_index = random.randrange(-len(anchor_list), -1)
                            start = -(ring_len - full_ring.index(anchor)) +1 #random_index + 1
                            for l in to_add_list: options_list[l] += f'piv{i}_1='+str(full_ring[start]) +' '
                            for j in range(1, ring_len-2):
                                for l in to_add_list: options_list[l] += f'int{i}_{j}='+str(full_ring[start+j])+' '
                            for l in to_add_list: options_list[l] += f'piv{i}_3='+str(full_ring[start+ring_len-2])+' '

                            tails = ripp_len_dict[ripp] - ring_len
                            tail_list = [x for x in range(1, ripp_len_dict[ripp]+1) if x not in full_ring]
                            for j in range(tails):
                                for l in to_add_list: options_list[l] += f'tail{i}_{j+1}='+str(tail_list[j])+' '

                            random_index = random.randint(start +1, start + ring_len - 3)
                            for l in to_add_list: options_list[l] += f'piv{i}_2='+str(full_ring[random_index])+' '
                    else:
                        ring_loop = []
                        for j in range(3):
                            ring_loop.append(max_ring - j)
                            ring_loop.append(min_ring + j)
                        for k in range(6):
                            if options_list == []:
                                for n in range(6): #for the six middle pivots to choose from
                                    options_list.append('')
                            to_add_list = list(range(k, len(options_list), 6))
                            for l in to_add_list: options_list[l] += f'piv{i}_1='+str(max_ring-3)+' '
                            for l in to_add_list: options_list[l] += f'piv{i}_3='+str(min_ring+3)+' '
                            for l in to_add_list: options_list[l] += f'int{i}_1='+str(max_ring-2)+' '
                            for l in to_add_list: options_list[l] += f'int{i}_2='+str(max_ring-1)+' '
                            for l in to_add_list: options_list[l] += f'int{i}_3='+str(max_ring)+' '
                            for l in to_add_list: options_list[l] += f'int{i}_4='+str(min_ring)+' '
                            for l in to_add_list: options_list[l] += f'int{i}_5='+str(min_ring+1)+' '
                            for l in to_add_list: options_list[l] += f'int{i}_6='+str(min_ring+2)+' '
                            for l in to_add_list: options_list[l] += f'piv{i}_2='+str(ring_loop[k])+' '

                            tails = ripp_len_dict[ripp] - ring_len
                            tail_list = [x for x in range(1, ripp_len_dict[ripp]+1) if x not in full_ring]
                            for j in range(tails):
                                for l in to_add_list: options_list[l] += f'tail{i}_{j+1}='+str(tail_list[j])+' '
                #add option to choose the xml file
                suffix = ''
                if randomize:
                    suffix += '_rand'
                if individual:
                    suffix += '_ind'
                if link_rlx:
                    suffix += '_lr'
                if add_linker:
                    suffix += '_rl'
                for item in options_list:
                    with open(f'options/{ripp}_{perm_count}{suffix}.options', 'w') as options_out:
                        base = ''
                        base += f'_{ripp}'
                        for ring in perm:
                            base += '_'+ring.split('_')[1] #This is the [0] in ripp topo dict #ring.split('_')[0]
                        base += suffix
                        xml = 'ripp_pred_xml' + base +'.xml'
                        ripp_upper = ripp.upper()
                        options_fin = options_str
                        options_fin = options_fin.replace('${XML}', xml).replace('6VLJ', ripp_upper).replace('$PREFIX', f'{ripp_upper}/{perm_count}_init_')
                        options_fin = options_fin.replace('$PDB', f'{renamed_path}{ripp_upper}_1_rename.pdb')
                        options_fin += item
                        options_out.write(options_fin)
                        perm_count += 1

relax_header = '''<ROSETTASCRIPTS>
    <SCOREFXNS>
    	<ScoreFunction name="r15" weights="ref2015"/>
    	<ScoreFunction name="r15_cst" weights="ref2015_cst"/>
    	<ScoreFunction name="r15_cart" weights="ref2015_cart_cst"/>
    </SCOREFXNS>
    <RESIDUE_SELECTORS>
        Index name="peptide" resnums="%%start%%-%%end%%"/>
        <Chain name="peptide" chains="A"/>
'''

relax_phase1_1 = '''    </RESIDUE_SELECTORS>
    <TASKOPERATIONS>
        <IncludeCurrent name="incl_curr" />
    </TASKOPERATIONS>
    <SIMPLE_METRICS>
        <RMSDMetric name="calc_rmsd_full" residue_selector="peptide" rmsd_type="rmsd_all_heavy" use_native="true" super="true"/>
        <RMSDMetric name="calc_rmsd_bb" residue_selector="peptide" rmsd_type="rmsd_protein_bb_heavy_including_O" use_native="true" super="true"/>
        <TotalEnergyMetric name="cst_eng" scorefxn="r15_cst" scoretype="atom_pair_constraint" />
        <TotalEnergyMetric name="eng_no_cst" scorefxn="r15" />
        <SasaMetric name="SASA" />
'''

relax_metrics = '''        <DihedralDistanceMetric name="ring1_dih" use_native="true" residue_selector="ring1"/>
        <TotalEnergyMetric name="ring1_eng" scorefxn="r15_cst" scoretype="total_score" residue_selector="ring1"/>
        <PeptideInternalHbondsMetric name="ring1_hbond" residue_selector="ring1"/>
        <TotalEnergyMetric name="ring1_nmr" scorefxn="r15_cst" scoretype="atom_pair_constraint" residue_selector="ring1"/>
        <RMSDMetric name="ring1_rmsd_full" residue_selector="ring1" rmsd_type="rmsd_all_heavy" use_native="true" super="true"/>
        <RMSDMetric name="ring1_rmsd_bb" residue_selector="ring1" rmsd_type="rmsd_protein_bb_heavy_including_O" use_native="true" super="true"/>
'''

relax_phase1_2 = '''    </SIMPLE_METRICS>
    <MOVERS>
'''

bb_cst = '''        <AddConstraints name="add_coord_cst">
                <CoordinateConstraintGenerator name="coord_cst" residue_selector="peptide" ca_only="true" sd="0.5"/> #0.5 is default
        </AddConstraints>
'''

full_cst = '''        <AddConstraints name="add_coord_cst">
                <CoordinateConstraintGenerator name="coord_cst" residue_selector="peptide" ca_only="false" sidechain="true" sd="0.5"/> #0.5 is default
        </AddConstraints>
'''

relax_phase2 = '''        <FastRelax name="init_frlx" repeats="5" scorefxn="r15_cst" min_type="dfpmin" task_operations="incl_curr" />
        RemoveConstraints name="rm_coord_cst" constraint_generators="coord_cst"/>
        <ConstraintSetMover name="add_nmr_cst" add_constraints="true" cst_file="8CWX_noe.cst"/>
        <ClearConstraintsMover name="rm_cst"/>
        <MinMover name="cart_min" cartesian="true" bondlength="0" bondangle="1" bb="false" omega="false" chi="true" scorefxn="r15_cart" />

    </MOVERS>
    <APPLY_TO_POSE>
    </APPLY_TO_POSE>
    <PROTOCOLS>
'''

relax_footer = '''        <Add metrics="calc_rmsd_full" labels="rmsd_full"/>
        <Add metrics="calc_rmsd_bb" labels="rmsd_bb"/>
        <Add metrics="SASA" labels="SASA"/>
        Add metrics="run_metrics" labels="run_metrics"/>
        <Add metrics="eng_no_cst" labels="eng_no_cst" />
    </PROTOCOLS>
    <OUTPUT scorefxn="r15_cst" />
</ROSETTASCRIPTS>'''

relax_options = '''-parser:protocol xml_name
-in:file:native file_loc/8CWX_1_rename.pdb
-nstruct 20
-out:file:scorefile 8CWX_free_relax.sc
-use_input_sc true
-in:file:fullatom
-ignore_zero_occupancy false
-linmem_ig 10
-out:pdb_gz
-out:prefix 8CWX_native/free_
-in:detect_disulf false
-in:file:extra_res_path /home/tydingcw/Documents/EGFR_antibodies/RiPP_design/NCAA/uff_params/final_params/
'''

def make_relax_scripts(ripp_topology_dict, mode, renamed_path):
    if mode not in ['bb', 'full', 'free']:
        print('Invalid mode, must be one of: bb, full, free')
        exit(1)

    for ripp in ripp_topology_dict.keys():
        num_rings = len(ripp_topology_dict[ripp][0])
        temp_xml = relax_header


        for i in range(num_rings):
            pos1 = min(int(ripp_topology_dict[ripp][0][i][0]), int(ripp_topology_dict[ripp][0][i][1]))
            pos2 = max(int(ripp_topology_dict[ripp][0][i][0]), int(ripp_topology_dict[ripp][0][i][1]))
            temp_xml += f'        <Index name="lan{i+1}" resnums="{pos1},{pos2}"/>\n'
            temp_xml += f'        <Index name="ring{i+1}" resnums="{pos1}-{pos2}"/>\n'

        temp_xml += relax_phase1_1

        for i in range(1, num_rings+1):
           temp_xml += relax_metrics.replace('ring1', f'ring{i}')

        temp_xml += relax_phase1_2 #simple metrics close

        for i in range(1, num_rings + 1):
            temp_xml += crosslink.replace('lan1', f'lan{i}').replace('cyclize1', f'cyclize{i}')

        #chose constraints based on mode
        if mode == 'bb':
            temp_xml += bb_cst
        elif mode =='full':
            temp_xml += full_cst

        temp_xml += relax_phase2.replace('8CWX', ripp)

        for i in range(1, num_rings+1):
            temp_xml += f'        <Add mover="cyclize{i}" />\n'

        if mode != 'free':
            temp_xml += f'        <Add mover="add_coord_cst" />\n'

        temp_xml += '        <Add mover="init_frlx" />\n'
        temp_xml += '        <Add mover="rm_cst" />\n'
        temp_xml += '        <Add mover="add_nmr_cst" />\n'
        temp_xml += '        <Add metrics="cst_eng" labels="nmr_cst"/>\n'

        #add ring metrics for nmr
        ring_nmr_string = ''
        for i in range(1, num_rings + 1):
            ring_nmr_string += f'ring{i}_nmr,'
        ring_nmr_string = ring_nmr_string[:-1]
        temp_xml += f'        <Add metrics="{ring_nmr_string}" labels="{ring_nmr_string}"/>\n'

        temp_xml += '        <Add mover="rm_cst" />\n'
        for i in range(1, num_rings+1):
            temp_xml += f'        <Add mover="cyclize{i}" />\n'
        temp_xml += '        <Add mover="cart_min" />\n'

        #add ring metrics
        ring_eng_string = ''
        for i in range(1, num_rings + 1):
            ring_eng_string += f'ring{i}_dih,ring{i}_eng,ring{i}_hbond,ring{i}_rmsd_full,ring{i}_rmsd_bb,'
        ring_eng_string = ring_eng_string[:-1]
        temp_xml += f'        <Add metrics="{ring_eng_string}" labels="{ring_eng_string}"/>\n'

        temp_xml += relax_footer

        with open(f'relax_{ripp}_{mode}.xml', 'w') as xml_out:
            xml_out.write(temp_xml)

        with (open(f'{ripp}_{mode}.options', 'w') as opt_out):
            opt_str = relax_options.replace('8CWX', ripp).replace('file_loc', renamed_path)
            opt_str = opt_str.replace('xml_name', f'relax_{ripp}_{mode}.xml').replace('free', mode)
            opt_out.write(opt_str)


if __name__ == "__main__":

    #parser = argparse.ArgumentParser(description = 'take in alignment and connect files and prepare for quick rosetta scan')
    #parser.add_argument('-i', '--input', help = 'input alignments', required = True)
    #parser.add_argument('-c', '--chain', help = 'chain of interest', default='A')
    #parser.add_argument('-p', '--prefix', help = 'prefix for output', default='')
    #argument=parser.parse_args()

    ripp_len_dict = {}
    ripp_len_dict['1AJ1'] = 19 #- Jung
    ripp_len_dict['6VHJ'] = 16
    ripp_len_dict['6PQG'] = 19 #file caused error - chainbreaker
    ripp_len_dict['6VLJ'] = 19
    ripp_len_dict['7JU9'] = 20 # has DBU
    ripp_len_dict['7JVF'] = 21 #has DBU

    renamed_path = 'renamed_pdb/'

    topo_dict = make_ripp_topo(ripp_len_dict, renamed_path, True)
    #make_ripp_xml(topo_dict, ripp_len_dict, False, False) #initial options
    rand = True #radomize input rama
    ind = True #predict the rings sequentially
    make_ripp_xml(topo_dict, ripp_len_dict, rand, ind)
    ##make_ripp_xml(topo_dict, ripp_len_dict, link_rlx=True) # not prefered
    make_ripp_options(topo_dict, renamed_path, rand, ind)

    # for the relax scripts, uncomment
    #for mode in ['bb', 'full', 'free']:
    #    make_relax_scripts(topo_dict, mode, renamed_path)
