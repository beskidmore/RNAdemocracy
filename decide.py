### Used to make decisions between RNA secondary structure prediction 
# programs                                                             ###
### Author: Benjamin Skidmore                                          ###
### Date: 10/9/2017                                                    ###

char = 'A', 'U', 'C', 'G', 'T', 'a', 'u', 'c', 'g', 't', 'N', 'n'
char2 = '.', '(', ')'
 

### SECTION 1 Structure identification                                 ###
### This section of the code is used to identify and count how many 
# types of each structure is found through ###
### the various first level predictions.                               ###


###This section is used to parse through and identify full stems in 
# individual predictions ###
###Four ways to identify: FFTL, SLF, IZ, FSLF ###
###Structures are written to all_stems.txt and stem_file_x ###
##Used later to replace Closeness dictonary keys with exactness keys##

def catch(st, ch):
    return [i for i, ltr in enumerate(st) if ltr == ch]
	
i = 0
j = 0
forward =0
backward =0
counter_close = 0
check_open = ''
loop_cutoff = 0
counter_open =0
counter_open2 =0
fequalblst = set()
stem_title_only = False
izlst = set()

with open('form_matrix') as windows:
    for l in windows:
        if '1>' in l and stem_title_only == False:
            seq_name = l
            i = i+1
            secondI= i
            stem_file = open('all_stems.txt', 'w')
            stem_title_only = True
            #print('NEXT SEQUENCE:')
            continue
        elif l.startswith('1U') or l.startswith('1G') or 
        l.startswith('1C') or l.startswith('1A'):
            SEQ = l[1:]
            #print(len(l))
            print(len(SEQ))
            continue
        else:
            j = j+1
            stem_file2 = open('stem_file_%i.txt' % (j), 'w')
            count_all = 0
            counter_open = 0
            backward = 0
            forward = 0
            counter_close = 0
            NEWCOUNT = 0
            NEWCOUNT_b = 0
            check_close = ''
            yesmam = l[1:]
            #print('\n'+ "NEXT PREDICTION:"+ '\n' + yesmam)
            stem_file.write(('\n'+ "NEXT PREDICTION:"+ '\n' + yesmam))
            stem_file2.write(('\n'+ "NEXT PREDICTION:"+ '\n' + yesmam))
            
            NEWCOUNT_sq = 0
            NEWCOUNT_b_sq = 0
            counter_open_sq = 0
            counter_close_sq = 0
            check_close_sq = ''

            NEWCOUNT_pt = 0
            NEWCOUNT_b_pt = 0
            counter_open_pt = 0
            counter_close_pt = 0
            check_close_pt = ''
            
            NEWCOUNT_ar = 0
            NEWCOUNT_b_ar = 0
            counter_open_ar = 0
            counter_close_ar = 0
            check_close_ar = ''
            fmlist = list(l.strip('1'))
            for xyz in yesmam:
                count_all +=1
                ##set condition 0##
                if xyz =='(':
                    check_open = True
                    counter_close = 0
                    forward +=1
                    front_holder = list(catch(yesmam, '(')) 
                    counter_open +=1
                    NewCount_i = list(catch(yesmam, '('))
                    open2_check = False
                    #print(NEWCOUNT, NewCount_i)
                    if check_close == False:
                        counter_open2 +=1
                        open2_check = True
                        #print(counter_open2)
                    if NEWCOUNT == NEWCOUNT_b:
                        ## Register the newcount i ##
                        NEWc_i = NEWCOUNT
                    NEWCOUNT +=1
                    
                elif xyz =='[':
                    check_open_sq = True
                    counter_close = 0
                    forward +=1
                    front_holder_sq = list(catch(yesmam, '[')) 
                    counter_open +=1
                    NewCount_i_sq = list(catch(yesmam, '['))
                    open2_check_sq = False
                    #print(NEWCOUNT, NewCount_i)
                    if check_close_sq == False:
                        counter_open2 +=1
                        open2_check_sq = True
                        #print(counter_open2)
                    if NEWCOUNT_sq == NEWCOUNT_b_sq:
                        ## Register the newcount i_sq ##
                        NEWc_i_sq = NEWCOUNT_sq
                    NEWCOUNT_sq +=1
                    
                elif xyz =='{':
                    check_open_pt = True
                    counter_close = 0
                    forward +=1
                    front_holder_pt = list(catch(yesmam, '{')) 
                    counter_open +=1
                    NewCount_i_pt = list(catch(yesmam, '{'))
                    open2_check_pt = False
                    #print(NEWCOUNT, NewCount_i)
                    if check_close_pt == False:
                        counter_open2 +=1
                        open2_check_pt = True
                        #print(counter_open2)
                    if NEWCOUNT_pt == NEWCOUNT_b_pt:
                        ## Register the newcount i_pt ##
                        NEWc_i_pt = NEWCOUNT_pt
                    NEWCOUNT_pt +=1
   
                elif xyz =='<':
                    check_open_ar = True
                    counter_close = 0
                    forward +=1
                    front_holder_ar = list(catch(yesmam, '<')) 
                    counter_open +=1
                    NewCount_i_ar = list(catch(yesmam, '<'))
                    open2_check_ar = False
                    #print(NEWCOUNT, NewCount_i)
                    if check_close_ar == False:
                        counter_open2 +=1
                        open2_check_ar = True
                    if NEWCOUNT_ar == NEWCOUNT_b_ar:
                        ## Register the NEWc_i_ar ##
                        NEWc_i_ar = NEWCOUNT_ar
                    NEWCOUNT_ar +=1
                    
                elif xyz == ')':
                    check_open = False
                    check_close = True
                    counter_close +=1
                    backward +=1
                    NEWCOUNT_b +=1
                    back_holder = list(catch(yesmam, ')'))
                    NewCount_j = list(catch(yesmam, ')'))
                    
                    ## From first to last ##
                    if counter_open == backward:
                        front_2 = min(front_holder[:counter_open])
                        back_2 = max(back_holder[:backward])
                        #print("f = b", str(front_2) + '-' + str(back_2))
                        fequalblst.add(str(front_2) + '-' + str(back_2))
                        #print(yesmam[front_2:(back_2+1)])
                        stem_file.write('\n'+ str(front_2) + '-' + str(back_2) +'\n')
                        stem_file.write(yesmam[front_2:(back_2+1)] + '\n')
                        stem_file2.write('\n'+ str(front_2) + '-' + str(back_2) +'\n')
                        stem_file2.write(yesmam[front_2:(back_2+1)] + '\n')
                        
                    ## Stem loop finder ##
                    if counter_open2 == counter_close and 
                    open2_check == True:
                        front_3 = max(front_holder[:(forward - counter_open2)])
                        back_3 = max(back_holder[:(backward+1)])
                        op_cl_count = 0
                        if front_3 >= back_3:
                            continue
                        else:
                            for op_clos in yesmam[front_3:back_3+1]:
                                if op_clos == '(':
                                    op_cl_count += 1
                                elif op_clos == ')':
                                    op_cl_count -= 1
                        if op_cl_count == 0:
                            #print("SLF", str(front_3) + '-' + str(back_3))
                            #print(yesmam[front_3:back_3+1] + '\n')
                            stem_file.write(str(front_3) + '-' + str(back_3) + '\n')
                            stem_file.write((yesmam[front_3:back_3+1] + '\n'))
                            counter_open2 = 0
                            stem_file2.write(str(front_3) + '-' + str(back_3) + '\n')
                            stem_file2.write((yesmam[front_3:back_3+1] + '\n'))
                            
                    ## Intermediate zeros (zeroed out in the middle of 
                    # sequence), ONLY OK FOR PSUDOKNOTS ##
                    if NEWCOUNT == NEWCOUNT_b:
                        #print(NewCount_i)
                        front_4 = NewCount_i[NEWc_i]
                        back_4 = NewCount_j[NEWCOUNT_b-1]
                        ## Do we reach 0? If so use that number where 
                        # first i != j if not keep going until we 
                        # find a j to make that i 0 
                        #print("NEWCOUNT", str(front_4) + '-' + str(back_4))
                        #print(yesmam[front_4:back_4+1])
                        stem_file.write(str(front_4) + '-' + str(back_4)+ '\n')
                        stem_file.write(yesmam[front_4:back_4+1]+ '\n')
                        stem_file2.write(str(front_4) + '-' + str(back_4)+ '\n')
                        stem_file2.write(yesmam[front_4:back_4+1]+ '\n')
                        izlst.add(str(front_4) + '-' + str(back_4))
                elif xyz == ']':
                    check_open_sq = False
                    check_close_sq = True
                    counter_close +=1
                    backward +=1
                    NEWCOUNT_b_sq +=1
                    back_holder_sq = list(catch(yesmam, ']'))
                    NewCount_j_sq = list(catch(yesmam, ']'))
                    
                    ## From first to last ##
                    if counter_open == backward:
                        front_2 = min(front_holder[:counter_open])
                        back_2 = max(back_holder[:backward])
                        #print("f = b", str(front_2) + '-' + str(back_2))
                        fequalblst.add(str(front_2) + '-' + str(back_2))
                        #print(yesmam[front_2:(back_2+1)])
                        stem_file.write('\n'+ str(front_2) + '-' + str(back_2) +'\n')
                        stem_file.write(yesmam[front_2:(back_2+1)] + '\n')
                        stem_file2.write('\n'+ str(front_2) + '-' + str(back_2) +'\n')
                        stem_file2.write(yesmam[front_2:(back_2+1)] + '\n')
                        
                    ## Stem loop finder ##
                    if counter_open2 == counter_close_sq and open2_check_sq == True:
                        #stem_file.write(str(open2_check)+str(counter_open2) + 'WOWOOW'+'\n')
                        front_3 = max(front_holder[:(forward - counter_open2)])
                        back_3 = max(back_holder[:(backward+1)])
                        op_cl_count = 0
                        if front_3 >= back_3:
                            continue
                        else:
                            for op_clos in yesmam[front_3:back_3+1]:
                                if op_clos == '[':
                                    op_cl_count += 1
                                elif op_clos == ']':
                                    op_cl_count -= 1
                        if op_cl_count == 0:                     
                            #print("SLF", str(front_3) + '-' + str(back_3))
                            #print(yesmam[front_3:back_3+1] + '\n')
                            stem_file.write(str(front_3) + '-' + str(back_3) + '\n')
                            stem_file.write((yesmam[front_3:back_3+1] + '\n'))
                            counter_open2 = 0
                            stem_file2.write(str(front_3) + '-' + str(back_3) + '\n')
                            stem_file2.write((yesmam[front_3:back_3+1] + '\n'))
                            
                    ## Intermediate zeros (zeroed out in the middle of 
                    # sequence), ONLY OK FOR PSUDOKNOTS ##
                    if NEWCOUNT_sq == NEWCOUNT_b_sq:
                        front_4 = NewCount_i_sq[NEWc_i_sq]
                        back_4 = NewCount_j_sq[NEWCOUNT_b_sq-1]
                        ## Do we reach 0? If so use that number where 
                        # first i != j if not keep going until we
                        #  find a j to make that i = 0 ## 
                        #print("NEWCOUNT", str(front_4) + '-' + str(back_4))
                        #print(yesmam[front_4:back_4+1])
                        pseudo_inside = 0 
                        for f4b4 in range(front_4,back_4+1):
                            if fmlist[f4b4] == '(' or fmlist[f4b4] == '{' or fmlist[f4b4] == '<' :
                                pseudo_inside += 1
                            elif fmlist[f4b4] == ')' or fmlist[f4b4] == '}' or fmlist[f4b4] == '>' :
                                pseudo_inside -= 1
                        n = 0
                        pseudo_outside_L = 0
                        pseudo_outside_R = 0
                        pi_po = False
                        while pseudo_inside >= 1 and pi_po == False:
                            ## Check right of the pseudoknot for
                            #  knot characters ##
                            n = n+1
                            fmlist[back_4+n]
                            if fmlist[back_4+n] == '(':
                                pseudo_outside_R = pseudo_outside_R+1
                            elif fmlist[back_4+n] == ')':
                                pseudo_outside_R = pseudo_outside_R-1
                            elif fmlist[back_4+n] == '[':
                                pseudo_outside_R = pseudo_outside_R+1
                            elif fmlist[back_4+n] == ']':
                                pseudo_outside_R = pseudo_outside_R-1
                            elif fmlist[back_4+n] == '{':
                                pseudo_outside_R = pseudo_outside_R+1
                            elif fmlist[back_4+n] == '}':
                                pseudo_outside_R = pseudo_outside_R-1
                            elif fmlist[back_4+n] == '<':
                                pseudo_outside_R = pseudo_outside_R+1
                            elif fmlist[back_4+n] == '>':
                                pseudo_outside_R = pseudo_outside_R-1
                            if pseudo_inside + pseudo_outside_R == 0 :
                                #print(yesmam, back_4, n, front_4)
                                #print(yesmam[front_4:back_4+n])
                                pi_po = True
                                stem_file2.write(str(front_4) + '-' + str(back_4+n)+ '\n')
                                stem_file2.write(yesmam[front_4:back_4+n+1] +'\n')
                                izlst.add(str(front_4) + '-' + str(back_4+n))
                        while pseudo_inside <= 1 and pi_po == False:
                            ## Check left of the pseudoknot for
                            #  knot characters ##
                            n = n-1
                            if fmlist[front_4+n] == '(':
                                pseudo_outside_L = pseudo_outside_L+1
                            elif fmlist[front_4+n] == ')':
                                pseudo_outside_L = pseudo_outside_L-1
                            elif fmlist[front_4+n] == '[':
                                pseudo_outside_L = pseudo_outside_L+1
                            elif fmlist[front_4+n] == ']':
                                pseudo_outside_L = pseudo_outside_L-1
                            elif fmlist[front_4+n] == '{':
                                pseudo_outside_L = pseudo_outside_L+1
                            elif fmlist[front_4+n] == '}':
                                pseudo_outside_L = pseudo_outside_L-1
                            elif fmlist[front_4+n] == '<':
                                pseudo_outside_L = pseudo_outside_L+1
                            elif fmlist[front_4+n] == '>':
                                pseudo_outside_L = pseudo_outside_L-1
                            if pseudo_inside + pseudo_outside_L == 0 :
                                #print(front_4, n, front_4+n, back_4)
                                #print(yesmam[front_4+n:back_4+1])
                                pi_po = True
                                stem_file2.write(str(front_4+n) + '-' + str(back_4)+ '\n')
                                stem_file2.write(yesmam[front_4+n:back_4+1] +'\n')
                                izlst.add(str(front_4+n) + '-' + str(back_4+n))
                ## FULL STEM LOOP FINDER ##
                if check_close == True and counter_close == 0:
                   stem_file.write('\n')
                   stem_file2.write('\n')
                   loop_cutoff = forward - backward
                   if loop_cutoff == 1:
                       front = max(front_holder[:loop_cutoff])
                       back = max(back_holder[:backward])
                       #print("FSLF", str(front) + '-' + str(back))
                       #print(yesmam[front:(back+1)] + '\n')
                       stem_file.write(str(front) + '-' + str(back) + '\n')
                       stem_file.write((yesmam[front:(back+1)] + '\n' +'\n'))
                       stem_file2.write(str(front) + '-' + str(back) + '\n')
                       stem_file2.write((yesmam[front:(back+1)] + '\n' +'\n'))
                   check_close = False
                   counter_open2 = 0
                   
        stem_file2.close()
stem_file.close()
stem_file2.close()
windows.close()

## Look into all_stem_x to see if there are duplicates of stems being 
# counted being counted from a single program ##
def stem_clean_set(seq):
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if ( x in seen or seen_add(x))]

def stem_clean(stem_file):
    buff = []
    for nex in stem_file:
        if nex.startswith("NEXT"):
            if buff: yield buff
            buff = [nex]
        elif nex.startswith('\n'):
            continue
        else:
            buff.append(nex)
    yield buff
    
duplicate =[]
'''
with open('all_stems.txt', 'r') as duplicate_tester:
    for head in stem_clean(duplicate_tester):
        stem_clean_set(head[1:])
        if stem_clean_set(head[1:]) != []:
            duplicate = duplicate + stem_clean_set(head[1:])'''
			
stem_file_reader = open('all_stems.txt', 'r')
for head in stem_clean(stem_file_reader):
    stem_clean_set(head[1:])
    if stem_clean_set(head[1:]) != []:
	    duplicate = duplicate + stem_clean_set(head[1:])

subtractionlist = []
for nu in duplicate:
    if not nu.startswith('(' or '.'):
        subtractionlist.append(nu.strip())

## HOLD ONTO SUBTRACTONLIST UNTIL WE START ADDING AND REMOVING
#  SCORES FROM DICT 2 ##
print('sublist',subtractionlist)
'''duplicate_tester.close()'''
stem_file_reader.close()

## Prepare data for input to dictonaries ##
searchable = []
secondsearchable = []

### Must reopen stem_file_reader here or error ###
stem_file_reader = open('all_stems.txt', 'r')
for odo in stem_file_reader:
    odo = odo.rstrip()
    ## Exactness-Stem loops exactly as they are, dot-bracket ##
    if odo.startswith('(') or odo.startswith('.'):
        searchable.append(odo)
    ## Closeness-Open to Close, stems or start to finish, 
    # just numbers ##
    elif '-' in odo:
        secondsearchable.append(odo)
stem_file_reader.close()
kirasearch = []
for kira in secondsearchable:
    if kira not in kirasearch:
        kirasearch.append(kira)

## Reference Dictionaries, score_dict holds scores and
#  will be updated as needed ##
#Structural Dictonary#
make_dict = {bable: searchable.count(bable) for bable in searchable} 
#Numerical Dictonary#
score_dict = {bable: secondsearchable.count(bable) for bable in secondsearchable} 

## Remove subtractlist from score_dict as the first edits ##
for closeness, iterations in score_dict.items():
    for subt in subtractionlist:
        if subt == closeness:
            score_dict[subt] -= 1


### SECTION 2 Identification of overlaps                               ###
### Prepare the data for comparison to itself. Must identify overlaps that
#  are potentially complementary or competitive. ###
### Start getting the 'outer' data ready for insert into models, prepare 
# for permutations of this data to be made        ###

## Terminal ends +/-2 search, makes 2 dictionaries of lists deepdict holds
#  opening basepairs +/-2 ##
## spacedict holds closing basepairs +/-2 as values and the i-j as key ##
deeplist = []
spacelist = []
deepdict = {}
spacedict = {}
for prophets in kirasearch:
    frontend = prophets.split('-')[0]
    frontendlist = []
    frontendlist.append(frontend)
    frontendlist.append(str(int(frontend)+1))
    frontendlist.append(str(int(frontend)+2))
    frontendlist.append(str(int(frontend)-1))
    frontendlist.append(str(int(frontend)-2))
    deepdict[prophets] = frontendlist
    
    backend = prophets.split('-')[1]
    backendlist = []
    backendlist.append(backend)
    backendlist.append(str(int(backend)+1))
    backendlist.append(str(int(backend)+2))
    backendlist.append(str(int(backend)-1))
    backendlist.append(str(int(backend)-2))
    spacedict[prophets] = backendlist
#print('kira', kirasearch)

these_overlap = {}
overset = set()
newoverlap = 0
## Look at things in frontend lists and see if they fit overlap parameters
for k, v in deepdict.items():
    mark1 = set()
    newOLback = 0
    gotone = 0
    for bk in v:
        for k2, v2 in deepdict.items():
            if bk in v2 and k != k2:
                mark1.add(k2)
                gotone = True
    ## Look at things in backend lists and see if they fit in overlap
    #  parameters
    if gotone == True:
        otherv = spacedict[k]
        for other_bk in otherv:
            for mk1 in mark1:
                otherv2 = spacedict[mk1]     
                if other_bk in otherv2 and k != mk1:
                    #print(k, mark1)
                    newoverlap += 1
                    newOLback += 1
                    gotantoherone = True
                    ###can use as a parameter of how well overlaps fit 
                    # together!! the more counts the more they fit ###
                    these_overlap[str(newoverlap) + '.' + str(newOLback)]= [k,mk1]
                    ### make dict with one copy of each overlap ###
                    overset.add(tuple([k,mk1])) 
                    
                elif other_bk not in otherv2:
                    gotanotherone3 = False
                    gotone = False

#print(these_overlap)
#print(overset)


## Remove duplicates from "these_overlap" ##
## these_overlap2 becomes the dict for reference and updating score ##
these_overlap2 = {}
over2 = []
overset = list(overset)
for overs in overset:
    if overs[::-1] in overset:
        overset.remove(overs[::-1])
ovcount = 0
for overs in overset:
    ovcount +=1
    these_overlap2[ovcount] = overs
#print(these_overlap2)
print('score_dict', score_dict)


### Adjust weight of specific input prediction program ###
""" Must give stem_file_x.txt identifier that lines up with correct 
prediction program
Same order as form_matrix:

CONTRAfold = 1
Centroidfold = 2
MFE = 3
MEA = 4
IPknots = 5
Pkiss = 6
"""
###### Give weight to Centroid (second program in form_matrix) +.5 
# in score_dict ####
"""
stem_file_2 = open('stem_file_2.txt', 'r')
for sf2 in stem_file_2:
    for k, v in score_dict.items():
        if k == sf2.strip('\n'):
            score_dict[k] = v + .5
    #        print('OKOKOKO', k)
stem_file_2.close()
#print(score_dict)
"""

## Find out if an area has more than 2 near perfect overlaps ##
## sometimes A overlaps B and C but B will not overlap C! ##
multinearovers = set()
howmany_mno = 0
for k, v in these_overlap2.items():
    for k2, v2 in these_overlap2.items():
        for over in v:
            if over in v2 and k != k2:
                multinearovers.add(over)
#print(multinearovers)

## Make sure things only overlap with close neighbors here ##
## Only grabbing things that are multioverlaps ###
for k,v in these_overlap2.items():
    sd_mno_count = 0
    mno_remover = set()
    for iv in v:
        if iv in multinearovers:
            sd_mno_count +=1
            caught_one = True
            mno_remover.add(iv)
    if sd_mno_count == 2:
        continue
    elif sd_mno_count == 1 and caught_one == True:
        for mnoR in mno_remover:
            multinearovers.remove(mnoR)
for x in multinearovers:
    howmany_mno +=1

## Find if i'-j' (B) overlaps i-j (A) and i''-j'' (C) or just i-j.
#  Do the same for i''-j'' ##
ovset = set()
mno_dict = {}
mno_count = 0
print("multinearovers", multinearovers)
#print(these_overlap2)
for mno in multinearovers:
    proc_count = 0 
    for mno2 in multinearovers:
        proc_count +=1
        hit_T = False
        hit_T2 = False
        hit_Tx2 = False
        ## Check the overlap of i in multinearovers ##
        for imno in range(int(mno.split('-')[0])-2, int(mno.split('-')[0])+2):
            if imno in range(int(mno2.split('-')[0])-2, int(mno2.split('-')[0])+2) and mno2 != mno:
                hit_T = True
                break
        if hit_T == True:
            ## If i overlaps check the overlap of j in multinearovers ##
            for imno in range(int(mno.split('-')[1])-2, int(mno.split('-')[1])+2):
                if imno in range(int(mno2.split('-')[1])-2,int(mno2.split('-')[1])+2):
                    hit_T2 = True
                    break
        if hit_T == True and hit_T2 == True:
            hit_Tx2 = True
        if hit_Tx2 == True:
            ## Both i and j overlap with the same structure: i'-j' ##
            already_here = False
            for k, v in mno_dict.items():
                for iv in v:
                    iv_len = int(iv.split('-')[1]) - int(iv.split('-')[0])
                    mno_len = int(mno.split('-')[1]) - int(mno.split('-')[0])
                    if int(mno.split('-')[0]) in range(int(iv.split('-')[0])-10, int(iv.split('-')[0])+10) or int(mno.split('-')[1]) in range(int(iv.split('-')[1])-10, int(iv.split('-')[1])+10):
                        ## The difference in size between "near overlaps" 
                        # cannot exceed 90% ##
                        if abs(mno_len - iv_len) > .9*iv_len:
                            continue
                        #if mno_len - iv_len > .1*mno_len or mno_len - iv_len > .1*iv_len or iv_len - mno_len > .1*mno_len or iv_len - mno_len > .1*iv_len:
                         #   continue
                        else:
                            nothere = True
                            for k2, v2 in mno_dict.items():
                                if mno in v2:
                                    nothere = False
                                    break
                            if nothere == True:
                                for imno in range(int(mno.split('-')[0]), int(mno.split('-')[1])+1):
                                    if imno in range(int(iv.split('-')[0]), int(iv.split('-')[1])+1):
                                        ## There is another overlap!! ##
                                        already_here = True
                                        new_k = k
                                        break
            if already_here == True:
                print('got added', mno, new_k)
                mno_dict[new_k].add(mno)
                break
            else:
                ## If the mno_dict is empty, because the loops just started we reach here and build the dict ##
                ## Dict holds one of the items in the multioverlap so the others can be compared ##
                stop_these = False
                for m, n in mno_dict.items():
                    if mno in n or mno2 in n:
                        stop_these = True
                if stop_these == True:
                    break
                else:
                    mno_count +=1
                    mno_dict[mno_count] = set()
                    mno_dict[mno_count].add(mno)
                    mno_dict[mno_count].add(mno2)
            break
    if hit_Tx2 == False and proc_count == howmany_mno:
        #print(mno, mno_count)
        mno_count +=1
        mno_dict[mno_count] = set()
        mno_dict[mno_count].add(mno)
        #mno_dict[mno_count].add(mno2)
print('\n')
### WATCH FOR SPECIAL CASES!!! ###
print(mno_dict)

## special_near_overs gets near overlaps from mno_count ...
## Captures weird near overlaps, B/C that don't overlap. Not seen fully 
# in mno_dict but found in these_overlap2 ##
outer_lst = set()
special_near_overs = set()
will_remove_mno = set()
for k, v in mno_dict.items():
    if len(v) == 1:
        will_remove_mno.add(k)
        for iv in v:
            score_iv = score_dict[iv]
            for k2, v2 in these_overlap2.items():
                if iv == v2[0]:
                    special_near_overs.add(v2[1])
                elif iv == v2[1]:
                    special_near_overs.add(v2[0])
print('special_near_overs', special_near_overs)

special_near_overs_score = 0
sno_dict = {}
special_near_overs_OMG = False
for sno in special_near_overs:
    if score_dict[sno] > score_iv:
        ## Means we have to figure out which one of the overlaps should 
        # get the score of the others and be used ##
        ## Sets the conditional to True leading to the decision maker 
        # between B/C when their scores are higher than A and they don't 
        # near overlap ##
        print('OH GOD WHAT HAVE I DONE, we are in for some chop')
        omg_var_count +=1
        sno_dict[score_dict[sno]].append(sno)
        special_near_overs_OMG = True
    else:
        ## Go with the 'middle one' aka the one that matches two while 
        # those two don't match each other ##
        ## Add the scores of B/C together in preparation to add to A ##
        special_near_overs_score += score_dict[sno]
        
#### This is unlikely to happen but is possible! ####
if special_near_overs_OMG == True:
    equals = False
    really_OMG = False
    for x, y in sno_dict.items():
        if len(y) > 1:
            if len(special_near_overs) -1 == omg_var_count:
                ## Multi equal more than one. Will break out of this loop,
                #  some structure shouldn't be in sno ##
                equals = True
            else:
                ## We must find the highest ##
                really_OMG = True
    if really_OMG == True:
        wew_score = 0
        ## Multi level high overlaps with different scores that dont all
        #  near overlap. wew. Getting crazy now. ##
        if len(sno_dict[max(x)]) > 1:
            sno_len = {}
            snoLengths = set()
            for g, h in sno_dict.items():
                if g == sno_dict[max(x)]:
                    ## Tempted to put these in the place of the segment 
                    # that was near overlapped with both and is in mno_dict ##
                    for ih in h:
                        sno_len[ih] = (ih.split('-')[1] - ih.split('-')[0])
                        snoLengths.add(ih.split('-')[1] - ih.split('-')[0])
            for f, j in sno_len.items():
                if j == max(snoLengths):
                    score_goes_up = f
                    for k, v in mno_dict.items():
                        if len(v) == 1:
                            for iv in v:
                                wew_score += score_dict[iv]
                                score_dict[iv] = 0
                    for k, v in sno_dict.items():
                        for iv in v:
                            wew_score += score_dict[iv]
            score_dict[score_goes_up] += wew_score
        ## Only one of the potential overlaps B/C has a highscore! Eazy 
        # peazy. Give scores from A to and loser to that ##
        elif len(sno_dict[max(x)]) == 1:
            for sno in special_near_overs:
                if sno != sno_dict[max(x)]:
                    wew_score += score_dict[sno]
                    score_dict[sno] = 0
            for k, v in mno_dict.items():
                if len(v) == 1:
                    for iv in v:
                        wew_score += score_dict[iv]
                        score_dict[iv] = 0
            score_dict[sno_dict[max(x)]] = wew_score
            outer_lst.add(sno_dict[max(x)])
## Go with middle structure: A ##
elif special_near_overs_OMG == False:
    for k, v in mno_dict.items():
        if len(v) == 1:
            for iv in v:
                ## Ok have score, now need to get sequence input ##
                ## Change the score of A so that it 
                # incorporates B/C scores##
                score_dict[iv] += special_near_overs_score
    ## Reduce score of B/C to 0 ##
    for sno in special_near_overs:
        score_dict[sno] = 0
    ## Add A to outer_lst for comparison with other overlaps ##
    for k, v in mno_dict.items():
        if len(v) == 1:
            for iv in v:
                outer_lst.add(iv)
for x in will_remove_mno:
    del mno_dict[x]
print(mno_dict)        

now_mno_dict = {}
now_mno_count = 0
for k, v in mno_dict.items():
    now_mno_count += 1
    now_mno_dict[now_mno_count] = v

## Make manager_dict, a triple tuple dictionary to hold information 
# about the near overlaps ##
manager_dict = {}
for k, v in now_mno_dict.items():
    howmany_bycase = 0
    lengths_mno = []
    for iv in v:
        howmany_bycase += 1
        lengths_mno.append(int(iv.split('-')[1]) - int(iv.split('-')[0]))
        manager_dict[iv] = (k, score_dict[iv], int(iv.split('-')[1]) - int(iv.split('-')[0]))

##  manager dict = (overlap #/unique identifier for these overlaps,
#  score, length) ##
## Add scores together, give to highest score (consensus winner) or 
# give to longest, remove from these_overlap2 and change scores 
# in score_dict##
proceed = 1
mno_dict2 = {}
uniques = 0
for x, y in mno_dict.items():
    uniques += 1
unique_count = 0 
for x, y in manager_dict.items():
    unique_count +=1


## This checks between items in manager_dict for larger score, 
# then size if scores are equal ##
## Essentially works as a double check for the THIS WILL NEVER HAPPEN, 
# but with manager_dict which is sourced from mno_dict. 
# instead of sno_dict ##
removeable = set()
rem_t = False
run_again = True
while unique_count != uniques :
    if removeable != set():
        proceed +=1 
    for rem in removeable:
        del manager_dict[rem]
        print('manager_dict', rem, manager_dict)
        unique_count -= 1
    removeable = set()    
    proc_seer = 0 
    for k, v in manager_dict.items():
        if proceed != v[0]:
                proc_seer += 1
        for k2, v2 in manager_dict.items():
            if v[0] == proceed and v2[0] == proceed and k != k2:
                print(k, v, k2, v2)
                if v[1] == v2[1]:
                    if v[2] > v2[2]:
                        ## need to find 3rd ##
                        continue
                    elif v[2] < v2[2]:
                        for g, h in now_mno_dict.items():
                            if g == v[0]:
                                rem_t = True
                                removeable.add(k)                                
                        ## need to find 3rd ##
                        continue
                    elif v[2] == v2[2]:
                        ## need to find 3rd ##
                        print('shouldnt be possible')
                        print(manager_dict)
                        only2 = 0 
                        for mdk, mdv in manager_dict.items():
                            if v[0] == mdv[0]:
                                only2 += 1
                        only2 -= len(removeable)
                        if only2 == 2:
                            removeable.add(k)
                            rem_t = True
                            print(mno_dict)
                elif v[1] > v2[1]:
                    continue
                elif v[1] < v2[1]:
                    for g, h in now_mno_dict.items():
                        if g == v[0]:
                            rem_t = True
                            removeable.add(k)
            elif proc_seer == len(manager_dict):
                for man_k, man_v in manager_dict.items():
                    if man_v == v2:
                        ## need to find 3rd ##
                        print('shouldnt be possible2')
                        print(manager_dict)
                        only2 = 0 
                        for mdk, mdv in manager_dict.items():
                            if man_v[0] == mdv[0]:
                                only2 += 1
                        only2 -= len(removeable)
                        if only2 == 2:
                            removeable.add(k)
                            rem_t = True
                            print(mno_dict)

## Remove multi near overlaps from these_overlap2 in prep for
#  JUST NEAR OVERLAPS ##
print(these_overlap2)
for k, v in manager_dict.items():
    changer = 0 
    for k2, v2 in mno_dict.items():
        for x in v2:
            if v[0] == k2 and k != x:
                #print(type(score_dict[x]), score_dict[x])
                changer += score_dict[x]
                subber = score_dict[x]
                score_dict[x] = subber - subber
    initial = score_dict[k]
    score_dict[k] = changer + initial
rem_these_overlap2 = set()
for k2, v2 in mno_dict.items():
    for x in v2:
        for s,d  in these_overlap2.items():
            if x in d:
                rem_these_overlap2.add(s)
for s in rem_these_overlap2:
    del these_overlap2[s]
print(these_overlap2)



## JUST NEAR OVERLAPS HERE ##
## Find if one of the overlaps has a higher score already ##
## 2 Near overlaps NOT multinear overlaps.
#  Those were removed from these_overlap2 ##
for marker, overlapers in these_overlap2.items():
    first_one = score_dict[overlapers[0]] 
    second_one = score_dict[overlapers[1]]
    print(overlapers[0], overlapers[1])
    if first_one > second_one and second_one == 1:
        print(first_one, second_one)
        print('F>S')
        old_second = score_dict[overlapers[1]]
        score_dict[overlapers[0]] += int(old_second)
        score_dict[overlapers[1]] -= int(old_second)
    elif second_one > first_one and first_one == 1:
        print(first_one, second_one)
        print('S>F')
        old_first = score_dict[overlapers[0]]
        score_dict[overlapers[1]] += int(old_first)
        score_dict[overlapers[0]] -= int(old_first)
    elif first_one == second_one:
        print(first_one, second_one)
        print('F=S')


## Remove anything in our dictionary below a value of 1, if 0 remove ##        
score_dict = {key: value for key, value in score_dict.items() if value > 0}

print('score_dict:', score_dict)


## Gernerator Functon ##
def by_dukat(du_check):
    buffer = []
    for odo in du_check:
        dukat = score_dict.get(odo.strip())
        ## Is it in dict2 ##
        if dukat is not None:
            yield buffer 
            buffer = [odo.strip('\n')]
        ## add in everthing that wasn't in the dict2 ##
        else:
            if 'NEXT' in odo:
                continue
            buffer.append(odo.strip('\n'))
    yield buffer


from itertools import chain
stem_lengths = []
new_odo = []
diff_list = []
final_dict = {}
ditch_first_count = 0
## Want to make new_odo which acts as a holding sequence of just 'h'
#  the same lenght as the sequence ##
## Makes final_dict which holds the actual dot-bracket structures found in
#  score_dict up to this point, will be used as a reference ##
with open('all_stems.txt', 'r') as stem_file_reader2:
    for stem_r in by_dukat(stem_file_reader2):        
        ## use ditch_first_count to clear things before the first 
        # dukat iteration ##
        ditch_first_count +=1
        ## indexing from by_dukat lists, want to remove duplicates!! ## 
        what_stem = stem_r[0]
        stem_actual = stem_r[:2]
        stem_actual.pop(0)
        full_length = stem_r[:2]
        if ditch_first_count == 1:
            new_odo_length = len(stem_actual.pop())
            new_odo = ['h']* new_odo_length
        elif ditch_first_count > 1:
            diff = int(what_stem.split('-')[1]) - int(what_stem.split('-')[0])
            diff_list.append(diff)
        ## need to force final_dict into diff_list order before 
        # adding to new_odo ## 
        if what_stem in score_dict:
            final_dict[what_stem] = stem_actual
    print('FD', final_dict)
    #print(diff_list) 
            
print('\n', '\n', 'AFTER NEAR OVERLAP SCORE CHANGE' ,'\n', '\n')


### SECTION 3 "outer" permutations                                     ###
### The rough start of this section. This is where we begin to make 
# permutations of 'outer' structures              ###
###                                                                    ###


CO = False
lst_all = []
## make lst_all from things in score_dict ##
for rg, score in score_dict.items():
    lst_all.append(rg)
lst_all = sorted(lst_all)
print(lst_all)
innerlst = set()
samelst = set()
## makes outer_lst from lst_all. i < i' and j' < j ##
for outer in lst_all:
    for rg, score in score_dict.items():
            ## rg inside of i-j add to outer_lst
        if int(rg.split('-')[0]) > int(outer.split('-')[0]) and int(rg.split('-')[1]) < int(outer.split('-')[1]):
            outer_lst.add(outer)
            ## Make a list of all part inner things too ##
        elif int(rg.split('-')[0]) == int(outer.split('-')[0]) and int(rg.split('-')[1]) < int(outer.split('-')[1]):
            samelst.add(rg)
        elif int(rg.split('-')[0]) > int(outer.split('-')[0]) and int(rg.split('-')[1]) == int(outer.split('-')[1]):
            samelst.add(rg)
## Make sure things in samelst are not in outerlst
for x in samelst:
    if x in outer_lst:
        outer_lst.remove(x)


## Move things loger than 90% of the length of the seq. to a new list.
#  These will still get their own models ##
## done to allow substructure overlaps to show up and non-overlaping 
# things to appear in the same models ##
## Stuff this long usually dominates the model, forces the creation of 
# a new model, or forces other outer structures to look like inners ##
longouters = []
longremove = []
for long in outer_lst:
    if (int(long.split('-')[1]) - int(long.split('-')[0])) > len(SEQ)*.9:
        longouters.append(long)
        longremove.append(long)

## Acutually removes >90% seq segments from outerlst ##
for longouter in longremove:
    if longouter in outer_lst:
        outer_lst.remove(longouter)

## Makes overlap_lst, which structures get added to if they are competing
#  for the same sequence space (overlaps) ##
## and one has a larger score than the other. The larger score gets added
#  to overlap_lst ##
## This works as the first level of narrowing down the competition ##
overlap_lst = set()
stopper = set()
print('outer_lst', outer_lst)
for ot in outer_lst:
    outer_area_f = int(ot.split('-')[0])
    outer_area_b = int(ot.split('-')[1])
    for rg, score in score_dict.items():
        if int(rg.split('-')[0]) in range(outer_area_f, outer_area_b+1) or int(rg.split('-')[1]) in range(outer_area_f, outer_area_b+1):
            if score_dict[rg] > score_dict[ot]:
                overlap_lst.add(rg)
                ## Adds to stopper if score are equal AND overlaps termini
                #  are within 4 nts of each other. Should have been
                #  resolved in near overlaps ##
            elif score_dict[rg] == score_dict[ot] and int(rg.split('-')[0]) in range(outer_area_f -2, outer_area_f +2) and int(rg.split('-')[1]) in range(outer_area_b -2, outer_area_b +2) and ot != rg:
                stopper.add(rg)
print('Overlap_lst', overlap_lst)

## A check back to lst_all/score_dict for things
#  in the overlap_lst
## This is used to make sure only one of the sides of the stem is
#  inside the outter stem. If the whole thing is inside it will be ##
## considered an 'inner', and used later for insertion into the outter
#  stems, removed from overlap_lst ##
for outer in lst_all:
    for rg, score in score_dict.items():
        if score >= 1 :
            if rg != outer and int(rg.split('-')[0]) >= int(outer.split('-')[0]) and int(rg.split('-')[1]) <= int(outer.split('-')[1]) and int(outer.split('-')[1])-int(outer.split('-')[0]) <= len(SEQ)-10 :
                if rg in overlap_lst:
                    overlap_lst.remove(rg)

## Compare scores of things in outer_lst, Go with higher score, and also
#  make sure they are not the same things ##
removal_set = set()
for ol in outer_lst:
    for ol2 in outer_lst:
        ## Check overlap
        for i in range(int(ol.split('-')[0]), int(ol.split('-')[1])+1):
            if i in range(int(ol2.split('-')[0]), int(ol2.split('-')[1])+1) and ol != ol2:
                ## Check score
                if ol in score_dict:
                    if score_dict[ol] > score_dict[ol2]:
                        removal_set.add(ol2)
                        #continue
                    elif score_dict[ol] == score_dict[ol2]:
                        continue
                        #removal_set.add(ol)

## Get rid of things in removal_set that are lower than
#  their overlaps in outer_lst ##                 
for k in removal_set:
    if k in outer_lst:
        outer_lst.remove(k)


NoRivals = outer_lst
NoRivals= list(NoRivals)



## Order the outer list largest to smallest ##

## Convert from set to list ##
outer_lst = list(outer_lst)
## Finds the lengths of things in outer_lst and stores in a reference:
#  size_dict ##
size_dict = {}
for k in outer_lst:
    g = int(k.split('-')[1]) - int(k.split('-')[0])
    size_dict[k] = g

## Orders the outer_lst by length of structures ##
import operator
temp_outer_lst =[]
for k,g in size_dict.items():
    vv = max(size_dict.items(), key=operator.itemgetter(1))[0]
    temp_outer_lst.append(vv)
    size_dict[vv] = 0
outer_lst = temp_outer_lst

## Want to remove thing found in overlap_lst from outer_lst if they
#  are stil there. Just a double check ##
print("OUTERLST", outer_lst)
about_to_remove = set()
for outer in outer_lst:
    for blocker in overlap_lst:       
        for i in range(int(blocker.split('-')[0]), int(blocker.split('-')[1])+1):
            if i in range(int(outer.split('-')[0]), int(outer.split('-')[1])+1):
                about_to_remove.add(blocker)
         
for x in about_to_remove:
    if x in outer_lst:
        outer_lst.remove(x)
        
print('\n',"Ordered OUTER_LST", outer_lst, '\n')

## start with highest scores and things that overlap with them ##
## Generate dicts with overlaps. The key is one i-j with the 
# value the other i-j's
outer_dict = {}
outer_dict_keep = {}
out_lst_count = 0
for ol in outer_lst:
    out_lst_count += 1
    outer_dict[ol] = set()
    outer_dict_keep[ol] = set()
    for ol2 in outer_lst:
        for i in range(int(ol.split('-')[0]),int(ol.split('-')[1])+1):
            if i in range(int(ol2.split('-')[0]),int(ol2.split('-')[1])+1) and ol != ol2:
                outer_dict[ol].add(ol2)
                outer_dict_keep[ol].add(ol2)
                break

## Order the Dicts into things that overlap the most to things that
#  overlap the least##
from collections import OrderedDict
outer_dict_sort = OrderedDict()
sort_OD_lst = []
outer_dict_keep = OrderedDict()
sort_ODK_lst = []
## Generate lists which will be used to make the ordered sorted dicts ##
for k in sorted(outer_dict, key = lambda k: len(outer_dict[k]), reverse = True):
    sort_OD_lst.append(k)
for k in sorted(outer_dict_keep, key = lambda k: len(outer_dict_keep[k]), reverse = True):
    sort_ODK_lst.append(k)

for sod in sort_OD_lst:
    for k, v in outer_dict.items():
        if sod == k:
            outer_dict_sort[sod] = v
for sodk in sort_ODK_lst:
    for k, v in outer_dict_keep.items():
        if sodk == k:
            outer_dict_keep[sodk] = v

## Look through outer_dict_sort for overlaps with multiple overlaps #

mark_for_multi_overlap = set()
for k, v in outer_dict_sort.items():
    for k2, v2 in outer_dict_sort.items():
        ## number of overlaps more than one?
        if k in v2 and len(v2) >= 2:
            #print('found a multioverlap')
            ## stuff in this multioverlap is in the outer_dict_sort
            #  dictionary. add the items that have been identified to
            #  overlap to mfmo ##
            for k3 in v2:
                if k3 in outer_dict_sort.keys():
                    if k in outer_dict_sort[k3]:
                        mark_for_multi_overlap.add(k)

print('mfmo', mark_for_multi_overlap)


### Start handeling MFMO more directly creating permutations. ###


## Remove things in mark_for_multi_overlap from outer_dict_sort ##
for mfmo in mark_for_multi_overlap:
    remover_k = set()
    for k, v in outer_dict_sort.items():
        kdel_now = False
        if mfmo == k:
            remover_k.add(mfmo)
            kdel_now = True
        if mfmo in v and kdel_now == False:
            new_ods = set()
            for u in v:
                ## For things that are not in mfmo make these the new
                #  values of outer_dict_sort ##
                if u != mfmo:
                    new_ods.add(u)
            outer_dict_sort[k] = new_ods
    for rk in remover_k:
        del outer_dict_sort[rk]
#print(ov_dict_sort)
        
## Make a list of lowest level overlaps as tuples, these will be
#  permuted into all options ##
now_outer_lst = []
overlap_tup = []
low_level_no_overlap = []
for k, v in outer_dict_sort.items():
    now_outer_lst.append(k)
    ## Overlap at the lowest level 
    # (1 vs. 1 with same score at same location) ##
    if v != set():
        for x in v:
            overlap_tup.append((k, x))
    ## Doesn't overlap at the lowest level ##
    elif v == set():
        low_level_no_overlap.append(k)
        
## Remove duplicate things that were just in reversed 
# order in overlap_tup ##
overlapt_up = list(overlap_tup)
if overlap_tup != []:
    overlap_tup.pop()
for w in overlap_tup:
    if tuple(reversed(w)) in overlap_tup:
        overlap_tup.remove(w)
for w in overlap_tup:
    if tuple(reversed(w)) in overlap_tup:
        overlap_tup.remove(w)
#print(overlap_tup)


## Clean up overlap_tup. If there is an unusual overlap 
# we can clean it here ##
for x in overlap_tup:
    for y in overlap_tup:
        for x2 in x:
            if x2 in y and x != y:
                overlap_tup.remove(y)
                overlap_tup.remove(x)


## Create all permutations of things in the lowest level of overlaps ##
import itertools
OLL =[]
for x in overlap_tup:
    OLL.append(set(x))
    
## Here we make all lower level permutations of overlaps ##
perms = [set((y) for y in x) for x in itertools.product(*OLL)]
   
## Add to our permutations things that would not be overlapping in 
# that same lower level ##
llno = list(low_level_no_overlap)
for pem in perms:
    for NR in llno:
        pem.add(NR)

#print('low_level_no_overlap', low_level_no_overlap)
print('PERMS', perms)
print('outer_dict_keep', outer_dict_keep)

## Indentification of permutaions
over_allthis = set()
new_perm_dict = {}
new_perm_dict2 = {}
new_perm_dict3 = {}
new_perm_LST = []
new_perm_LST_2 = []
already_used_perms = set()
never_past = False
for k, v in outer_dict_keep.items():
    iv_count = 0
    if k in mark_for_multi_overlap:
        ## Things here overlap every item in outer_lst. 
        # Is it's own permutation ##
        if len(v)+1 == out_lst_count:
            perms.append([k])
            over_allthis.add(k)
        ## Something else overlaps a lot but doesn't cover the whole 
        # sequence, must create a permutation that contains options ##
        elif len(v)+1 != out_lst_count:
            unusual_over = 0
            unusual_lst = []
            uul_len = 0
            ## Finds overlap where A goes over B and C but B doesn't
            #  go over C. In the lowest level ##
            for k2, v2 in outer_dict_keep.items():
                if k in v2 and k2 not in over_allthis:
                    unusual_lst.append(k2)
                    unusual_over += 1
            ## Make permutaion with just this k, 
            # then one with just things in v 
            # (assuming they only overlap with k) ##
            ## More difficult case here. Multiple non-overarching 
            # overlaps in outer_dict_keep and mark_for_multi_overlap.
            #  Now in unusual_lst ##
            how_long_uul = len(unusual_lst)
            if unusual_over > 1:
                never_past = True
                uul_count = 0
                odd_uul_out_count = 0
                seen_uul = set()
                for uul in unusual_lst:
                    ## Make sure we've seen all. 
                    # uul_count must == how_long_uul ##
                    uul_count += 1
                    uul_len = len(outer_dict_keep[uul]) - len(over_allthis)                    
                    lone_uul = []
                ## This means we have a 'B' or 'C' with multi overlaps
                #  one being 'A', the other is not 
                # corresponding 'b' or 'c'. RARE ##
                ## A in B & C & D, C not in B or D, but B, C, or D 
                # could have another overlap A does not have if there
                #  is no over_allthis AND/OR over_allthis captured
                #  every structure ## 
                if uul_len > 1 and uul_count == how_long_uul:
                    print('A in B & C & D, C not in B or D, but B, C, or D could have another overlap A does not have but over_allthis captured')
                    #print(k, v)
                if uul_len == 1 and uul_count == how_long_uul:
                    #print('we only get here if i unusual_over > 1', k, v)
                    pem = list(outer_lst)
                    pem2 = list(outer_lst)
                    ## Make permutations using 1. k or "A"; 2. 'B' & 'C';
                    #  This 'may' require more than just 2 permutations if
                    #  other things in outer_dict_keep ##
                    for item in v:
                        if item in pem:
                            ## A ## 
                            pem.remove(item)
                    ## B & C ## 
                    pem2.remove(k)
                    ## Get rid of overarching overlaps if present.
                    #  These will throw off permutes ##
                    for item in over_allthis:
                        if item in pem2:
                            pem2.remove(item)
                    """ THIS SHOULD REPLACE OWN_PERM/REMOVAL_LST """
                    for item in over_allthis:
                        if item in pem:
                            pem.remove(item)
                    ## Loop through outer_dict_keep again for production
                    #  of new_perm_dict ##
                    for k3, v3 in outer_dict_keep.items():
                        if k3 not in over_allthis and k3 not in v and k3 != k and k3 not in new_perm_dict.keys():
                            if (len(outer_dict_keep[k3]) - len(over_allthis)) >= 1:
                                if new_perm_dict == {}:
                                    ## If we get here we must duplicate
                                    #  pem once for every time/type that we get here ##                                
                                    for iv3 in v3:
                                        if iv3 not in over_allthis:
                                            new_perm_dict[k3] = list(pem)
                                            new_perm_dict[k3].remove(k3)
                                            new_perm_dict[iv3] = list(pem)
                                            new_perm_dict[iv3].remove(iv3)
                                else:
                                    ## Must duplciate new_perm_dict
                                    #  for new types then 
                                    # alter how they look##
                                    new_perm_LST = []
                                    ## Make enough duplicates for 
                                    # every key in outer_dict_keep
                                    #  execpt for those in over_allthis ##
                                    for count_ov_dict_duplicate in range(0, int(len(outer_dict_keep[k3]) - len(over_allthis))+1):
                                        old_unique_types = 0
                                        new_unique_types = 0 
                                        for npd_k, npd_v in new_perm_dict.items():
                                            new_perm_LST.append(npd_v)
                                            old_unique_types += 1
                                    ## Iterates thru unique types and 
                                    # removes iv3, the other, other 
                                    # overlap type...## 
                                    # NOT Weird types...I think
                                    for iv3 in v3:
                                        if iv3 not in over_allthis:
                                            new_unique_types +=old_unique_types
                                            for X_IV in range((0+new_unique_types),(old_unique_types+new_unique_types)):
                                                if iv3 in new_perm_LST[int(X_IV)]:
                                                    new_perm_LST[int(X_IV)].remove(iv3)
                    #print('after_pem', new_perm_dict)

                    ### DEAL WITH TYPE 2, B & C ###
                    ## First if B & C (& others) don't all overlap 
                    # with each other. make permutations where 
                    # non overlaps included ##
                    ## This complements the 'A' type permute created in
                    #  new_perm_dict/new_perm_LST ## 
                    for uul in unusual_lst:
                        for ul_k, ul_v in outer_dict_keep.items():
                            if uul in ul_v and ul_k in unusual_lst and uul not in seen_uul:
                                seen_uul.add(uul)
                                odd_uul_out_count += 1                    
                    if odd_uul_out_count != uul_count and odd_uul_out_count != 0 :
                        no_Cs_perm_set = set()
                        for uul2 in unusual_lst:
                            for uul2_v in outer_dict_keep[uul2]:
                                if uul2_v in unusual_lst and uul2_v not in over_allthis and uul2_v != k :
                                    no_Cs_perm_set.add(uul2_v)
                                    no_Cs_perm_set.add(uul2)
                        for NO_C in no_Cs_perm_set:
                            new_perm_dict2[NO_C] = list(pem2)
                            new_perm_dict2[NO_C].remove(NO_C)
                        ## Identifies things not in new_perm_dict2 but 
                        # still in outer_dict_keep, same approach done
                        #  to type 'A' above ##
                        new_perm_dict3 = {}
                        for k3, v3 in outer_dict_keep.items():
                            new_perm_dict3 = {}
                            if k3 not in over_allthis and k3 not in v and k3 != k and k3 not in new_perm_dict2.keys():
                                if (len(outer_dict_keep[k3]) - len(over_allthis)) >= 1:
                                    if new_perm_dict2 != {}:
                                        ## If we get here we must duplicate
                                        #  pem. once for every time/type 
                                        # that we get here. ##                         
                                        for weird_k, weird_v in new_perm_dict2.items():
                                            for iv3 in v3:
                                                if iv3 not in over_allthis and iv3 not in already_used_perms:
                                                    npd3_key = weird_k + k3
                                                    npd3_key2 = weird_k + iv3
                                                    new_perm_dict3[str(npd3_key)] = list(weird_v)
                                                    new_perm_dict3[str(npd3_key)].remove(k3)
                                                    new_perm_dict3[str(npd3_key2)] = list(weird_v)
                                                    new_perm_dict3[str(npd3_key2)].remove(iv3)
                                                    already_used_perms.add(iv3)
                                                    already_used_perms.add(k3)
                            ## Must replace new_perm_dict2 because it still
                            #  holds overlapping items! ##
                            if new_perm_dict3 != {}:
                                new_perm_dict2 = dict(new_perm_dict3)
                        print('new_perm_dict2 non overlaping B/C ', new_perm_dict2, new_perm_dict3)

                    ## Second if all B & C overlap ##
                    elif odd_uul_out_count == uul_count and odd_uul_out_count != 0:
                        for k3, v3 in outer_dict_keep.items():
                            #print(k3, over_allthis, v, k)
                            if k3 not in over_allthis and k3 not in v and k3 != k and k3 not in new_perm_dict2.keys():
                                if (len(outer_dict_keep[k3]) - len(over_allthis)) >= 1:
                                    if new_perm_dict2 == {}:
                                        ## If we get here we must 
                                        # duplicate pem. once for every 
                                        # time/type that we get here. ##                              
                                        for iv3 in v3:
                                            if iv3 not in over_allthis:
                                                new_perm_dict2[k3] = list(pem2)
                                                new_perm_dict2[k3].remove(k3)
                                                new_perm_dict2[iv3] = list(pem2)
                                                new_perm_dict2[iv3].remove(iv3)
                                    else:
                                        ## Must duplciate new_perm_dict 
                                        # for new types then alter how 
                                        # they look ##
                                        new_perm_LST_2 = []
                                        for count_outer_dict_duplicate in range(0, int(len(outer_dict_keep[k3]) - len(over_allthis))+1):
                                            old_unique_types = 0
                                            new_unique_types = 0 
                                            for npd_k, npd_v in new_perm_dict.items():
                                                new_perm_LST_2.append(npd_v)
                                                old_unique_types += 1
                                        ## Iterates thru unique types and 
                                        # removes iv3, the other, other 
                                        # overlap type...## 
                                        for iv3 in v3:
                                            if iv3 not in over_allthis:
                                                new_unique_types +=  old_unique_types
                                                for X_IV in range((0+new_unique_types),(old_unique_types+new_unique_types)):
                                                    if iv3 in new_perm_LST_2[int(X_IV)]:
                                                        new_perm_LST_2[int(X_IV)].remove(iv3)


                    print('after_pem 2', new_perm_dict2)
            ## There is only one overlap in the structure,
            #  not including overarching. For mark_for_multi_overlap
            #  to be true MUST have 'Z'/overarching structure ##
            ## Easier case here. Only one thing in mark_for_multi_overlap
            #  and outer_dict_keep that is non-overarching ##
            elif unusual_over == 1 and never_past == False:
                #print('yes we do come here also', k, unusual_lst)
                outer_lst_copy = list(outer_lst)
                out_c_count = 0
                for item in v:
                    if item in outer_lst_copy:
                        outer_lst_copy.remove(item)
                ## May never reach this point. It is built for very
                #  complex sequences where 2 sections A & B have various
                # levels of overlaps that require a permutation to 
                # build in a model there ##
                for p in outer_lst_copy:
                    out_c_count += 1
                    if p != k:
                        olt_count = 0
                        for olt in overlap_tup:
                            olt_count += 1
                            ## Works but only if lower overlap with 
                            # only one other thing. Catches outer_lst_copy
                            #  one vs one overlaps ##
                            if p in olt:
                                checkedall = False
                                if p == olt[1]:
                                    p_overlap = olt[0]
                                    ## Tried overlaptups and there are 
                                    # matches. make iterations of 
                                    # permutations then add to perms ##
                                    if p_overlap in outer_lst_copy:
                                        olc2 = outer_lst_copy.remove(p_overlap)
                                        perms.append(olc2)
                                elif p == olt[0]:
                                    p_overlap = olt[1]
                                    ## Tried overlaptups and there are 
                                    # matches. make iterations of 
                                    # permutations then add to perms ##
                                    if p_overlap in outer_lst_copy:
                                        olc2 = outer_lst_copy.remove(p_overlap)
                                        perms.append(olc2)
                            ## Tried all overlaptups but there are no 
                            # matches. we add in pem to perms ##
                            elif p not in olt and olt_count == len(overlap_tup) and out_c_count == len(outer_lst_copy):
                                #checkedall = True
                                perms.append(pem)
                        if overlaptup == []:
                            perms.append(outer_lst_copy)
                            break


#for k,v in new_perm_dict.items():
    #perms.append(v)
if new_perm_dict3 == {}:
    for k,v in new_perm_dict2.items():
        perms.append(v)
elif new_perm_dict3 != {}:
    for k, v in new_perm_dict2.items():
        perms.append(v)
if new_perm_LST == []:
    for k,v in new_perm_dict.items():
        perms.append(v)
for iv in new_perm_LST:
    perms.append(iv)
for iv in new_perm_LST_2:
    perms.append(iv)
                        
while None in perms:
    perms.remove(None)
            
### These are all mixed outer permutations ### the lower level perms 
# only are in '{}'. Stuff including at least one upper level perm in '[]'.
#print('PERMS', perms)

## Look for duplicates. remove until one of each type ##
remove_perm = []
[remove_perm.append(p) for p in perms if p not in remove_perm]    
perms = list(remove_perm)


### SECTION 4 Initial model production                                 ###
### This is where we do our first actual model production and insertion of
#  outer structures into these.                                        ###
###                                                                    ###

## Does actual insertion into holder dict. Makes dictionary that holds scores of perms in them for later use ##
outLST = set()
next_model = 0
outLST_MAT = {}
model_score_dict = {}
model_dict = OrderedDict()
for model in perms:
    next_model += 1
    outLST_matcher = 'outLST_%i' % (next_model)
    outLST_MAT[outLST_matcher] = set()
    build_odo = list(new_odo[:])
    ## Insert items one at a time into the empty sequence, remove 
    # sequence space that it occupies ##
    for seg in model:
        ## outLST_MAT keeps a tab of all the segments added to a model ##
        outLST_MAT[outLST_matcher].add(seg)
        outLST.add(seg)
        ## Actual insertion of segments into model ##
        jarvis = list(chain.from_iterable(final_dict[seg]))
        build_odo[int(seg.split('-')[0])] = jarvis
        del build_odo[int(seg.split('-')[0])+1:int(seg.split('-')[0])+len(jarvis)]
        ## Flatten the lists out ##
        build_odo = list(chain.from_iterable(build_odo))
    model_dict['new_model%i' % (next_model)] = list(chain.from_iterable(build_odo))
    ## Make scores for the outer inserts ##
    modelscore = 0 
    for mdl in model:
        for k, v in score_dict.items():
            if mdl == k:
                modelscore += v
    model_score_dict['score_model%i' % (next_model)] = modelscore
#print(model_dict)
    
for i, j in model_dict.items():
    j = list(j)
    j = ''.join(list(chain.from_iterable(j)))

### OK AT THIS POINT WE HAVE INSERTED OUTTER THINGS IN MOST POSSIBLE 
# ITERATIONS, GOT A SCORE OF THESE THINGS, AND MARKED WHAT WAS 
# INSERTED WHERE. ###
### NEXT WE INSERT INNER THINGS, NR, AND MID PIECES, 
# THEN WE SHOULD BE DONE!!! ###

print('\n', 'AFTER OUTER INSERTION:', '\n')

## This checks to make sure all open bases equal all close bases ##
for i, j in model_dict.items():
    spot =0 
    jfcount =0
    jbcount =0
    for x in j:
        if x == '(' or x == '[' or x == '{' or x == '<':
            jfcount +=1
        elif x == ')' or x == ']' or x == '}' or x == '>':
            jbcount +=1
    print(jfcount, jbcount)
    if jfcount != jbcount:
        print(i, '<--- messed up')
        print(''.join(list(chain.from_iterable(j))))
        continue
    elif jfcount == jbcount:
        print(i, 'all good boss man')
        #print(''.join(list(chain.from_iterable(j))), len(j))
        continue

#print(outLST)
    



### SECTION 5  "Inner" segemnt insertion                               ###
### Prepare and insert the inner segments to the models, if/when we can###
###                                                                    ###

## Ok we have model_dict with the outer things inserted in various models.
## Then we move on to inner inserserts which shouldn't have to be changed.

## Make sorted list of tuples, representative of a dict ##
score_dict_x = sorted(score_dict.items(), key=operator.itemgetter(1), reverse = True)
        
seen_inners = set()
inner_space = set()
finner = []
binner =[]
fblst=[]

"""CHANGED FROM LST_ALL TO outLST WATCH FOR RESULTS!!!!!!!!!!!!!!!!!!!"""

for inner in outLST:
    #print(inner)
    finner.append(inner.split('-')[0])
    binner.append(inner.split('-')[1])
fblst = zip(finner,binner)

keep_in_lst = set()
midpieces = set()
for f, b in fblst:
    inner_hits_pseudo = False        
    for rg, score in score_dict_x:
        ## Check to see if the inners overlap with part of a pseudoknot ## 
        if rg in izlst and int(f) > int(rg.split('-')[0] and int(rg.split('-')[1])) > int(b):
            sivraj = list(chain.from_iterable(final_dict[rg]))
            for iz in range( (int(f) - int(rg.split('-')[0])) , (int(f) - int(rg.split('-')[0])) + (int(b) - int(f)) ):
                print(iz)
                if sivraj[iz] == '(' or '[' or'<' or '{' or ')' or ']' or '>' or '}':
                    inner_hits_pseudo = True
                else:
                    inner_hits_pseudo = False			
                
        in_lst = set()
		## Highest score first ##
        if score == 6 and (int(f) < int(rg.split('-')[0]) and (int(rg.split('-')[1])) < int(b)) and inner_hits_pseudo == False:
            if f+'-'+b in outLST:
                outseg = f+'-'+b
                if outseg in score_dict.keys():
                    outscore = score_dict[outseg]
                    if outscore > score:
                        continue
                    ## Insert into model_dict ##
                    elif outscore <= score:
                        sivraj = list(chain.from_iterable(final_dict[rg]))
                        for i, j in model_dict.items():
                            j = list(j)
                            j[int(rg.split('-')[0])] = sivraj
                            j[int(rg.split('-')[0]) + 1 : int(rg.split('-')[1]) + 1] = '0'*(int(rg.split('-')[1]) - int(rg.split('-')[0]))
                            in_lst.add(rg)
                            keep_in_lst.add(rg)
                    for fou in in_lst:
                        for fs in range(int(fou.split('-')[0]), int(fou.split('-')[1]) + 1):
                                inner_space.add(fs)
                                continue
                        ######print(score_dict_x)
                        ######indexTuple(score_dict_x, score, rg)
        ## Two options, 5/4 cannot occupy the same space so it's ok. ##
        elif score == 5 or score == 4 and (int(f) < int(rg.split('-')[0]) and (int(rg.split('-')[1])) < int(b)) and inner_hits_pseudo == False:
            if f+'-'+b in outLST:
                outseg = f+'-'+b
                if outseg in score_dict.keys():
                    outscore = score_dict[outseg]
                    if outscore > score:
                        continue
                    elif outscore <= score:
                        keep_in_lst.add(rg)
                        in_lst.add(rg)
                        for fou in sorted(in_lst):
                            ## If rgs overlap with things already in 
                            # inner_space then continue to get 
                            # out of the loop##
                            if any((True for r in range(int(rg.split('-')[0]), int(rg.split('-')[1])+1) if r in inner_space)):
                                continue
                            ## Insert into model_dict ##
                            else:
                                print('inserted','3', rg)
                                sivraj = list(chain.from_iterable(final_dict[rg]))
                                for i, j in model_dict.items():
                                    ## KEEP ##
                                    j = list(j)
                                    j[int(rg.split('-')[0])] = sivraj
                                    j[int(rg.split('-')[0])+1 : int(rg.split('-')[1])+1] = '0'*(int(rg.split('-')[1])-int(rg.split('-')[0]))
                            for fs in range(int(fou.split('-')[0]), int(fou.split('-')[1])+1):
                                inner_space.add(fs)
                                continue
                            in_lst.remove(fou)
        elif score == 3 and (int(f) < int(rg.split('-')[0]) and (int(rg.split('-')[1])) < int(b)) and inner_hits_pseudo == False:
            if f+'-'+b in outLST:
                outseg = f+'-'+b
                if outseg in score_dict.keys():
                    outscore = score_dict[outseg]
                    if outscore > score:
                        continue
                    elif outscore <= score:
                        keep_in_lst.add(rg)
                        in_lst.add(rg)
                        for fou in sorted(in_lst):
                            ## If rgs overlap with things already in
                            #  inner_space then continue to get 
                            # out of the loop##
                            if any((True for r in range(int(rg.split('-')[0]), int(rg.split('-')[1])+1) if r in inner_space)):
                                    continue
                            ## Insert into model_dict ##
                            else:
                                sivraj = list(chain.from_iterable(final_dict[rg]))
                                for i, j in model_dict.items():
                                    j = list(j)
                                    j[int(rg.split('-')[0])] = sivraj
                                    j[int(rg.split('-')[0])+1 : int(rg.split('-')[1])+1] = '0'*(int(rg.split('-')[1])-int(rg.split('-')[0]))
                            for fs in range(int(fou.split('-')[0]), int(fou.split('-')[1])+1):
                                inner_space.add(fs)
                                continue
                            in_lst.remove(fou)
        ## Two options, 2 will always go first because we ordered
        #  the tuple score_dict_x ##
        elif score == 2 or score == 1 and (int(f) < int(rg.split('-')[0]) and (int(rg.split('-')[1])) < int(b)) and inner_hits_pseudo == False:
            ## Check that the inner loop wasn't graded lower overall
            #  than the outter, reject if it was ##
            if f+'-'+b in outLST:
                outseg = f+'-'+b
                if outseg in score_dict.keys():
                    outscore = score_dict[outseg]
                    if outscore > score:
                        continue
                    elif outscore <= score:
                        keep_in_lst.add(rg)
                        for fou in sorted(in_lst):
                                ## If rgs overlap with things already
                                #  in inner_space then continue to get 
                                # out of the loop##
                                if any((True for r in range(int(rg.split('-')[0]), int(rg.split('-')[1])+1) if r in inner_space)):
                                    continue
                                ## Insert into model_dict ##
                                else:
                                    sivraj = list(chain.from_iterable(final_dict[rg]))
                                    for i, j in model_dict.items():
                                        j = list(j)
                                        j[int(rg.split('-')[0])] = sivraj
                                        j[int(rg.split('-')[0])+1 : int(rg.split('-')[1])+1] = '0'*(int(rg.split('-')[1])-int(rg.split('-')[0]))
                                for fs in range(int(fou.split('-')[0]), int(fou.split('-')[1])+1):
                                    inner_space.add(fs)
                                    continue
                                in_lst.remove(fou)
                                
## CLEAR OUT EMPTY SPACE '0'S AND FLATTEN OUT INNER LISTS ##    
count_ad = 0
for i, j in model_dict.items():
    count_ad +=1
    j = list(chain.from_iterable(j))
    while '0' in j:
        del j[j.index('0')]
    model_dict[i] = j

print('\n', 'AFTER INNER INSERTIONS:', '\n') 

## Check for any errors in inserters. make sure
#  all ( have corresponding ) ##
for i, j in model_dict.items():
    spot =0 
    jfcount =0
    jbcount =0
    for x in j:
        if x == '(' or x == '[' or x == '{' or x == '<':
            jfcount +=1
        elif x == ')' or x == ']' or x == '}' or x == '>':
            jbcount +=1
    print(jfcount, jbcount)
    if jfcount != jbcount:
        print(i, '<--- messed up')
        print(''.join(list(chain.from_iterable(j))))
        continue
    elif jfcount == jbcount:
        print(i, 'all good boss man')
        #print(''.join(list(chain.from_iterable(j))), len(j))
        continue
    
### SECTION 6 Remainder insertion                                      ###
### Now we look for any structures remaining that may fit into our 
# models without interference. It is unlikely to find                  ###
### too many but hey! Insert these structures if possible              ###


## Counts up the score of outer things that have alredy been 
# inserted to model_dict ## 
consensus_count_dict = {}
for k, v in model_dict.items():
    consensus_count_dict[k] = 0
for rg, sc in score_dict.items():
    for olm, mal in outLST_MAT.items():
        if rg in mal:
            for cck, ccv in consensus_count_dict.items():
                if olm.split('_')[1] == cck.split('model')[1]:
                    consensus_count_dict[cck] += sc

print("consensus_count_dict before mid-pieces added", consensus_count_dict, '\n', '\n')            


if consensus_count_dict == {}:
    consensus_count_dict['new_model1'] = 0
            

## INSERT NON-INNER/OUTTER STRANDS THAT HAVE CONSENSUS SCORES ##

## For each model produce a dict that includes all things CAN'T 
# be inserted into it. stopper_dicts ##
stopper_dict = {}
for k, v in model_dict.items():
    stopper_dict['stopper_' + str(k.split('model')[1])] = set()

## Captures items from score_dict in stopper_dict if they are not in 
# outer or inner things but they overlap with outerLST_MAT things ##
didwereachhere = False
for rg, score in score_dict.items():
    dont_proc = 0
    if rg not in outLST and rg not in keep_in_lst:
        #print(rg, 'NOT in either outter or inner')
        for ke, val in outLST_MAT.items():
            for real_val in val:
                for i in range(int(rg.split('-')[0]), int(rg.split('-')[1])+1):
                    if i in range(int(real_val.split('-')[0]), int(real_val.split('-')[1])+1):                        
                        didwereachhere = True
                        for k, v in model_dict.items():
                            if str(ke.split('_')[1]) == str(k.split('model')[1]) :
                                stopper_dict['stopper_'+str(k.split('model')[1])].add(rg)
                                break
    elif rg in outLST or rg in keep_in_lst:
        stopper.add(rg)

## We must have didwereachhere to prevent strange cases. There are 
# sometimes weird cases that do overlap that have not appeared yet 
# because they have no underneath structures ##
## Add things to midpieces if they havn't been used yet and are not 
# rivaled in their location in the sequence ##
if didwereachhere == False:
    for rg, score in score_dict.items():
        if rg not in outLST and rg not in keep_in_lst and rg in NoRivals:
            midpieces.add(rg)
            
## Generates total_in_stop_dict for use in next segment ##
total_in_stop_dict =0
for y, l in stopper_dict.items():
    total_in_stop_dict +=1

## if equal score in same location, go with longer ##
all_data = []
for rg, score in score_dict.items():
    mid_add = True
    THIS_RG = 0
    mark_these = set()
    rg_count = 0
    stop_key_count = 0
    for stop_k, stop_rg in stopper_dict.items():
        stop_key_count += 1
        if rg in stop_rg:
            rg_count +=1
        elif rg not in stop_rg:
            mark_these.add(stop_k.split('_')[1])
        ## Every stopper_dict has this RG ##
        if rg_count == total_in_stop_dict:
            THIS_RG = True
        if rg_count != total_in_stop_dict and stop_key_count == total_in_stop_dict and rg_count > 0:
            rg_count = 0
            stop_key_count = 0
            ## Find which model the rg doesn't overlap with. We want to 
            # insert it there. Already identified by mark_these ##            
            for k, v in model_dict.items():
                for num in mark_these:
                    if num == k.split('model')[1]:
                        ## Create list of tripple tuples to add later, 
                        # (i, j, k) where i = model identifier, j = rg, 
                        # k = score ##
                        all_d = (num, rg, score)
                        all_data.append(all_d)
        elif total_in_stop_dict > 0 and stop_key_count == total_in_stop_dict and rg_count == 0:
            THIS_RG = False
    ## rg is found in every stopper_dict. This exists as a double check ## 
    if THIS_RG == True:
        THIS_RG = 0
        if rg not in outLST and rg not in keep_in_lst:
            for ol in outLST:
                for i in range(int(rg.split('-')[0]), int(rg.split('-')[1])+1):
                    if i in range(int(ol.split('-')[0]), int(ol.split('-')[1])+1):
                        stopper.add(rg)
                        break
            for kl in keep_in_lst:
                ## If rg is not in outLST, keep_in_lst or range of stuff 
                # of stopper, and has score >= 3 add to midpieces ##
                for xx in range(int(rg.split('-')[0]),int(rg.split('-')[1])+1):
                    if xx in range(int(kl.split('-')[0]), int(kl.split('-')[1])+1):
                        stopper.add(rg)
                        continue
                    else:
                        if score >= 3 and rg not in stopper and (int(rg.split('-')[1])-int(rg.split('-')[0])) < (len(SEQ)-10):
                            midpieces.add(rg)
            ## Add to midpieces if no inner pieces exist and its not in 
            # outer pieces if score >= 3 ## 
            if keep_in_lst == set():
                if score >= 3 and rg not in stopper and (int(rg.split('-')[1])-int(rg.split('-')[0])) < (len(SEQ)-10):
                    midpieces.add(rg)
    elif THIS_RG == False:
        if score >= 3 and rg not in stopper and (int(rg.split('-')[1])-int(rg.split('-')[0])) < (len(SEQ)-10):
            for k, v in model_dict.items():
                for iv in range(int(rg.split('-')[0]), int(rg.split('-')[1])+1):
                    if v[iv] != 'h':
                        mid_add = False
            if mid_add == True:
                midpieces.add(rg)
        continue
    


## Want to really make sure that thigns in midpieces don't interfere 
# with outLST stuff, (from outer permutations) ## 
for mp in midpieces:
    for g in outLST:
        for i in range(int(mp.split('-')[0]), int(mp.split('-')[1])+1):
            if i in range(int(g.split('-')[0]), int(g.split('-')[1])+1):
                stopper.add(mp)


## Want to see how many potential inserts from all_data we have for every 
# model, count_ad is a count of the models ##
reqtup = {}
for xer in range(0, int(count_ad)+1):
    howmany = 0
    for x, y, z in all_data:
        if str(xer) == str(x):
            howmany += 1
            reqtup[x] = howmany
            
## This differentiates the leftovers into 3 types ##
new_one =0
effector = set()
ineffector = set()
maybe_effector = set()
for xer in range(0, int(count_ad)+1):
    new_one += 1
    same_one = 0
    hold_round= {}
    hold_round_tup = []
    for x, y, z in all_data:
        ## Find a unique model identifier that matches xer, this allows 
        # for differentiation between (y,z) tuples ##
        if str(x) == str(xer):
            same_one += 1
            hold_round[str(same_one)] = y
            hold_round_tup.append((y, z))
            ## Look through reqtup, which holds all_data model ident. 
            # as keys and numbers of these as values ##
            for typeer, howmany in reqtup.items():
                ## When all of one model type have been seen 
                # (same_one == howmany) and the type is still the 
                # same (typeer == x) ##
                if same_one == howmany and str(typeer) == str(x):
                    ## This process makes sure that the items in 
                    # all_dict for a specific model don't overlap 
                    # with each other ##
                    for h, g in hold_round.items():
                        checked_all = 0
                        count_all = 0
                        for m, n in hold_round.items():
                            count_all += 1
                            first_score = 0
                            second_score = 0
                            ## Check for overlap. If they do overlap,
                            #  check score, then length ##
                            for i in range(int(g.split('-')[0]), int(g.split('-')[1])+1):
                                if i in range(int(n.split('-')[0]), int(n.split('-')[1])+1) and h != m:
                                    ## Check scores, go with higher ##
                                    for rang, scor in hold_round_tup:
                                        if g == rang:
                                            first_score = scor
                                        if n == rang:
                                            second_score = scor
                                    if first_score > second_score:
                                        effector.add((x, g))
                                        ineffector.add((x, n))
                                    elif first_score < second_score:
                                        effector.add((x,n))
                                        ineffector.add((x, g))
                                    ## If scores equal check length ##
                                    elif first_score == second_score:
                                        ## Check length, go with longer ##
                                        if (int(g.split('-')[1]) - int(g.split('-')[0])) > (int(n.split('-')[1]) - int(n.split('-')[0])):
                                            effector.add((x, g))
                                            ineffector.add((x,n))
                                        elif (int(g.split('-')[1]) - int(g.split('-')[0])) < (int(n.split('-')[1]) - int(n.split('-')[0])):
                                            effector.add((x, n))
                                            ineffector.add((x,g))
                                        ## length equal, just go with 
                                        # first one.. 
                                        # This is very unlikely ## 
                                        else:
                                            print("TOSS UP......AKA YOU HAVE A CRAZY SEQUENCE!!!")
                                            effector.add((x, g))
                                            ineffector.add((x,n))
                                ## If no overlap add to maybe_effector 
                                # for potential additon to models at a 
                                # later point ##
                                elif i not in range(int(n.split('-')[0]), int(n.split('-')[1])+1) and h !=m:
                                    checked_all += 1
                                    if checked_all == count_all:
                                        maybe_effector.add((x, g))


for ie in ineffector:
    if ie in maybe_effector:
        maybe_effector.remove(ie)
    if ie in effector:
        effector.remove(ie)
for E in effector:
    if E in maybe_effector:
        maybe_effector.remove(E)
        
## Things in effector are overlaping, probably in other sets also. ## 
## Can be the same structure but but one is missing the 
# longer end or segemnt ## 

## Add in the actual structures from effector. Add up score in 
# consensus_count_dict ##
for mark, rg in effector:
    for k, v in model_dict.items():
        if k.split('model')[1] == mark:
            chain_maker = list(chain.from_iterable(final_dict[rg]))
            v[int(rg.split('-')[0])] = chain_maker
            v[int(rg.split('-')[0])+1 : int(rg.split('-')[1])+1] = '0'*(int(rg.split('-')[1])+1-int(rg.split('-')[0]))
            v = v.remove(v[int(rg.split('-')[0])+1])         
    for k, score_v in score_dict.items():
        if rg == k:
            for cck, ccv in consensus_count_dict.items():
                if cck.split('model')[1] == mark:
                    consensus_count_dict[cck] += score_v

for i, j in model_dict.items():
    while '0' in j:
        del j[j.index('0')]
    model_dict[i] = j
           
moreleftover_data = set()
for ME in maybe_effector:
    if ME not in moreleftover_data:
        moreleftover_data.add(ME)
        
if moreleftover_data != maybe_effector:
    print("WHAT IN TARNATION!")

## If nothing got added in effector or ineffector but data is 
# in all_data then add this stuff to moreleftover_data ###
if effector == set() and ineffector == set() and all_data != set():
    print(' Nothing got identified in all_data but something is there')
    for x, y, z in all_data:
        moreleftover_data.add((x, y))

## If stuff in moreleftover_data or midpieces overlaps with stuff already
#  in the models it should be added to clean_mlod/clean_midpieces ##
clean_mlod = set()
clean_midpieces = set()
for k, v in model_dict.items():
    v = (list(chain.from_iterable(v)))
    model_dict[k] = v
    for mark, rg in moreleftover_data:
        if mark == k.split('model')[1]:
            for i in range(int(rg.split('-')[0]), int(rg.split('-')[1])+1):
                if v[i] == '(' or v[i] == ')' or v[i] == '.' or 
                v[i] == '[' or v[i] == ']' or v[i] == '{' or 
                v[i] == '}' or v[i] == '<' or v[i] == '>' :
                    clean_mlod.add((mark, rg))
    for mp in midpieces:
        for i in range(int(mp.split('-')[0]), int(mp.split('-')[1])+1):
            if v[i] == '(' or v[i] == ')' or v[i] == '.' or v[i] == '[' or
             v[i] == ']' or v[i] == '{' or v[i] == '}' or
              v[i] == '<' or v[i] == '>' :
                clean_midpieces.add(mp)
                       
## Remove cleam_mlod stuff from moreleftove_data ##
for cmlod in clean_mlod:
    if cmlod in moreleftover_data:
        moreleftover_data.remove(cmlod)

## Remove clean_midpieces from midpieces ##
for cmp in clean_midpieces:
    if cmp in midpieces:
        midpieces.remove(cmp)

## Add moreleftover_data to the models, update scores in 
# consensus_count_dict ##
for k, v in model_dict.items():
    knumb = str(k.split('model')[1])
    for mark, rg in moreleftover_data:
        print(mark)
        if str(mark) == knumb:         
            chain_maker = list(chain.from_iterable(final_dict[rg]))
            v[int(rg.split('-')[0])] = chain_maker
            v[int(rg.split('-')[0])+1 : int(rg.split('-')[1])+1] = '0'*(int(rg.split('-')[1])+1-(int(rg.split('-')[0])))

        for k, score_v in score_dict.items():
            if rg == k:
                for cck, ccv in consensus_count_dict.items():
                    if cck.split('model')[1] == mark:
                        consensus_count_dict[cck] += score_v
   
for i, j in model_dict.items():
    j = list(''.join(list(chain.from_iterable(j))))
    while '0' in j:
        del j[j.index('0')]
    model_dict[i] = j

print('\n', 'AFTER EFFECTOR & MAYBE_EFFECTOR:', '\n')
    
## Checks open to close bases after moreleftover_data and 
# effector are added ##
for i, j in model_dict.items():
    spot =0 
    jfcount =0
    jbcount =0
    for x in j:
        if x == '(' or x == '[' or x == '{' or x == '<':
            jfcount +=1
        elif x == ')' or x == ']' or x == '}' or x == '>':
            jbcount +=1
    print(jfcount, jbcount)
    if jfcount != jbcount:
        print(i, '<--- messed up')
        print(''.join(list(chain.from_iterable(j))))
        continue
    elif jfcount == jbcount:
        print(i, 'all good boss man')
        #print(''.join(list(chain.from_iterable(j))), len(j))
        continue

    
for st in stopper:
    if st in midpieces:
        midpieces.remove(st)

## Sort midpieces ##
#print('\n', 'MIDPIECES', midpieces, '\n')
midpieces = sorted(list(midpieces))

## IF NO OUTER/INNER PIECES EXIST MUST MAKE A DICT TO INSERT MIDS INTO ##
if model_dict == OrderedDict():
    model_dict['new_model1'] = list('h'*int(int(len(SEQ))-int(1)))

for i, j in model_dict.items():
    j = list(''.join(list(chain.from_iterable(j))))
    while '0' in j:
        del j[j.index('0')]
    model_dict[i] = j
    
## Insert midpieces ##
for k, v in model_dict.items():
    for mp in midpieces:        
        chain_maker = list(chain.from_iterable(final_dict[mp]))
        v[int(mp.split('-')[0])] = chain_maker
        del v[int(mp.split('-')[0])+1:int(mp.split('-')[0])+len(chain_maker)]
        v = list(chain.from_iterable(v))
        ## Add scores of midpieces to consensus_count_dict ## 
        for i, j in score_dict.items():
            if mp == k:
                for cck, ccv in consensus_count_dict.items():
                    consensus_count_dict[cck] += j            

                
## Flatten out the models and remove 0 if present ##
for i, j in model_dict.items():
    j = list(''.join(list(chain.from_iterable(j))))
    while '0' in j:
        del j[j.index('0')]
    model_dict[i] = j


## MAKE SURE model_dicts ARE NAMED CORRECTLY ##
al_o_count = 0
n_model_dict = OrderedDict()
n_consensus_count_dict = OrderedDict()
changed = False
must_do = False
for k, v in model_dict.items():
    al_o_count +=1
    if k != ('new_model%i' % (al_o_count)):
        must_do = True
    elif k == ('new_model%i' % (al_o_count)):
        continue
    
if must_do == True:
    dumblist = []
    ccdumblist = []
    for i in range(1, al_o_count+1):
        dumblist.append(('new_model%i' % i))
        ccdumblist.append(('new_model%i' % i))
    for x, y in model_dict.items():
        for l in dumblist:   
            n_model_dict[l] = y
            dumblist.pop(0)
    for g, h in consensus_count_dict.items():
        for l in ccdumblist:
            n_consensus_count_dict[l] = h
            ccdumblist.pop(0)
    consensus_count_dict = n_consensus_count_dict
    model_dict = n_model_dict



### SECTION 7 Prepare for Round_2                                      ###
### Double check the open vs closed base pairs, write in required 
# format for second level prediction programs                          ###
### close files                                                        ###



## HERE WE TEST EACH TO BE SURE THAT WE HAVE EQUAL NUMBERS OF OPEN AND 
# CLOSING BASE PAIRS after everything 
# (first check after midpieces added) ##
## Prepare for the second level of prediction. 
sec_count =0
ccd_count =0 

sec_run = open('sec_run_MFE.txt', 'w')
sec_run2 = open('sec_run_cent.txt', 'w')
catchall = open('sec_catchall', 'a+')
h_count_file = open('h_counts.txt' , 'a+')
for i, j in model_dict.items():
    ccd_count +=1
    spot =0 
    jfcount =0
    jbcount =0
    for x in j:
        if x == '(':
            jfcount +=1
        elif x == ')' :
            jbcount +=1
            continue
    if jfcount != jbcount:
        print(i, '<--- messed up')
        print(''.join(list(chain.from_iterable(j))))
        continue
    elif jfcount == jbcount:
        h_count = 0 
        print(i, 'all good boss man')
        print(''.join(list(chain.from_iterable(j))), len(j))
        sec_run.write('\n'+'>new_'+str(seq_name.strip('1>').strip('\n'))+'_'+i+'\n')
        sec_run.write(SEQ)
        sec_run2.write('\n'+'>new_'+str(seq_name.strip('1>').strip('\n'))+'_'+i+'\n')
        sec_run2.write(SEQ)
        catchall.write('\n'+str(seq_name.strip('\n')))
        catchall.write('\n'+'>'+i+'\n')
        catchall.write(SEQ)
        for ob in j:
            spot +=1
            if ob in ( '(', ')'):
                catchall.write(ob)
                sec_run.write(ob)
                sec_run2.write(ob)
            elif ob in ('.'):
                catchall.write('.')
                sec_run.write('x')
                sec_run2.write('.')
            ## h is here to account for not b, o, or c in the structure.##
            elif ob in ('h'):
                catchall.write('h')
                h_count += 1
                sec_run.write('.')
                sec_run2.write('?')
        for k, v in consensus_count_dict.items():
            if k == 'new_model%i' % (ccd_count):
                h_count_file.write(str(seq_name.strip('\n'))+'    '+str(k)+'    '+str(h_count)+'/'+str(len(SEQ))+'\n')
                print(k,'score:', str(v),'\n')
                sec_run.write('\n'+str(v))
                sec_run2.write('\n'+str(v))
                catchall.write('\n'+str(v))
        continue


catchall.close()
h_count_file.close()
sec_run.close()
sec_run2.close()
