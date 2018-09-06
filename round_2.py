### Uses constraints from decide.py to predict a final model for users ###
### Author: Benjamin Skidmore                                          ###
### Date: 10/9/2017                                                    ###

import os
from subprocess import Popen, STDOUT, call, check_output
import subprocess
from collections import OrderedDict

count =0
char = 'A', 'U', 'C', 'G', 'T', 'a', 'u', 'c', 'g', 't', 'N', 'n'
seqstuff = '.', '(', 'x', '?'
inters = '1', '2', '3', '4', '5', '6', '7', '8', '9', '0'

final_MFE_results = open('final_MFE_results.txt', 'w')
mlst = ["None"]*count
secrun1 = open('sec_run_MFE.txt', 'r')
after_merg = open('sec_run_ready.txt', 'w')

titlelst = []
fsints = []
scins = []
ccd_dict = {}
otherccd = []
ccdset = set()
fTF = False
sTF = False
srcount = 0
for sr in secrun1:
    if sr.startswith('>'):
        srcount +=1
        titlelst.append((srcount, sr))
    elif sr.startswith(char):
        fsints.append((srcount, sr))
        fTF = True
    elif sr.startswith(seqstuff):
        scins.append((srcount, sr))
        secins = sr
        sTF = True
    elif sr.startswith(inters):
        ccd_dict[srcount] = sr.strip('\n')
        otherccd.append(sr.strip('\n'))
        ccdset.add(sr.strip('\n'))

print(ccd_dict)
occdict = OrderedDict()
howmany = 0
cscount = 0

## CREATE NEW DICT LIKE CCD_DICT BUT IN CORECT ORDER ##
for cs in ccdset:
    ## this gives the total number of types ##
    cscount += 1
othercount =0
for otc in otherccd:
    othercount +=1
## prevents potentail problem where we have more in the set but these 
# are not of the max type ##
## Tests to see if there is a predicted model that is the same in 
# sec_run_1 ##
if othercount == cscount:
    ## then do max, remove, repeat ##
    while otherccd != str(0):
        for k, v in ccd_dict.items():
            if v == max(otherccd):
                occdict[k] = v
                otherccd.remove(v)
            if otherccd == []:
                otherccd = str(0)
## There are identical models ##
elif othercount != cscount:
    maxc = 0
    comedown = False
    holdthese = 0
    while otherccd != str(0):
        for k, v in ccd_dict.items():
            if otherccd == []:
                otherccd = str(0)
            elif v == max(otherccd) and k not in occdict.keys():
                holdthese = v
                maxc +=1
                if maxc > 1:
                    comedown = True
                    occdict[k] = v
                    maxc -= 1
                    otherccd.remove(v)
            elif holdthese != max(otherccd):
                maxc = 0
                comedown = False
                
## Make sec_run_ready.txt from occdict stuff ##
for x, y in occdict.items():
    for j, k in titlelst:
        if j == x:
            after_merg.write(k)
    for j, k in fsints:
        if j == x:
            after_merg.write(k)
    for j, k in scins:
        if j == x:
            after_merg.write(k)
            after_merg.write(str(y)+'\n')
after_merg.close()
after_merg2 = open('sec_run_ready.txt', 'r')

firTF = False
secTF = False
og_mod_numb = []
for am in after_merg2:
    if am.startswith('>'):
        og_mod_numb = am
    elif am.startswith(char):
        firstins = am
        firTF = True
    elif am.startswith(seqstuff):
        secins = am
        secTF = True
    elif am.startswith(inters):
        ccd_per = am        
        if firTF == True and secTF == True:
            firstins = firstins.encode()
            secins = secins.encode()
            resulta = Popen(['RNAfold -C - - enforceConstraint'], shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=STDOUT)
            resultb = resulta.communicate(firstins+secins)[0]
            firTF = False
            secTF = False
            final_MFE_results.write(og_mod_numb.strip('\n') + '_MFE' + '\n')
            final_MFE_results.write(str(resultb.decode()))
            final_MFE_results.write(ccd_per)

print('Done with RNAfold')            
after_merg3 = open('sec_run_cent.txt', 'r')
final_cent_results = open('final_cent_results.txt', 'w')
mod_numb = []
seqTF = False
ssTF = False
for am in after_merg3:
    if am.startswith('>'):
        mod_numb = am
    if am.startswith(char):
        afm3_seq = am
        seqTF = True
    elif am.startswith(seqstuff):
        afm3_ss = am
        ssTF = True
    elif am.startswith(inters):
        afm3_int = am        
        mr2=["None"]*1
        if seqTF == True and ssTF == True:            
            sec_cent = open('sec_cent.fasta', 'w')
            sec_cent.write(mod_numb.strip('\n')+'_Cent'+'\n')
            sec_cent.write(afm3_seq)
            sec_cent.write(afm3_ss)
            sec_cent.close()
            resulta2 = str(check_output(['centroid_fold', '-C', 'sec_cent.fasta'], universal_newlines=True))
            mr2[0] = resulta2
            final_cent_results.write(resulta2)
            final_cent_results.write(afm3_int)
            seqTF = False
            ssTF = False
final_cent_results.close()                        
final_MFE_results.close()
sec_cent.close()
os.remove('sec_cent.fasta')
print('Done with centroid_fold')

## checks to see if the results are all the same or different ##
fmr2 = open('final_MFE_results.txt', 'r')
fcr_2 = open('final_cent_results.txt', 'r')
check_diff = []
MFE_titles = []
MFE_inters = []
use_seq1 = False
fr_count = 0 
for fr in fmr2:
    if fr.startswith('>'):
        MFE_titles.append(fr)
    elif fr.startswith(char):
        MFE_seq = fr
        use_seq1 = True
        continue
    elif fr.startswith(inters):
        MFE_inters.append(fr)
    elif fr.startswith(seqstuff):
        fr_count += 1
        check_diff.append(str(fr.split(' ')[0]))
    elif fr in check_diff:
        print('already here!')
        MFE_titles.pop(-1)
        MFE_inters.pop(-1)
        
## Parse final_centroid_results, if the model aready exists (rare)
#  it will give a 'aready here' ##
cent_titles = []
cent_inters = []
for fr in fcr_2:
    if fr.startswith('>'):
        cent_titles.append(fr)
    elif fr.startswith(char):
        continue
    elif fr.startswith(inters):
        cent_inters.append(fr)
    elif fr.startswith(seqstuff):
        fr_count += 1
        check_diff.append(str(fr.split(' ')[0]))
    elif fr in check_diff:
        print('already here!')
        cent_titles.pop(-1)
        cent_inters.pop(-1)

## Write all the models to the 'final_results.txt' file ##        
Titles = MFE_titles + cent_titles
Weights = MFE_inters + cent_inters

## Grab title ##
fftitle = MFE_titles[0].split('>new_')[1]
print(fftitle)
fftitle2 = fftitle.split('.FASTA')[0]

finalfinal = open('final_results.txt', 'w')
keep_final = open('final_results_keeper', 'a')
finalfinal.write(MFE_seq)
keep_final.write('\n'+'\n'+'>'+fftitle2+'.fasta'+'\n')
keep_final.write(MFE_seq)
for i in range(0,fr_count):
    finalfinal.write(Titles[i])
    finalfinal.write(check_diff[i]+'\n')
    finalfinal.write(Weights[i])
    keep_final.write(Titles[i])
    keep_final.write(check_diff[i]+'\n')
    keep_final.write(Weights[i])
finalfinal.close()
keep_final.close()
os.remove('rna.ps')
