### Initiate preditions of secondary structure on RNAs using established 
# methods                            ###
### Author: Benjamin Skidmore                                          ###
### Date: 10/14/2017                                                   ###

from subprocess import *
import os
import sys
import time

""" must provide the DATAPATH to RNAstructure for MFE and MEA to work 
properly. This will be in the users install location for these 
programs """
# export DATAPATH=programs/RNAstructure/data_tables/

## Set the time to start. measure time to run entire predict.py ##
start_time = time.time()
					
### Prep files for Running Programs ###
makeUpper = open("job_file.txt", 'r')
makeUpperWrite = open("job_file_fixed.txt", 'w')
for s in makeUpper:
    print(s)
    s = str(s).upper()
    makeUpperWrite.write(s)
makeUpper.close()
makeUpperWrite.close()

record_iter = open('job_file_fixed.txt', 'r')
for ln in record_iter:
    if ln.startswith('>'):
        filename = "1.fasta" 
        newfile = "MFE_1.ct" 
        newfile2 = "ME_1.ct"
        nowfasta = open("now_1.fasta", 'w')
        nowfasta2 = open("now2_1.fasta", 'w')
        fasta_open = open(filename, "w")
        seqfile = open('1.seq', 'w')
        fasta_open.write(ln)
    else:
        fasta_open.write(ln)
        seqfile.write(ln.strip('\n'))
fasta_open.close()
seqfile.close() 

### Begin to run programs ###
filename = "1.fasta"
char = 'A', 'U', 'C', 'G', 'T', 'a', 'u', 'c', 'g', 't', 'N', 'n'
char2 = '.', '(', '[', '<', '{'
merge = open('merge.txt', 'w')
makego = [1]
## CONTRAfold ##
handle1 = open('result1.txt','w')
mlist1 = ['None']*1
for x in makego:
    result1 = str(check_output(['contrafold', 'predict', filename], universal_newlines = True))
    print('contra_done')
    mlist1[0] = result1
for i in mlist1[0]:
    handle1.write(i)
handle1.close()
## Make merge file with all the outputs ## 
handle1_re = open('result1.txt', 'r')
for i in handle1_re:
    if i.startswith('>'):
        merge.write(str(1) + i)
    elif i.startswith(char):
        merge.write(str(1) + i)
    elif i.startswith(char2):
        merge.write(str(1) + i)
handle1_re.close()

## Centroidfold ##
handle2 = open('result2.txt','w')
centkeeper = open('centkeeper.txt', 'a')
mlist2 = []
for x in makego:
    result2 = str(check_output(['centroid_fold', filename], universal_newlines = True))
    print('Centroid_done')
    mlist2 = result2
for i in mlist2:
    handle2.write(i)
    centkeeper.write(i)
handle2.close()
centkeeper.close()
## Add to merge file ## 
handle2_re = open('result2.txt','r')
for i in handle2_re:
    if i.startswith(char2):
        merge.write(str(1) + i + '\n')
handle2_re.close()

## MFE / RNAstructure Fold##
handle3 = open('result3.txt', 'w')
newfile = "MFE_1.ct"
nowfasta = 'now_1.fasta'
grabthat = open('now_1.fasta', 'r')
check_output(['Fold', filename, newfile ], universal_newlines = True)
check_output(['ct2dot', newfile, '1', nowfasta], universal_newlines = True)
print('MFEdone')
for i in grabthat:
    handle3.write(i)
grabthat.close()    
handle3.close()
os.remove(newfile)
os.remove(nowfasta)

## MEA / RNAstructure MaxExpect ##
handle4 = open('result4.txt', 'w')
newfile2 = "ME_1.ct"
nowfasta2 = 'now2_1.fasta'
grabthat2 = open('now2_1.fasta', 'r')
check_output(['MaxExpect', filename, newfile2, '--sequence'], universal_newlines = True)
check_output(['ct2dot', newfile2, '1', nowfasta2], universal_newlines = True) 
print('MEdone')
for i in grabthat2:
    handle4.write(i)
grabthat2.close()
handle4.close()
os.remove(newfile2)
os.remove(nowfasta2)

## Pseudoknot predictors ##
'''
## IpKnots ##
handle5 = open('result5.txt', 'w')
mlist5 = []
for x in makego:
    result5 = str(check_output(['ipknot', filename], universal_newlines = True))
    print('Ipknots done')
    mlist5 = result5
for i in mlist5:    
    handle5.write(i)
handle5.close()

## pKiss ##
handle6 = open('result6.txt', 'w')
result6 = (check_output(['./../programs/pseudoknots/fold-grammars/Misc/Applications/pKiss/x86_64-redhat-linux/pKiss_mfe', '-f', '1.seq'], universal_newlines=True))
print('pKiss_done')
handle6.write(result6)
handle6.close()

##clean pKiss output file##
clean_r6 = open('result6.txt', 'r')
r6_1 = open('result6_1.txt', 'w')
r6_1_catch = []
start_capture = False
capture_count = 0 
for i in clean_r6:
    for cc in i:
        r6_1_catch.append(cc)
        if cc == ',' :
            start_capture = True
        elif cc == ' ' and start_capture == True:
            capture_count +=1
            if capture_count == 1:
                 r6_1_catch = []
            elif capture_count == 2:
                 r6_1.write(str(''.join(r6_1_catch))+'\n')
                 r6_1_catch = []
            elif capture_count > 2:
                 continue
os.remove('result6.txt')
r6_1.close()
'''
## Merge the rest of the files ##
### COMBINE THE OTHER FILES HERE ###

## MFE ##
reopenMFE = open('result3.txt', 'r')
for lne in reopenMFE:
    if lne.startswith(char2):
        merge.write(str(1) + lne +'\n')
reopenMFE.close()

## ME ##
reopenME = open('result4.txt', 'r')
countMEpred = 0 
for lne in reopenME:
    if lne.startswith(char2):
        merge.write(str(1) + lne +'\n')
reopenME.close()
'''
## IpKnots ##
reopenIPK = open('result5.txt', 'r')
for lne 'in reopenIPK:
    if lne.startswith(char2):
        merge.write(str(1) + lne +'\n')
reopenIPK.close()

## pKiss ##
reopenHK = open('result6_1.txt', 'r')
for lne in reopenHK:
    if lne.startswith(char2):
        merge.write(str(1) + lne + '\n')
reopenHK.close()
'''
merge.close()

form_matrix=open('form_matrix', 'w')
f = open('merge.txt', 'r')
pred_list = []
pl_count = 0
char_1 = '1A', '1U', '1C', '1G', '1T', '1a', '1u', '1c', '1g', '1t', '1N', '1n'
char2_1 = '1.', '1(', '1[', '1<', '1{'


for line in f:
    if line.startswith('1'):
        form_matrix.write(line.split(' ')[0].strip('\n') + '\n')
    if line.startswith('1>'):
        title = line
    elif line.startswith(char_1):
        seq = line
    elif line.startswith(char2_1):
        pred_list.append(line.split(' ')[0])
        pl_count += 1
f.close()
os.remove('merge.txt')
os.remove('1.seq')
os.remove('1.fasta')
actual = open('actual.txt', 'r')
real_stru = []
for lin in actual:
        if lin.startswith('R'):
            real_stru.append(lin.split('R')[1])
actual.close()

form_matrix.close()

pred_t = open('time_to_predict.txt', 'w')
pred_t.write("%s" % (time.time() - start_time))
pred_t.close
