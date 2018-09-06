import os
import subprocess, shutil
from itertools import chain
actual = open('actual.txt', 'r')
a = []
b = []
alist = ''
i = 0
intcount =0
ccd_per = []
bases = 'G', 'g', 'A', 'a', 'C', 'c', 'U', 'u', 'T', 't', 'N', 'n'
case = '.', '(', ')', '[', ']', '{', '}', '<'
inters = '1', '2', '3', '4', '5', '6', '7', '8', '9', '0'
model_progs = []
minereal = open('mine_real.txt', 'w')
for line in actual:
    print(line)
    i += 1
    if line.startswith('R'):
        a = line.replace('R', '')
        alist = alist.strip('\n')+a
        line_ln = len(line)-1 
        score = line_ln
        print(a.strip('\n'))        
    elif line.startswith(bases):
        seq = line
    elif line.startswith(case):
        b.append(line)
    elif line.startswith(inters):
        intcount +=1
        ccd_per.append((str(intcount), line.strip('\n')))
    elif line.startswith(">") and not line.startswith('>new'):
        seqname = line
    elif line.startswith('>new'):
        model_progs.append(line)
print(model_progs)
print( ccd_per)
newb = []
newlinechar_count = 0
do_zero = True
for l in b:
    if ',' in l:
        if do_zero == True:
            newb.append(l.split(',')[0])
        do_zero = False
        newlinechar_count +=1
        newb.append(l.split(',')[newlinechar_count])
       
print('newb', newb)

minereal.write('>actual')
minereal.write('\n')
minereal.write(seq)
minereal.write(alist)
minereal.close()
subprocess.call(['dot2ct', 'mine_real.txt', 'mine_act.txt'])

count = 0
masterscore = open('masterscore.txt', 'a')

if newb == []:
    for l in b:
        #print(l)
        count +=1
        mine = open('mine_%i' %(count), 'w')
        mine.write('>prediction_%i' %(count))
        mine.write('\n')
        mine.write(seq)
        mine.write(l)
        mine.close()
        minename = 'mine_%i' %(count)
        mine2name = 'mine2_%i' %(count)
        subprocess.call(['dot2ct', minename, mine2name])
        scoreout = 'score_%i' %(count)
        subprocess.call(['scorer', mine2name, 'mine_act.txt', scoreout])
        #subprocess.call(['./programs/RNAstructure/exe/scorer', mine2name, 'mine_act.txt', scoreout])
        masterscore.write(seqname)
        masterscore.write(model_progs[count-1])
        for fn, sn in ccd_per:
            if str(count) == fn:
                masterscore.write(sn+'\n')
        shutil.copyfileobj(open(scoreout), masterscore)
        os.remove(mine2name)
        os.remove(minename)
        os.remove(scoreout)
 
elif newb != []:
    for l in newb:
        l = l.strip(' ')
        l = l.strip("'")
        count +=1
        mine = open('mine_%i' %(count), 'w')
        mine.write('>model_%i' %(count))
        mine.write('\n')
        mine.write(seq)
        mine.write(l)
        mine.close()
        minename = 'mine_%i' %(count)
        mine2name = 'mine2_%i' %(count)
        subprocess.call(['dot2ct', minename, mine2name])
        scoreout = 'score_%i' %(count)
        subprocess.call(['scorer', mine2name, 'mine_act.txt', scoreout])
        #subprocess.call(['./programs/RNAstructure/exe/scorer', mine2name, 'mine_act.txt', scoreout])
        masterscore.write(seqname)
        masterscore.write(model_progs[count])
        shutil.copyfileobj(open(scoreout), masterscore)
        
        os.remove(mine2name)
        os.remove(minename)
        os.remove(scoreout)
actual.close()

#os.remove('actual.txt')
os.remove('mine_real.txt')
os.remove('mine_act.txt')
masterscore.close()
