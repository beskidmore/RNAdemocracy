1  ### Master file for forming consensus models of RNA secondary structure 
2  # predictions from established methods                                 ###
3  ### Author: Benjamin Skidmore                                          ###
4  ### Date: 10/14/2017                                                   ###
5  ### license: GNU GPL 3.0-or-later                                      ###
6  import subprocess, os
7  import time
8  start_time = time.time()
9  
10 char = 'a', 'g', 'c', 'u', 't', 'n', 'A', 'G', 'C', 'U', 'T', 'N'               
11 char2 ='(', '.', '[', ')', ']', '{', '}', '<', '>' 
12 inters = '1', '2', '3', '4', '5', '6', '7', '8', '9', '0'
13 
14 database = open("split1", 'r')
15 
16 #os.remove('time_to_predict.txt')
17 #os.remove('total_time_run.txt')
18 #os.remove('masterscore1')
19 #os.remove('masterscore2')
20 #os.remove('masterscore3')
21 #os.remove('masterscore4')
22 
23 newtitle = False
24 newSEQ = False
25 parsing = False
26 actualTEST = False
27 ender = False
28 parsingSEQ = False
29 parseACT = False
30 NOWRUN = False
31 realLST = []
32 #final_results = open('final_results.txt', 'a')
33 ## Makes input.txt a clean .fasta runable file. Initates the run by calling
34 ##  The working scripts ##
35 for line in database:
36     if line.startswith('# File'):
37         job_file = open('job_file.txt', 'w')
38         actual = open('actual.txt', 'w')
39         title = line.strip('# File ')
40         job_file.write('>' + title)
41         actual.write('>' + title)
42         half_actual = open('actual_half.txt', 'w')
43         half_actual.write('>' + title)
44         #final_results.write('>'+title)
45         #final_results.close()
46     elif line.startswith(char):
47         newSEQ = True
48         parsingSEQ = True
49     elif line.startswith(char2):
50         actualTEST = True
51         parsingSEQ = False
52         parseACT = True
53         actual.write('R' + line)
54         realLST.append('R' + line)
55         half_actual.write('R' + line)
56     elif line.startswith('\n') and parseACT == True:
57         ender = True
58         parseACT = False
59         NOWRUN = True
60         job_file.close()
61     if parsingSEQ == True:
62         line = list(line)
63         ch_count = 0 
64         for ch in line:
65             ch_count += 1
66             if ch not in char and ch is not '\n' :
67                 point = ch_count
68                 print('point', point)
69                 del line[point-1]
70                 line.insert(point-1, 'N')
71         lines = ''.join(line)
72         job_file.write(lines)
73         newtitle = False
74         newSeq = False
75     if NOWRUN == True:
76         start_time2 = time.time()
77         actual.close()
78         ln_count = 0
79         ln_lst = []
80         mm_count = 0
81         ## RUN SCRIPTS HERE ##
82         subprocess.call(['python3.5', 'predict.py', '&'])
83         subprocess.call(['python3.5', 'decide.py', '&'])
84         subprocess.call(['python3.5', 'round_2.py', '&'])
85         
86 		## Check runtime without testing ##
87         run_time = open('time_no_test.txt', 'a')
88         run_time.write("%s" % (time.time() - start_time2)+ '\n')
89         run_time.close()
90 		
91         form_matrix = open('form_matrix', 'r')
92 		#### USED ONLY FOR TESTING ####
93 		## Prepares the individual established predictions for testing ##
94         for ln in form_matrix:
95             if ln.startswith('1.') or ln.startswith('1(') or 
96             ln.startswith('1[') or ln.startswith('1{') or 
97             ln.startswith('1<'):
98                 ln_count +=1
99                 allopenall = ('actual_'+'{0}'.format(ln_count)+'.txt')
100                if allopenall in os.walk('.'):
101                    os.remove(allopenall)
102                nln = ln.strip('1')
103                ln_lst.append(nln)
104            elif ln.startswith('1>'):
105                continue
106            else:
107                realSEQ = ln.strip('1')
108        for mm in ln_lst:
109            mm_count += 1
110            with open(('actual_'+'{0}'.format(mm_count)+'.txt'), 'w') as openall:
111                openall.write('>'+title)
112                for rL in realLST:
113                    openall.write(rL)
114                openall.write(realSEQ)
115                openall.write(mm)
116            openall.close()
117        realLST = []
118        form_matrix.close()
119		## Testing each established methods prediction against the 
120        # reference structure ##             
121        subprocess.call(['python3.5', 'testing_1.py', '&'])
122        subprocess.call(['python3.5', 'testing_2.py', '&'])
123        subprocess.call(['python3.5', 'testing_3.py', '&'])
124        subprocess.call(['python3.5', 'testing_4.py', '&'])
125        #subprocess.call(['python3.5', 'testing_5.py', '&'])
126        #subprocess.call(['python3.5', 'testing_6.py', '&'])
127        ## OUR PROGRAMS FINAL OUTPUT ##
128        final_results = open('final_results.txt', 'r')
129        actual = open('actual.txt', 'a')
130
131        ## Prepare to test against RNAdemocracy results ##
132        for ln in final_results:
133            actual.write(ln)
134            if ln.startswith(char):
135                half_actual.write(ln)
136        actual.close()
137        final_results.close()
138		## test_half is meant to test the predictons before the second
139        # round, so only the consensus structures are tested without the
140        # filled in space haveing any assigne structure ##
141        sec_rr = open('sec_run_ready.txt', 'r')
142        for ln in sec_rr:
143            if ln.startswith(char):
144                continue
145            else:
146                half_actual.write(ln)
147        half_actual.close()
148        sec_rr.close()
149        subprocess.call(['python3.5', 'testing.py', '&'])
150        subprocess.call(['python3.5', 'testing_half.py', '&'])
151        #os.remove('actual.txt')
152        with_testing_t = open('time_with_test.txt', 'a')
153        with_testing_t.write("%s" % (time.time() - start_time2)+ '\n') 
154        with_testing_t.close()
155        NOWRUN = False
156        
157database.close()
158actual.close()
159form_matrix.close()
160final_results.close()
161job_file.close()
162
163total_t = open('time_run_total.txt', 'w')
164total_t.write("%s" % (time.time() - start_time)+ '\n')
165total_t.close
166print("--- %s seconds ---" % (time.time() - start_time))
167
