import os 
import sys
import re

#python analyzer 

#read all the directories

#for each directory read paml files
#paml files for mafft and true results
#aggregate in a single file 

def analysePaml(dist):
    p = re.compile(r'(\d+)')
    mafft_res = []
    true_res = []
    #search for all the paml output files
    for pamlentry in os.listdir():
        if pamlentry.startswith('paml'):
            pamlfile = open(pamlentry,'r')
            contents = pamlfile.read()
            pamlfile.close()
            values = extractValues(contents);
            values.append(dist)
            #m = p.search(pamlentry)
            #if m:
            #    values.append(m.group(1))
            if pamlentry.find('mafft') != -1:
                #print ('mafft')
                values.append('mafft')
                mafft_res.append(values)
            if pamlentry.find('true') != -1:
                #print ('true')
                values.append('true')
                true_res.append(values)
            #print(values)
            #for val in values:
            outfile.write(values[0] + '\t' + values[1] + '\t' +values[2])
            outfile.write('\n')
            
def extractValues(strc):
    p1 = re.compile(r'tree length =\s+([0-9.]+)')
    p2 = re.compile(r'Rate parameters:.*')
    p3 = re.compile(r'([0-9.]+)')
    m1 =  p1.search(strc)
    m2 =  p2.search(strc)
    retlist = []
    if m1:
        retlist.append(m1.group(1))
    #if m2:
        #print('matched!')
        #for gr in p3.findall(m2.group()):
            #retlist.append(gr)
    return retlist

def process_hmmout(filename):
    p1 = re.compile('Analysis step ([0-9.]+)')
    p2 = re.compile('([0-9.]+)')
    curr_dist = 'error'
    for line in filename:
        m1 = p1.search(line)
        m2 = p2.search(line)
        if m1:
            curr_dist = m1.group(1)
        elif m2:
            matches = p2.findall(line)
            outfile.write(matches[5] + '\t' + curr_dist + '\thmm\n')

true_ext = 'true'
mafft_ext = 'mafft'
rates_match = 'Rate parameters:'
tree_match = 'tree length ='

if len(sys.argv) != 3:
    print('Usage: script.py seqence_length indel_rate', file=sys.stderr)
    raise SystemExit(1)

seq_len = int(sys.argv[1])
rate = round(float(sys.argv[2]),2)

outfile = open('out' + str(seq_len) + '_' + sys.argv[2] + '.txt', 'w')
hmmfile = open('hmm_'+  str(seq_len) + '_' + sys.argv[2] +'_16','r')

outfile.write('infer\treal\ttype\n')

for entry in os.listdir():
    if os.path.isdir(entry):
        p = re.compile(r'(\d+)_([0-9.]+).*GTR_([0-9.]+)')
        m = p.search(entry)
        if m:
            length = m.group(1)
            if int(length) != seq_len:
                continue
            indel_r = m.group(2)
            if rate != round(float(indel_r),2):
                continue
            distance = m.group(3)
            print(distance)
            print(length)
            print(indel_r)
            os.chdir(entry)
            analysePaml(distance)
            os.chdir('..')

process_hmmout(hmmfile)

outfile.close()
hmmfile.close()
