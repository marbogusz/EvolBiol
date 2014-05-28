import subprocess
import sys
import os
import shutil
import time
import random
from threading import Thread

def transformToClean(filename):
    ifl = open(filename, 'r')
    new_name = filename.replace('fas','hmm') 
    ofl = open(new_name ,'w')
    for line in ifl:
        if line.find('>') == -1:
            ofl.write(line.rstrip() + '\n')
    ifl.close()
    ofl.close()
    return new_name


def callHMM(executable, filename, descriptor):
    #print ("Calling " + executable + "on " + filename) 
    params = []
    params.append(executable)
    params += hmm_base_params
    params.append(filename)
    params += hmm_gtr_params
    params += hmm_indel_params
    params += hmm_misc_params

    #print(params)

    subprocess.call(params,stdout=descriptor, stderr=subprocess.DEVNULL) 
    #vals = [1,2,3,4,5,6]
    #time.sleep(random.choice(vals))
    #print(filename + ' dummy')


def alignBatch(count):
    for i in range(count):
        curr_file  = indelible_output + '_' + str(i+1) + '.fas' 
        mafft_file = 'mafft_' + str(i+1) + '.fas'
        mfd = open(mafft_file,'w')
        subprocess.call([mafft_exec, curr_file],stdout=mfd,stderr=subprocess.DEVNULL)
        mfd.close()

def runHMMbatch(count,filed):
    threads = []
    for i in range(count):
        clean_name = transformToClean(indelible_output + '_' + str(i+1) + '.fas')
        t = Thread(target=callHMM, args=(hmm_binary, clean_name,filed,))
        threads.append(t)
        t.start()
    for th in threads:
        th.join()

def runPaml(count, paml_template):
    for i in range(count):
        paml_true_fc = paml_template.format(seqfile=indelible_output + '_TRUE_' + str(i+1) + '.fas', outfile='paml_true_' + str(i+1), modelno=modelno);
        paml_mafft_fc = paml_template.format(seqfile='mafft_' + str(i+1) + '.fas', outfile='paml_mafft_' + str(i+1), modelno=modelno);
        paml_viterbi_fc = paml_template.format(seqfile='viterbi_' + str(i+1) + '.fas', outfile='paml_viterbi_' + str(i+1), modelno=modelno);
        pt = open('baseml_true.ctl','w')
        pm = open('baseml_mafft.ctl','w')
        pv = open('baseml_viterbi.ctl','w')
        pt.write(paml_true_fc);
        pm.write(paml_mafft_fc);
        pm.write(paml_viterbi_fc);
        pt.close()
        pm.close()
        pv.close()
        subprocess.call(['baseml', 'baseml_true.ctl'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        subprocess.call(['baseml', 'baseml_mafft.ctl'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        subprocess.call(['baseml', 'baseml_viterbi.ctl'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

def simulate(s,r):
    ifile = open(file_prefix+model_suffix+file_suffix, 'r')
    ictl_template = ifile.read()
    ifile.close()

    current_dist = distance_step

    while s > 0:
        print ("Simulation step {}".format(round(current_dist,1)))
        ofile = open('control.txt', 'w');
        ofile.write(ictl_template.format(distance=round(current_dist/2,2), rate=round(indel_rate,2), output=indelible_output, length=seq_len, replicates=r))
        ofile.close()
        #mkdir GTR + indelible + distance 
        current_dir = str(seq_len) + '_' + str(round(indel_rate,2)) + '_Indelible_' + str(r) + '_' + model_suffix + '_' + str(round(current_dist,1))
        os.mkdir(current_dir)
        shutil.copy('control.txt',current_dir) 
        shutil.copy('pair.trees',current_dir) 
        os.chdir(current_dir)
        #execute indelible
        subprocess.call('indelible',stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        #alignBatch(replicates)
        #runHMMbatch(replicates)
        #runPaml(replicates,pml_template)
        #go back
        os.chdir('..');
        current_dist += distance_step
        s -= 1



def analyze(s,r):
    pfl = open('baseml.paml', 'r');
    outfile = open('hmm_' + str(seq_len) + '_' + str(round(indel_rate,2))+ '_' + str(r),'w',1)
    pml_template = pfl.read()
    pfl.close()

    current_dist = distance_step

    while s > 0:
        print("Analysis step {}".format(round(current_dist,1)))
        outfile.write("Analysis step {} \n".format(round(current_dist,1)))
        outfile.flush()
        current_dir =str(seq_len) + '_' + str(round(indel_rate,2)) + '_Indelible_' + str(r) + '_' + model_suffix + '_' + str(round(current_dist,1))
        shutil.copy('pair.trees',current_dir) 
        os.chdir(current_dir)
        #execute indelible
        runHMMbatch(r,outfile)
        alignBatch(r)
        runPaml(r,pml_template)
        #go back
        os.chdir('..');
        current_dist += distance_step
        s -= 1
    outfile.close()

#definitions

if len(sys.argv) != 4:
    print('Usage: script.py seqence_length replicates indel_rate', file=sys.stderr)
    raise SystemExit(1)


#indicates how many distance steps to generate
steps = 15
#sequence length for indelible
seq_len = int(sys.argv[1])
#number of replicates for indelible
replicates = int(sys.argv[2])
#indel rate for indelible
indel_rate  = float(sys.argv[3])
#distance step 
distance_step = 0.1

paml_binary = 'baseml'
indelible_binary = 'indelible'
hmm_binary = 'HMMest'

hmm_base_params = ['-F','--rev','--in']
hmm_gtr_params = ['--param_rev','1.39','0.2','0.22','0.3','0.25']
hmm_indel_params = ['-i', '0.05', '0.5'] 
hmm_misc_params = ['-b', '1', '-o', '0', '--ov']

file_prefix = 'control'
file_suffix = '.indelible'

gtr_suffix = 'GTR'
hky_suffix = 'HKY85'

model_suffix = gtr_suffix
modelno = '7';

indelible_output = 'idlbl';

mafft_exec = 'mafft-linsi'

#Read the template control file for indelilble

#open output control file

#main control

print("HMM analysis for {} steps with {} replicates. Step size : {}".format(steps,replicates,distance_step))

simulate(steps,replicates);
analyze(steps,replicates);

#def processStep(distance):    #phylogenetic distance
#    pass
#
#def generateFile(template, distance):    #return filename
#    pass

