import subprocess
import sys
import os
import shutil
import time
import random
import re
from threading import Thread

class HmmGenerator:
    def __init__(self,sequenceLength, replicatesNo, indelRate, modelName, alpha):
        #sequence length for indelible
        self.seq_len = int(sequenceLength)
        #number of replicates for indelible
        self.replicates = int(replicatesNo)
        #indel rate for indelible
        self.indel_rate  = float(indelRate)
        #model name 
        self.model = modelName
        #indicates how many distance steps to generate
        self.steps = 15
        #distance step 
        self.alpha = alpha
        self.distance_step = 0.1
        self.paml_binary = 'baseml'
        if (self.model == 'LG'):
            self.paml_binary = 'codeml'
        self.indelible_binary = 'indelible'
        self.hmm_binary = 'HMMest'
        self.gamma=''
        if (alpha):
            self.gamma = 'gamma'

        self.paml_binary
        self.hmm_base_params = ['-F','--in']
        self.hmm_gtr_params = ['--param_rev','1.39','0.2','0.22','0.3','0.25']
        self.hmm_hky_params = ['--param_hky','2']
        self.hmm_indel_params = ['-i', str(round(self.indel_rate,2)), '0.5'] 
        self.hmm_misc_params = ['-b', '1', '-o', '0', '--bf', '10']
        self.hmm_alpha_params = ['--rateCat', '5', '--initAlpha', '0.5']
        self.file_prefix = 'control'
        self.file_suffix = '.indelible'
        self.gtr_suffix = 'GTR'
        self.hky_suffix = 'HKY85'
        self.lg_suffix = 'LG'
        self.model_suffix = self.model
        self.modelno = '7';
        self.indelible_output = 'idlbl';
        self.mafft_exec = 'mafft-linsi'

        print("HMM analysis for {} steps with {} replicates. Step size : {}".format(self.steps,self.replicates,self.distance_step))

    def run(self):
        self.simulate(self.steps,self.replicates,self.model);
        self.analyze(self.steps,self.replicates,self.model);



    def transformToClean(self, filename):
        ifl = open(filename, 'r')
        new_name = filename.replace('fas','hmm') 
        ofl = open(new_name ,'w')
        for line in ifl:
            if line.find('>') == -1:
                ofl.write(line.rstrip() + '\n')
        ifl.close()
        ofl.close()
        return new_name
    
    def callHMM(self, executable, filename, descriptor,full):
        #print ("Calling " + executable + "on " + filename) 
        params = []
        params.append(executable)
        params += self.hmm_base_params
        params.append(filename)
        if (not full):
            params += self.hmm_indel_params
        params += self.hmm_misc_params
        if (self.model == 'GTR'):
            params.append('--rev')
            if(not full):
                params += self.hmm_gtr_params
        if (self.model == 'HKY'):
            params.append('--hky')
            if (not full):
                params += self.hmm_hky_params
        if (self.model == 'LG'):
            params.append('--lg')
        if (self.alpha):
            params +=self.hmm_alpha_params
        if (full):
            params.append('--ov')
    
        #print(params)
    
        subprocess.call(params,stdout=descriptor, stderr=subprocess.DEVNULL) 
        #vals = [1,2,3,4,5,6]
        #time.sleep(random.choice(vals))
        #print(filename + ' dummy')
    
    
    def alignBatch(self, count):
        for i in range(count):
            curr_file  = self.indelible_output + '_' + str(i+1) + '.fas' 
            mafft_file = 'mafft_' + str(i+1) + '.fas'
            mfd = open(mafft_file,'w')
            subprocess.call([self.mafft_exec, curr_file],stdout=mfd,stderr=subprocess.DEVNULL)
            mfd.close()
    
    def runHMMbatch(self, count,filed,full):
        threads = []
        for i in range(count):
            clean_name = self.transformToClean(self.indelible_output + '_' + str(i+1) + '.fas')
            t = Thread(target=self.callHMM, args=(self.hmm_binary, clean_name,filed,full,))
            threads.append(t)
            t.start()
        for th in threads:
            th.join()
    
    def runPaml(self, count, paml_template):
        for i in range(count):
            paml_true_fc = paml_template.format(seqfile=self.indelible_output + '_TRUE_' + str(i+1) + '.fas', outfile='paml_true_' + str(i+1), modelno=self.modelno);
            paml_mafft_fc = paml_template.format(seqfile='mafft_' + str(i+1) + '.fas', outfile='paml_mafft_' + str(i+1), modelno=self.modelno);
            paml_viterbi_fc = paml_template.format(seqfile='viterbi_' + self.indelible_output + '_' + str(i+1) + '.fas', outfile='paml_viterbi_' + str(i+1), modelno=self.modelno);
            pt = open(self.paml_binary + '_true.ctl','w')
            pm = open(self.paml_binary + '_mafft.ctl','w')
            pv = open(self.paml_binary + '_viterbi.ctl','w')
            pt.write(paml_true_fc);
            pm.write(paml_mafft_fc);
            pv.write(paml_viterbi_fc);
            pt.close()
            pm.close()
            pv.close()
            subprocess.call([self.paml_binary, self.paml_binary + '_true.ctl'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
            subprocess.call([self.paml_binary, self.paml_binary + '_mafft.ctl'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
            subprocess.call([self.paml_binary, self.paml_binary + '_viterbi.ctl'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    
    def simulate(self, s,r,modelname):
    
        ifile = open(self.file_prefix+modelname+self.gamma+self.file_suffix, 'r')
        ictl_template = ifile.read()
        ifile.close()
    
        current_dist = self.distance_step
    
        while s > 0:
            print ("Simulation step {}".format(round(current_dist,1)))
            ofile = open('control.txt', 'w');
            ofile.write(ictl_template.format(distance=round(current_dist/2,2), rate=round(self.indel_rate,2), output=self.indelible_output, length=self.seq_len, replicates=r))
            ofile.close()
            #mkdir GTR + indelible + distance 
            current_dir = str(self.seq_len) + '_' + str(round(self.indel_rate,2)) + '_Indelible_' + str(r) + '_' + self.model_suffix + '_' + str(round(current_dist,1)) + self.gamma
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
            current_dist += self.distance_step
            s -= 1
    
    
    def analyze(self, s,r,model):
        pfl = open(self.paml_binary + self.gamma + '.paml', 'r');
        outfile_all = open(self.model + '_' + 'hmm_' + str(self.seq_len) + '_' + str(round(self.indel_rate,2))+ '_' + str(r) +'all' + self.gamma,'w',1)
        outfile_ltd = open(self.model + '_' + 'hmm_' + str(self.seq_len) + '_' + str(round(self.indel_rate,2))+ '_' + str(r) +'ltd' + self.gamma,'w',1)
        pml_template = pfl.read()
        pfl.close()
    
        current_dist = self.distance_step
    
        while s > 0:
            print("Calculation step {}".format(round(current_dist,1)))
            outfile_all.write("Analysis step {} \n".format(round(current_dist,1)))
            outfile_ltd.write("Analysis step {} \n".format(round(current_dist,1)))
            outfile_all.flush()
            outfile_ltd.flush()
            current_dir =str(self.seq_len) + '_' + str(round(self.indel_rate,2)) + '_Indelible_' + str(r) + '_' + self.model_suffix + '_' + str(round(current_dist,1)) + self.gamma 
            shutil.copy('pair.trees',current_dir) 
            os.chdir(current_dir)
            #execute indelible
            self.runHMMbatch(r,outfile_all,True)
            self.runMMbatch(r,outfile_ltd,False)
            self.alignBatch(r)
            self.runPaml(r,pml_template)
            #go back
            os.chdir('..');
            current_dist += self.distance_step
            s -= 1
        outfile_all.close()
        outfile_ltd.close()
    

class HmmAnalyzer:

    def __init__(self,sequenceLength, replicatesNo, indelRate, modelName,alpha):
        self. true_ext = 'true'
        self.mafft_ext = 'mafft'
        self.rates_match = 'Rate parameters:'
        self.tree_match = 'tree length ='

        self.gamma =''
        if (alpha):
            self.gamma='gamma'
        self.seq_len = int(sequenceLength)
        self.rate = round(float(indelRate),2)
        self.stype = modelName
        self.replicates = replicatesNo

        self.outfilename = modelName + '_out' + str(self.seq_len) + '_' + indelRate + '_' + replicatesNo  + self.gamma + '.txt'

        self.outfile = open(self.outfilename, 'w')
        self.hmmfile_all = open(modelName + '_hmm_'+  str(self.seq_len) + '_' + indelRate +'_' + replicatesNo + 'all' + self.gamma,'r')
        self.hmmfile_ltd = open(modelName + '_hmm_'+  str(self.seq_len) + '_' + indelRate +'_' + replicatesNo + 'ltd' + self.gamma,'r')
        
        self.idict = {'HKY':1,'GTR':5,'LG':0,'AA':0}

    
    
    def analysePaml(self, dist):
        p = re.compile(r'(\d+)')
        mafft_res = []
        true_res = []
        viterbi_res = []
        #search for all the paml output files
        for pamlentry in os.listdir():
            if pamlentry.startswith('paml'):
                pamlfile = open(pamlentry,'r')
                contents = pamlfile.read()
                pamlfile.close()
                values = self.extractValues(contents);
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
                    values.append('True')
                    true_res.append(values)
                if pamlentry.find('viterbi') != -1:
                    #print ('true')
                    values.append('viterbi')
                    viterbi_res.append(values)
                #print(values)
                #for val in values:
                self.outfile.write(values[0] + '\t' + values[1] + '\t' +values[2])
                self.outfile.write('\n')
                
    def extractValues(self, strc):
        p1 = re.compile(r'tree length =\s+([0-9.]+)')
    
        m1 =  p1.search(strc)
        retlist = []
        if m1:
            retlist.append(m1.group(1))
        #FIXME - so far, only return divergence times
            
        #if (stype == 'GTR' or stype == 'HKY'):
        #    p2 = re.compile(r'Rate parameters:.*')
        #    p3 = re.compile(r'([0-9.]+)')
        #    m2 =  p2.search(strc)
        #    if m2:
        #        for gr in p3.findall(m2.group()):
        #            retlist.append(gr)
        return retlist
    
    def process_hmmout(self, filename,desc):
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
                self.outfile.write(matches[self.idict[self.stype]] + '\t' + curr_dist + '\t' + desc +'\n')
    
    def run(self):
    
        self.outfile.write('infer\treal\ttype\n')

        for entry in os.listdir():
            if os.path.isdir(entry):
                p = re.compile(r'(\d+)_([0-9.]+).*' + self.stype +'_([0-9.]+)' + self.gamma)
                m = p.search(entry)
                if m:
                    length = m.group(1)
                    if int(length) != self.seq_len:
                        continue
                    indel_r = m.group(2)
                    if self.rate != round(float(indel_r),2):
                        continue
                    distance = m.group(3)
                    print('Processing : ' + entry)
                    os.chdir(entry)
                    self.analysePaml(distance)
                    os.chdir('..')
        
    
        self.process_hmmout(self.hmmfile_all,'hmmAll')
        self.process_hmmout(self.hmmfile_ltd,'hmmLtd')
        
        self.outfile.close()
        self.hmmfile_all.close()
        self.hmmfile_ltd.close()
    
        pngfile_m = self.outfilename.replace('.txt','_m.png'); 
        pngfile_s = self.outfilename.replace('.txt','_s.png'); 

        subprocess.call(['Rscript', 'graph_m.R',self.outfilename, pngfile_m, pngfile_m ],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        subprocess.call(['Rscript', 'graph_s.R',self.outfilename, pngfile_s, pngfile_s ],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

