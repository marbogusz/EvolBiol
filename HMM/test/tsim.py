#!/usr/bin/env python2

import random
import dendropy
from dendropy import treesim
import sys
import os
import subprocess
import shutil
import re

#print(tree.as_newick_string())
#print(tree.as_string('newick'))

class TreeGenerator:
    def __init__(self):
        self.fullTaxonSet = []
        for i in range(100):
            self.fullTaxonSet.append('S'+str(i))

    def rescaleTree(self, tree, factor):
        nds = [nd for nd in tree.postorder_node_iter()]
        for i, n in enumerate(nds):
            n.edge_length = random.gauss(n.edge_length*factor, n.edge_length*0.1*factor) 

    def getTree(self, size, birthParam):
        tree = treesim.birth_death(birth_rate=birthParam, death_rate=0, taxon_set=dendropy.TaxonSet(self.fullTaxonSet[0:size]))
        tree.deroot()
        #randomize slightly
        self.rescaleTree(tree,1.0)
        return tree

class HmmDistanceGenerator:
    def __init__(self,numberTaxa, sequenceLength, replicatesNo, modelName):
        #sequence length for indelible
        self.seq_len = int(sequenceLength)
        #number of replicates for indelible
        self.replicates = int(replicatesNo)
        #model name 
        self.model = modelName
        
        self.taxaNo = numberTaxa
        
        self.treegen = TreeGenerator()
        
        #indicates how many birth rate steps to generate
        self.steps = 15
        self.paml_binary = 'baseml'
        if (self.model == 'LG'):
            self.paml_binary = 'codeml'
        self.indelible_binary = 'indelible'
        self.hmm_binary = 'HMMestF1'
        self.hmm_base_params = ['-F','--in']
        self.hmm_misc_params = ['-b', '1', '-o', '0', '--bf', '20']
        self.hmm_alpha_params = ['--rateCat', '5', '--initAlpha', '0.75']
        self.hmm_alpha_est_params = ['--estimateAlpha', '1']
        self.file_prefix = 'control'
        self.file_suffix = '.indelible'
        self.gtr_suffix = 'GTR'
        self.hky_suffix = 'HKY85'
        self.lg_suffix = 'LG'
        self.model_suffix = self.model
        self.modelno = '7';
        self.indelible_output = 'idlbl';
        self.mafft_exec = 'mafft-linsi'
        self.phyml_exec = 'phyml'
        self.prank_exec = 'prank'
        self.muscle_exec = 'muscle'

        self.logfile = open('debug_' + str(self.taxaNo) + '.log', 'w')
        
        self.gamma = 'gamma'

        print("HMM analysis for {} steps with {} replicates.".format(self.steps,self.replicates))

    def getAlpha(self):
        return random.uniform(0.1,4.0)
      
    def getLambda(self):
        return random.gauss(0.03,0.02)
      
    def getEpsilon(self):
        return random.uniform(0.25,0.75)

    def alignMafft(self, fid):
        curr_file  = self.indelible_output + '_' + str(fid+1) + '.fas' 
        mafft_file = 'mafft_' + str(fid+1) + '.fas'
        mfd = open(mafft_file,'w')
        subprocess.call([self.mafft_exec, curr_file],stdout=mfd,stderr=self.logfile)
        mfd.close()

    def alignMuscle(self, fid):
        curr_file  = self.indelible_output + '_' + str(fid+1) + '.fas' 
        muscle_file = 'muscle_' + str(fid+1) + '.fas'
        subprocess.call([self.muscle_exec, '-in', curr_file, '-out', muscle_file],stdout=self.logfile,stderr=self.logfile)


    def alignPrank(self, fid):
        curr_file  = self.indelible_output + '_' + str(fid+1) + '.fas' 
        prank_file = 'prank_' + str(fid+1)
        subprocess.call([self.prank_exec, '-d='+curr_file, '-o='+prank_file],stdout=self.logfile,stderr=self.logfile)
    
    def alignBatch(self, count):
        for i in range(count):
            self.alignMafft(i)
            self.alignMuscle(i)
            #self.alignPrank(i)

    def runHMMbatch(self, count,filed,full):
        threads = []
        for i in range(count):
            clean_name = self.indelible_output + '_' + str(i+1) + '.fas'
            t = Thread(target=self.callHMM, args=(self.hmm_binary, clean_name,filed,full,))
            threads.append(t)
            t.start()
        for th in threads:
            th.join()

    def run(self):
        self.simulate(self.steps,self.replicates,self.model);
        self.analyze(self.steps,self.replicates,self.model);

    def simulate(self, s,r,modelname):
    
        ifile = open(self.file_prefix+modelname+self.gamma+self.file_suffix, 'r')
        ictl_template = ifile.read()
        ifile.close()
    
    
        while s > 0:
        
            birth_rate = 0.1 * s
        
            tree = self.treegen.getTree(self.taxaNo,birth_rate)
            ofile = open('control.txt', 'w');

            p1 = re.compile('(\(.*\)).*')
            m1 = p1.search(tree.as_newick_string() )
            
            newicktext = m1.group(1) 
        
            templateDict = {'alpha' : round(self.getAlpha(),3), 'epsilon' : round(self.getEpsilon(),3), 'lambda' : round(self.getLambda(),3), 'newick' : newicktext, 'output' : self.indelible_output, 'length' : self.seq_len, 'replicates' : r}

            ofile.write(ictl_template.format(**templateDict))
            ofile.close()
            #mkdir GTR + indelible + distance 
            current_dir = str(self.seq_len) + '_' + str(round(birth_rate,1)) + '_Indelible_' + str(r) + '_' + self.model_suffix
            os.mkdir(current_dir)
            shutil.copy('control.txt',current_dir) 
            #shutil.copy('star.trees',current_dir) 
            os.chdir(current_dir)
            #execute indelible
            subprocess.call('indelible',stdout=self.logfile,stderr=self.logfile)
            #alignBatch(replicates)
            #runHMMbatch(replicates)
            #runPaml(replicates,pml_template)
            #go back
            os.chdir('..');
            s -= 1
    
    
    def analyze(self, s,r,model):
      #pfl = open(self.paml_binary + self.gamma + '.paml', 'r');
      #if (not onlyPaml):
      #    outfile_all = open(self.model + '_' + 'hmm_' + str(self.seq_len) + '_' + str(round(self.indel_rate,2))+ '_' + str(r) +'all' + self.gamma,'w',1)
      #    outfile_ltd = open(self.model + '_' + 'hmm_' + str(self.seq_len) + '_' + str(round(self.indel_rate,2))+ '_' + str(r) +'ltd' + self.gamma,'w',1)
      #pml_template = pfl.read()
      #pfl.close()
      
        while s > 0:
            birth_rate = 0.1 * s
            print("Calculation step {}".format(round(birth_rate,1)))
            current_dir = str(self.seq_len) + '_' + str(round(birth_rate,1)) + '_Indelible_' + str(r) + '_' + self.model_suffix
            #shutil.copy('star.trees',current_dir) 
            os.chdir(current_dir)
            #execute indelible
            #self.runHMMbatch(r,outfile_all,True)
            #self.runHMMbatch(r,outfile_ltd,False)
            self.alignBatch(r)
            #self.runPaml(r,pml_template)
            #go back
            os.chdir('..');
            s -= 1
      #if (not onlyPaml):
      #    outfile_all.close()
      #    outfile_ltd.close()

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print('Usage: ' + sys.argv[0] + ' numberTaxa seqence_length replicate_no model')
        raise SystemExit(1)

    runner = HmmDistanceGenerator(int(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4])
    runner.run()




