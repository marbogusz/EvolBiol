#!/usr/bin/env python2

import random
import dendropy
from dendropy import treesim

#print(tree.as_newick_string())
#print(tree.as_string('newick'))

class TreeGenerator:
    def __init__(self, maxNumber=20):
        self.fullTaxonSet = []
        for i in range(maxNumber):
            self.fullTaxonSet.append('S'+str(i))

    def rescaleTree(self, tree, factor):
        nds = [nd for nd in tree.postorder_node_iter()]
        for i, n in enumerate(nds):
            n.edge_length = random.gauss(n.edge_length*factor, n.edge_length*0.1*factor) 

    def getTree(self, size, minTotalLength=100):
        tree = treesim.pop_gen_tree(taxon_set=dendropy.TaxonSet(self.fullTaxonSet[0:size]))
        tree.deroot()
        if tree.length() < minTotalLength :
            self.rescaleTree(tree,minTotalLength/tree.length())
        return tree

class HmmDistanceGenerator:
    def __init__(self,sequenceLength, replicatesNo, modelName):
        #sequence length for indelible
        self.seq_len = int(sequenceLength)
        #number of replicates for indelible
        self.replicates = int(replicatesNo)
        #model name 
        self.model = modelName
        
        self.treegen = TreeGenerator()
        
        #indicates how many distance steps to generate
        self.steps = 15
        self.paml_binary = 'baseml'
        if (self.model == 'LG'):
            self.paml_binary = 'codeml'
        self.indelible_binary = 'indelible'
        self.hmm_binary = 'HMMestF1'
        self.hmm_base_params = ['-F','--in']
        self.hmm_indel_params = ['-i', str(round(self.indel_rate,2)), '0.5'] 
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
        
        self.gamma = 'gamma'

        print("HMM analysis for {} steps with {} replicates. Step size : {}".format(self.steps,self.replicates,self.distance_step))

    def getAlpha(self):
	return random.uniform(0.1,4.0)
    def getLambda(self):
	return random.gauss(0.03,0.02)
    def getEpsilon(self)
	return rendom.uniform(0.25,0.75)

    def run(self):
        self.simulate(self.steps,self.replicates,self.model);
        self.analyze(self.steps,self.replicates,self.model,False);

    def simulate(self, s,r,modelname):
	
	    ifile = open(self.file_prefix+modelname+self.gamma+self.file_suffix, 'r')
	    ictl_template = ifile.read()
	    ifile.close()
	
	
	    while s > 0:
		
		tree = self.treegen.getTree(5,5)
		ofile = open('control.txt', 'w');
		
		templateDict = {'alpha' : round(self.getAlpha(),3), 'epsilon' : round(self.getEpsilon(),3), 'lambda' : round(self.getLambda(),3), {newick} : tree.as_newick_string(), 'output' : self.indelible_output, 'length' : self.seq_len, 'replicates' : r}
		newvals = list(self.proportions);
		
		for i in range(len(newvals)):
			newvals[i] = round(newvals[i] * current_dist/20.0,3);
		distancesDict = dict(zip(self.keys, newvals))
		templateDict.update(distancesDict)

		ofile.write(ictl_template.format(**templateDict))
		ofile.close()
		#mkdir GTR + indelible + distance 
		current_dir = str(self.seq_len) + '_' + str(round(self.indel_rate,2)) + '_Indelible_' + str(r) + '_' + self.model_suffix + '_' + str(round(current_dist,1)) + self.gamma
		os.mkdir(current_dir)
		shutil.copy('control.txt',current_dir) 
		shutil.copy('star.trees',current_dir) 
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
    
    
      def analyze(self, s,r,model,onlyPaml):
	  pfl = open(self.paml_binary + self.gamma + '.paml', 'r');
	  if (not onlyPaml):
	      outfile_all = open(self.model + '_' + 'hmm_' + str(self.seq_len) + '_' + str(round(self.indel_rate,2))+ '_' + str(r) +'all' + self.gamma,'w',1)
	      outfile_ltd = open(self.model + '_' + 'hmm_' + str(self.seq_len) + '_' + str(round(self.indel_rate,2))+ '_' + str(r) +'ltd' + self.gamma,'w',1)
	  pml_template = pfl.read()
	  pfl.close()
      
	  current_dist = self.distance_step
      
	  while s > 0:
	      print("Calculation step {}".format(round(current_dist,1)))
	      if (not onlyPaml):
		  outfile_all.write("Analysis step {} \n".format(round(current_dist,1)))
		  outfile_ltd.write("Analysis step {} \n".format(round(current_dist,1)))
		  outfile_all.flush()
		  outfile_ltd.flush()
	      current_dir =str(self.seq_len) + '_' + str(round(self.indel_rate,2)) + '_Indelible_' + str(r) + '_' + self.model_suffix + '_' + str(round(current_dist,1)) + self.gamma 
	      shutil.copy('star.trees',current_dir) 
	      os.chdir(current_dir)
	      #execute indelible
	      if (not onlyPaml):
		  self.runHMMbatch(r,outfile_all,True)
		  self.runHMMbatch(r,outfile_ltd,False)
	      self.alignBatch(r)
	      self.runPaml(r,pml_template)
	      #go back
	      os.chdir('..');
	      current_dist += self.distance_step
	      s -= 1
	  if (not onlyPaml):
	      outfile_all.close()
	      outfile_ltd.close()



