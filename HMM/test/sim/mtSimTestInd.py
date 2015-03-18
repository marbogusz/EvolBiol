#!/usr/bin/env python2

import random
import dendropy
from dendropy import Tree, treesim
import sys
import os
import subprocess
import shutil
import re
import math
from threading import Thread
from multiprocessing import Process, Manager

########################################################################
#  Script version 18.03.2015 - final
#######################################################################


class HmmDistanceGenerator:
    def __init__(self,numberTaxa, sequenceLength, replicatesNo, modelName, ttpe, itype):
        #sequence length for indelible
        self.seq_len = int(sequenceLength)
        #number of replicates for indelible
        self.replicates = int(replicatesNo)
        #model name 
        self.model = modelName
        
        self.taxaNo = numberTaxa
        
        self.treegen = TreeGenerator()

        self.idist = itype
        
        self.t_c = 'chained'
        self.t_r = 'random'
        self.t_b = 'balanced'
        self.treeType = self.t_r

        if ttpe == 'C':
            self.treeType = self.t_c
        if ttpe == 'B':
            self.treeType = self.t_b

        self.cores = 4;
        self.window = 4; 
        self.wins = []

        self.sim_distances = [0.25, 0.5, 1, 2, 4, 8, 16]

        for d in self.sim_distances:
            i =  self.replicates - 1 
            while i >= (self.window-1):
                j = i-self.window+1
                self.wins.append((d,i-self.window+1, i))
                i -= self.window
            if i != -1:
                self.wins.append((d,0,i))

        self.man = Manager()
        self.distances = self.man.list(self.wins)

        self.indelible_binary = 'indelible'
        self.hmm_binary = 'HMMtreeRC1'
        self.hmm_base_params = ['--lD', '-F','--in']
        self.hmm_misc_params = ['-b', '1', '-o', '0', '--bf', '20']
        self.hmm_alpha_params = ['--initAlpha']
        self.hmm_alpha_est_params = ['--estimateAlpha', '1', '--rateCat', '5']
        self.file_prefix = 'control'
        self.hmmtreefile = '.hmm.tree'

        self.hmmExtTrees = []
        self.hmmExtBins = [] 
        
        self.original_treefile = 'tree.sim'


# TODO - a config file would be cool 
        self.file_suffix = '.indelible'
        self.gtr_suffix = 'GTR'
        self.hky_suffix = 'HKY85'
        self.lg_suffix = 'LG'
        self.model_suffix = self.model
        self.modelno = '7';
        self.indelible_output = 'idlbl';
        self.mafft_exec = 'mafft-linsi'
        self.muscle_exec = 'muscle'
        self.prank_exec = 'prank'
        self.clustal_exec = 'clustal'
        self.raxml_bin = 'raxmlHPC'
        self.raxml_prefix = 'RAxML_bestTree.'
        self.raxml_GTR_params = ['-m', 'GTRGAMMA', '-p', '12345']
        self.raxml_LG_params = ['-m', 'PROTGAMMALG', '-p', '12345']
        self.raxml_HKY_params = ['-m', 'GTRGAMMA', '-p', '12345']

        self.raxml_params = self.raxml_GTR_params
        if (self.model == 'LG'):
            self.raxml_params = self.raxml_LG_params

        self.raxml_params.append('-s')

        self.logfile = open('sim_' + str(self.taxaNo) + '_taxa_' + self.model_suffix + '_' + str(self.seq_len) + '_'  + str(self.replicates) + '_' + self.treeType + '_' + self.idist +  '.log', 'w')
        self.gamma = 'gamma'

        print("HMM analysis for {} discrete distanceswith {} replicates.".format(len(self.sim_distances) ,self.replicates))

        if(self.idist == 'L'):
            print('Using low indel rates')
        elif(self.idist == 'M'):
            print('Using medium indel rates')
        else:
            print('Using high indel rates')

        
    def run(self):
        self.simulate(self.replicates,self.model);
        self.calculate_v2(self.replicates,self.model);
        self.analyseOutput(self.replicates,self.model);
        
    def simulate(self,r,modelname):
    
        ifile = open(self.file_prefix+modelname+self.gamma+self.file_suffix, 'r')
        ictl_template = ifile.read()
        ifile.close()
        index = 0;

    
        for th in self.sim_distances:
            
            treeHeight = round(th,2)
        
            current_dir = str(self.taxaNo) + '_taxa_' + self.model_suffix + '_' + str(self.seq_len) + '_' + str(treeHeight) + '_' + self.treeType +  '_' + self.idist +'_' + str(r) 
            if os.path.exists(current_dir):
                print('Simulation direcory ' + current_dir + ' exists. Skipping\n') 
                continue

            os.mkdir(current_dir)

            for rpl in range(r):
                
                if self.treeType == self.t_r:
                    tree = self.treegen.getTreeByHeight(self.taxaNo,treeHeight)
                elif self.treeType == self.t_b:
                    tree = self.treegen.getBalancedTreeByHeight(self.taxaNo,treeHeight)
                else:
                    tree = self.treegen.getChainedTreeByHeight(self.taxaNo,treeHeight)

                ofile = open('control.txt', 'w');

                p1 = re.compile('(\(.*\)).*')
                m1 = p1.search(tree.as_newick_string() )
            
                newicktext = m1.group(1) 
                outidl =  self.indelible_output + '_' + str(rpl+1)

                alpha = self.getAlpha()
                epsilon = self.getEpsilon()

                if(self.idist == 'L'):
                    lmbda = self.getLambdaUL()
                elif(self.idist == 'M'):
                    lmbda = self.getLambdaUM()
                else:
                    lmbda = self.getLambdaUH()

        
                if (modelname == 'LG'):
                    templateDict = {'alpha' : alpha, 'epsilon' : epsilon, 'lambda' : lmbda, 'newick' : newicktext, 'output' : outidl , 'length' : self.seq_len, 'replicates' : 1}
                elif (modelname == 'GTR'):
                    rates = self.getRevRates()
                    pis = self.getNucleotideFrequencies()
                    templateDict = {'a' : rates[0], 'b' : rates[1], 'c' : rates[2],'d' : rates[3], 'e' : rates[4], 'pi1': pis[0], 'pi2': pis[1], 'pi3': pis[2], 'pi4': pis[3], 'alpha' : alpha, 'epsilon' : epsilon, 'lambda' : lmbda, 'newick' : newicktext, 'output' : outidl, 'length' : self.seq_len, 'replicates' : 1}


                ofile.write(ictl_template.format(**templateDict))
                ofile.close()
                shutil.copy('control.txt',current_dir) 
                os.chdir(current_dir)
                shutil.copy('control.txt','control_'+str(rpl+1)+'.txt') 
                tfile = open(self.original_treefile+'_'+str(rpl+1),'w')
                tfile.write(newicktext)
                #execute indelible
                subprocess.call('indelible',stdout=self.logfile,stderr=self.logfile)
                os.chdir('..')

    def calculate_v2(self,r,model):

        threads = []
        for j in range(self.cores):
            t = Process(target=self.calcThread_v2, args=(j,r,model,))
            threads.append(t)
            t.start()
            
        for th in threads:
            th.join()
        self.logfile.close()

    def calcThread_v2(self, idx,  r, model):
        print('*******Calculation thread ' + str(idx) + ' starts'); 
        while True:
            try:
                th = self.distances.pop()
            except:
                break
        
            os.chdir(sys.path[0])
            treeHeight = round(th[0],2)
            rnge = (th[1],th[2])
            print('####Thread ' + str(idx) +   " calculation of tree height {}".format(treeHeight))
            current_dir = str(self.taxaNo) + '_taxa_' + self.model_suffix + '_' + str(self.seq_len) + '_' + str(treeHeight) + '_' + self.treeType +  '_' + self.idist +'_' + str(r) 
            os.chdir(current_dir)
            print('--Thread ' + str(idx) + ' tree ' + str(treeHeight) + ' HMM')  
            self.runHMMbatch(rnge)
            print('--Thread ' + str(idx) + ' tree ' + str(treeHeight) + ' Alignment')  
            self.alignBatch(rnge)
            print('--Thread ' + str(idx) + ' tree ' + str(treeHeight) + ' RAxML')  
            self.runRaxml(rnge)
        print('******Thread ' + str(idx) + ' is done')
    
    def runRaxml(self, rnge):
	for i in range(rnge[0], rnge[1]+1):
            print('Rax repl ' + str(i+1))
	    true_fc = self.indelible_output + '_' + str(i+1) +  '_TRUE_1.fas'
            mafft_fc = 'mafft_' + str(i+1) + '.fas'
            muscle_fc = 'muscle_' + str(i+1) + '.fas'
            clustal_fc = 'clustal_' + str(i+1) + '.fas'
            prank_fc = 'prank_' + str(i+1) + '.fas.best.fas'

            params_true = self.raxml_params[:]
            params_mafft = self.raxml_params[:]
            params_muscle = self.raxml_params[:]
            params_clustal = self.raxml_params[:]
            params_prank = self.raxml_params[:]

            params_true += [true_fc,'-n', 'true'+str(i+1)]
            params_mafft += [mafft_fc,'-n', 'mafft'+str(i+1)]
            params_muscle += [muscle_fc,'-n', 'muscle'+str(i+1)]
            params_clustal += [clustal_fc,'-n', 'clustal'+str(i+1)]
            params_prank += [prank_fc,'-n', 'prank'+str(i+1)]

            #self.callRaxml(([self.raxml_bin]+params_muscle))
            self.callRaxml(([self.raxml_bin]+params_true))
            #self.callRaxml(([self.raxml_bin]+params_mafft))
            self.callRaxml(([self.raxml_bin]+params_clustal))
            self.callRaxml(([self.raxml_bin]+params_prank))

    def callRaxml(self, params):
        if os.path.isfile(self.raxml_prefix+params[-1]):
            return
        subprocess.call(params,stdout=self.logfile,stderr=self.logfile)
        #print "Raxml", params[-1]

    def runHMMbatch(self, rnge):
        i = rnge[0]
        while i <= rnge[1]:
            clean_name = self.indelible_output + '_' + str(i+1) + '_1' + '.fas'
            i +=1
            print('HMM repl ' + str(i))
            if not os.path.isfile(clean_name + '.hmm.tree'):
                self.callHMM(self.hmm_binary, clean_name,'control_'+str(i)+'.txt')
            for binary in self.hmmExtBins:
                self.callHMM(binary, clean_name,'control_'+str(i)+'.txt')

    def callHMM(self, executable, filename, controlfile):
	#print ("Calling " + executable + "on " + filename) 
        ind = []
        sub = []
        alp = self.parseControl(controlfile, ind, sub)
        params = []
        params.append(executable)
        params += self.hmm_base_params
        params.append(filename)
        #params.append('-i')
        #params += ind
        if (self.model == 'GTR'):
            params.append('--rev')
            #params.append('--param_rev')
            #params += sub
        if (self.model == 'HKY'):
            params.append('--hky')
        if (self.model == 'LG'):
            params.append('--lg')
        #params +=self.hmm_alpha_params
        #params.append(alp)
        params += self.hmm_alpha_est_params
        #print(params)
        subprocess.call(params,stdout=self.logfile) 

      
    def alignBatch(self, rnge):
        for i in range(rnge[0], rnge[1]+1):
            #self.alignMafft(i)
            #self.alignMuscle(i)
            print('algnmt repl ' + str(i))
            self.alignPrank(i)
            self.alignClustal(i)
    
                
    def alignMafft(self, fid):
        curr_file  = self.indelible_output + '_' + str(fid+1) + '_1.fas' 
        mafft_file = 'mafft_' + str(fid+1) + '.fas'
        if os.path.isfile(mafft_file):
            return
        mfd = open(mafft_file,'w')
        subprocess.call([self.mafft_exec, curr_file],stdout=mfd,stderr=self.logfile)
        mfd.close()

    def alignMuscle(self, fid):
        curr_file  = self.indelible_output + '_' + str(fid+1) +'_1.fas' 
        muscle_file = 'muscle_' + str(fid+1) + '.fas'
        if os.path.isfile(muscle_file):
            return
        subprocess.call([self.muscle_exec, '-in', curr_file, '-out', muscle_file],stdout=self.logfile,stderr=self.logfile)


    def alignPrank(self, fid):
        curr_file  = self.indelible_output + '_' + str(fid+1) +'_1.fas' 
        prank_file = 'prank_' + str(fid+1)+ '.fas'
        if os.path.isfile(prank_file+'.best.fas'):
            return
        subprocess.call([self.prank_exec, '-d='+curr_file, '-o='+prank_file],stdout=self.logfile,stderr=self.logfile)
    
    def alignClustal(self, fid):
        curr_file  = self.indelible_output + '_' + str(fid+1) +'_1.fas' 
        clustal_file = 'clustal_' + str(fid+1)+ '.fas'
        if os.path.isfile(clustal_file):
            return
        subprocess.call([self.clustal_exec, '-i', curr_file, '-o', clustal_file],stdout=self.logfile,stderr=self.logfile)

    def analyseOutput(self, r, model):
        filepref = str(self.taxaNo) + '_taxa_' + self.model_suffix + '_' + str(self.seq_len) + '_'  + self.treeType + '_' + self.idist + '_'+ str(r) 
        
        #Robinson Foulds results
        rfname = 'RF_' + filepref + '.txt'
        tdname = 'TD_' + filepref + '.txt'

        resultsRF = open(rfname,'w')
        resultsRF.write('TreeHeight\tDistance\talgorithm\n');
        #Total Tree distance results
        #resultsTD = open(tdname,'w')
        #resultsTD.write('TreeHeight\tRealDistance\tInferredDistance\talgorithm\n');

        for th in self.sim_distances:
            treeHeight = round(th,2)
            print("Analysis step {}".format(treeHeight))
            current_dir = str(self.taxaNo) + '_taxa_' + self.model_suffix + '_' + str(self.seq_len) + '_' + str(treeHeight) + '_' + self.treeType + '_' + self.idist + '_' + str(r) 
            os.chdir(current_dir)
            #create trees based on results
            rmft = []
            rtru = []
            rmus = []
            rclu = []
            rpra = []
            hmmt = []

            

            for i in range(r):
                if not os.path.isfile(self.indelible_output + '_' + str(i+1) + '_1.fas'):
                    print('Failed simulation for file ' + str(i+1))
                    continue
                #print('analysing replicate ' + str(i+1))
                reftree = dendropy.Tree.get_from_stream(open(self.original_treefile+'_'+ str(i+1), 'rU'), "newick", tree_offset=0)
                #true rax
                try:
                    rtru.append(dendropy.Tree.get_from_stream(open(self.raxml_prefix + 'true'+str(i+1), 'rU'), "newick", tree_offset=0))
                    self.writeRF(resultsRF, treeHeight, reftree.symmetric_difference(rtru[-1]), 'True+RAxML')
                except:
                    print('True tree file read error ' + str(i+1))
                #mafft rax
                try:
                    rmft.append(dendropy.Tree.get_from_stream(open(self.raxml_prefix + 'mafft'+str(i+1), 'rU'), "newick", tree_offset=0))
                    self.writeRF(resultsRF, treeHeight, reftree.symmetric_difference(rmft[-1]), 'MAFFT+RAxML')
                except:
                    print('Mafft tree file read error ' + str(i+1))
                #muscle rax
                try:
                    rmus.append(dendropy.Tree.get_from_stream(open(self.raxml_prefix + 'muscle'+str(i+1), 'rU'), "newick", tree_offset=0))
                    self.writeRF(resultsRF, treeHeight, reftree.symmetric_difference(rmus[-1]), 'MUSCLE+RAxML')
                except:
                    print('Muscle tree file read error ' + str(i+1))
                #clustal rax
                try:
                    rclu.append(dendropy.Tree.get_from_stream(open(self.raxml_prefix + 'clustal'+str(i+1), 'rU'), "newick", tree_offset=0))
                    self.writeRF(resultsRF, treeHeight, reftree.symmetric_difference(rclu[-1]), 'Clustal+RAxML')
                except:
                    print('Clustal tree file read error ' + str(i+1))
                #prank rax
                try:
                    rpra.append(dendropy.Tree.get_from_stream(open(self.raxml_prefix + 'prank'+str(i+1), 'rU'), "newick", tree_offset=0))
                    self.writeRF(resultsRF, treeHeight, reftree.symmetric_difference(rpra[-1]), 'Prank+RAxML')
                except:
                    print('Clustal tree file read error ' + str(i+1))
                #hmm
                try:
                    hmmt.append(dendropy.Tree.get_from_stream(open(self.indelible_output + '_' + str(i+1) + '_1.fas' + self.hmmtreefile, 'rU'), "newick", tree_offset=0))
                    self.writeRF(resultsRF, treeHeight, reftree.symmetric_difference(hmmt[-1]), 'alignment-free_HMM')

                    for trf in self.hmmExtTrees:
                        res = dendropy.Tree.get_from_stream(open(self.indelible_output + '_' + str(i+1) + '_1.fas' + trf, 'rU'), "newick", tree_offset=0)
                        self.writeRF(resultsRF, treeHeight, reftree.symmetric_difference(res), 'HMM' + trf)
                except:
                    print('Hmm tree read error: ' +  self.indelible_output + '_' + str(i+1) + '_1.fas')

            os.chdir('..');
        self.logfile.close()

        resultsRF.close();
        #resultsTD.close();


        subprocess.call(['Rscript', 'plotViolin.R',rfname])
        #subprocess.call(['Rscript', 'graph_TDs.R',tdname, pngfile_d, pngfile_d ])
   
########################         Simulation parameter genaration    ###############################
###################################################################################################
    def getAlpha(self):
        return round(random.uniform(0.1,4.0),3)
      
    def getLambdaUL(self):
        ret = round(random.uniform(0.005,0.03333),3)
        while ret  <= 0:
            ret = round(random.uniform(0.005,0.03333),3)
        return ret

    def getLambdaUM(self):
        ret = round(random.uniform(0.03333,0.06666),3)
        while ret  <= 0:
            ret = round(random.uniform(0.03333,0.06666),3)
        return ret
      
    def getLambdaUH(self):
        ret = round(random.uniform(0.06666,0.1),3)
        while ret  <= 0:
            ret = round(random.uniform(0.06666,0.1),3)
        return ret
      
      
    def getLambda(self):
        ret = round(random.gauss(0.03,0.02),3)
        while ret  <= 0:
            ret = round(random.gauss(0.03,0.02),3)
        return ret
      
    def getEpsilon(self):
        return round(random.uniform(0.15,0.75),3)
    
    def getNucleotideFrequencies(self):
        f1=f2=f3=f4=-1.0
        while f1  <= 0:
            f1 = round(random.gauss(0.25,0.05),3)
        while f2  <= 0 or f1+f2 > 1.0:
            f2 = round(random.gauss(0.25,0.05),3)
        while f3  <= 0 or f1+f2+f3 > 0.99:
            f3 = round(random.gauss(0.25,0.05),3)
        f4 = round(1.0 -f1 -f2 -f3,3);
        return [f1,f2,f3,f4]

    def getRevRates(self):
        a=b=c=d=e=f=-1.0;
        while a  <= 0:
            a = round(random.gauss(0.9,0.3),3)
        while b  <= 0:
            b = round(random.gauss(0.9,0.3),3)
        while c  <= 0:
            c = round(random.gauss(0.9,0.3),3)
        while d  <= 0:
            d = round(random.gauss(0.9,0.3),3)
        while e  <= 0:
            e = round(random.gauss(0.9,0.3),3)
        return [a,b,c,d,e]

###################    Parse indelible control file #################################################3
    def parseControl(self, filename, ind, sub):
        ctl = open(filename,'r')
        if ctl:
            contents = ctl.read()
        else:
            print('cannot open ' + filename)

        #get GTR section and numbers
        pSub = re.compile(r'\s+\[submodel\]\s+GTR\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+).*')
        #get inndels
        pLam = re.compile(r'\s+\[insertrate\]\s+([0-9.]+).*')
        pEps = re.compile(r'\s+\[indelmodel\]\s+NB\s+([0-9.]+).*')

        #get alpha 
        pAlp = re.compile(r'\s+\[rates\]\s+0\s+([0-9.]+).*');

        if self.model == 'GTR':
            mg = pSub.search(contents)
            sub.append(mg.group(1))
            sub.append(mg.group(2))
            sub.append(mg.group(3))
            sub.append(mg.group(4))
            sub.append(mg.group(5))

        ml = pLam.search(contents)
        me = pEps.search(contents)
    
        if(ml and me):
            ind.append(ml.group(1))
            ind.append(me.group(1))
        ma = pAlp.search(contents)
        if ma:
            alp = ma.group(1)
        return alp


#######################   Write to analysis file ###################################################    
    def writeRF(self,fd,bd,dist,which):
        fd.write(str(bd)+'\t' +str(dist) + '\t' + which + '\n')

    def writeTD(self, fd,bd,dist1, dist2,which):
        fd.write(str(bd)+'\t' +str(dist1) + '\t' +str(dist2) + '\t' + which + '\n')




####################################################################################################    
####################################################################################################
####################################################################################################
class TreeGenerator:
    def __init__(self):
        self.taxCtr = 0
        self.fullTaxonSet = []
        self.bal_newick =''
        for i in range(256):
            self.fullTaxonSet.append('S'+str(i))
    def get_birth_rate_from_expected_height(self, ntips, expected_height):
        tip_sum = 0
        for i in range(2, ntips + 1):
            tip_sum += float(ntips) / i
        return tip_sum / (expected_height * ntips)

    def rescaleTree(self, tree, factor):
        nds = [nd for nd in tree.postorder_node_iter()]
        for i, n in enumerate(nds):
            elen = n.edge_length;
            if elen < 0.001:
                elen = 0.001
            n.edge_length = round(elen,3) 
            #n.edge_length = random.gauss(n.edge_length*factor, n.edge_length*0.1*factor) 

    def getHeight(self,tree):
        node_prev = None
        node = tree.leaf_nodes()[0]
    
        fd  = node.distance_from_root()

        #print(fd)

        while node.parent_node is not None:
            node = node.parent_node
            node_prev = node
        rd = node_prev.distance_from_root()

        #print(rd)

        return (fd-rd)

    def getBalancedTreeByHeight(self, size, th):
        clades  = int(math.log(size,2)) 
        unit = (th*1.0) / clades
        self.taxCtr = 0;

        self.bal_newick = ''
        self.buildBalancedString(1,clades,unit)
        self.bal_newick += ';'
        #print(self.bal_newick)
        tree = Tree.get_from_string(self.bal_newick,"newick")
        #print(tree)
        tree.deroot()
        return tree

    def buildBalancedString(self, lvl,max_lvl,unit):
        self.bal_newick += '('
        if lvl == max_lvl:
            self.bal_newick += self.fullTaxonSet[self.taxCtr] + ':' + str(round(unit,3))
            self.taxCtr += 1
        else: 
            self.buildBalancedString(lvl+1, max_lvl, unit)
        self.bal_newick += ','
        if lvl == max_lvl:
            self.bal_newick += self.fullTaxonSet[self.taxCtr] + ':' + str(round(unit,3))
            self.taxCtr += 1
        else: 
            self.buildBalancedString(lvl+1, max_lvl, unit)
        self.bal_newick += ')' + ':' + str(round(unit,3))



    def getChainedTreeByHeight(self, size, th):
        unit = (th*1.0) / size
        #print('size: ' + str(size) + ' height ' + str(th) + ' unit len ' + str(unit)) 
        newick = ''
        for i in range(size-1):
            newick += '('
        newick += self.fullTaxonSet[0] + ':' + str(round(unit*2,3)) + ',' + self.fullTaxonSet[1] + ':' + str(round(unit*2,3)) + ')'
        for i in range(2,size):
            newick += str(round(unit,3)) + ',' + self.fullTaxonSet[i] + ':' +  str(round(unit*(i+1),3)) + ')'
        newick += ';'

        tree = Tree.get_from_string(newick,"newick")
        #print(tree)
        tree.deroot()
        return tree

    def getTreeByHeight(self, size, th):

        tree = self.getTree(size,self.get_birth_rate_from_expected_height(size,th))

        while abs(self.getHeight(tree) - th) > (th * 0.05) :
            tree = self.getTree(size,self.get_birth_rate_from_expected_height(size,th))

        tree.deroot()
        return tree

    def getTree(self, size, birthParam):
        tree = treesim.birth_death(birth_rate=birthParam, death_rate=0, taxon_set=dendropy.TaxonSet(self.fullTaxonSet[0:size]))
        #tree.deroot()
        #print(tree)
        #randomize slightly
        self.rescaleTree(tree,1.0)
        return tree



    
    
if __name__ == "__main__":
    if len(sys.argv) != 7:
        print('Usage: ' + sys.argv[0] + ' numberTaxa seqence_length replicate_no model{GTR|HKY|LG} tree_type{B|C|R Indel_dist{L|M|H}')
        raise SystemExit(1)

    runner = HmmDistanceGenerator(int(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    runner.run()




