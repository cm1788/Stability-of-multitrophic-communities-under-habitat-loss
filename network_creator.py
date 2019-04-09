'''
Created on Aug 3, 2011

@author: miguel
'''

import math
from random import Random, randrange, choice, sample
import itertools

from scipy.stats import beta, uniform
import networkx as nx

from configure import SPECIES, CONNECTANCE, HERBIVORY , PRIMARY_PRODS
from configure import PROB_NICHE_MODEL, MIN_NESTEDNESS, MAX_NESTEDNESS, MUTUALISM, MIN_MUTUALISTS
from configure import POSSIBLE_INVADERS, TOP_PREDATORS

from web import Network
import utilities as utls

class NetworkCreator():
    '''
    classdocs
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self.rnd = Random()
        self.rnd_uniform = uniform()
        self.net = Network()
    
#        self.invaders = []
    
    def reset_state(self):
        self.rnd = Random()
        self.rnd_uniform = uniform()
        self.net.clear()
    
    def create_niche_model_network(self, S, C, prob=False):
        '''
        prob = False; whether to create the network following the probabilistic niche model
        Niche model
        
        This is an implementation of the niche model
        '''
        ns = self.rnd_uniform.rvs(size=S+POSSIBLE_INVADERS)
        bet = (1/(2*C)) - 1
        self.beta_dist = beta(1,bet)
        
        producers = int(math.floor(S*PRIMARY_PRODS))
        current_producers = 0
        
        herbivores = int(math.floor(S*HERBIVORY))
        current_herbivores = 0 
        
        predators = int(math.floor(S*TOP_PREDATORS))
        current_predators = 0    
        while self.net.size() == 0 or not nx.is_connected(self.net.to_undirected()) or current_herbivores < herbivores or current_producers < producers or current_predators < predators: #or max(self.net.get_shortest_chain_length().values()) < 3:
            self.net.clear()
#            self.invaders = []
            
            basal = None
            smallest_n = None
            
            #here we obtain the fundamental niche values for each one of the species in the network
            # n, c, r
            for i in xrange(1,S+1):
                self.net.add_node(i)
                self.net.node[i]['n'] = float(ns[i-1]) #self.rnd_uniform.rvs()
                self.net.node[i]['r'] = self.beta_dist.rvs() * self.net.node[i]['n']
                self.net.node[i]['c'] = self.rnd.uniform((self.net.node[i]['r']/2), min(self.net.node[i]['n'], (1-(self.net.node[i]['r']/2))))
                #the original value of c (commented bit) as originally presented in the Nature 2000 paper
                #was changed after reading the niche model specification presented in the JAE 2008 paper
                #the new specification ensures that the r of all the species always lie within the niche interval [0,1]
                #self.net.node[i]['c'] = self.rnd.uniform(self.net.node[i]['r']/2, self.net.node[i]['n'])
                self.net.node[i]['invader'] = False
                
                if smallest_n == None or self.net.node[i]['n'] < smallest_n:
                    smallest_n = self.net.node[i]['n']
                    basal = i
            
            #if INVASION:
#            for i in xrange(S+1, (S+1)+(POSSIBLE_INVADERS)):
##                invader = dict()
##                invader['n'] = self.rnd_uniform.rvs()
##                invader['r'] = self.beta_dist.rvs() * self.net.node[i]['n']
##                invader['c'] = self.rnd.uniform(self.net.node[i]['r']/2, self.net.node[i]['n'])
##                invader['invader'] = True
#                
#                self.net.add_node(i)
#                self.net.node[i]['n'] = float(ns[i-1]) #self.rnd_uniform.rvs()
#                self.net.node[i]['r'] = self.beta_dist.rvs() * self.net.node[i]['n']
#                self.net.node[i]['c'] = self.rnd.uniform((self.net.node[i]['r']/2), min(self.net.node[i]['n'], (1-(self.net.node[i]['r']/2))))
#                self.net.node[i]['invader'] = False
                
#                self.invaders.append(invader)
                          
            self.net.node[basal]['r'] = 0.0    
            self._create_links_based_on_fundamental_niche()
            
            current_herbivores = 0
            current_producers = 0
            current_predators = 0
            tls = self.net.get_trophic_levels()
            for k in tls.keys():
                if tls[k] == 1:
                    current_herbivores += 1 
                elif tls[k] == 0:
                    current_producers += 1
                elif tls[k] == 3:
                    current_predators += 1
            
        return self.net
            
    def _create_links_based_on_fundamental_niche(self, nodes_to_link=None):
        #based on the fundamental niche values obtained above we construct the network by adding
        #the corresponding links according to the species' niche and feeding range
        
        if nodes_to_link == None:
            nodes_to_link = set(self.net.nodes())
        
        for i in self.net.nodes():
            if self.net.node[i]['r'] == 0.0:
                continue
            
            r_lower_bound = self.net.node[i]['c'] - (self.net.node[i]['r']/2)
            r_upper_bound = self.net.node[i]['c'] + (self.net.node[i]['r']/2) 
            for j in self.net.nodes():
                if i not in nodes_to_link and j not in nodes_to_link:
                    continue
                
                if self.net.node[j]['n'] >= r_lower_bound and self.net.node[j]['n'] <= r_upper_bound:  
                    self.net.add_edge(j,i)
        
        #for disconnected or duplicated nodes
        disc_nodes = set()        
        self_loops = set(self.net.selfloop_edges())
        
        for i in nodes_to_link:
            disconnected = False
            
            if self.net.degree(i) == 0:
                disconnected = True
            elif self.net.in_degree(i) == 1 and (i,i) in self_loops:
#                print 'producer with selfloop'
                if self.net.out_degree(i) > 1:
                    self.net.remove_edge(i,i)
                else:
                    disconnected = True
                    attrs = self.net.node[i]
                    self.net.remove_node(i)
                    self.net.add_node(i, attr_dict=attrs)
            else:
                i_succs = set(self.net.successors(i))
                i_predecs = set(self.net.predecessors(i))
                for j in self.net.nodes():
                    j_succs = self.net.successors(j)
                    j_predecs = self.net.predecessors(j)
                    
                    if i_succs == j_succs and i_predecs == j_predecs:
                        disconnected = True
                        attrs = self.net.node[i]
                        self.net.remove_node(i)
                        self.net.add_node(i, attr_dict=attrs)
                        print 'duplicated node'
                        break
                    
            #we reassign the fundamental niche values to nodes that are disconnected or duplicated        
            if disconnected:
                self.net.node[i]['n'] = self.rnd_uniform.rvs()
                self.net.node[i]['r'] = self.beta_dist.rvs() * self.net.node[i]['n']
                self.net.node[i]['c'] = self.rnd.uniform(self.net.node[i]['r']/2, self.net.node[i]['n'])
                disc_nodes.add(i)
                
        
        if len(disc_nodes) > 0:
            self._create_links_based_on_fundamental_niche(disc_nodes)
        
        return
        

def obtain_interactions_network():
    if MIN_NESTEDNESS > 0.0:
        return obtain_interactions_network_with_nestedness()
      
    nc = NetworkCreator()
    net = nc.create_niche_model_network(SPECIES, CONNECTANCE, PROB_NICHE_MODEL) 
    b, producers = net.basal()
    
    herbs = set()
    tls = net.get_trophic_levels()
    for n in net.nodes():
        if tls[n] == 1:
            herbs.add(n)
    
    
    if MUTUALISM != 0.0:
#        failed = True
#        while failed:
#            n_mut_producers = int(math.ceil(len(producers)*MUTUALISM))
#            n_mutualists = int(math.ceil(len(herbs)*MUTUALISM))
#            
#            print 'herbs', len(herbs)
#            print 'muts n', n_mutualists
#            print 'mut prods n', n_mut_producers
#
#            prods_sets = set(itertools.combinations(producers, n_mut_producers))
#            mut_set = set()
#            for s in prods_sets:
#                mut_set.clear()
#                for prod in s:
#                    successors = herbs & set(net.successors(prod))
#                    mut_set |= successors
#                
#                mut_prod_set = set()
#                for m in mut_set:
#                    tmp = set()
#                    tmp.add(m)
#                    preds = set(net.predecessors(m)) - tmp
#                    mut_prod_set |= preds            
#
#                if len(mut_set) == n_mutualists: # and len(mut_prod_set) == n_mut_producers:
#                    failed = False
#                    print 'muts', mut_set
#                    break
#            
#            if failed:
#                nc.reset_state()
#                net = nc.create_niche_model_network(SPECIES, CONNECTANCE, PROB_NICHE_MODEL)
#                
#                b, producers = net.basal()
#                herbs = set()
#                tls = net.get_trophic_levels()
#                for n in net.nodes():
#                    if tls[n] == 1:
#                        herbs.add(n)
        
        print 'about to obtain combinations'
        
        n_mutualists = int(math.ceil(len(herbs)*MUTUALISM))
        #muts_sets = set(itertools.combinations(herbs, n_mutualists))
        
        print 'already obtained combinations'
        
        #mut_set = muts_sets.pop()
        
        mut_set = sample(herbs, n_mutualists)
        
        mut_prod_set = set()
        for m in mut_set:
            tmp = set()
            tmp.add(m)
            preds = set(net.predecessors(m)) - tmp
            mut_prod_set |= preds
        
    else:
        mut_prod_set = set()
        mut_set = set()
    
    ##here, we eliminate cannibalistic loops on herbivore nodes
    self_loops = set(net.selfloop_edges())
    print self_loops
    for h in herbs:
        if (h,h) in self_loops:
            net.remove_edge(h,h)
    
    
    print 'herbs', herbs
    print 'muts', mut_set
    print 'mut_prods', mut_prod_set       
    for n in net.nodes():
        if n in mut_prod_set:
            net.node[n]['mut_prod'] = True
        else:
            net.node[n]['mut_prod'] = False
        
        if n in mut_set:
            net.node[n]['mut'] = True
        else:
            net.node[n]['mut'] = False
            
    return net
    
    
def obtain_interactions_network_with_nestedness():
    nc = NetworkCreator()
    net = nc.create_niche_model_network(SPECIES, CONNECTANCE, PROB_NICHE_MODEL)
    
    #we pick the invaders randomly from each trophic level
    tls = net.get_trophic_levels()
    invaders = 0
    
    trophic_ls = set(tls.values())
    sps_in_tls = dict.fromkeys(trophic_ls, False)
    for sp in tls.keys():
        if sps_in_tls[tls[sp]] == False:
            sps_in_tls[tls[sp]] = [sp]
        else:
            sps_in_tls[tls[sp]].append(sp)
        
    temp_tl = trophic_ls.copy()
    #this is the network that is going to be considered for obtaining the mutualistic subnetwork (without the invaders)
    net_temp = net.copy()
    while invaders < POSSIBLE_INVADERS:
        invader_found = False
        current_tl = temp_tl.pop()
        candidates = sps_in_tls[current_tl]
        
        nc.rnd.shuffle(candidates)
        for inv_id in candidates:
            net_c = net.copy()
            net_c.remove_node(inv_id)
            net_c.remove_edges_from(net_c.selfloop_edges())
            
            if 0 in net_c.degree().values():
                continue
            else:        
                if not net.node[inv_id]['invader']:
                    net.node[inv_id]['invader'] = True
                    invaders += 1
                    invader_found = True
                    net_temp.remove_node(inv_id)
                    print 'invader =', inv_id
                    print 'trophic level:', tls[inv_id]
                    break
        
        if not invader_found:
            print 'could not find a suitable invader from trophic level:', current_tl
                
        if len(temp_tl) == 0:
            temp_tl = trophic_ls.copy()
    
    
    
    if MUTUALISM != 0.0:  
        while True:
            max_nodf = 0.0
            nodf = 0.0
            mutualistic_matrix, prods, herbs = get_mutualistic_matrix(net_temp)
            ordered_matrix, prods, herbs = utls.order_mutualistic_matrix(mutualistic_matrix, prods, herbs)
            number_of_mutualists = int(math.ceil(len(herbs)*MUTUALISM))
            
            if number_of_mutualists < MIN_MUTUALISTS:
                number_of_mutualists = MIN_MUTUALISTS
            
            while len(ordered_matrix) != number_of_mutualists and len(ordered_matrix[0]) > 2:
                ordered_matrix.pop()
                herbs.pop()
                
                rows = len(ordered_matrix)
                cols = len(ordered_matrix[0])
                to_remove = []
                for i in range(cols):
                    sumval = 0
                    for j in range(rows):
                        if ordered_matrix[j][i] == 1:
                            sumval += 1
                            break
                    if sumval == 0:
                        to_remove.append(i)
            
                for i in sorted(to_remove,reverse=True):
                    prods.pop(i)
                    for j in ordered_matrix:
                        j.pop(i)
                               
                ordered_matrix, prods, herbs = utls.order_mutualistic_matrix(ordered_matrix, prods, herbs)
                nodf = utls.calculate_nodf(ordered_matrix)
                if nodf > max_nodf:
                    max_nodf = nodf
    
            #nodf = utls.calculate_nodf(ordered_matrix)
            if nodf >= MIN_NESTEDNESS and nodf <= MAX_NESTEDNESS:
                print 'final nestedness = ', nodf
                print 'max nestedness = ', max_nodf
                break
            else:
                nc.reset_state()
                net = nc.create_niche_model_network(SPECIES, CONNECTANCE, PROB_NICHE_MODEL)
    
                #we pick the invaders randomly from each trophic level
                tls = net.get_trophic_levels()
                invaders = 0
                
                trophic_ls = set(tls.values())
                sps_in_tls = dict.fromkeys(trophic_ls, False)
                for sp in tls.keys():
                    if sps_in_tls[tls[sp]] == False:
                        sps_in_tls[tls[sp]] = [sp]
                    else:
                        sps_in_tls[tls[sp]].append(sp)
                    
                temp_tl = trophic_ls.copy()
                #this is the network that is going to be considered for obtaining the mutualistic subnetwork (without the invaders)
                net_temp = net.copy()
                while invaders < POSSIBLE_INVADERS:
                    invader_found = False
                    current_tl = temp_tl.pop()
                    candidates = sps_in_tls[current_tl]
                    
                    nc.rnd.shuffle(candidates)
                    for inv_id in candidates:
                        net_c = net.copy()
                        net_c.remove_node(inv_id)
                        net_c.remove_edges_from(net_c.selfloop_edges())
                        
                        if 0 in net_c.degree().values():
                            continue
                        else:        
                            if not net.node[inv_id]['invader']:
                                net.node[inv_id]['invader'] = True
                                invaders += 1
                                invader_found = True
                                net_temp.remove_node(inv_id)
                                print 'invader =', inv_id
                                print 'trophic level:', tls[inv_id]
                                break
                    
                    if not invader_found:
                        print 'could not find a suitable invader from trophic level:', current_tl
                            
                    if len(temp_tl) == 0:
                        temp_tl = trophic_ls.copy()
                
    else:
        prods = set()
        herbs = set()
    
    print 'connectance = ', net.connectance()
    print 'fraction of nodes in loops and number of cycles = ', net.fraction_in_loops()        


    ##here, we eliminate cannibalistic loops on herbivore nodes
    herbivores = sps_in_tls[1]
    self_loops = set(net.selfloop_edges())
    for h in herbivores:
        if (h,h) in self_loops:
            net.remove_edge(h,h)

    #here, we include the information about mutualists and their partners into the network as nodes' attributes
    for n in net.nodes():
        if n in prods:
            net.node[n]['mut_prod'] = True
        else:
            net.node[n]['mut_prod'] = False
        
        if n in herbs:
            net.node[n]['mut'] = True
        else:
            net.node[n]['mut'] = False
        
    return net

def get_mutualistic_matrix(net):
    '''
    obtains the mutualistic matrix of herbivores consuming primary producers from net
    where the producers are represented by the columns and the herbivores or mutualists
    by the rows
    '''
    high_predators = set()
    intermediate_predators = set()
    herbivores = set()
    producers = set()

    loops = False
    loop_edges = net.selfloop_edges()
    if len(loop_edges) > 0:
        net.remove_edges_from(loop_edges)
        loops = True
        
    for n in net.nodes():
        if net.in_degree(n) == 0:
            producers.add(n)
    
    print producers  
    
    for n in net.nodes():
        if n not in producers:
            n_predecs = set(net.predecessors(n))
            if n_predecs <= producers:
                herbivores.add(n)
            else:
                if net.out_degree(n) == 0:
                    high_predators.add(n)
                else:
                    intermediate_predators.add(n)
    if loops:
        net.add_edges_from(loop_edges)
    
    prods = []
    for n in producers:
        consumers = set(net.successors(n))
        if len(herbivores & consumers) != 0:
            prods.append(n)
        
    prods.sort()
    herbs = sorted(herbivores)
   
    print herbs
    
    mutualistic_matrix = []
    producers_edges = net.edges(producers)
    
    for i in range(len(herbs)):
        row = []
        for j in range(len(prods)):
            if (prods[j],herbs[i]) in producers_edges:
                row.append(1)
            else:
                row.append(0)
        mutualistic_matrix.append(row)
    
        
    return mutualistic_matrix, prods, herbs      
            