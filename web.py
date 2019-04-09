
"""
:mod:`~web` presents the implementation of a Network data structure, which adapts
the structure of a complex network for the particular use of this object in the
field of Ecological Networks.
"""

import os
import datetime

import networkx as nx
import math
#import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from operator import itemgetter
import random
import subprocess

import csv

from cluster import *

class Network(nx.DiGraph):
    """
    The Network class implements the necessary methods for representing the network of interaction between species
    in an ecological network. 
    
    It conveniently provides a series of calculations and metrics useful for the analysis and study of this type of
    networks.
    
    In this representation of an ecological network, the nodes possess two common attributes:
    
        =========    ===================================================================================================================================
        Attribute    Description
        =========    ===================================================================================================================================
        biomass      the mean biomass of an individual of the species represented by the node
        group        the major group the species belongs to ('Bird', 'BirdPrey', 'Reptile', 'Mammal', 'MammalCarn', 'Fish', 'Amphibian', 'Invertebrate')
        =========    ===================================================================================================================================
    
    This class inherits some of its behaviour from the parent class DiGraph from the NetworkX library for complex
    networks analyses (http://networkx.lanl.gov/reference/classes.digraph.html) 
    """
    
    added = False
    #this structure keeps the distance of the longest path to each node (after calling the longest_path_length method)
    lentgh_to = None

    def connectance(self):
        """
        Calculates the connectance of the Network, i.e. the fraction of realised links
        of all the possible within it if all the nodes were connected to one another. ::
        
            C = L/S^2
            
        where L is the number of links and S the number of species (or nodes) in the network respectively
        
        A float value of the network's connectance is returned.
        
        .. seealso::
    
            :meth:`~linkage_density`
                For the calculation of mean number of links per node
        """
        if self.size() == 0:
            return 0.0
        
        self.bipartite = True
        nodes_in = set()
        nodes_out = set()
        
        for n in self.nodes():
            in_deg = self.in_degree(n)
            out_deg = self.out_degree(n)
            if in_deg > 0 and out_deg > 0:
                self.bipartite = False
                break
        
            if in_deg > 0:
                nodes_in.add(n)
            if out_deg > 0:
                nodes_out.add(n)
        
        if self.bipartite and nodes_in.isdisjoint(nodes_out):
            self.bipartite = True
        else:
            self.bipartite = False
        
        if self.bipartite:
            c = float(self.size())/float(len(nodes_in)*len(nodes_out))
        else:
            if self.number_of_selfloops() == 0:
                c = (float(self.size())/float(self.number_of_nodes()*(self.number_of_nodes()-1)))
            else:
                c = (float(self.size())/float(self.number_of_nodes()**2))
        return c


    def linkage_density(self):
        """
        This method calculates the linkage density (mean number of links per node) of the network. ::
        
            linkage_density = L/S
        
        where L is the number of links and S the number of species (or nodes) in the network respectively
        
        It returns the float representation of the linkage density in the network.
        """
        if self.size() == 0 or self.number_of_nodes() == 0:
            return 0.0
        else:
            return (float(self.size())/float(self.number_of_nodes()))
    
    def get_nodes_by_species_group(self, group):
        """
        Gets a list of the species (or nodes) in the network belonging to the group ``group``    
    
        *group*
          a string from the set::
              
              'Bird', 'BirdPrey', 'Reptile', 'Mammal', 'MammalCarn', 'Fish', 'Amphibian', 'Invertebrate'
          
        An array with the id of the nodes belonging the to the group ``group`` is returned
        """
        nodes_group = []
        
        nodes = self.nodes(data=True)
        for n in nodes:
            atts = n[1]
            try:
                if atts['group'] == group:
                    nodes_group.append(n[0])
            except:
                pass
        
        return nodes_group
    
    
    def top_predators(self):
        """
        This method calculates the fraction of top predator species (species that are not preyed upon by any other)
        present in the food web. i.e. the fraction of nodes in the network possessing only incoming links.    
    
        It returns a float representing the fraction of this type of nodes and a set containing the top predator
        species/nodes with only incoming connections.
    
        .. seealso::
        
            :meth:`basal` 
                for the calculation of the basal species
            
            :meth:`intermediate` 
                for the calculation of the intermediate species
        
        """
        cannibals = self.nodes_with_selfloops()
        top_species = set()
        preds = 0
        if self.number_of_nodes() == 0:
            return preds, top_species
        
        for n in self.nodes():
            if self.out_degree(n) == 0 or (self.out_degree(n) == 1 and n in cannibals):
                preds += 1
                top_species.add(n)
        
        return float(preds)/float(self.number_of_nodes()), top_species
    
    def basal(self, heterotrophs=False):
        """
        Calculates the fraction of basal species (species that do not prey on any other)
        present in the food web. i.e. the fraction of nodes in the network possessing only outgoing links.    
        
        *heterotrophs = False*
            a boolean value specifying whether to take into account heterotroph species as basal ones. This flag allows
            for the elimination of basal nodes that should not be considered for other network parameters calculations.
    
        It returns a float representing the fraction of basal nodes and a set containing the basal
        species (nodes with only outgoing connections).
    
        .. seealso::
    
            :meth:`top_predators`
                for the calculation of the top predator species
            
            :meth:`intermediate` 
                for the calculation of the intermediate species
        """
        basal_species = set()
        basal = 0
        if self.number_of_nodes() == 0:
            return basal, basal_species
        
        if heterotrophs:
            basals = set()
            for n in self.nodes():
                if self.in_degree(n) == 0:
                    basals.add(n) 
            
            temp_net = self.copy()
            temp_net.remove_nodes_from(basals)
            
            print basals
            
            for n in temp_net.nodes():
                if temp_net.in_degree(n) == 0:
                    basal += 1
                    basal_species.add(n)
            
            basal = float(basal)/float(temp_net.number_of_nodes())
            
        else:
            for n in self.nodes():
                if self.in_degree(n) == 0:
                    basal += 1
                    basal_species.add(n)
            
            basal = float(basal)/float(self.number_of_nodes())
            
        return basal, basal_species
    
    def intermediate(self, heterotrophs=False):
        """
        Calculates the fraction of intermediate species (species that are preyed upon and at the same time prey on other species in the network)
        present in the food web. i.e. the fraction of nodes in the network possessing both incoming and outgoing links.    
        
        *heterotrophs = True*
            a boolean value specifying whether to take into account heterotroph species as basal ones. This flag allows
            for the elimination of basal nodes that should not be considered for other network paramenters calculations.
    
        It returns a float representing the fraction of intermediate nodes and a set containing the intermediate
        species (nodes with both incoming and outgoing connections).
    
        .. seealso::
    
            :meth:`top_predators` 
                for the calculation of the top predator species
            
            :meth:`basal` 
                for the calculation of the basal species
        """
        cannibals = self.nodes_with_selfloops()
        inter_species = set()
        inter = 0
        if self.number_of_nodes() == 0:
            return inter, inter_species
        
        if heterotrophs:
            basals = set()
            for n in self.nodes():
                if self.in_degree(n) == 0:
                    basals.add(n) 
            
            temp_net = self.copy()
            temp_net.remove_nodes_from(basals)
            
            for n in temp_net.nodes(): 
                if temp_net.in_degree(n) > 0 and (temp_net.out_degree(n) > 1 or (temp_net.out_degree(n) == 1 and n not in cannibals)):
                    inter += 1
                    inter_species.add(n)
            
            inter = float(inter)/float(temp_net.number_of_nodes())
            
        else:
            for n in self.nodes():
                if self.in_degree(n) > 0 and (self.out_degree(n) > 1 or (self.out_degree(n) == 1 and n not in cannibals)):
                    inter += 1
                    inter_species.add(n)
            
            inter = float(inter)/float(self.number_of_nodes())
        
        return inter, inter_species
    
    def cannibalism(self):
        """
        This method obtains the fraction of cannibalism in the network.
        
        It returns a float representing the fraction of nodes that possess self loops in the network.
        """
        if self.number_of_nodes() == 0 or self.number_of_selfloops() == 0:
            return 0
        else:
            return float(self.number_of_selfloops())/float(self.number_of_nodes())
    
    def fraction_in_loops(self):
        fraction = 0.0
        number = 0
        if self.size() == 0 or self.number_of_nodes() == 0:
            return fraction, number
        
        has_loops = False
        if self.number_of_selfloops() > 0:
            loops = self.selfloop_edges(data=True)
            self.remove_edges_from(loops)
            has_loops = True
        
	#updated to fix conflict between networkx versions
	# nx1.7 returns cycles as list, nx1.9 returns a generator
	if nx.__version__=='1.7': 
            cycles = nx.algorithms.cycles.simple_cycles(self)
	else:
	    cycles = list(nx.algorithms.cycles.simple_cycles(self))

        nodes_in_cycles = set()
        number = len(cycles)
        if number > 0:
            for c in cycles:
                nodes_in_cycles |= set(c) 
                
            fraction = float(len(nodes_in_cycles))/float(self.number_of_nodes())
        
        if has_loops:
            self.add_edges_from(loops)
             
        return fraction, number  
        
    def longest_path_length(self, draw_path=False):
        """
        Obtains the longest path within the network and consequently, its length. This method offers the possibility to
        draw the path using matplotlib and the NetworkX library capabilities to draw graphs. If there are more than one
        path with the maximum length the path drawn is randomly chosen from the set of possible longest paths.    
        
        *draw_path = True*
            a boolean flag specifying whether to display the graphical representation of the longest path obtained for the network.
    
        It returns the length of the longest path in the network.
        
        .. seealso::
    
            :meth:`mean_path_length` 
                for retrieving the mean length of the shortest paths in the network
        """
        if self.number_of_nodes() == 0:
            return 0 
        
        has_loops = False
        if self.number_of_selfloops() > 0:
            loops = self.selfloop_edges(data=True)
            self.remove_edges_from(loops)
            has_loops = True
        
        has_cycles = False
	#updated to fix conflict between networkx versions
	# nx1.7 returns cycles as list, nx1.9 returns a generator
	if nx.__version__=='1.7': 
            cycles = nx.algorithms.cycles.simple_cycles(self)
	else:
	    cycles = list(nx.algorithms.cycles.simple_cycles(self))

        if len(cycles) > 0:
            has_cycles = True
            removed_edges_cycles = []
            for c in cycles:
                removed_edges_cycles.append((c[0], c[1]))
                
            self.remove_edges_from(removed_edges_cycles)
        
        #we make this structure an attribute of the class to be able to infer trophic levels
        #based on the longest path to each node
        self.length_to = dict.fromkeys(self.nodes(), 0)
        longest_paths = dict.fromkeys(self.nodes())
        
        top_order_nodes = nx.topological_sort(self) 
        for v in top_order_nodes:
            for e in self.out_edges(v):
                if self.length_to[e[1]] <= self.length_to[e[0]] + 1:
                    self.length_to[e[1]] = self.length_to[e[0]] + 1
                    longest_paths[e[1]] = e
        
        max_d = max(self.length_to.values())
        
        #from here on we calculate how many longest paths there are and we trace back one of them
        #to draw it as a network
        if draw_path:
            qty = 0
            for k in self.length_to.keys():
                if self.length_to[k] == max_d:
                    n = k
                    qty += 1
            
            print 'there are', qty, 'paths of length =', max_d 
            
            longest_path = []
            longest_path.append(longest_paths[n])
            next = longest_paths[n][0]
            fin = False
            while not fin:
                if self.length_to[next] == 0:
                    fin = True
                    break
                longest_path.append(longest_paths[next])
                next = longest_paths[next][0]
            
            print longest_path
            
            longest_path_net = nx.DiGraph()
            longest_path_net.add_edges_from(longest_path)
            
            layout = nx.graphviz_layout(longest_path_net, prog="dot", args='-Gnodesep=.07, -Granksep=.1, -Grankdir=BT')
            #layout = nx.spring_layout(net)
            
            #fig = plt.figure()
            #network_plot = fig.add_subplot(111)
            #nx.draw_networkx(longest_path_net, layout, ax=network_plot)    
            
        
        if has_loops:
            self.add_edges_from(loops)
        
        if has_cycles:
            self.add_edges_from(removed_edges_cycles)
        
        return max_d
    
    def mean_path_length(self):
        """
        This method obtains the mean length of the shortest paths from each node to any other node in the network.
        Taking advantage of the :meth:`~networkx.algorithms.shortest_paths.generic.shortest_path_length` algorithm
        for the calculation of the shortest paths from any node in the network, this function obtains an average of
        these lengths.
    
        It returns the average value of the shortest paths from any node to each other in the network.
    
        .. seealso::
    
            :meth:`longest_path_length` 
                for obtaining the longest path in the network
        """
        if self.number_of_nodes() == 0:
            return 0
        else:
            shortest_paths = nx.algorithms.shortest_paths.generic.shortest_path_length(self)
            
            sum_d = 0
            n_paths = 0
            
            for src in shortest_paths.keys():
                for tgt in shortest_paths[src].keys():
                    length = shortest_paths[src][tgt] 
                    if length != 0:
                        sum_d += length
                        n_paths += 1
            if n_paths == 0:
                return 0
            else:
                return float(sum_d)/float(n_paths)
    
    def complexity(self):
        """
        Calculates the complexity index SC (Number of nodes * Connectance): of a network.
        This index has been traditionally linked to the stability of food webs.  
      
        It returns the value of the complexity index number_of_nodes*connectance
        """
        return self.number_of_nodes()*self.connectance()
    
    def dynamic_complexity(self):
        """
        Calculates the dynamic complexity index <i>*sqrt(SC) of a network, as proposed by May.
        This index has been traditionally linked to the stability of food webs.
        
        NOTE: It assumes link weights calculated using a proxy measure for the May's per capita interaction
        strength; and stored in the network as an attribute called 'weight'
      
        It returns the value of the complexity index <interaction_strength>*sqrt(number_of_nodes*connectance)
        """
        
        if self.size() == 0:
            return 0.0
        
        average_i = 0.0
        for u,v,atts in self.edges(data=True):
            if not atts.has_key('weight'):
                return None
            else:
                average_i += atts['weight']
        
        average_i = float(average_i)/self.size()
        
        return average_i*math.sqrt(self.number_of_nodes()*self.connectance())
    
    def components(self):
        """
        This method obtains the number of connected components in the network by calling the method
        :meth:`~networkx.algorithms.components.connected.number_connected_components` from the NetworkX library.
        
        It returns the number of connected components in the network.
        """
        if self.number_of_nodes() == 0:
            return 0
        
        temp_net = self.to_undirected()
        return nx.algorithms.components.connected.number_connected_components(temp_net)

     
    def omnivory(self):
        """
        This method calculates the degree of omnivory in the food web (i.e. the fraction of species that feed in two
        or more trophic levels)
        
            NOTE: this method must be called after calling the method :meth:`longest_path_length`, since
            it relies on the information contained in the data structure length_to in order to arrange species in trophic positions
        
        It returns the fraction of omnivore species within the network and a set containing the species that fulfil this 
        characteristic (i.e. the feed on at least two different trophic levels).
        """
        omnivores = 0
        omni_sps = set()
        if self.number_of_nodes() == 0:
            return omnivores, omni_sps
        
        for n in self.nodes():
            position = self.length_to[n]
            predecessors = self.predecessors(n)
            level_current_prey = -1
            for pred in predecessors:
                if n != pred and position != self.length_to[pred]:
                    if level_current_prey != -1 and level_current_prey != self.length_to[pred]:
                        omnivores += 1
                        omni_sps.add(n)
                        break
                    else:
                        level_current_prey = self.length_to[pred]
                                             
        return float(omnivores)/float(self.number_of_nodes()), omni_sps

    
    def find_all_paths(self, start, end, path=[]):
        """
        This method finds all the existing paths between two given nodes.
        
        call signature::
    
            find_all_paths(start, end, path=[])
        
        *start*
            the node from which to start the search (the source node)
            
        *end*
            the sink or target node, that is, the end node in each one of the paths found
    
        *path*
            since this method is recursive it needs this optional additional argument to pass the
            route obtained so far for a particular path
    
        An array containing all the paths between the *start* and *end* nodes is returned.
        """
        path = path + [start]
        if start == end:
            return [path]
        if not start in self.nodes():
            return []
        paths = []
        for node in self.successors(start):
            if node not in path:
                newpaths = self.find_all_paths(node, end, path)
                for newpath in newpaths:
                    paths.append(newpath)
        return paths

    def find_trophic_positions(self):
        """
        This method calculates the trophic position for each node/species in the network.
            
        It returns:
            1.- a dictionary, where the keys are the node ids, containing the trophic position for each node
            
            2.- a dictionary, where the keys are again the node ids, containing the number of paths existing 
            in the network from all of the basal species to that node
            
            3.- the average length of all the paths that exist between the basal species of the network and 
            each one of the other nodes
            
        .. seealso::
    
            :meth:`find_all_paths` 
                for consulting how the paths between any pair of nodes are found 
        """
        trophic_positions = dict.fromkeys(self.nodes(), 0.0)
        number_of_paths = dict.fromkeys(self.nodes(), 0)
        overall_mean_length = 0
        overall_number_of_paths = 0
        self.path_length_variance = 0.0
        self.path_length_sd = 0.0 
        self.total_number_paths = 0
        
        #we remove the self loops in nodes in order to prevent infinite recursion
        has_loops = False
        if self.number_of_selfloops() > 0:
            loops = self.selfloop_edges(data=True)
            self.remove_edges_from(loops)
            has_loops = True
        
        if self.number_of_nodes() == 0 or self.size() == 0:
            return trophic_positions, number_of_paths, overall_mean_length
        
        #if there are cycles we have to break them before the topological sort
        has_cycles = False

	#updated to fix conflict between networkx versions
	# nx1.7 returns cycles as list, nx1.9 returns a generator
	if nx.__version__=='1.7': 
            cycles = nx.algorithms.cycles.simple_cycles(self)
	else:
	    cycles = list(nx.algorithms.cycles.simple_cycles(self))

        if len(cycles) > 0:
            has_cycles = True
            removed_edges_cycles = []
            for c in cycles:
                removed_edges_cycles.append((c[0], c[1]))
                
            self.remove_edges_from(removed_edges_cycles)
        
        #obtain the topological order to make more efficient the calculation
        top_order_nodes = nx.topological_sort(self) 
        paths = dict()
        
        #we obtain the basal species, from which all the paths will be constructed
        basals = set()
        for n in self.nodes():
            if self.in_degree(n) == 0 and self.out_degree(n) > 0:
                basals.add(n)
                trophic_positions[n] = 1.0
                number_of_paths[n] = 0 
                paths[n] = [[n]]
        
        #finding all the paths to all of the nodes from the base
        for n in top_order_nodes:
            if n in basals:
                continue
            paths[n] = []
            
            predecs = self.predecessors(n)
            
            for predecessor in predecs:
                if len(paths[predecessor]) == 0:
                    all_pred_paths = []
                    for b in basals:
                        pred_paths = self.find_all_paths(b,predecessor)
                        all_pred_paths += pred_paths
                    paths[predecessor] = all_pred_paths
                    
                pred_paths = paths[predecessor]
                    
                for p in pred_paths:
                    p_new = p + [n]
                    paths[n].append(p_new)
        
            mean_length = 0
            number_of_paths[n] = len(paths[n])
            if number_of_paths[n] == 0:
                trophic_positions[n] = 0
                continue
                
            for path in paths[n]:
                mean_length += len(path)-1
                
            overall_mean_length += mean_length
            overall_number_of_paths += number_of_paths[n]
            
            trophic_positions[n] = (float(mean_length)/float(number_of_paths[n])) + 1
            
        if overall_number_of_paths == 0:
            overall_mean_length = 0
        else:
            overall_mean_length = float(overall_mean_length)/float(overall_number_of_paths)
        
        self.total_number_paths = overall_number_of_paths
        
        for sp in paths.keys():
            if sp in basals:
                continue
            for path in paths[sp]:
                self.path_length_variance += ((len(path)-1) - overall_mean_length)**2
        
        if self.total_number_paths == 0:
            self.path_length_variance = 0.0
        else:
            self.path_length_variance = self.path_length_variance/self.total_number_paths
        self.path_length_sd = math.sqrt(self.path_length_variance)
        
        if has_loops:
            self.add_edges_from(loops)
        
        if has_cycles:
            self.add_edges_from(removed_edges_cycles)
        
        return trophic_positions, number_of_paths, overall_mean_length

    def get_path_length_feats(self):
        if self.total_number_paths == 0:
            path_no = 0
        else:
            path_no = math.log10(self.total_number_paths)
             
        return self.path_length_variance, self.path_length_sd, path_no 
        

    def get_cummulative_in_degree_dist(self):
        """
        Obtains the cummulative in-degree distribution for the network, i.e. the cummulative distribution of
        the frequencies with which nodes with a given in-degree occur in the network
            
        It returns an array with the numbers of possible in-links of the nodes in the network and another array
        with the number of nodes possessing that number of in-links (the arrays are correlated in respect to 
        the order of their values i.e. value1 in the first array corresponds to the number of in-links that 
        the number of nodes given by value1 in the second array, possess).  
            
        .. seealso::
    
            :meth:`get_cummulative_out_degree_dist` 
                for the cummulative out-degree distribution
                
            :meth:`get_cummulative_degree_dist` 
                for the cummulative degree distribution
        """
        inds = self.in_degree().values()
        ins = set(inds)
        ins.discard(0)
        set_in = sorted(ins, reverse=True)
        
        #here we calculate the in-degree cumulative degree distribution
        l_in = []
        cum_in = []
        current = 0
        for deg in set_in:
            current += inds.count(deg)
            l_in.append(deg)
            cum_in.append(current)
        
        return l_in, cum_in
    
    def get_cummulative_out_degree_dist(self):
        """
        Obtains the cummulative out-degree distribution for the network, i.e. the cummulative distribution of
        the frequencies with which nodes with a given out-degree occur in the network
            
        It returns an array with the numbers of possible out-links of the nodes in the network and another array
        with the number of nodes possessing that number of out-links (the arrays are correlated in respect to 
        the order of their values i.e. value1 in the first array corresponds to the number of out-links that 
        the number of nodes given by value1 in the second array, possess).  
            
        .. seealso::
    
            :meth:`get_cummulative_in_degree_dist` 
                for the cummulative in-degree distribution
                
            :meth:`get_cummulative_degree_dist` 
                for the cummulative degree distribution
        """
        outds = self.out_degree().values()
        outs = set(outds)
        outs.discard(0)
        set_out = sorted(outs, reverse=True)
        
        #here we calculate the out-degree cumulative degree distribution
        l_out = []
        cum_out = []
        current = 0
        for deg in set_out:
            current += outds.count(deg)
            l_out.append(deg)
            cum_out.append(current)
        
        return l_out, cum_out
        
    def get_cummulative_degree_dist(self):
        """
        Obtains the cummulative degree distribution for the network, i.e. the cummulative distribution of
        the frequencies with which nodes with a given degree occur in the network
            
        It returns an array with the numbers of possible links the nodes in the network might possess and another array
        with the number of nodes possessing that quantity of links (the arrays are correlated in respect to 
        the order of their values i.e. value1 in the first array corresponds to the number of links that 
        the number of nodes given by value1 in the second array, possess).  
            
        .. seealso::
            
            :meth:`get_cummulative_in_degree_dist` 
                for the cummulative in-degree distribution
            
            :meth:`get_cummulative_out_degree_dist` 
                for the cummulative out-degree distribution
        """
        totalds = self.degree().values()
        
        #here we calculate the general cumulative degree distribution
        l_total = []
        cum_total = []
        
        max_deg = max(totalds)
        degrees = sorted(range(1,max_deg+1), reverse=True)
        current = 0
        for deg in degrees:
            current += totalds.count(deg)
            l_total.append(deg)
            cum_total.append(current)
        
        return l_total, cum_total
    

    def generate_random_graph(self):
        """
        This method generates a random graph the the same number of nodes and edge probability (connectance) of the current
        network using the method :meth:`~networkx.generators.random_graphs.erdos_renyi_graph` from the NetworkX library.
        
        It returns a random graph with the same characteristics of number of nodes and edge probability of this network.
        """
        nodes = self.number_of_nodes()
        edge_probability = self.connectance()
        
        return nx.erdos_renyi_graph(nodes, edge_probability, directed=True)
    
    def degree_of_compartmentalization(self):
        """
        Returns the measure of compartmentalization proposed by Pimm and Lawton in ('Are food webs divided into compartments?'. JAE, 49, pp. 879-898. 1980)
        """
        mean_c = 0
        if self.order() == 0 or self.size() == 0 or self.order() == 1:
            return mean_c
            
        vertices = self.nodes()
        for i in range(len(vertices)):
            node_i = vertices[i]
            neighs_i = set(self.successors(node_i)) | set(self.predecessors(node_i))
            for j in range(len(vertices)):
                if j != i:
                    node_j = vertices[j]
                    neighs_j = set(self.successors(node_j)) | set(self.predecessors(node_j))
                    
                    if len(neighs_i | neighs_j) != 0:
                        mean_c += len(neighs_i & neighs_j)/len(neighs_i | neighs_j)
        
        mean_c = float(mean_c)/float(self.order()*(self.order()-1))
                
        return mean_c
    
    def chi_square_deg_freq(self, distribution='poisson'):
        """
        This method obtains the value of the chi square statistic and its corresponding p-value for the cummulative
        degree distribution when compared with another given distribution in the same value range.
        For now, only comparison against a poisson distribution is implemented
        
        *distribution = 'poisson'*
            a string specifying the type of distribution to be compared against the network's distribution of degrees
    
        Returns the chi square and correspnding p-value 
        """
        chi_sq = 0.0
        p_val = 0.0
        if self.order() == 0 or self.size() == 0:
            return chi_sq, p_val
        
        links, distro = self.get_cummulative_degree_dist()  
        control_distro = []
        if distribution == 'poisson':
            poisson = np.random.poisson(self.linkage_density(), sum(distro))
            poisson = list(poisson)
            
            for i in links:
                times = poisson.count(i)
                control_distro.append(times)
            
            #this line is to avoid errors due to zeros in the vector of poisson values
            a1 = np.ma.masked_equal(control_distro, 0)
            
        chi_sq, p_val = stats.chisquare(distro, a1)
    
        return chi_sq, p_val
    
    def robustness(self, criterion='conn', ordering='asc', interval=10, max_removed=60, cumulative=False, weighted=False, beta=0.0, ext_threshold=0.01):
        """
        This method executes extinction experiments over the network in order to determine its response
        to disturbances exerted on the structure and the how its configuration is affected by these
        disturbances.
        
        call signature::
    
          robustness(criterion='conn', ordering='asc', interval=10, max_removed=60, cumulative=False, weighted=False, beta=0.0, ext_threshold=0.01)
    
        This method delegates the tasks of finding the robustness indexes for the intervals and maximum values specified to other two methods:
            `meth`:`robustness_dynamic`
            `meth`:`robustness_static`
        
        Optional keyword arguments:
    
            =============    ===============================================================================================================
            Keyword          Description
            =============    ===============================================================================================================
            criterion        a string specifying the extinction criterion to be used (can take any value from: 'conn'= primary extinctions 
                             by connectance, 'mass'= primary extinctions by the species' biomasses, and 'random'= random primary extinctions 
            ordering         a string specifying the order to be applied on the nodes according to the criterion given by *criterion*.
                             Possible values: 'asc' for ascending and 'desc' for descending. This argument does not apply when the selected
                             criterion is 'random'
            interval         an integer specifying the fraction of the interval (over 100) applied in each interation of the extinction
                             experiments
            max_removed      an integer that determines the maximum fraction of species to be removed from the network for the extinction
                             experiments to stop
            cumulative       a boolean value stating whether the extinction scenario is cumulative (i.e. all the primary extinctions occur
                             before calculating the amount of secondary ones)
            weighted         a boolean flag that when set to True enables the calculations of the robustness indexes in a dynamic way,
                             i.e. considering the links' weights for the calculation of the secondary extinctions
            beta             a float value that is used in the dynamic experiments of robustness to establish the fraction of link weight
                             that a predator is able to redirect towards other prey when one of its preys gets extinct
            ext_threshold    a float value used in the dynamic calculations of robustness to state the fraction of summed weights allowed
                             for a predator to loose before getting extinct (threshold to extinction over the sum of the weights)  
            =============    ===============================================================================================================

        Returns a dictionary with the robustness indexes for each fraction of primary extinctions calcualted as fraction = fraction + interval
        and initial fraction = interval (e.g. [10, 20, 30, 40, 50, 60]) where these fractions are the dictionary keys 
        
        .. seealso::
    
            :meth:`~robustness_dynamic`
                for the calculation of robustness in a dynamical way (looking at the links' weights)
            
            :meth:`~robustness_static`
                for the calculation of the static robustness (based solely on the network structure)
        """
        if criterion == 'conn':
            order = self.get_ordered_degrees(ordering=ordering, type='out')
            
        elif criterion == 'mass':
            order = self.get_ordered_biomasses(ordering=ordering)
            
        elif criterion == 'random':
            order = self.nodes()
            random.shuffle(order)
        
        elif criterion == 'trophic_position':
            order = self.get_ordered_trophic_positions(ordering=ordering)
        
        if weighted:
            return self.robustness_dynamic(order, criterion, interval, max_removed, beta, ext_threshold)
        else:
            return self.robustness_static(order, criterion, interval, max_removed, cumulative)
    
    def robustness_dynamic(self, order, criterion, interval, max_removed, beta, ext_threshold):
        """
        This method assumes the tasks delegated by :meth:`robustness` when the robustness experiments to be performed over
        the network are dynamic, in the sense that the secondary extinctions are calculated based on the sum of the weights
        of the interactions between species rather than on the topological isolation of the nodes.
        
        Although it assumes part of the job of :meth:`robustness` this method can also be called using the following
        signature and arguments.
        
        call signature::
    
          robustness_dynamic(order, criterion, interval, max_removed, beta, ext_threshold)
    
        Arguments:
    
            =============    ===============================================================================================================
            Keyword          Description
            =============    ===============================================================================================================
            order            an array containing the ids of the network's nodes ordered in a particular fashion. This is the order according
                             to which the nodes are going to be removed from the network 
            criterion        a string that usually takes one of the values available for the same argument in :meth:`robustness` but in this
                             case only if it is equal to *random* provides information to the logic of the method
            interval         an integer specifying the fraction of the interval (over 100) applied in each interation of the extinction
                             experiments
            max_removed      an integer that determines the maximum fraction of species to be removed from the network for the extinction
                             experiments to stop
            beta             a float value that is used to establish the fraction of link weight that a predator is able to redirect towards
                             other prey when one of its preys gets extinct
            ext_threshold    a float value used to specify the fraction of summed weights allowed for a predator to loose before getting 
                             extinct (threshold to extinction over the sum of the weights)  
            =============    ===============================================================================================================

        Returns a dictionary with the robustness indexes for each fraction of primary extinctions calcualted as fraction = fraction + interval
        and initial fraction = interval (e.g. [10, 20, 30, 40, 50, 60]) where these fractions are the dictionary keys 
        
        .. seealso::
    
            :meth:`~robustness`
                for the method that delegates the robustness calculations to this, acting as a generic method for static and dynamic robustness
                calculations
            
            :meth:`~robustness_static`
                for the calculation of the static robustness (based solely on the network structure)
        """
        robustness = dict()
        
        #we obtain the basal species because they are not going to get extinct as a cascade effect
        basal = set()
        in_degrees = self.in_degree()
        for sp in in_degrees.keys():
            if in_degrees[sp] == 0:
                basal.add(sp)
        
        ratios = []
        curr_rat = interval
    
        while curr_rat <= max_removed:
            ratios.append(curr_rat)
            curr_rat += interval
        
        for f in ratios:
            net_temp = self.copy()
            net_temp.remove_edges_from(net_temp.selfloop_edges(data=True))
            net_temp.obtain_interactions_strengths()
            
            to_remove = math.ceil(float(self.order())*(float(f)/100))
            removed = 0
            total_loss = 0
            original_n = self.order()
            
            for i in xrange(0,len(order)):        
                if total_loss >= (original_n - 1):
                    break
                
                if criterion == 'random':
                    sp = order[i]
                else:
                    sp, val = order[i]
          
                if sp in net_temp.nodes():
                    if beta != 0.0:
                        predators = net_temp.successors(sp)
                        for pred in predators:
                            old_weight = net_temp[sp][pred]['weight']
                            links_predator = net_temp.in_edges(pred, data=True)
                            avail_preys = len(links_predator)-1
                            
                            if avail_preys == 0:
                                continue
                            
                            delta_weight = (beta*old_weight)/avail_preys
                            
                            for x,y,atts in links_predator:
                                if x == sp:
                                    continue
                                atts['weight'] += delta_weight
                        
                    net_temp.remove_node(sp)
                    removed += 1
                    total_loss += 1
                
                    #secondary extinctions
                    net_temp.update_total_weights(type='in')
                    for v in net_temp.nodes():
                        if v in basal:
                            continue
                        current_w = net_temp.node[v]['current_in_weight']
                        original_w = net_temp.node[v]['original_in_weight']
                        threshold_w =  original_w - (original_w*ext_threshold)
                        if current_w < threshold_w:
                            net_temp.remove_node(v)
                            total_loss += 1
                
                    #stopping condition, if we have removed all we were supposed to remove
                    if removed >= to_remove:
                        value = float(total_loss)/float(original_n)
                        robustness[f] = value
                        break
        
        return robustness
        
        
    def robustness_static(self, order, criterion, interval, max_removed, cumulative):
        """
        This method assumes the tasks delegated by :meth:`robustness` when the robustness experiments to be performed over
        the network are static, in the sense that the secondary extinctions are calculated solely based on the topological
        isolation of nodes that occurs as a consequence of an intentional removal of other species in the network.
        
        Although it assumes part of the job of :meth:`robustness` this method can also be called using the following
        signature and arguments.
        
        call signature::
    
          robustness_dynamic(order, criterion, interval, max_removed, cumulative)
    
        Arguments:
    
            =============    ===============================================================================================================
            Keyword          Description
            =============    ===============================================================================================================
            order            an array containing the ids of the network's nodes ordered in a particular fashion. This is the order according
                             to which the nodes are going to be removed from the network 
            criterion        a string that usually takes one of the values available for the same argument in :meth:`robustness` but in this
                             case only if it is equal to *random* provides information to the logic of the method
            interval         an integer specifying the fraction of the interval (over 100) applied in each interation of the extinction
                             experiments
            max_removed      an integer that determines the maximum fraction of species to be removed from the network for the extinction
                             experiments to stop
            cumulative       a boolean value stating whether the extinction scenario is cumulative (i.e. all the primary extinctions occur
                             before calculating the amount of secondary ones)
            =============    ===============================================================================================================

        Returns a dictionary with the robustness indexes for each fraction of primary extinctions calcualted as fraction = fraction + interval
        and initial fraction = interval (e.g. [10, 20, 30, 40, 50, 60]) where these fractions are the dictionary keys 
        
        .. seealso::
    
            :meth:`~robustness`
                for the method that delegates the robustness calculations to this, acting as a generic method for static and dynamic robustness
                calculations
            
            :meth:`~robustness_dynamic`
                for the calculation of robustness in a dynamical way (looking at the links' weights)
        """
        robustness = dict()
        #we obtain the basal species because they are not going to get extinct as a cascade effect
        basal = set()
        in_degrees = self.in_degree()
        for sp in in_degrees.keys():
            if in_degrees[sp] == 0:
                basal.add(sp)
        
        if cumulative:
            net_temp = self.copy()
            net_temp.remove_edges_from(net_temp.selfloop_edges(data=True))
            to_remove = math.ceil(float(self.order())*(float(interval)/100))
            removed = 0
            total_removed = 0
            total_loss = 0
            fraction_removed = interval
            original_n = self.order()
            
            for i in xrange(0,len(order)):
                if total_loss >= (original_n - 1):
                    break
                
                if criterion == 'random':
                    sp = order[i]
                else:
                    sp, val = order[i]
                
                if sp in net_temp.nodes():
                    net_temp.remove_node(sp)
                    removed += 1
                    total_loss += 1
                    
                if removed >= to_remove:
                    #secondary extinctions
                    new_order = net_temp.get_ordered_degrees(type='in')
                    for n, d in new_order:
                        if n in basal:
                            continue
                        if d > 0:
                            break
                        net_temp.remove_node(n)
                        total_loss += 1
                        
                    total_removed += removed
                    removed = 0
                    key = str(total_removed)+'/'+str(original_n)
                    val = str(total_loss)+'/'+str(original_n)
                    k = float(total_removed)/float(original_n)
                    value = float(total_loss)/float(original_n)
                    robustness[fraction_removed] = value
                    if fraction_removed >= max_removed:
                        break 
                    fraction_removed += interval
        else:
            ratios = []
            curr_rat = interval
        
            while curr_rat <= max_removed:
                ratios.append(curr_rat)
                curr_rat += interval
            
            for f in ratios:
                net_temp = self.copy()
                net_temp.remove_edges_from(net_temp.selfloop_edges(data=True))
                to_remove = math.ceil(float(self.order())*(float(f)/100))
                removed = 0
                total_loss = 0
                original_n = self.order()
                
                for i in xrange(0,len(order)):        
                    if total_loss >= (original_n - 1):
                        break
                    
                    if criterion == 'random':
                        sp = order[i]
                    else:
                        sp, val = order[i]
                    
                    if sp in net_temp.nodes():
                        net_temp.remove_node(sp)
                        removed += 1
                        total_loss += 1
                        
                    #secondary extinctions
                    new_order = net_temp.get_ordered_degrees(type='in')
                    for n, d in new_order:
                        if n in basal:
                            continue
                        if d > 0:
                            break
                        net_temp.remove_node(n)
                        total_loss += 1
                    
                    if removed >= to_remove:
                        value = float(total_loss)/float(original_n)
                        robustness[f] = value
                        break
                        
        
        return robustness
        
    def get_ordered_degrees(self, ordering='asc', type='all'):
        """
        Returns an array containing the ids of the nodes in the network ordered according to a particular type of node degree and in
        the specified order.
        
        call signature::
    
          get_ordered_degrees(ordering='asc', type='all')
    
        Optional keyword arguments:
    
            =============    =============================================================================================================
            Keyword          Description
            =============    =============================================================================================================
            ordering         a string specifying the order to be applied on the nodes. Possible values: 'asc' for ascending and 'desc' for 
                             descending.
            type             a string that determines the type of node degree to be considered for ordering the nodes. Possible values:
                             *all*, *in*, and *out*; for degree, in-degree, and out-degree respectively.
            =============    =============================================================================================================
        
        **Examples:** ::
        
            get_ordered_degrees(ordering='asc', type='all')
        
        returns an array of the nodes in ascending (from lowest to highest) order according to their degree ::
        
            get_ordered_degrees(ordering='desc', type='out')
        
        returns an array of the nodes in descending (from highest to lowest) order according to their out-degree (number of outgoing links)    
        """
        if type == 'all':
            degs = self.degree()
        elif type == 'out':
            degs = self.out_degree()
        elif type == 'in':
            degs = self.in_degree()
        else:
            return []
            
        if ordering == 'asc':
            return sorted(degs.items(), key=itemgetter(1))
        elif ordering == 'desc':
            return sorted(degs.items(), key=itemgetter(1), reverse=True)
        else:
            return []
    
    def get_ordered_biomasses(self, ordering='asc'):
        """
        Returns an array containing the ids of the nodes in the network ordered according to their biomass. The type of order is specified
        by the argument *ordering*
        
        call signature::
    
          get_ordered_biomasses(ordering='asc')
    
        Optional keyword arguments:
    
            =============    ================================================================================
            Keyword          Description
            =============    ================================================================================
            ordering         a string specifying the order to be applied on the nodes. Possible values: 'asc' 
                             for ascending and 'desc' for descending.
            =============    ================================================================================
        
        **Example:** ::
        
            get_ordered_biomasses(ordering='asc')
        
        returns an array of the nodes in ascending (from lowest to highest) order according to their biomass    
        """
        biomasses = dict()
        for n in self.nodes():
            mass = self.node[n]['biomass']
            if mass == ' ':
                biomasses[n] = 0.0
            else:
                biomasses[n] = float(mass)
        
        if ordering == 'asc':
            return sorted(biomasses.items(), key=itemgetter(1))
        elif ordering == 'desc':
            return sorted(biomasses.items(), key=itemgetter(1), reverse=True)
        else:
            return []
    
    def get_ordered_weights(self):
        """
        Returns an array containing the ids of the nodes in the network in ascending order according to the sum of the weights of their links to all
        of its predecessors/preys.
        """
        weights = dict()
        for n in self.nodes():
            w = 0.0
            predecessors = self.predecessors(n)
            if len(predecessors) == 0:
                weights[n] = w
                continue
            
            for v in predecessors:
                w += self[v][n]['weight']
            
            weights[n] = w
        
        return sorted(weights.items(), key=itemgetter(1))
    
    def get_ordered_trophic_positions(self, ordering='asc'):
        """
        Returns an array containing the ids of the nodes in the network ordered according to the trophic positions
        of the species in the food web represented by each node, either in ascending or descending order
        
        call signature::
    
          get_ordered_trophic_positions(ordering='asc')
    
        Optional keyword arguments:
    
            ========    ================================================================================
            Keyword     Description
            ========    ================================================================================
            ordering    a string specifying the order to be applied on the nodes. Possible values: 'asc' 
                        for ascending and 'desc' for descending.
            ========    ================================================================================
        
        .. seealso::
    
            :meth:`~find_trophic_positions`
                For finding the trophic positions of the nodes
        """
        tps, nops, m_length = self.find_trophic_positions()
        
        if ordering == 'asc':
            return sorted(tps.items(), key=itemgetter(1))
        elif ordering == 'desc':
            return sorted(tps.items(), key=itemgetter(1), reverse=True)
        else:
            return []
    
    
    def mean_body_mass_ratio(self):
        """
        This method calculates the average of the predator:prey body mass ratio for all the pairs of nodes that are connected in the network.
        Returns: The logarithm in base 10 of the mean
        """
        mean_ratio = 0
        edges_count = 0
        
        for prey, predator in self.edges():
            mass_prey = self.node[prey]['biomass']
            mass_predator = self.node[predator]['biomass']
            if mass_prey == ' ' or mass_prey == '0.0' or mass_prey == None or mass_predator == ' ' or mass_predator == '0.0' or mass_predator == None:
                continue
        
            ratio = float(mass_predator)/float(mass_prey)
            self[prey][predator]['mass_ratio'] = ratio
            mean_ratio += ratio
            edges_count += 1
        
        if edges_count == 0:
            return 0.0
        
        return math.log10(mean_ratio/edges_count)
    
    def get_min_max_biomasses(self):
        """
        This method returns three float values that report the minimum, maximum, and average biomass of all the species present in the 
        ecological network respectively.
        """
        min_biomass = None
        max_biomass = None
        mean_biomass = 0.0
        species_count = 0
        
        for n in self.nodes():
            current_mass = self.node[n]['biomass']
            if current_mass == ' ' or current_mass == '0.0' or current_mass == None:
                continue
            
            current_mass = float(current_mass)
            mean_biomass += current_mass
            species_count += 1
            
            if max_biomass == None or current_mass > max_biomass:
                max_biomass = current_mass
            
            if min_biomass == None or current_mass < min_biomass:
                min_biomass = current_mass
        
        if species_count > 0:
            mean_biomass = mean_biomass/species_count
        
        return min_biomass, max_biomass, mean_biomass
    
    def obtain_interactions_strengths(self, scaling='3/4', normalise=False):
        """
        This method calculates the interaction strengths for each of the links present between any two species in the network. It is designed
        to support the calculation of these strengths based on different criteria, even though for now it only supports their calculation
        based on the 3/4 scaling law of predator/prey biomass ratio: i.e. (predator biomass)^0.75 / (prey biomass). The value for each
        strength is stored in the network as an edge attribute for each link.
        
        call signature::
    
          obtain_interactions_strengths(scaling='3/4', normalise=False)
    
        Optional keyword arguments:
    
            =============    ===============================================================================================================
            Keyword          Description
            =============    ===============================================================================================================
            scaling          a string symbolising the type of scaling to be considered for the calculations of the strengths of interactions
                             between species in the network
            normalise        a boolean flag used to determine whether to normalise the values obtained for each interaction based on the 
                             lowest negative value which absolute value is then added to each of the other calculated values to avoid
                             having negative interaction strengths facilitating thus the data analysis
            =============    ===============================================================================================================
        
            NOTE: the 3/4 scaling is calculated using the log10 of the term in order to correct for the bias of the fractions' differences
        
        .. seealso::
    
            :meth:`update_total_weights`
                for the calculation of the sum of all the weights for each node
        """
        if scaling == '3/4':
            if normalise:
                larger_negative = 0.0
            
            for prey, predator, atts in self.edges(data=True):
                mass_prey = self.node[prey]['biomass']
                mass_predator = self.node[predator]['biomass']
                if mass_prey == ' ' or mass_prey == '0.0' or mass_predator == ' ' or mass_predator == '0.0':
                    atts['weight'] = 0.0
                    continue
                
                value = math.log10((float(mass_predator)**(0.75))/(float(mass_prey)))
                atts['weight'] = value
                
                if normalise:
                    if value < larger_negative:
                        larger_negative = value
            
            if normalise:
                larger_negative = math.ceil(math.fabs(larger_negative))
                for prey, predator, atts in self.edges(data=True):
                    if atts['weight'] != 0.0:
                        atts['weight'] += larger_negative
            
            self.update_total_weights(update_originals=True)
    
    
    def update_total_weights(self, update_originals=False, type='both'):
        """
        This method obtains the sum of the interaction strengths of in-edges and out-edges for each node and stores these values as node's
        attributes. It stores these values duplicated as *original* and *current* weight, the latter to be modified when any of the strengths
        changes, facilitating in this way to keep track of the changes suffered by each species and its interactions.
        
        call signature::
    
            update_total_weights(update_originals=False, type='both')
    
        Optional keyword arguments:
    
            ================     =======================================================================================================
            Keyword              Description
            ================     =======================================================================================================
            update_originals     a boolean value stating whether to update the original sum of the weights. Ideally the original sum of 
                                 weights should be kept constant throughout the analyses, changing the current weights sum and comparing
                                 with the former to observe the changes.
            type                 a string stating the type of connections to be considered for the sum calculations. It can take one of
                                 the following values: 'in' for incoming connections, 'out' for outgoing connections, or 'both' for
                                 calculating the sum of both types of connections for each node
            ================     =======================================================================================================
        
        .. seealso::
    
            :meth:`obtain_interactions_strengths`
                for the calculation of the interaction strenghts
        """ 
        for v in self.nodes():
            if type == 'in' or type == 'both':
                total_in_weight = 0.0
                preds = self.predecessors(v)
                for p in preds:
                    total_in_weight += self[p][v]['weight']
                self.node[v]['current_in_weight'] = total_in_weight
                
                if update_originals: 
                    self.node[v]['original_in_weight'] = total_in_weight
            
            if type == 'out' or type == 'both':
                total_out_weight = 0.0
                sucs = self.successors(v)
                for s in sucs:
                    total_out_weight += self[v][s]['weight']
                self.node[v]['current_out_weight'] = total_out_weight
                
                if update_originals:    
                    self.node[v]['original_out_weight'] = total_out_weight
                
            
    def get_copy_removed_nodes(self, removed_nodes=[]):
        """
        Returns a copy of the network when the nodes in its only argument *removed_nodes* are removed from it. By default this argument,
        which is optional, is an empty array, in which case this method simply returns a copy of the original network.
        """ 
        new_net = self.copy()
        if len(removed_nodes) > 0:
            for n in removed_nodes:
                try:
                    new_net.remove_node(n)
                except:
                    #print 'The node '+ n +' does not exist'
                    pass
            for n in new_net.nodes():
                if new_net.degree(n) == 0:
                    new_net.remove_node(n)
            
        return new_net
    
    def get_invaded_network(self, introduced_group=None, invasive_name='Invasive sp.', biomass_introduced=None, generalism=0.5, predators=0.5, ext_threshold=0.01, beta=0.5):
        """
        This method is used to obtain a copy of the network when it is invaded by a characterised species. The introduced species is
        defined by the optional arguments that are received by this method. The algorithm implemented by this method obtains a set
        of likely interactions that the introduced species could possess based on its similarity to other species already in the
        network that belong to the same species group as the introduced one and on the degree of generalism of the invasor. In addition
        to determining its structural position within the network in the manner described above, this method also calculates the interaction
        strengths of the newly introduced species with its interactions partners and adjusts the weights of the other links the second
        species possess based on a criterion similar to that used in the robustness to extinctions experiments, in which a fraction
        of the new weight is either added or removed to the weights of the other links depending on whether the links are to predators
        or preys of the partner species being considered. If when rebalancing the interactions strenghts of any species in this way
        the new sum of its weights exceeds or goes below (for out and in links respectively) a given threshold, extinctions may occur, which
        implies that the new, invaded network could lack one or more of its original species.
        
        call signature::
    
          get_invaded_network(introduced_group=None, invasive_name='Invasive sp.', biomass_introduced=None, generalism=0.5, ext_threshold=0.01, beta=0.5)
    
        Optional keyword arguments:
    
            ==================    ===============================================================================================================
            Keyword               Description
            ==================    ===============================================================================================================
            introduced_group      a string specifying the species group the introduced species belongs to. It must take values within the range
                                  of available groups: 'Bird', 'BirdPrey', 'Reptile', 'Mammal', 'MammalCarn', 'Fish', 'Amphibian', 'Invertebrate'
                                  Its default value is set to *None*, in which case the method returns nothing. 
            invasive_name         a string used to name the species/node to be introduced into the network
            biomass_introduced    a float number stating the biomass of the introduced species, which will be used to calculate its interations
                                  strenghts. Its default value is set to None, in which case the biomass of the introduced species is
                                  taken as the mean of the species already present in the network that belong to the same group as the invader
            generalism            a float number between 0 and 1 that determines the fraction of incoming links (preys) possessed by all of the 
                                  species in the invader's group that it will possess. This is a measure of how generalist we want the introduced
                                  species to be.
            predators             a float number between 0 and 1 that determines the fraction of outgoing links (predators) possessed by all of 
                                  the species in the invader's group that it will possess. This is a measure of how vulnerable to predation we 
                                  want the introduced species to be.
            ext_threshold         a float value between 0 and 1 used to state the fraction of summed weights allowed for a predator to loose 
                                  or for a prey to get before getting extinct (threshold to extinction over the sum of the weights) when the 
                                  new species has been introduced and the sum of weights of the other species calculated
            beta                  a float value between 0 and 1 that determines the fraction of the weight of the link between the newly
                                  introduced species and each of its prey that the other predators of that predator are going to loose because
                                  of the effect of the new species and the fraction of weight of the link between each of those predators and 
                                  each of their preys that the other preys are going to gain as a consequence of the pressure of the new predator 
            ==================    ===============================================================================================================

        Returns a new network with the introduced species and possibly some of the former species extinct due to the proccesses described above.
        """
        if introduced_group == None:
            return False
        
        new_net = self.copy()
        self.obtain_interactions_strengths(normalise=True)
        nodes_group = new_net.get_nodes_by_species_group(introduced_group)
        
        if len(nodes_group) == 0:
            #print 'There are no nodes in the network belonging to the group "',introduced_group, '" - No network returned'
            message = str('there are no nodes in the network belonging to the group ') + str(introduced_group) + str(' - No network returned')
            raise NotInvadableNetwork(message, None)
            return False
        
        average_in_weight = 0.0
        average_out_weight = 0.0
        speciesin_count = 0
        speciesout_count = 0
        
        for n in nodes_group:
            if self.node[n]['original_in_weight'] > 0.0:
                average_in_weight += self.node[n]['original_in_weight']
                speciesin_count += 1
                
            if self.node[n]['original_out_weight'] > 0.0:
                average_out_weight += self.node[n]['original_out_weight']
                speciesout_count += 1
                
        if speciesin_count == 0:
            average_in_weight = 0.0
        else:
            average_in_weight = average_in_weight/speciesin_count
        
        if speciesout_count == 0:
            average_out_weight = 0.0
        else:
            average_out_weight = average_out_weight/speciesout_count
        
        
        if biomass_introduced == None or biomass_introduced == 0.0:
            biomass = 0.0
            species_count = 0
            for n in nodes_group:
                biomass_node = float(new_net.node[n]['biomass'])
                if biomass_node > 0.0:
                    biomass += biomass_node
                    species_count += 1
                    
            biomass_introduced = biomass/species_count
        
        attr_new = {'group':introduced_group, 'biomass':biomass_introduced}    
        new_net.add_node(invasive_name, attr_dict=attr_new)
        
        succs = []
        preds = []
        
        for n in nodes_group:
            succs += new_net.successors(n)
            preds += new_net.predecessors(n)
             
        succs_set = set(succs)
        preds_set = set(preds)
        
        succs = []
        for n in succs_set:
            succs.append(n)
            
        random.shuffle(succs)
        number_of_succs = len(succs)
        succs_links = int(math.ceil(float(number_of_succs)*predators))
        added = 0
        for n in succs:
            if new_net.node[n]['biomass'] > biomass_introduced:
                new_net.add_edge(invasive_name, n)
                added += 1
                if added >= succs_links:
                    break
                
#        for n in succs:
#            if new_net.node[n]['biomass'] > biomass_introduced:
#                new_net.add_edge(invasive_name, n)
#        
        preds = []
        for n in preds_set:
            if float(new_net.node[n]['biomass']) == 0.0:
                continue 
            preds.append(n)
            
        random.shuffle(preds)
        number_of_preds = len(preds)
        pred_links = int(math.ceil(float(number_of_preds)*generalism))
        added = 0
        invasive_succs = set(new_net.successors(invasive_name))
        
        for n in preds:
            n_biomass = float(new_net.node[n]['biomass']) 
            if n not in invasive_succs and n_biomass < biomass_introduced:
                new_net.add_edge(n, invasive_name)
                added += 1
                if added >= pred_links:
                    break
        
        if new_net.degree(invasive_name) == 0:
            message = str('invasive species has no links')
            raise NotInvadableNetwork(message, None)
            return False
        
        if number_of_preds > 0 and new_net.in_degree(invasive_name) == 0:
            print 'invasive species has no prey'
            message = str('invasive species has no prey')
            raise NotInvadableNetwork(message, None)
            return False
        
        #after adding all the new links we calculate the links' weights for all the network...
        new_net.obtain_interactions_strengths(normalise=True)
        
        if new_net.node[invasive_name]['original_in_weight'] < average_in_weight:
            #print 'invasion failed --- invasive in weight: ',new_net.node[invasive_name]['original_in_weight'], ', average in weight: ', average_in_weight, 'not enough resources'
            message = str('invasion failed --- invasive in weight: ') + str(new_net.node[invasive_name]['original_in_weight']) + str(' average in weight: ') + str(average_in_weight) + str(' not enough resources')
            raise NotInvadableNetwork(message, new_net)
            return False
        
        if new_net.node[invasive_name]['original_out_weight'] > average_out_weight:
            #print 'invasion failed --- invasive out weight: ',new_net.node[invasive_name]['original_out_weight'], ' average out weight: ', average_out_weight, 'predation pressure too large'
            message = str('invasion failed --- invasive out weight: ') + str(new_net.node[invasive_name]['original_out_weight']) + str(', average out weight: ') + str(average_out_weight) + str(' predation pressure too large')
            raise NotInvadableNetwork(message, new_net)
            return False
        
        #...and then apply the relevant transformations to the introduced species links' weights
        
        net_post = new_net.copy()
        
        preys = net_post.predecessors(invasive_name)
        for v in preys:
            predators = net_post.successors(v)
            gamma = net_post[v][invasive_name]['weight']
            n_predators = len(predators)-1
            delta_weights_predators = ((1/n_predators)*beta*gamma)
            for n in predators:
                if n != invasive_name:
                    net_post[v][n]['weight'] = net_post[v][n]['weight'] - delta_weights_predators 
                    
                    second_level_preys = net_post.predecessors(n)
                    n_preys = len(second_level_preys)
                    delta_weights_second_preys = ((1/n_preys)*beta*delta_weights_predators)
                    for p in second_level_preys:
                        net_post[p][n]['weight'] = net_post[p][n]['weight'] + delta_weights_second_preys
        
        net_post.update_total_weights()
        #here we verify whether the pressure over the preys of the introduced species has reached a certain threshold
        #after which they get extinct. Similarly, we check whether the amount of prey is enough to sustain
        #each predator species
        for v in net_post.nodes():
            if v != invasive_name:
                if net_post.has_node(v):
                    current_w = net_post.node[v]['current_out_weight']
                    original_w = self.node[v]['original_out_weight']
                    threshold_w =  original_w + (original_w*ext_threshold)
                    if current_w > threshold_w:
                        net_post.remove_node(v)
                
                if net_post.has_node(v):  
                    current_w = net_post.node[v]['current_in_weight']
                    original_w = self.node[v]['original_in_weight']
                    threshold_w =  original_w - (original_w*ext_threshold)
                    if current_w < threshold_w:
                        net_post.remove_node(v)
        
        degrees = net_post.degree()
        for v in net_post.nodes():
            if degrees[v] == 0:
                net_post.remove_node(v)
        
        return new_net, net_post
    
    
    def get_links_fractions_between_levels(self, basal=None, intermediate=None, top=None):
        if basal == None:
            b, basal = self.basal(heterotrophs=False)
        if intermediate == None:
            i, intermediate = self.intermediate(heterotrophs=False)
        if top == None:
            t, top = self.top_predators()
            
        top_links = set(self.in_edges(top))
        intermediate_links = set(self.in_edges(intermediate)) & set(self.out_edges(intermediate))
        basal_links = set(self.out_edges(basal))
        
        t_b_links = 0
        t_i_links = 0
        i_i_links = 0
        i_b_links = 0
        
        for prey, pred in top_links:
            if prey in intermediate:
                t_i_links += 1
            if prey in basal:
                t_b_links += 1

        for prey, predator in intermediate_links:
            if prey in intermediate and predator in intermediate:
                i_i_links += 1
        
        for prey, predator in basal_links:
            if predator in intermediate:
                i_b_links += 1
        
        fractions = dict()
        if self.size() == 0:
            fractions['tb'] = 0
            fractions['ti'] = 0
            fractions['ii'] = 0
            fractions['ib'] = 0
        else:
            fractions['tb'] = float(t_b_links)/float(self.size())
            fractions['ti'] = float(t_i_links)/float(self.size())
            fractions['ii'] = float(i_i_links)/float(self.size())
            fractions['ib'] = float(i_b_links)/float(self.size())
        
        return fractions
            
                  
    def generality_vulnerability_sd(self):
        if self.size() == 0 or self.number_of_nodes() == 0:
            return 0.0, 0.0
        
        mean_gen = 0.0 
        mean_vul = 0.0
        generalities = []
        vulnerabilities = []
        
        ld = self.linkage_density()
        for n in self.nodes():
            gn = float(self.in_degree(n))/ld
            vn = float(self.out_degree(n))/ld
            
            mean_gen += gn
            mean_vul += vn
            
            generalities.append(gn)
            vulnerabilities.append(vn)
            
        mean_gen = mean_gen/self.number_of_nodes()
        mean_vul = mean_vul/self.number_of_nodes()
        
        gen_variance = 0.0
        for g in generalities:
            gen_variance += (g-mean_gen)**2
            
        gen_variance = gen_variance/len(generalities)
            
        vul_variance = 0.0
        for v in vulnerabilities:
            vul_variance += (v-mean_vul)**2
        vul_variance = vul_variance/len(vulnerabilities)
        
        gen_sd = math.sqrt(gen_variance)
        vul_sd = math.sqrt(vul_variance)    
        
        return gen_sd, vul_sd    
        
    
    def maximum_similarity(self):
        max_sim = 0.0
        if self.size() == 0 or self.number_of_nodes() < 2:
            return max_sim
        
        for i in self.nodes():
            predators_i = set(self.successors(i))
            preys_i = set(self.predecessors(i))
            max_sim_ij = None
            
            for j in self.nodes():
                if j == i:
                    continue
                predators_j = set(self.successors(j))
                preys_j = set(self.predecessors(j))
                
                if float((len(predators_i | predators_j) + len(preys_i | preys_j))) == 0.0:
                    sim_ij = 0.0
                else:
                    sim_ij = float((len(predators_i & predators_j) + len(preys_i & preys_j))) / float((len(predators_i | predators_j) + len(preys_i | preys_j))) 
                
                if max_sim_ij == None or sim_ij > max_sim_ij:
                    max_sim_ij = sim_ij
                
            max_sim += max_sim_ij
            
        max_sim = max_sim/self.number_of_nodes()
        
        return max_sim
        
    def modularity_rgraph(self, seed, iter_factor, cooling_factor, randoms):
        self.modularity_of_randomizations = 0.0
        self.sd_mod_of_randomizations = 0.0
        
        if self.number_of_nodes < 2 or self.size() == self.number_of_selfloops():
            return 0.0, 0
        
        tmp_dir = '../temp_'+str(datetime.datetime.now()).replace(' ','')+'/'
        os.mkdir(tmp_dir)
        file_out = open(tmp_dir+'temp_out.out', 'w')
        try:
            retcode = subprocess.call(['netcarto_cl'], stdout=file_out, stderr=subprocess.STDOUT)
        except OSError:
            print 'Netcarto is not installed or is not available on the PATH in this computer.'
            print 'Cannot obtain the modularity using the rgraph library.'
            print 'Check your rgraph installation.'
            return 0.0, 0
        
        dict_nodes_numbers = dict()
        dict_numbers_nodes = dict()
        nodes = self.nodes()
        for idx in range(len(nodes)):
            dict_nodes_numbers[nodes[idx]] = idx
            dict_numbers_nodes[idx] = nodes[idx]
        
        file_net_temp = open(tmp_dir+'temp_net.adj', 'w')
        
        for u,v in self.edges():
            if u != v:
                file_net_temp.write(str(dict_nodes_numbers[u])+' '+str(dict_nodes_numbers[v])+'\n')
        
        file_net_temp.close()
        
        args = ['netcarto_cl', 'temp_net.adj', str(seed), '-1', str(iter_factor), str(cooling_factor), str(randoms)]
        subprocess.call(args, stdout=file_out, stderr=subprocess.STDOUT, cwd=tmp_dir)
        file_out.close()
        
        file_modules = open(tmp_dir+'modules.dat', 'r')
        modules = []
        for line in file_modules:
            a,b,mods = line.partition('---')
            
            if mods != '':
                modules.append(mods.split())
            else:
                start = line.index('=')
                modularity = line[start+1:].strip()
        
        file_modules.close()
                    
        self.number_of_modules = len(modules)
        self.modularity = modularity
        
        #we assign each node to its corresponding module
        current_module = 1
        for module in modules:            
            for n in module:
                self.node[dict_numbers_nodes[int(n)]]['module'] = current_module
        
            current_module += 1
        
        #we also find the role of each node and add this information to the network
        file_roles = open(tmp_dir+'roles.dat', 'r')
        roles = dict()
        for line in file_roles:
            a,b,nodes = line.partition('---')            
            role = a[0]
            roles[role] = nodes.split()
            
        file_roles.close()
        
        for r in roles.keys():
            for n in roles[r]:
                self.node[dict_numbers_nodes[int(n)]]['role'] = r
        
        if randoms > 0:
            file_randoms = open(tmp_dir+'randomized_mod.dat', 'r')
        
            for line in file_randoms:
                if line.startswith('#'):
                    continue
                
                values = line.split()
                self.modularity_of_randomizations = values[1]
                self.sd_mod_of_randomizations = values[2]
                    
            file_randoms.close()
        
        return self.modularity, self.number_of_modules
        
        
        
    def robustness_to_removal(self, to_remove, beta=0.5, ext_threshold=0.5):                
        #we obtain the basal species because they are not going to get extinct as a cascade effect
        basal = set()
        in_degrees = self.in_degree()
        for sp in in_degrees.keys():
            if in_degrees[sp] == 0:
                basal.add(sp)
        
        net_removed = self.copy()
        net_removed.obtain_interactions_strengths()
        
        for r in to_remove:        
            if r in net_removed.nodes():
                if beta != 0.0:
                    predators = net_removed.successors(r)
                    for pred in predators:
                        old_weight = net_removed[r][pred]['weight']
                        links_predator = net_removed.in_edges(pred, data=True)
                        avail_preys = len(links_predator)-1
                        
                        if avail_preys == 0:
                            continue
                        
                        delta_weight = (beta*old_weight)/avail_preys
                        
                        for x,y,atts in links_predator:
                            if x == r:
                                continue
                            atts['weight'] += delta_weight
                    
                net_removed.remove_node(r)
            
                #secondary extinctions
                net_removed.update_total_weights(type='in')
                for v in net_removed.nodes():
                    if net_removed.degree(v) == 0:
                        net_removed.remove_node(v)
                        continue
                    
                    if v in basal:
                        continue
                    
                    current_w = net_removed.node[v]['current_in_weight']
                    original_w = net_removed.node[v]['original_in_weight']
                    threshold_w = original_w - (original_w*ext_threshold)
                    if current_w < threshold_w:
                        net_removed.remove_node(v)
                      
        loss = self.number_of_nodes() - net_removed.number_of_nodes()
        robustness = float(loss)/float(self.number_of_nodes())
                
        return net_removed, robustness
    
    
    def get_shortest_chain_length(self):
        scl = dict.fromkeys(self.nodes(), None)
        b, b_species = self.basal()
        
        for sp in b_species:
            scl[sp] = 0
        
        for sp_b in b_species:
            short_paths = nx.single_source_shortest_path_length(self, sp_b)
            for sp in short_paths.keys():
                if scl[sp] == None or short_paths[sp] < scl[sp]:
                    scl[sp] = short_paths[sp]
        
        return scl
        
    def get_trophic_levels(self):
        
        has_loops = False
        if self.number_of_selfloops() > 0:
            loops = self.selfloop_edges(data=True)
            self.remove_edges_from(loops)
            has_loops = True

        tls = dict.fromkeys(self.nodes(), None)
        producers = set()
        
        for n in self.nodes():
            if self.in_degree(n) == 0:
                producers.add(n)
                tls[n] = 0
            
        for n in self.nodes():
            if n not in producers:
                n_predecs = set(self.predecessors(n))
                if n_predecs <= producers:
                    tls[n] = 1
                else:
                    if self.out_degree(n) == 0:
                        tls[n] = 3
                    else:
                        tls[n] = 2
        
        if has_loops:
            self.add_edges_from(loops)
 
        return tls
    
    def get_aggregated_network(self, sim_threshold):
        net_temp = self.copy()
        
        groups = dict()
        for n,atts in self.nodes(data=True):
            if not groups.has_key(atts['group']):
                groups[atts['group']] = []
            
            groups[atts['group']].append(n)
            
        
        clusters_all = []
        
        for g in groups.keys():
            cl = HierarchicalClustering(groups[g], self.nodes_similarity)
            clusters = cl.getlevel(sim_threshold)
            
            if type(clusters[0]) is str:
                clusters = [clusters]
                
            clusters_all += clusters
            
        new_net = Network()
        for cluster_no in range(len(clusters_all)):
            new_net.add_node(cluster_no)
            
        for i in new_net.nodes():
            current_cluster = clusters_all[i]
            for sp in current_cluster:
                predators = set(self.successors(sp))
                prey = set(self.predecessors(sp))
                
                for j in new_net.nodes():
                    if len(predators & set(clusters_all[j])) > 0:
                        new_net.add_edge(i,j)
                    
                    if len(prey & set(clusters_all[j])) > 0:
                        new_net.add_edge(j,i)
                        
        
        for n,atts in new_net.nodes(data=True):
            atts['group'] = self.node[clusters_all[n][0]]['group']
            atts['biomass'] = self.node[clusters_all[n][0]]['biomass']
            
        return new_net    
            
            
            
    
    def nodes_similarity(self, a, b):
        predators_a = set(self.successors(a))
        preys_a = set(self.predecessors(a))
            
        predators_b = set(self.successors(b))
        preys_b = set(self.predecessors(b))
            
        sim_ab = float((len(predators_a & predators_b) + len(preys_a & preys_b))) / float((len(predators_a | predators_b) + len(preys_a | preys_b))) 
        
        return 1-sim_ab
    
    
    def write_network_files(self, files_prefix):
        #nx.write_edgelist(self, files_prefix+'.edgelist')
        
        nx.write_dot(self,files_prefix+'.dot')
        
        nodes_dict = dict()
        nodes = self.nodes()
        
        biomass = False
        if self.node[nodes[0]].has_key('biomass'):
            headers = ['sp_number','species_name', 'body_mass']
    
            out_file_nodes = open(files_prefix+'.nodes', 'w')
            out_nodes = csv.DictWriter(out_file_nodes, headers, delimiter=' ')
        
            out_nodes.writeheader()
            out_row_nodes = dict()
            
            biomass = True
        
        for i in xrange(len(nodes)):
            nodes_dict[nodes[i]] = i
            
            if biomass:
                out_row_nodes['sp_number'] = i
                out_row_nodes['species_name'] = nodes[i]
                out_row_nodes['body_mass'] = self.node[nodes[i]]['biomass']
                
                out_nodes.writerow(out_row_nodes)
        
        if biomass:
            out_file_nodes.close()
        
        headers = ['prey', 'predator']
    
        out_file_links = open(files_prefix+'.tro', 'w')
        out_links = csv.DictWriter(out_file_links, headers, delimiter=' ')
    
        out_links.writeheader()
        
        out_row_links = dict()
        
        for u,v in self.edges():
            out_row_links['prey'] = nodes_dict[u] 
            out_row_links['predator'] = nodes_dict[v]
            
            out_links.writerow(out_row_links)
        
        out_file_links.close()
    
class NotInvadableNetwork(Exception):
    def __init__(self, value, network):
        self.value = value
        self.net = network
    
    def __str__(self):
        return repr(self.value)
