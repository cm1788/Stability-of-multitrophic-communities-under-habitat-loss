
"""
	File name: ecosystem.py
	Author: Miguel Lurgi Rivera
	Date created: 03/08/2011

	Created in 2011 as part of my PhD dissertation The assembly and disassembly of ecological networks in a changing world
	Submitted to obtain the PhD degree to the Autonomous University of Barcelona
	Part of this work was funded by Microsoft Research

"""

from random import Random
from math import floor, log
from copy import copy
import types
import threading

#from pylab import *
import math
import numpy as np
import scipy.spatial
import random 

from web import Network
from individual import Individual

from configure import ROWS, COLUMNS, MAX_HABITATS, HABITATS, OCCUPIED_CELLS, LOST_HABITAT
from configure import MOVE_RADIUS, MOVE_RADIUS_MUTUALISTS, CAPTURE_PROB, REPRODUCTION_RATE
from configure import INVADER_NUMBER, DOUBLE_HERBIVORY, MATING_SPATIAL_RATIO, INMIGRATION
from configure import HABITAT_LOSS_TYPE, DISPERSAL_KERNEL, SPATIAL_VARIATION, P_HL_CORR


__author__ = ["Miguel Lurgi", "Chris McWilliams"]
__copyright__ = "Copyright 2019"
__credits__ = ["Miguel Lurgi", "Chris McWilliams"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Miguel Lurgi"
__email__ = "miguel.lurgi@swansea.ac.uk"
__status__ = "Production"

class Ecosystem():
    
    def __init__(self, network, drawing=False):
        self.rnd = Random()
        
        self.mutualists = set()
        self.mutualistic_producers = set()
        self.potential_invaders = []
        for n in network.nodes(data=True):
            if network.node[n[0]]['invader']:
                self.potential_invaders.append(n)
                network.remove_node(n[0])
            else:
                if network.node[n[0]]['mut']:
                    self.mutualists.add(n[0])
                if network.node[n[0]]['mut_prod']:
                    self.mutualistic_producers.add(n[0])
            
        self.net = network
        self.world = []
        self.species = self.net.nodes()
        
        self.habitats = range(1,HABITATS+1)
        self.habitats_dist = []    
        self.occupied_habitats = set()
        self.realised_net = Network()
        
        self.drawing = drawing
        
        #if self.drawing:
            #ion()
            #self.fig = figure()
            #self.sp_dist_fig = self.fig.add_subplot(111)
            #self.sp_dist_plot = None
            
            #self.fig_v = figure()
            #self.sp_dist_fig_v = self.fig_v.add_subplot(111)
            #self.sp_dist_plot_v = None            
            
        self.populations = dict.fromkeys(self.species, 0)
        self.species_scl = self.net.get_trophic_levels()
        
        self.new_inds_reproduction = dict.fromkeys(self.species, 0)
        self.new_inds_inmigration = dict.fromkeys(self.species, 0)
        self.dead_individuals = dict.fromkeys(self.species, 0)
        
        self.species_removed = []
        
        if SPATIAL_VARIATION:
            self.centroids = dict.fromkeys(self.species, 0)
            self.areas = dict.fromkeys(self.species, 0)
        
        
    def initialise_world(self, homogeneous=False):
        self._assign_species_habitats(type=2)
        
        b, producers = self.net.basal()
        
        if self.drawing:    
            Z = []
            Z_v = []
            
        if homogeneous:
            original_set = set(copy(self.species))
            original_set -= producers
            ids = list(original_set)     
            
            if DOUBLE_HERBIVORY:
                for n in self.species_scl.keys():
                    if self.species_scl[n] == 1:
                        ids.append(n)
                        
            original_pool = copy(ids)
            prod_list = list(producers)
                        
            for i in range(ROWS):
                row = []
                dist_row = []
                dist_row_v = [] 
                for j in range(COLUMNS):
                    new_cell = Cell()
                    
                    if self.rnd.random() <= OCCUPIED_CELLS:
                        species_id = self.rnd.choice(ids)
                        
                        #this is done to ensure that all species in the pool are sampled
                        #are represented in the world by at least one individual
                        ids.remove(species_id)
                        if len(ids) == 0:
                            ids.extend(original_pool)
                            
                        ind = Individual(species_id)
                        
                        if species_id in self.mutualists:
                            ind.become_mutualist()
                            dist_row.append(max(self.species_scl.values())+1)
                        else:
                            dist_row.append(self.species_scl[species_id])
                         
                        new_cell.visitor = ind
                        new_cell.habitat = self.rnd.choice(sorted(self.net.node[species_id]['habitats']))
                        
                        prod_sp_id = self.rnd.choice(prod_list)
                        
                        while new_cell.habitat not in self.net.node[prod_sp_id]['habitats']:
                            prod_sp_id = self.rnd.choice(prod_list)
                        
                        ind_prod = Individual(prod_sp_id)
                        if prod_sp_id in self.mutualistic_producers:
                            ind_prod.become_mutualistic_producer()
                        new_cell.inhabitant = ind_prod
                    else:
                        #new_cell.habitat = self.rnd.choice(self.habitats)
                        prod_sp_id = self.rnd.choice(prod_list)
                        ind_prod = Individual(prod_sp_id)
                        if prod_sp_id in self.mutualistic_producers:
                            ind_prod.become_mutualistic_producer()
                        new_cell.inhabitant = ind_prod
                        new_cell.habitat = self.rnd.choice(sorted(self.net.node[prod_sp_id]['habitats']))
                        dist_row.append(-1)
                    
                    row.append(new_cell)
                    dist_row_v.append(-1)
                    
                self.world.append(row)
                
                if self.drawing:
                    Z.append(dist_row)
                    Z_v.append(dist_row_v)
                
            self._get_species_per_habitat()
        else:
            self.create_continous_habitats()
            dict_sp_habitats = self._get_species_per_habitat()
            for i in range(ROWS):
                dist_row = []
                dist_row_v = []         
                for j in range(COLUMNS):
                    cell = self.world[i][j]
                    if len(dict_sp_habitats[cell.habitat]) == 0:
                        continue
                    
                    #we populate the first layer of the world with primary producers
                    prod_sp_id = self.rnd.choice(dict_sp_habitats[cell.habitat])
                    while prod_sp_id not in producers:
                        prod_sp_id = self.rnd.choice(dict_sp_habitats[cell.habitat])
                    
                    ind_prod = Individual(prod_sp_id)
                    if prod_sp_id in self.mutualistic_producers:
                            ind_prod.become_mutualistic_producer()
                    cell.inhabitant = ind_prod
                    
                    #we then populate a fraction of the second layer of the world with species other than primary producers
                    if self.rnd.random() <= OCCUPIED_CELLS:                            
                        species_id = self.rnd.choice(dict_sp_habitats[cell.habitat])
                        while species_id in producers:
                            species_id = self.rnd.choice(dict_sp_habitats[cell.habitat])
                        ind = Individual(species_id)
                        
                        if species_id in self.mutualists:
                            ind.become_mutualist()
                            dist_row.append(max(self.species_scl.values())+1)
                        else:
                            dist_row.append(self.species_scl[species_id])
                        cell.visitor = ind
                    else:
                        dist_row.append(-1)
                    
                    dist_row_v.append(-1)
                
                if self.drawing:    
                    Z.append(dist_row)
                    Z_v.append(dist_row_v)
        
        #if self.drawing:   
            #z = array(Z)
            #self.sp_dist_plot = self.sp_dist_fig.pcolor(z, cmap='gist_rainbow', edgecolors='k', linewidths=0)
            
            #z_v = array(Z_v)
            #self.sp_dist_plot_v = self.sp_dist_fig_v.pcolor(z_v, cmap='gist_rainbow', edgecolors='k', linewidths=0)
            
            #self.fig.canvas.draw()
            #self.fig_v.canvas.draw()
        
    def _get_species_per_habitat(self):
        sp_habitats = dict.fromkeys(self.habitats)
        for h in sp_habitats.keys():
            sp_habitats[h] = []
        
        for n,atts in self.net.nodes(data=True):
            for h in atts['habitats']:
                sp_habitats[h].append(n)
        
        self.occupied_habitats.clear()
        for h in sp_habitats.keys():
            if len(sp_habitats[h]) > 0:
                self.occupied_habitats.add(h) 
        
        return sp_habitats
    
    def create_continous_habitats(self):
        self._init_empty_cells()
        cells_per_habitat = floor(float(ROWS*COLUMNS)/len(self.habitats))
        for h in self.habitats:
            cell_count = 0
            x_coord = self.rnd.randint(0,ROWS-1)
            y_coord = self.rnd.randint(0,COLUMNS-1)
            
            while self.world[x_coord][y_coord].habitat != None:
                x_coord = self.rnd.randint(0,ROWS-1)
                y_coord = self.rnd.randint(0,COLUMNS-1)
            
            self.world[x_coord][y_coord].habitat = h
            
            cell_count = 1
            x_count = 1
            y_count = 0
            
            x_next = 1
            y_next = 1
            
            x_offset = -1
            y_offset = 1
            while cell_count < cells_per_habitat:
                if x_count > 0:
                    x_coord = (x_coord+x_offset)%ROWS
                    x_count -= 1
                    
                    if x_count == 0:
                        x_next += 1
                        if x_offset == 1:
                            x_offset = -1
                        elif x_offset == -1:
                            x_offset = 1
                        
                        y_count = y_next
                        
                elif y_count > 0:                
                    y_coord = (y_coord+y_offset)%COLUMNS
                    y_count -= 1
                    
                    if y_count == 0:
                        y_next += 1
                        if y_offset == 1:
                            y_offset = -1
                        elif y_offset == -1:
                            y_offset = 1
                        
                        x_count = x_next
                
                if self.world[x_coord][y_coord].habitat == None:
                    self.world[x_coord][y_coord].habitat = h
                    cell_count += 1
    
#this piece of code shows a plot displaying the initial arrangement of habitats
#in the ecosystem after initialising it using the algorithm above    
#        Z = []
#        for i in range(ROWS):
#            row = []
#            for j in range(COLUMNS):
#                if self.world[i][j].habitat == None:
#                    print 'This cell does not have habitat', i, j 
#                row.append(self.world[i][j].habitat)
#            Z.append(array(row))
#        
#        z = array(Z)
#
#        new_fig = figure()
#        new_plot = new_fig.add_subplot(111)
#        new_plot.pcolor(z, edgecolors='k', linewidths=1)
#        draw()
        
        #show()
    
    
    def apply_habitat_loss(self, type=HABITAT_LOSS_TYPE):
        
	cells_to_loose = floor((ROWS*COLUMNS)*LOST_HABITAT)
	lost_habitat = 0
	
	if type==1:

		x_coord = self.rnd.randint(0,ROWS-1)
		y_coord = self.rnd.randint(0,COLUMNS-1)
		    
		while not self.world[x_coord][y_coord].habitat in self.occupied_habitats:
		    x_coord = self.rnd.randint(0,ROWS-1)
		    y_coord = self.rnd.randint(0,COLUMNS-1)
		
		self.world[x_coord][y_coord].habitat = lost_habitat
		self.world[x_coord][y_coord].inhabitant = None
		self.world[x_coord][y_coord].visitor = None
		    
		cells_lost = 1
		x_count = 1
		y_count = 0
		
		x_next = 1
		y_next = 1
		
		x_offset = -1
		y_offset = 1
		while cells_lost < cells_to_loose:
		    if x_count > 0:
		        x_coord = (x_coord+x_offset)%ROWS
		        x_count -= 1
		        
		        if x_count == 0:
		            x_next += 1
		            if x_offset == 1:
		                x_offset = -1
		            elif x_offset == -1:
		                x_offset = 1
		            
		            y_count = y_next
		            
		    elif y_count > 0:                
		        y_coord = (y_coord+y_offset)%COLUMNS
		        y_count -= 1
		        
		        if y_count == 0:
		            y_next += 1
		            if y_offset == 1:
		                y_offset = -1
		            elif y_offset == -1:
		                y_offset = 1
		            
		            x_count = x_next
		    
		    self.world[x_coord][y_coord].habitat = lost_habitat
		    self.world[x_coord][y_coord].inhabitant = None
		    self.world[x_coord][y_coord].visitor = None
		    cells_lost += 1

	if type==2:

		cells_lost = 0
		#cells_deleted = []
		while cells_lost < cells_to_loose:

			x_coord = self.rnd.randint(0,ROWS-1)
			y_coord = self.rnd.randint(0,COLUMNS-1)
			    
			while not self.world[x_coord][y_coord].habitat in self.occupied_habitats or self.world[x_coord][y_coord].habitat==lost_habitat:#(x_coord,y_coord) in cells_deleted:
			    x_coord = self.rnd.randint(0,ROWS-1)
			    y_coord = self.rnd.randint(0,COLUMNS-1)
		
			self.world[x_coord][y_coord].habitat = lost_habitat
			self.world[x_coord][y_coord].inhabitant = None
			self.world[x_coord][y_coord].visitor = None

			cells_lost += 1
			#cells_deleted.append((x_coord,y_coord))
	if type==3:

		cells_lost = 0
		while cells_lost < cells_to_loose:

			x_coord = self.rnd.randint(0,ROWS-1)
			y_coord = self.rnd.randint(0,COLUMNS-1)
			    
			while not self.world[x_coord][y_coord].habitat in self.occupied_habitats or self.world[x_coord][y_coord].habitat==lost_habitat:
			    x_coord = self.rnd.randint(0,ROWS-1)
			    y_coord = self.rnd.randint(0,COLUMNS-1)
		
			self.world[x_coord][y_coord].habitat = lost_habitat
			self.world[x_coord][y_coord].inhabitant = None
			self.world[x_coord][y_coord].visitor = None

			cells_lost += 1
			x_count = 1
			y_count = 0
		
			x_next = 1
			y_next = 1
		
			x_offset = -1
			y_offset = 1
			
			while random.random() < P_HL_CORR and cells_lost < cells_to_loose:
				## do contiguous:
			
			        if x_count > 0:
		        		x_coord = (x_coord+x_offset)%ROWS
				        x_count -= 1
		        
			    	    	if x_count == 0:
				            x_next += 1
		           		    if x_offset == 1:
				                x_offset = -1
				            elif x_offset == -1:
				                x_offset = 1
		            
				            y_count = y_next
		            
			        elif y_count > 0:                
		        		y_coord = (y_coord+y_offset)%COLUMNS
				        y_count -= 1
		        
				        if y_count == 0:
				            y_next += 1
				            if y_offset == 1:
				                y_offset = -1
				            elif y_offset == -1:
				                y_offset = 1
		            
				            x_count = x_next
		    
			        self.world[x_coord][y_coord].habitat = lost_habitat
		    		self.world[x_coord][y_coord].inhabitant = None
				self.world[x_coord][y_coord].visitor = None
				cells_lost += 1
		    

#the following piece of codes displays a figure showing what happens to the ecosystem
#habitats after the habitat loss event implement using the algorithm above
#        Z = []
#        for i in range(ROWS):
#            row = []
#            for j in range(COLUMNS):
#                if self.world[i][j].habitat == None:
#                    print 'This cell does not have habitat', i, j 
#                row.append(self.world[i][j].habitat)
#            Z.append(array(row))
#        
#        z = array(Z)
#        
#    
#        new_fig = plt.figure()
#        new_plot = new_fig.add_subplot(111)
#        new_plot.pcolor(z, edgecolors='k', linewidths=1)
#        plt.draw()
#        #show()
    
    
       
    #def draw_species_distribution(self):
    #    draw_thread = threading.Thread(target=self.threaded_drawing_sp_dist, args=(self.net, self.sp_dist_plot, self.world, self.species_scl))
    #    draw_thread.start()
        #draw_thread.join()
        
    #    draw_thread2 = threading.Thread(target=self.threaded_drawing_sp_dist, args=(self.net, self.sp_dist_plot_v, self.world, self.species_scl, True))
    #    draw_thread2.start()
        
    #    self.fig.canvas.draw()
    #    self.fig_v.canvas.draw()
    
    #def threaded_drawing_sp_dist(self, net, plot, world, scl_rank, visitor=False):
        #Z = []
        #for i in range(ROWS):
        #    row = []
        #    for j in range(COLUMNS):
        #        if visitor:
        #            if world[i][j].habitat == 0:
        #                row.append(-2)
        #            elif world[i][j].visitor == None:
        #                row.append(-1)
        #            else:
        #                if net.node[world[i][j].visitor.species_id].has_key('invader') and net.node[world[i][j].visitor.species_id]['invader']:
        #                    row.append(max(scl_rank.values()) + 1)
        #                elif world[i][j].visitor.species_id in self.mutualists:
        #                    row.append(max(scl_rank.values()) + 1)
        #                else:
        #                    row.append(scl_rank[world[i][j].visitor.species_id])
        #        else:
        #            if world[i][j].habitat == 0:
        #                row.append(-2)
        #            elif world[i][j].inhabitant == None:
        #                row.append(-1)
        #            else:
        #                if net.node[world[i][j].inhabitant.species_id].has_key('invader') and net.node[world[i][j].inhabitant.species_id]['invader']:
        #                    row.append(max(scl_rank.values()) + 1)
        #                elif world[i][j].inhabitant.species_id in self.mutualists:
        #                    row.append(max(scl_rank.values()) + 1)
        #                elif world[i][j].inhabitant.species_id in self.mutualistic_producers:
        #                    row.append(max(scl_rank.values()) + 2)
        #                else:
        #                    row.append(scl_rank[world[i][j].inhabitant.species_id])
        #    Z.append(array(row))
        
        #z = array(Z)
        
        #z[1][1] = -2
        #z[1][2] = -1
        #z[1][3] = 0
        #z[1][4] = 1
        #z[1][5] = 2
        #z[1][6] = 3
        #z[1][7] = 4
        #z[1][8] = 5
        
        
        
        #plot.set_array(z.ravel())
        #plot.autoscale()
        
        
    def _init_empty_cells(self):
        self.world = []
        for i in range(ROWS):
            row = []
            for j in range(COLUMNS):
                new_cell = Cell()
                row.append(new_cell)
            self.world.append(row)
    
    def inmigration(self, cell):
        inmigrant_sp = self.rnd.choice(self.species)
        while inmigrant_sp in self.species_removed:
            inmigrant_sp = self.rnd.choice(self.species)
        
        sp_habitats = self.net.node[inmigrant_sp]['habitats']
          
        #if the habitat available in the cell is not one of the species' do nothing
        if not cell.habitat in sp_habitats:
            return False
        
        inmigrant = Individual(inmigrant_sp)
        if inmigrant_sp in self.mutualists:
            inmigrant.become_mutualist()
        elif inmigrant_sp in self.mutualistic_producers:
            inmigrant.become_mutualistic_producer()
        
        if cell.inhabitant == None:
            cell.inhabitant = inmigrant
            self.new_inds_inmigration[inmigrant_sp] += 1
            return True
        else:
            if self.net.in_degree(inmigrant_sp) == 0:
                return False
            
            current_indiv = cell.inhabitant
            
            if current_indiv.species_id == inmigrant_sp:
                return False
            
            if (current_indiv.species_id, inmigrant_sp) in self.net.edges():
                if self.net.in_degree(current_indiv.species_id) == 0:
                    if cell.visitor == None:
                        if self.species_scl[inmigrant_sp] > 1:
                            inmigrant.eat(current_indiv, herbivorous=True, omnivore=True)
                        else:
                            inmigrant.eat(current_indiv, herbivorous=True)
                            
                        cell.visitor = inmigrant
                        self._update_realised_network(current_indiv, inmigrant)
                        
                        if inmigrant.mutualist:
                            inmigrant.set_mutualistic_partner(current_indiv)
                        
                        self.new_inds_inmigration[inmigrant_sp] += 1
                        return True
                else:
                    if self.rnd.random() < CAPTURE_PROB:
                        inmigrant.eat(current_indiv)
                        cell.inhabitant = inmigrant
                        self._update_realised_network(current_indiv, inmigrant)
                        
                        self.dead_individuals[current_indiv.species_id] += 1
                        self.new_inds_inmigration[inmigrant_sp] += 1
                        
                        return True
#            elif (current_indiv.species_id, second_indiv.species_id) in self.net.edges():
#                if self.rnd.random() < CAPTURE_PROB:
#                    second_indiv.eat(current_indiv)
#                    self._update_realised_network(current_indiv, second_indiv)
#                else:
#                    return
        if self.net.in_degree(inmigrant_sp) != 0:
            if cell.visitor == None:
                cell.visitor = inmigrant
                self.new_inds_inmigration[inmigrant_sp] += 1
                return True
            else:
                current_indiv = cell.visitor
                
                if current_indiv.species_id == inmigrant_sp:
                    return False
                
                if (current_indiv.species_id, inmigrant_sp) in self.net.edges() and self.rnd.random() < CAPTURE_PROB:
                    inmigrant.eat(current_indiv)
                    cell.visitor = inmigrant
                    self._update_realised_network(current_indiv, inmigrant)
                    
                    self.dead_individuals[current_indiv.species_id] += 1
                    self.new_inds_inmigration[inmigrant_sp] += 1
                    
                    return True
            
        
    def update_world(self):
        idx_col = self.rnd.randint(0,COLUMNS-1)
        init_idx_col = idx_col
        idx_row = self.rnd.randint(0,ROWS-1)
        init_idx_row = idx_row
        
        self.new_inds_reproduction = dict.fromkeys(self.species, 0)
        self.new_inds_inmigration = dict.fromkeys(self.species, 0)
        self.dead_individuals = dict.fromkeys(self.species, 0)
        
        row_count = 0
        for i in range(idx_row, ROWS):
            if row_count > 0:
                idx_col = 0
            for j in range(idx_col, COLUMNS):
                current_cell = self.world[i][j]
                
                if self.rnd.random() < INMIGRATION:
                    self.inmigration(current_cell)    
                
                if current_cell.inhabitant == None and current_cell.visitor == None:
                    continue
                else:
                    if current_cell.visitor != None:
                        if current_cell.inhabitant == None:
                            current_cell.inhabitant = current_cell.visitor
                            current_cell.visitor = None
                        else:
                            self._move_individual(i, j, True)
                        
                    self._move_individual(i, j)
            
            row_count += 1
        
        row_count = 0
        end_col = COLUMNS
        for i in range(0, idx_row+1):
            if row_count == init_idx_row:
                end_col = init_idx_col
            for j in range(0, end_col):
                current_cell = self.world[i][j]
                
                if self.rnd.random() < INMIGRATION:
                    self.inmigration(current_cell)
                
                if current_cell.inhabitant == None and current_cell.visitor == None:
                    continue
                else:
                    if current_cell.visitor != None:
                        if current_cell.inhabitant == None:
                            current_cell.inhabitant = current_cell.visitor
                            current_cell.visitor = None
                        else:
                            self._move_individual(i, j, True)
                        
                    self._move_individual(i, j)
            
            row_count += 1
                     
    
    def _move_individual(self, i, j, visitor=False):
        current_cell = self.world[i][j]
        
        if visitor:
            current_indiv = current_cell.visitor
        else:
            current_indiv = current_cell.inhabitant
        
        producer = False
        if self.net.in_degree(current_indiv.species_id) == 0:
            if visitor:
                print 'there is a producer in a visitor spot'
            producer = True
        
        if current_indiv.live(producer) == False:
            self.dead_individuals[current_indiv.species_id] += 1
            if visitor:
                current_cell.visitor = None
            else:
                current_cell.inhabitant = None
                if current_cell.visitor != None:
                    current_cell.inhabitant = current_cell.visitor
                    current_cell.visitor = None
            return
        
        #if the individual is a mutualist cool off its efficiency
        if current_indiv.mutualist and current_indiv.current_host != None:
            current_indiv.mutualistic_cool_off()
            
            #... and if it is in a cell with an empty space for producers... reproduce is mutualistic host
            if not visitor and current_cell.visitor == None and self.rnd.random() < current_indiv.mut_efficiency and current_indiv.current_host != None:
                current_cell.visitor = current_indiv
                
                current_cell.inhabitant = Individual(current_indiv.current_host)
                self.new_inds_reproduction[current_indiv.current_host] += 1
                
                current_cell.inhabitant.become_mutualistic_producer()
                current_indiv.reset_mutualistic_state()  
                
                return
        
        #this is to state whether the current individual is a plant that can auto reproduce
        auto_reproductive = False
        #if the individual is a primary producer it cannot move, so, continue...
        if producer:
            current_indiv.synthesis()   #... but it must feed
            if not current_indiv.mutualistic_producer:
                auto_reproductive = True
            else:
                return
        else:
            #here animals can reproduce. If the current individual do reproduces, then returns (it doesn't do anything else)
            if self.sexual_reproduction(current_indiv, i, j):
#                if visitor:
#                    print 'visitor reproducing', current_indiv.species_id
                return
            
        if current_indiv.mutualist:
            new_idx_x = self.rnd.randint(-MOVE_RADIUS_MUTUALISTS, MOVE_RADIUS_MUTUALISTS)
            new_idx_y = self.rnd.randint(-MOVE_RADIUS_MUTUALISTS, MOVE_RADIUS_MUTUALISTS)
        else:
            new_idx_x = self.rnd.randint(-MOVE_RADIUS, MOVE_RADIUS)
            new_idx_y = self.rnd.randint(-MOVE_RADIUS, MOVE_RADIUS)
            
        new_cell_x = (i+new_idx_x)%ROWS
        new_cell_y = (j+new_idx_y)%COLUMNS
        
        #if the cell is the same do nothing (stay)
        if new_cell_x == i and new_cell_y == j:
            return
        
        new_cell = self.world[new_cell_x][new_cell_y]
        sp_habitats = self.net.node[current_indiv.species_id]['habitats']
        
        #if the habitat available in the cell is not one of the individuals' do nothing
        if not new_cell.habitat in sp_habitats:
            return
        
        if new_cell.inhabitant == None and new_cell.visitor != None:
            new_cell.inhabitant = new_cell.visitor
            new_cell.visitor = None
        
        #this is where the reproduction of primary producers (plants) occur
        #if the current individual is a plant that can auto reproduce (wind dispersal) then...
        if auto_reproductive:
            if self.rnd.random() < REPRODUCTION_RATE:
                if new_cell.inhabitant == None:
                    if current_indiv.mutualist:
                        print 'I am a mutualist auto reproducing'
                    new_cell.inhabitant = Individual(current_indiv.species_id)
                    self.new_inds_reproduction[current_indiv.species_id] += 1
                elif new_cell.visitor == None and self.species_scl[new_cell.inhabitant.species_id] > 0:
                    new_cell.visitor = new_cell.inhabitant
                    new_cell.inhabitant = Individual(current_indiv.species_id)
                    self.new_inds_reproduction[current_indiv.species_id] += 1
                
            return
        
        if new_cell.inhabitant == None:
            #if the current individual is a mutualist moving to an empty cell it can (depending on its
            #mutualistic efficiency and the availability of space) facilitate the creation of a new
            #individual of its previous host
            if new_cell.visitor == None and current_indiv.mutualist and current_indiv.current_host != None and new_cell.habitat in self.net.node[current_indiv.current_host]['habitats'] and self.rnd.random() < current_indiv.mut_efficiency:
                new_cell.inhabitant = Individual(current_indiv.current_host)
                
                self.new_inds_reproduction[current_indiv.current_host] += 1
                
                new_cell.inhabitant.become_mutualistic_producer()
                current_indiv.reset_mutualistic_state()
                new_cell.visitor = current_indiv         
            else:
                new_cell.inhabitant = current_indiv
            
            if visitor:
                current_cell.visitor = None
            else:
                current_cell.inhabitant = None
        else:
            
            if new_cell.visitor != None:
                #no interactions involving plants are possible in the following case
                second_indiv = new_cell.visitor
                if (current_indiv.species_id, second_indiv.species_id) in self.net.edges():
                    if self.rnd.random() < CAPTURE_PROB:
                        #print 'this is a carnivorous link between prey', current_indiv.species_id, 'and predator', second_indiv.species_id, 'new cell visitor not empty'
                        second_indiv.eat(current_indiv)
                        self._update_realised_network(current_indiv, second_indiv)
                        self.dead_individuals[current_indiv.species_id] += 1
                        
                        if visitor:
                            current_cell.visitor = None
                        else:
                            current_cell.inhabitant = None
                            
                elif (second_indiv.species_id, current_indiv.species_id) in self.net.edges():
                    if self.rnd.random() < CAPTURE_PROB:
                        #print 'this is a carnivorous link between prey', second_indiv.species_id, 'and predator', current_indiv.species_id, 'the inhabitant eats the newcomer'
                        current_indiv.eat(second_indiv)
                        new_cell.visitor = None
                        new_cell.visitor = current_indiv
                        self._update_realised_network(second_indiv, current_indiv)
                        self.dead_individuals[second_indiv.species_id] += 1
                        if visitor:
                            current_cell.visitor = None
                        else:
                            current_cell.inhabitant = None

            else:
                second_indiv = new_cell.inhabitant
                
                #if the individuals belong to the same species they can either mate or,
                #if they share a connection, one can eat the other 
                if second_indiv.species_id == current_indiv.species_id:
                    if (second_indiv.species_id, current_indiv.species_id) in self.net.edges():
                        if self.rnd.random() < CAPTURE_PROB:
                            #print 'this is a cannibalistic link and hence a predator-prey interaction between', current_indiv.species_id
                            if self.rnd.random() < 0.5:
                                current_indiv.eat(second_indiv)
                                new_cell.inhabitant = current_indiv
                                
                                self._update_realised_network(second_indiv, current_indiv)
                                
                                self.dead_individuals[second_indiv.species_id] += 1
                            else:
                                second_indiv.eat(current_indiv)
                                
                                self._update_realised_network(current_indiv, second_indiv)
                                self.dead_individuals[current_indiv.species_id] += 1
                            if visitor:
                                current_cell.visitor = None
                            else:
                                current_cell.inhabitant = None
                
                elif (second_indiv.species_id, current_indiv.species_id) in self.net.edges():
                    #this is the only case in which the second individual may be a primary producer
                    
                    #if the second individual is a primary producer it remains alive, although losing some resource
                    #and if the current individual is a mutualist it will record the second individual information
                    if self.net.in_degree(second_indiv.species_id) == 0:
                        
                        #print 'this is an herbivorous link between producer', second_indiv.species_id, 'and herbivore', current_indiv.species_id
                        if new_cell.visitor == None:
                            if self.species_scl[current_indiv.species_id] > 1:
                                current_indiv.eat(second_indiv, herbivorous=True, omnivore=True)
                            else:
                                current_indiv.eat(second_indiv, herbivorous=True)
                                
                            new_cell.visitor = current_indiv
                            self._update_realised_network(second_indiv, current_indiv)
                            
                            if current_indiv.mutualist:
                                current_indiv.set_mutualistic_partner(second_indiv)
                        else:
                            return
                    else:
                        #print 'this is a carnivorous link between prey', second_indiv.species_id, 'and predator', current_indiv.species_id, '(the visitor is the predator)'
                        if self.rnd.random() < CAPTURE_PROB:
                            current_indiv.eat(second_indiv)
                            new_cell.inhabitant = current_indiv
                            self._update_realised_network(second_indiv, current_indiv)
                            self.dead_individuals[second_indiv.species_id] += 1
                        else:
                            return
                            
                    if visitor:
                        current_cell.visitor = None
                    else:
                        current_cell.inhabitant = None
                
                elif (current_indiv.species_id, second_indiv.species_id) in self.net.edges():
                    
                    if self.net.in_degree(current_indiv.species_id) == 0 and self.species_scl[second_indiv.species_id] > 1:
                        print 'omnivore predator eating without paying omnivore penalty'
                        
                    
                    if self.rnd.random() < CAPTURE_PROB:
                        
                        #print 'this is a carnivorous link between prey', current_indiv.species_id, 'and predator', second_indiv.species_id, 'the host is the predator'
                        
                        second_indiv.eat(current_indiv)
                        self._update_realised_network(current_indiv, second_indiv)
                        self.dead_individuals[current_indiv.species_id] += 1
                    else:
                        return
                    
                    if visitor:
                        current_cell.visitor = None
                    else:
                        current_cell.inhabitant = None
                
                else:
                    if new_cell.visitor == None and self.species_scl[current_indiv.species_id] != 0:
                        new_cell.visitor = current_indiv
                        if visitor:
                            current_cell.visitor = None
                        else:
                            current_cell.inhabitant = None
                        
               
                    
    
    def sexual_reproduction(self, individual, i, j):
        if not individual.ready_to_mate():
            return False
         
        x = (i-MATING_SPATIAL_RATIO)%ROWS
        start_y = (j-MATING_SPATIAL_RATIO)%COLUMNS
        
        mating_cell = None
        partner = None
        
        for x_offset in range(MATING_SPATIAL_RATIO*2):
            x = (x+1)%ROWS
            y = start_y
            for y_offset in range(MATING_SPATIAL_RATIO*2):
                y = (y+1)%COLUMNS
                
                if x == i and y == j:
                    continue
                
                temp_cell = self.world[x][y]
                
                if mating_cell == None and (temp_cell.inhabitant == None or temp_cell.visitor == None) and temp_cell.habitat in self.net.node[individual.species_id]['habitats']:
                    mating_cell = temp_cell
                
                if partner == None:
                    if temp_cell.inhabitant != None and temp_cell.inhabitant.species_id == individual.species_id and temp_cell.inhabitant.ready_to_mate():
                        partner = temp_cell.inhabitant
                    elif temp_cell.visitor != None and temp_cell.visitor.species_id == individual.species_id and temp_cell.visitor.ready_to_mate():
                        partner = temp_cell.visitor
#                        print 'visitor chosen for reproduction', partner.species_id
                        
                if mating_cell != None and partner != None:
                    break    
        
        if mating_cell != None and partner != None:
            newborn = Individual(individual.species_id)
            self.new_inds_reproduction[individual.species_id] += 1
            if individual.species_id in self.mutualists:
                newborn.become_mutualist()
            newborn.resource = (individual.reproduce() + partner.reproduce())*2
            
            if mating_cell.inhabitant == None:
                mating_cell.inhabitant = newborn
            else:
                mating_cell.visitor = newborn
            
            return True
        else:
            return False
        
        
                       
    def _create_habitats_distribution(self):
        largest_r = None
        smallest_r = None
        for n, atts in self.net.nodes(data=True):
            if largest_r == None or atts['r'] > largest_r:
                largest_r = atts['r']
            if smallest_r == None or atts['r'] < smallest_r:
                smallest_r = atts['r']

        self.habitats_dist = []
        self.habitats_dist.append(smallest_r)
        for i in range(1,MAX_HABITATS):
            self.habitats_dist.append( self.habitats_dist[i-1] + ((largest_r-self.habitats_dist[i-1])/2) )
        
        self.habitats_dist.append(largest_r)
        

    def _assign_species_habitats(self, type=1):
        self._create_habitats_distribution()
        if type == 1:
            for n in self.net.nodes():            
                current_r = self.net.node[n]['r']
                habitats_no = 1
                for i in range(len(self.habitats_dist)):
                    if current_r >= self.habitats_dist[i] and current_r <= self.habitats_dist[i+1]:
                        break
                    habitats_no += 1 
                
                self.rnd.shuffle(self.habitats)
                self.net.node[n]['habitats'] = set()
                for i in range(habitats_no):
                    self.net.node[n]['habitats'].add(self.habitats[i])
        elif type == 2:
            producers = set()
            for n in self.net.nodes():
                if self.net.in_degree(n) == 0:
                    current_r = self.net.node[n]['r']
                    habitats_no = 1
                    for i in range(len(self.habitats_dist)):
                        if current_r >= self.habitats_dist[i] and current_r <= self.habitats_dist[i+1]:
                            break
                        habitats_no += 1
                    
                    self.net.node[n]['habitats'] = set()
                    if habitats_no == 1:
                        self.net.node[n]['habitats'].add(self.rnd.choice(self.habitats))
                    else:    
                        self.rnd.shuffle(self.habitats)
                        for i in range(habitats_no):
                            self.net.node[n]['habitats'].add(self.habitats[i])
        
                    producers.add(n)
            
            
            nodes = set(self.net.nodes())
            nodes -= producers
            
            nodes_sorted = sorted(nodes, self.compare_nodes_niches)
            for n in nodes_sorted:
                if not self.net.node[n].has_key('habitats'):
                    self._assign_node_habitats(n)
            
         
    def _assign_node_habitats(self, n):
        self.net.node[n]['habitats'] = set()
        predecessors = sorted(self.net.predecessors(n), self.compare_nodes_niches)
        for pre in predecessors:
            if n == pre:
                continue
            
            if not self.net.node[pre].has_key('habitats'):
                self._assign_node_habitats(pre)

            self.net.node[n]['habitats'] |= self.net.node[pre]['habitats']
    
    
    def _update_realised_network(self, prey, predator):
        if not prey.species_id in self.realised_net.nodes():
            node_data = self.net.node[prey.species_id]
            self.realised_net.add_node(prey.species_id, attr_dict=node_data)
        
        if not predator.species_id in self.realised_net.nodes():
            node_data = self.net.node[predator.species_id]
            self.realised_net.add_node(predator.species_id, attr_dict=node_data)
    
        if not (prey.species_id, predator.species_id) in self.realised_net.edges():
            self.realised_net.add_edge(prey.species_id, predator.species_id)
            self.realised_net[prey.species_id][predator.species_id]['is'] = 1
        else:
            self.realised_net[prey.species_id][predator.species_id]['is'] += 1
    
    
    def clear_realised_network(self):
        self.realised_net.clear()
    
    
    def count_individuals(self, spatial=False):
        self.populations = dict.fromkeys(self.species, 0)
        if spatial:
            positions = dict.fromkeys(self.net.nodes(), None)
        
        for i in range(ROWS):
            for j in range(COLUMNS):
                cell = self.world[i][j]
                if cell.inhabitant is not None and type(cell.inhabitant) is not types.NoneType:
                    try:
                        self.populations[cell.inhabitant.species_id] += 1
                        
                        if spatial:
                            if not(positions[cell.inhabitant.species_id] is None):
                                positions[cell.inhabitant.species_id].append((i,j))
                            else:
                                positions[cell.inhabitant.species_id] = [(i,j)]
                        
                    except:
                        print 'inhabitant in cell ' + str(i) + ',' + str(j) + ' declared None but is:'
                        print type(cell.inhabitant)                    
                if cell.visitor is not None and type(cell.visitor) is not types.NoneType:                    
                    try:
                        self.populations[cell.visitor.species_id] += 1
                        
                        if spatial:
                            if not(positions[cell.visitor.species_id] is None):
                                positions[cell.visitor.species_id].append((i,j))
                            else:
                                positions[cell.visitor.species_id] = [(i,j)]
                    except:
                        print 'visitor in cell ' + str(i) + ',' + str(j) + ' declared None but is:'
                        print type(cell.visitor)      
        
        if spatial:
            for s in self.species:         
                if (not(positions[s] is None)) and np.size(positions[s], axis=0) > 0:
                    
                    self.centroids[s] = np.average(positions[s], axis=0)
                    
                    distances = scipy.spatial.distance.cdist(positions[s], [self.centroids[s]])
                    to_remove = math.floor(self.populations[s]*DISPERSAL_KERNEL)
                    
                    sorted_dists = np.sort(distances,axis=0)
                    length = np.size(distances)
                    length = int(length-to_remove)
                    
                    if length > 10:
                        radius = sorted_dists[length-1][0]
                        self.areas[s] = math.pi * (radius**2) 
                        #['frac_occ'] = species_autocorr[s]['area'] / (self.populations[s] - to_remove)
                
        return self.populations
    
    def compare_nodes_niches(self, a, b):
        n_a = self.net.node[a]['n'] 
        min_range_a = self.net.node[a]['c'] - (self.net.node[a]['r']/2)
        max_range_a = self.net.node[a]['c'] + (self.net.node[a]['r']/2)
        
        n_b = self.net.node[b]['n']
        min_range_b = self.net.node[b]['c'] - (self.net.node[b]['r']/2)
        max_range_b = self.net.node[b]['c'] + (self.net.node[b]['r']/2)
        
        
        if n_a >= min_range_b and n_a <= max_range_b and n_b >= min_range_a and n_b <= max_range_a:
            if min_range_a > min_range_b:
                return 1
            elif min_range_b > min_range_a:
                return -1
        elif n_a >= min_range_b and n_a <= max_range_b:
            return -1
        elif n_b >= min_range_a and n_b <= max_range_a:
            return 1
        else:
            if n_a > n_b:
                return 1
            elif n_b > n_a:
                return -1
            else:
                return 0
    
    def invade(self):
        self.rnd.shuffle(self.potential_invaders)
        count = len(self.potential_invaders)
        seen = 0
        
        inv_id = max(self.net.nodes()) + 1
        while seen < count:
            invader = self.potential_invaders[seen]
            self.net.add_node(inv_id, attr_dict=invader)
            
            for i in self.net.nodes():
                if self.net.node[i]['r'] == 0.0:
                    continue
                
                r_lower_bound = self.net.node[i]['c'] - (self.net.node[i]['r']/2)
                r_upper_bound = self.net.node[i]['c'] + (self.net.node[i]['r']/2) 
                for j in self.net.nodes():
                    if i != inv_id and j != inv_id:
                        continue
                    
                    if self.net.node[j]['n'] >= r_lower_bound and self.net.node[j]['n'] <= r_upper_bound:  
                        self.net.add_edge(j,i)
        
            
            disconnected = False
            
            if self.net.degree(inv_id) == 0:
                disconnected = True
            elif self.net.in_degree(inv_id) == 1 and inv_id in self.net.nodes_with_selfloops():
                disconnected = True
            else:
                i_succs = set(self.net.successors(inv_id)) - set([inv_id])
                i_predecs = set(self.net.predecessors(inv_id)) - set([inv_id])
                for j in self.net.nodes():
                    if j == inv_id:
                        continue
                    j_succs = set(self.net.successors(j)) - set([j])
                    j_predecs = set(self.net.predecessors(j)) - set([j])
                    
                    if i_succs == j_succs and i_predecs == j_predecs:
                        disconnected = True
                        break
            
            if not disconnected:
                self._assign_node_habitats(inv_id)
                for sp in self.net.successors(inv_id):
                    self.net.node[sp]['habitats'] |= self.net.node[inv_id]['habitats']
                break
            
            self.net.remove_node(inv_id)
            seen += 1
            
        if seen == count:
            print 'None of the possible invaders can invade this network due to lack of links'
            return
        
        invaders_count = 0
        for i in range(INVADER_NUMBER):
            individual = Individual(inv_id)
            inv_x, inv_y = self._get_random_empty_cell(self.net.node[inv_id]['habitats'])
            
            if inv_x == -1 and inv_y == -1:
                print 'no empty cell found'
                continue
            else:
                self.world[inv_x][inv_y].inhabitant = individual
            
            invaders_count += 1
        
        if invaders_count > 0: 
            self.species_scl = self.net.get_trophic_levels()
            print 'ecosystem invaded', invaders_count, 'individuals'
        else:
            print 'There are no empty cells for this species to invade the ecosystem' 
        
    
    def _get_random_empty_cell(self, habitats):
        count = 0
        x_coord = self.rnd.randint(0,ROWS-1)
        y_coord = self.rnd.randint(0,COLUMNS-1)
            
        while not self.world[x_coord][y_coord].habitat in habitats or self.world[x_coord][y_coord].inhabitant != None or count < ROWS*COLUMNS:
            x_coord = self.rnd.randint(0,ROWS-1)
            y_coord = self.rnd.randint(0,COLUMNS-1)
            count += 1
            
        if count == ROWS*COLUMNS:
            return -1,-1
        else:
            return x_coord, y_coord 
        
    def get_groups_counts(self, spatial=False):
        producer_sp = mut_prod_sp = herbivore_sp = mutualist_sp = primary_pred_sp = secondary_pred_sp = 0
        
        producer_count = mut_prod_count = herbivore_count = mutualist_count = primary_pred_count = secondary_pred_count = 0
        
        producer_new_reprod = mut_prod_new_reprod = herbivore_new_reprod = mutualist_new_reprod = primary_pred_new_reprod = secondary_pred_new_reprod = 0
        
        producer_new_inmig = mut_prod_new_inmig = herbivore_new_inmig = mutualist_new_inmig = primary_pred_new_inmig = secondary_pred_new_inmig = 0
        
        producer_dead = mut_prod_dead = herbivore_dead = mutualist_dead = primary_pred_dead = secondary_pred_dead = 0
        
        self.count_individuals(spatial)
        
        for k in self.species_scl.keys():
            if self.populations[k] > 0:
                if self.species_scl[k] == 0:
                    if k in self.mutualistic_producers:
                        mut_prod_count += self.populations[k]
                        mut_prod_sp += 1
                    else:
                        producer_count += self.populations[k]
                        producer_sp += 1
                elif self.species_scl[k] == 1:
                    if k in self.mutualists:
                        mutualist_count += self.populations[k]
                        mutualist_sp += 1
                    else:
                        herbivore_count += self.populations[k]
                        herbivore_sp += 1
                elif self.species_scl[k] == 2:
                    primary_pred_count += self.populations[k]
                    primary_pred_sp += 1
                elif self.species_scl[k] == 3:
                    secondary_pred_count += self.populations[k]
                    secondary_pred_sp += 1
            #couting the new individuals due to reproduction in each trophic level
            if self.new_inds_reproduction[k] > 0:
                if self.species_scl[k] == 0:
                    if k in self.mutualistic_producers:
                        mut_prod_new_reprod += self.new_inds_reproduction[k]
                    else:
                        producer_new_reprod += self.new_inds_reproduction[k]
                elif self.species_scl[k] == 1:
                    if k in self.mutualists:
                        mutualist_new_reprod += self.new_inds_reproduction[k]
                    else:
                        herbivore_new_reprod += self.new_inds_reproduction[k]
                elif self.species_scl[k] == 2:
                    primary_pred_new_reprod += self.new_inds_reproduction[k]
                elif self.species_scl[k] == 3:
                    secondary_pred_new_reprod += self.new_inds_reproduction[k]
            #couting the new individuals due to inmigration in each trophic level
            if self.new_inds_inmigration[k] > 0:
                if self.species_scl[k] == 0:
                    if k in self.mutualistic_producers:
                        mut_prod_new_inmig += self.new_inds_inmigration[k]
                    else:
                        producer_new_inmig += self.new_inds_inmigration[k]
                elif self.species_scl[k] == 1:
                    if k in self.mutualists:
                        mutualist_new_inmig += self.new_inds_inmigration[k]
                    else:
                        herbivore_new_inmig += self.new_inds_inmigration[k]
                elif self.species_scl[k] == 2:
                    primary_pred_new_inmig += self.new_inds_inmigration[k]
                elif self.species_scl[k] == 3:
                    secondary_pred_new_inmig += self.new_inds_inmigration[k]
            #counting dead individuals
            if self.dead_individuals[k] > 0:
                if self.species_scl[k] == 0:
                    if k in self.mutualistic_producers:
                        mut_prod_dead += self.dead_individuals[k]
                    else:
                        producer_dead += self.dead_individuals[k]
                elif self.species_scl[k] == 1:
                    if k in self.mutualists:
                        mutualist_dead += self.dead_individuals[k]
                    else:
                        herbivore_dead += self.dead_individuals[k]
                elif self.species_scl[k] == 2:
                    primary_pred_dead += self.dead_individuals[k]
                elif self.species_scl[k] == 3:
                    secondary_pred_dead += self.dead_individuals[k]


        total_sp = producer_sp + mut_prod_sp + herbivore_sp + mutualist_sp + primary_pred_sp + secondary_pred_sp
        total_count = producer_count + mut_prod_count + herbivore_count + mutualist_count + primary_pred_count + secondary_pred_count
        groups_counts = dict()
        groups_counts['total'] = (total_sp, total_count)
        groups_counts['prods'] = (producer_sp, producer_count)
        groups_counts['mut_prods'] = (mut_prod_sp, mut_prod_count)
        groups_counts['herbs'] = (herbivore_sp, herbivore_count)
        groups_counts['muts'] = (mutualist_sp, mutualist_count)
        groups_counts['prim_preds'] = (primary_pred_sp, primary_pred_count)
        groups_counts['second_preds'] = (secondary_pred_sp, secondary_pred_count)
        
        groups_counts['prods_new'] = (producer_new_reprod, producer_new_inmig)
        groups_counts['mut_prods_new'] = (mut_prod_new_reprod, mut_prod_new_inmig)
        groups_counts['herbs_new'] = (herbivore_new_reprod, herbivore_new_inmig)
        groups_counts['muts_new'] = (mutualist_new_reprod, mutualist_new_inmig)
        groups_counts['prim_preds_new'] = (primary_pred_new_reprod, primary_pred_new_inmig)
        groups_counts['second_preds_new'] = (secondary_pred_new_reprod, secondary_pred_new_inmig)
        
        groups_counts['prods_dead'] = producer_dead
        groups_counts['mut_prods_dead'] = mut_prod_dead
        groups_counts['herbs_dead'] = herbivore_dead
        groups_counts['muts_dead'] = mutualist_dead
        groups_counts['prim_preds_dead'] = primary_pred_dead
        groups_counts['second_preds_dead'] = secondary_pred_dead
                
        return groups_counts
        
    def get_shannon_biodiversity_index(self, groups_count=None):
        if groups_count == None:
            gc = self.get_groups_counts()
            #self.count_individuals()
        else:
            gc = groups_count
        
        bio_index = 0
        total_individuals = sum(self.populations.values())
        sps = 0
        
        prod_sps = gc['prods'][0] + gc['mut_prods'][0] 
        total_prods = gc['prods'][1] + gc['mut_prods'][1] 
        
        herb_sps = gc['herbs'][0] + gc['muts'][0]
        total_herbs = gc['herbs'][1] + gc['muts'][1]
        
        int_sps, total_ints = gc['prim_preds']
        top_sps, total_tops = gc['second_preds']
        
        bio_index_prod = 0
        bio_index_herb = 0
        bio_index_int = 0
        bio_index_top = 0
        
        for i in self.populations.keys():
            if self.populations[i] > 0:
                p_i = float(self.populations[i])/total_individuals
                bio_index += (p_i * math.log(p_i,2))
                sps += 1
                
                if self.species_scl[i] == 0:
                    p_i = float(self.populations[i])/total_prods
                    bio_index_prod += (p_i * math.log(p_i,2))
                if self.species_scl[i] == 1:
                    p_i = float(self.populations[i])/total_herbs
                    bio_index_herb += (p_i * math.log(p_i,2))
                if self.species_scl[i] == 2:
                    p_i = float(self.populations[i])/total_ints
                    bio_index_int += (p_i * math.log(p_i,2))
                if self.species_scl[i] == 3:
                    p_i = float(self.populations[i])/total_tops
                    bio_index_top += (p_i * math.log(p_i,2))
                
        bio_index = -1*bio_index
        bio_index_prod = -1*bio_index_prod
        bio_index_herb = -1*bio_index_herb
        bio_index_int = -1*bio_index_int
        bio_index_top = -1*bio_index_top
        
        indexes = dict()
        if sps <= 1:
            eq = 0.0
        else:
            eq = (bio_index/(math.log(sps,2))) 
        indexes['general'] =  (bio_index, eq)
        
        if prod_sps <= 1:
            eq = 0.0
        else:
            eq = (bio_index_prod/(math.log(prod_sps,2)))
        indexes['prods'] =  (bio_index_prod, eq)
        
        if herb_sps <= 1:
            eq = 0.0
        else:
            eq = (bio_index_herb/(math.log(herb_sps,2)))
        indexes['herbs'] =  (bio_index_herb, eq)
        
        if int_sps <= 1:
            eq = 0.0
        else:
            eq = (bio_index_int/(math.log(int_sps,2)))
        indexes['ints'] =  (bio_index_int, eq)
        
        if top_sps <= 1:
            eq = 0.0
        else:
            eq = (bio_index_top/(math.log(top_sps,2)))
        indexes['tops'] =  (bio_index_top, eq)
        
        return indexes 
    
    def extinguish_species(self, level, fraction):
        candidates = []
        if level == 'top':
            for sp in self.species_scl.keys():
                if self.species_scl[sp] == 3:
                    candidates.append(sp)
        
        number_of_extinctions = int(len(candidates)*fraction)
        
        self.rnd.shuffle(candidates)
        self.species_removed = []
        for i in range(number_of_extinctions):
            self.species_removed.append(candidates[i])
    
        print 'candidates', candidates, 'species to remove', self.species_removed
    
        for i in range(ROWS):
            for j in range(COLUMNS):
                cell = self.world[i][j]
                if cell.inhabitant != None and cell.inhabitant.species_id in self.species_removed:
                    cell.inhabitant = None
                    
                if cell.visitor != None and cell.visitor.species_id in self.species_removed:
                    cell.visitor = None
        
        return self.species_removed


    
    
    ###this method implements Moran's I and Geary's C indexes for spatial autocorrelation of each species in order to quantify 
    ###their clustering. This function relays on a faithful representation of population numbers on the world and hence the
    ###function 'count_individuals' has to be executed before it.
    def calculate_spatial_autocorrelation(self):
        
        species_autocorr = dict.fromkeys(self.net.nodes(), None)
        for s in species_autocorr.keys():
            species_autocorr[s] = dict(mean=self.populations[s]/(ROWS*COLUMNS),num=0, num_fn=0, num2=0, num2_fn=0, den=0, cur_Xi=0, morans_i=0, gearys_c=0, morans_i_fn=0, gearys_c_fn=0, empties=[], positions=[], mass_centre=(0,0), area=0, density=0)
        
#       weight_sum = 0
#       weight_sum_fn = 0
        for i in range(ROWS):
            for j in range(COLUMNS):
                cell = self.world[i][j]
                
                if cell.inhabitant is not None and type(cell.inhabitant) is not types.NoneType:
                    inhabitant = cell.inhabitant.species_id
                else:
                    inhabitant = -1
                    
                if cell.visitor is not None and type(cell.visitor) is not types.NoneType:
                    visitor = cell.visitor.species_id
                else:
                    visitor = -1
                
                if inhabitant != -1:
                    species_autocorr[inhabitant]['positions'].append((i,j))
                
                if visitor != -1:
                    species_autocorr[visitor]['positions'].append((i,j))
                
                for s in species_autocorr.keys():
                    if s == inhabitant or s == visitor:
                        Xi = 1
                    else:
                        Xi = 0
                        species_autocorr[s]['empties'].append((i,j))
                        
#                    species_autocorr[s]['cur_Xi'] = Xi
                    species_autocorr[s]['den'] += (Xi - species_autocorr[s]['mean'])**2
                
#                for l in range(ROWS):
#                    for n in range(COLUMNS):
#                    
#                        if i == l and j == n:
#                            continue
#                        
#                        #####this long IF statement is to verify whether cells are adjacent (i.e. they are connected)                        
#                        if ((i == ROWS-1 and l == 0) or (i == 0 and l == ROWS-1) or (j == 0 and n == COLUMNS-1) or (j == COLUMNS-1 and n == 0)):
#                            if (not ( ( (i,j) == (0,0) and (l,n) == (ROWS-1,COLUMNS-1) ) or ( (l,n) == (0,0) and (i,j) == (ROWS-1,COLUMNS-1) ) or ( (i,j) == (0,COLUMNS-1) and (l,n) == (ROWS-1,0) ) or ( (l,n) == (0,COLUMNS-1) and (i,j) == (ROWS-1,0) ) ) ) and ((abs( (i) - (l) ) % (ROWS-1) == 0 ) and (abs( (j) - (n) ) % (COLUMNS-1) == 0 )): 
#                                w_ij = 1
#                            else:
#                                w_ij = 0
#                        else:
#                            if ((abs( (i) - (l) ) % (ROWS-1) == 1 ) and (abs( (j) - (n) ) % (COLUMNS-1)  == 0 )) or  ((abs( (i) - (l) ) % (ROWS-1) == 0 ) and (abs( (j) - (n) ) % (COLUMNS-1)  == 1 )): 
#                                w_ij = 1
#                            else:
#                                w_ij = 0
#                        ######
#                        
#                        #### we also calculate neighbour status for full neighborhood
#                        if ((abs( (i) - (l) ) % (ROWS-1) <= 1 ) and (abs( (j) - (n) ) % (COLUMNS-1)  <= 1 )): 
#                            w_ij_fn = 1
#                        else:
#                            w_ij_fn = 0
#                        
#                        weight_sum += w_ij
#                        weight_sum_fn += w_ij_fn
#                        
#                        cell_b = self.world[i][j]
#                        
#                        if cell_b.inhabitant is not None and type(cell_b.inhabitant) is not types.NoneType:
#                            inhabitant_b = cell_b.inhabitant.species_id
#                        else:
#                            inhabitant_b = -1
#                    
#                        if cell_b.visitor is not None and type(cell_b.visitor) is not types.NoneType:
#                            visitor_b = cell_b.visitor.species_id
#                        else:
#                            visitor_b = -1
#                        
#                        for s in species_autocorr.keys():
#                            if s == inhabitant_b or s == visitor_b:
#                                Xj = 1
#                            else:
#                                Xj = 0
#                                
#                            species_autocorr[s]['num'] += w_ij * (species_autocorr[s]['cur_Xi'] - species_autocorr[s]['mean']) * (Xj - species_autocorr[s]['mean'])
#                            species_autocorr[s]['num2'] += w_ij * (species_autocorr[s]['cur_Xi'] - Xj)**2
#                            
#                            species_autocorr[s]['num_fn'] += w_ij_fn * (species_autocorr[s]['cur_Xi'] - species_autocorr[s]['mean']) * (Xj - species_autocorr[s]['mean'])
#                            species_autocorr[s]['num2_fn'] += w_ij_fn * (species_autocorr[s]['cur_Xi'] - Xj)**2
        

        #this part was re-coded to make the algorithm more efficient. When you know your cells are in a 2D grid connected to specific neighbour cells you don't have to 
        #visit every other cell...        
        weight_sum_fn = (ROWS*COLUMNS)*8
        for s in species_autocorr.keys():
            for (i,j) in species_autocorr[s]['positions']:
                neighbors = [ (i,(j+1)%(COLUMNS)), (i,(j-1)%(COLUMNS)), ((i+1)%(ROWS),j), ((i-1)%(ROWS), j), ((i+1)%(ROWS), (j+1)%(COLUMNS)), ((i-1)%(ROWS), (j+1)%(COLUMNS)), ((i+1)%(ROWS), (j-1)%(COLUMNS)), ((i-1)%(ROWS), (j-1)%(COLUMNS)) ]
                Xi = 1
                for n in neighbors:
                    if n in species_autocorr[s]['positions']:
                        Xj = 1
                    else:
                        Xj = 0
                    
                    species_autocorr[s]['num_fn'] += (Xi-species_autocorr[s]['mean']) * (Xj-species_autocorr[s]['mean'])
                    species_autocorr[s]['num2_fn'] += (Xi - Xj)**2
                    
            for (i,j) in species_autocorr[s]['empties']:
                neighbors = [ (i,(j+1)%(COLUMNS)), (i,(j-1)%(COLUMNS)), ((i+1)%(ROWS),j), ((i-1)%(ROWS), j), ((i+1)%(ROWS), (j+1)%(COLUMNS)), ((i-1)%(ROWS), (j+1)%(COLUMNS)), ((i+1)%(ROWS), (j-1)%(COLUMNS)), ((i-1)%(ROWS), (j-1)%(COLUMNS)) ]
                Xi = 0
                for n in neighbors:
                    if n in species_autocorr[s]['empties']:
                        Xj = 0
                    else:
                        Xj = 1
                    
                    species_autocorr[s]['num_fn'] += (Xi-species_autocorr[s]['mean']) * (Xj-species_autocorr[s]['mean'])
                    species_autocorr[s]['num2_fn'] += (Xi - Xj)**2
                                                
        for s in species_autocorr.keys():
            #species_autocorr[s]['morans_i'] = ( float(ROWS*COLUMNS) / weight_sum ) * ( species_autocorr[s]['num'] / species_autocorr[s]['den'] )
            #species_autocorr[s]['gearys_i'] =  ( float((ROWS*COLUMNS)-1) * species_autocorr[s]['num2'] ) / ( 2 * weight_sum * species_autocorr[s]['den'] )
            
            try:
                species_autocorr[s]['morans_i_fn'] = ( float(ROWS*COLUMNS) / weight_sum_fn ) * ( species_autocorr[s]['num_fn'] / species_autocorr[s]['den'] )
            except:
                species_autocorr[s]['morans_i_fn'] = 'N/A'
            
            try:
                species_autocorr[s]['gearys_c_fn'] = ( float((ROWS*COLUMNS)-1) * species_autocorr[s]['num2_fn'] ) / ( 2 * weight_sum_fn * species_autocorr[s]['den'] )
            except:
                species_autocorr[s]['gearys_c_fn'] = 'N/A'
            
            species_autocorr[s]['mass_centre'] = self.centroids[s]
            species_autocorr[s]['area'] = self.areas[s]
            
            if self.areas[s] != 0:
                species_autocorr[s]['density'] = self.populations[s] / species_autocorr[s]['area']
            
                    
        print 'finished calculating spatial correlation'
        
        return species_autocorr

    def calculate_spatial_state(self):

	inhabitant_matrix = np.zeros((ROWS, COLUMNS))
	visitor_matrix    = np.zeros((ROWS, COLUMNS))

        for i in range(ROWS):
            for j in range(COLUMNS):
                cell = self.world[i][j]
                
                if cell.inhabitant is not None and type(cell.inhabitant) is not types.NoneType:
                    inhabitant = cell.inhabitant.species_id
                else:
                    inhabitant = -1
                    
                if cell.visitor is not None and type(cell.visitor) is not types.NoneType:
                    visitor = cell.visitor.species_id
                else:
                    visitor = -1

		inhabitant_matrix[i,j] = inhabitant
		visitor_matrix[i,j] = visitor

	return (inhabitant_matrix, visitor_matrix)


class Cell():
    
    def __init__(self):
        self.habitat = None
        self.inhabitant = None
        self.visitor = None    
