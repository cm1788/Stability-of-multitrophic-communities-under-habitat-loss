

"""
	File name: main_cluster.py
	Author: Miguel Lurgi Rivera
	Date created: 03/08/2011

	Created in 2011 as part of my PhD dissertation The assembly and disassembly of ecological networks in a changing world
	Submitted to obtain the PhD degree to the Autonomous University of Barcelona
	Part of this work was funded by Microsoft Research

"""

import os
from datetime import datetime
import shutil
import networkx as nx
import csv
import math

from copy import copy

from web import Network
from network_creator import obtain_interactions_network
from ecosystem import Ecosystem

from utilities import get_out_row, get_eco_state_row, NetStats, EcosystemStats, write_spatial_analysis, write_spatial_state

from configure import ITERATIONS, HABITAT_LOSS, HABITAT_LOSS_ITER, INVASION, INVASION_ITER, NETWORK_RESET, SPATIAL_VARIATION
from configure import REFRESH_RATE, REMOVAL_LEVEL, REMOVAL_FRACTION, EXTINCTION_EVENT, TIME_WINDOW
from configure import SRC_NET_FILE, READ_FILE_NETWORK, NETWORK_RECORD, ITERATIONS_TO_RECORD, INT_STRENGTHS, RECORD_SPATIAL_VAR

__author__ = ["Miguel Lurgi", "Chris McWilliams"]
__copyright__ = "Copyright 2019"
__credits__ = ["Miguel Lurgi", "Chris McWilliams"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Miguel Lurgi"
__email__ = "miguel.lurgi@swansea.ac.uk"
__status__ = "Production"


if __name__ == '__main__':
    start_sim = datetime.now()
    
    ##### these modifications are performed so different replicates can be run on the cluster
    #job = '1' #os.environ['JOB_ID'];
    #task = '1' #os.environ['SGE_TASK_ID'];  
    
    job = os.environ['PBS_JOBID'];
    job = job[0:7]
    task = os.environ['PBS_ARRAYID'];
    
    output_dir = './' + job + '_' + task;

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if (os.path.isfile('./output/'+SRC_NET_FILE)):
        shutil.copy('./output/'+SRC_NET_FILE, output_dir)
    ##############################################
    
    header_names = ['iteration','S', 'L', 'L/S','C', 'T', 'B', 'I', 'Ca', 'Loop', 'NCycles', 'O', 'T-B', 'T-I', 'I-I', 'I-B', 'GenSD', 'VulSD', 'MxSim', 'MaxChainLength', 'MeanFoodChainLength', 'ChnSD', 'ChnNo', 'complexity', 'dynamic_complexity', 'components', 'cc', 'compartmentalisation', 'mean_tp', 'sd_tp', 'stable', 'mean_cv', 'nodf', 'h2', 'G_qi', 'V_qi', 'G_q', 'V_q', 'spatially_stable', 'mean_cv_centroid', 'mean_cv_area', 'mean_cv_density'];
    file_net = open(output_dir+'/output_network.csv', 'w')
    out = csv.DictWriter(file_net, header_names)
    
    ###this is for python < 2.7
#    headers_dict = dict()
#    for n in header_names:
#        headers_dict[n] = n
#        
#    out.writerow(headers_dict)
    
    ### for python >= 2.7 comment the above block and uncomment the following line 
    out.writeheader()
    
    header_names = ['iteration', 'total_sp', 'total_count', 'prod_sp', 'prod_count', 'mut_prod_sp', 'mut_prod_count', 'herb_sp', 'herb_count', 'mut_sp', 'mut_count', 'prim_pred_sp', 'prim_pred_count', 'sec_pred_sp', 'sec_pred_count', 'shannon_index', 'shannon_eq', 'shannon_index_prods', 'shannon_eq_prods', 'shannon_index_herbs', 'shannon_eq_herbs', 'shannon_index_interm', 'shannon_eq_interm', 'shannon_index_top', 'shannon_eq_top']
    file_eco = open(output_dir+'/output_ecosystem.csv', 'w')
    out_eco = csv.DictWriter(file_eco, header_names)
    
    ###this is for python < 2.7
#    headers_dict = dict()
#    for n in header_names:
#        headers_dict[n] = n
#        
#    out_eco.writerow(headers_dict)
    
    ### for python >= 2.7 comment the above block and uncomment the following line
    out_eco.writeheader()
    
    network_file = output_dir+'/'+SRC_NET_FILE
    if READ_FILE_NETWORK:
        graph = nx.read_graphml(network_file)
        net = Network(graph)
        
        print 'connectance = ', net.connectance()
        
        tls = net.get_trophic_levels()
        
        top, top_preds = net.top_predators()
        basal, basal_sps = net.basal()
        for u,v in net.edges():
            if u in basal_sps and v in top_preds and tls[v] == 3:
                net.remove_edge(u,v)
                
        print 'new connectance = ', net.connectance()
    else:
        net = obtain_interactions_network()
        net_to_save = net.copy()
        #nx.write_graphml(net_to_save, network_file)  ## new networkx doesn't like numpy floats
    
    ecosystem = Ecosystem(net, drawing=False)
    ecosystem.initialise_world(True)
    
    out_row = get_out_row(0, net, '', 0, '', '')
    out.writerow(out_row)
    
    out_row_eco = get_eco_state_row(0, ecosystem)
    out_eco.writerow(out_row_eco)
    
    
    series_counts = dict()
    if SPATIAL_VARIATION:
        centroids_counts = dict()
        areas_counts = dict()
    ##this structure holds the numbers of immigration, birth and dead of individuals
    ##for each species during the last ITERATIONS_TO_RECORD iterations
    cumulative_sps_stats = dict.fromkeys(net.nodes(), None)
    stats = ['immigrants', 'born', 'dead', 'tps']
    for sp in cumulative_sps_stats.keys():
        cumulative_sps_stats[sp] = dict.fromkeys(stats, 0)
    
    threshold_iter = math.ceil(ITERATIONS - (ITERATIONS*ITERATIONS_TO_RECORD))
    
    for i in range(1, ITERATIONS+1):
        print i
        ecosystem.update_world()
        #ecosystem.draw_species_distribution()
        
        if i >= threshold_iter:
            for sp in cumulative_sps_stats.keys():
                sp_stats = cumulative_sps_stats[sp]
                
                sp_stats['immigrants'] += ecosystem.new_inds_inmigration[sp]
                sp_stats['born'] += ecosystem.new_inds_reproduction[sp]
                sp_stats['dead'] += ecosystem.dead_individuals[sp] 
                                
#        if HABITAT_LOSS and i == HABITAT_LOSS_ITER:
#            net_temp = ecosystem.realised_net.copy()
##            layout_temp = nx.circular_layout(net_temp)
##            
##            fig_temp = figure()
##            network_p_temp = fig_temp.add_subplot(111)
##            nx.draw_networkx(net_temp, layout_temp, ax=network_p_temp)
##    
##            print 'S realised =', net_temp.order(), 'L realised =', net_temp.size(), 'C realised =', net_temp.connectance()
##            
#            ecosystem.clear_realised_network()
#            
#            out_row = get_out_row(i, net_temp)
#            out.writerow(out_row)
#            
#            ecosystem.apply_habitat_loss()
#            
#            #ecosystem.draw_species_distribution()
#        
#        if INVASION and i == INVASION_ITER:
#            net_temp = ecosystem.realised_net.copy()
#            ecosystem.invade(invaders)
#            
#            out_row = get_out_row(i, net_temp)
#            out.writerow(out_row)
#    
#        if i == (ITERATIONS - iteration_to_reset):
#            ecosystem.clear_realised_network()
            
        #eco_temp = copy(ecosystem)
        
        if SPATIAL_VARIATION and i%20 == 0:
            EcosystemStats(out_eco, ecosystem, i, series_counts, centroids_counts, areas_counts)
        else:
            EcosystemStats(out_eco, ecosystem, i, series_counts)
    
        if i%NETWORK_RECORD == 0 or i == ITERATIONS:
            net_temp = ecosystem.realised_net.copy()
            
            ##here we obtain the trophic position of each species so at the end we can calculate
            ##its mean and standard deviation for the species statistics
            tps, a, b = net_temp.find_trophic_positions()
            for sp in cumulative_sps_stats.keys():
                if cumulative_sps_stats[sp]['tps'] == 0:
                    cumulative_sps_stats[sp]['tps'] = []
                
                if tps.has_key(sp):
                    cumulative_sps_stats[sp]['tps'].append(tps[sp])
            
            if SPATIAL_VARIATION:
                NetStats(out, net_temp, i, NETWORK_RECORD, series_counts, INT_STRENGTHS, output_dir, centroids_counts, areas_counts)
            else:
                NetStats(out, net_temp, i, NETWORK_RECORD, series_counts, INT_STRENGTHS, output_dir)
            
            ecosystem.clear_realised_network()
        
        if HABITAT_LOSS and i == HABITAT_LOSS_ITER: 
            ecosystem.apply_habitat_loss()

        ##calculate spatial variation metrics
        if SPATIAL_VARIATION and (i%RECORD_SPATIAL_VAR == 0 or i == ITERATIONS):
            start = datetime.now()
            write_spatial_analysis(ecosystem, i, output_dir)
            write_spatial_state(ecosystem,i, output_dir)
            stop = datetime.now()
            elapsed = stop-start
            print elapsed
        
    file_net.close()
    file_eco.close()
    
    
    # here we write the output file for the species populations dynamics
    header_names = ['iteration']
    for sp in sorted(ecosystem.species):
        header_names.append(sp)
        
    file_populations = open(output_dir+'/output_populations.csv', 'w')
    out_populations = csv.DictWriter(file_populations, header_names)
    
    ###this is for python < 2.7
    headers_dict = dict()
    for n in header_names:
        headers_dict[n] = n
    
    out_populations.writerow(headers_dict)
    
    out_row_pops = dict()
    for iter in series_counts.keys():
        out_row_pops['iteration'] = iter
        for sp in sorted(series_counts[iter].keys()):
            out_row_pops[sp] = series_counts[iter][sp]
        
        out_populations.writerow(out_row_pops)
        
    file_populations.close()
    
    if SPATIAL_VARIATION:
        file_centroids = open(output_dir+'/output_centroids.csv', 'w')
        out_centroids = csv.DictWriter(file_centroids, header_names)
        
        out_centroids.writeheader()
        
        out_row_cents = dict()
        for iter in sorted(centroids_counts.keys()):
            out_row_cents['iteration'] = iter
            for sp in sorted(centroids_counts[iter].keys()):
                out_row_cents[sp] = centroids_counts[iter][sp]
            
            out_centroids.writerow(out_row_cents)
            
        file_centroids.close()
        
        file_areas = open(output_dir+'/output_areas.csv', 'w')
        out_areas = csv.DictWriter(file_areas, header_names)
        
        out_areas.writeheader()
        
        out_row_areas = dict()
        for iter in sorted(areas_counts.keys()):
            out_row_areas['iteration'] = iter
            for sp in sorted(areas_counts[iter].keys()):
                out_row_areas[sp] = areas_counts[iter][sp]
            
            out_areas.writerow(out_row_areas)
            
        file_areas.close()
    
    
    header_names = ['species', 'init_tl', 'final_tl', 'mutualist', 'mutualistic_producer', 'mean_tp', 'tp_sd', 'individuals', 'immigrants', 'born', 'dead']
    file_species = open(output_dir+'/output_species.csv', 'w')
    out_species = csv.DictWriter(file_species, header_names)
    
    out_species.writeheader()
    out_row_species = dict()
    
    for sp in sorted(cumulative_sps_stats.keys()):
        out_row_species['species'] = sp
        
        init_tls = net.get_trophic_levels()
        out_row_species['init_tl'] = init_tls[sp]
        
        final_tls = net_temp.get_trophic_levels()
        if final_tls.has_key(sp):
            out_row_species['final_tl'] = final_tls[sp]
        else:
            out_row_species['final_tl'] = 'N/A'
        
        out_row_species['mutualist'] = net.node[sp]['mut']
        out_row_species['mutualistic_producer'] = net.node[sp]['mut_prod']
        
        tps = cumulative_sps_stats[sp]['tps']
        if len(tps) == 0:
            mean_tps = 'N/A'
            sd_tps = 'N/A'
        else:
            mean_tps = sum(tps)/len(tps)
            sd_tps = 0.0
            for n in tps:
                sd_tps += (n-mean_tps)**2
        
            sd_tps = math.sqrt(sd_tps/len(tps))
        
        out_row_species['mean_tp'] = mean_tps
        out_row_species['tp_sd'] = sd_tps
        
        out_row_species['individuals'] = series_counts[ITERATIONS][sp]
        out_row_species['immigrants'] = cumulative_sps_stats[sp]['immigrants']
        out_row_species['born'] = cumulative_sps_stats[sp]['born']
        out_row_species['dead'] = cumulative_sps_stats[sp]['dead']
        
        out_species.writerow(out_row_species)
        
    file_species.close()
    
    stop_sim = datetime.now()
    elapsed_sim = stop_sim-start_sim
    print 'time for simulation' , elapsed_sim    
    
    
