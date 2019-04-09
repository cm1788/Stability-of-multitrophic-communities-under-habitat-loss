
from operator import itemgetter
import math
import csv
from copy import copy
import threading
import networkx as nx

import numpy
from scipy.sparse import coo_matrix

#import matplotlib.pyplot as plt

from web import Network



def get_mutualistic_matrix_from_network(net):
    '''
    obtains the mutualistic matrix of interactions between mutualits and their hosts
    provided that each of this set of species is identified by the tags 'mut' and 
    'mut_prod' respectively as an attribute of each node.
    
    It returns:
        1.- an array of arrays which represents the matrix of interactions (0 = no interaction,
            1 = interaction), between the hosts and the visitors by row.
        2.- an array of arrays which represents the matrix of weighted interactions between
            the hosts and the visitors by row. Weights are obtained from the attribute 'is'
            of each link in the network, which is supposed to carry the interaction strenghts.
        3.- a list with the ids of the hosts in the network
        4.- a list with the ids of the visitors in the network
        
    '''
    muts = [];
    hosts = [];
    
    for n in net.nodes():
        if net.node[n]['mut']:
            muts.append(n);
        if net.node[n]['mut_prod']:
            hosts.append(n);
 
    mutualistic_matrix = [];
    weighted_matrix = [];
    producers_edges = net.edges(hosts);
    
    for i in range(len(hosts)):
        row = [];
        weighted_row =[];
        for j in range(len(muts)):
            if (hosts[i], muts[j]) in producers_edges:
                row.append(1);
                if net[hosts[i]][muts[j]].has_key('is'):
                    weighted_row.append(net[hosts[i]][muts[j]]['is']);
                else:
                    weighted_row.append(0);
            else:
                row.append(0);
                weighted_row.append(0);
        
        mutualistic_matrix.append(row);
        weighted_matrix.append(weighted_row);
    
    print 'hosts :', hosts
    print 'muts: ', muts
    
    return mutualistic_matrix, weighted_matrix, hosts, muts


def compare_rows(a, b):
    sum_a = 0
    sum_b = 0
    length = len(a)
    
    for i in range(length):
        if a[i] == 1:
            sum_a += 1
        if b[i] == 1:
            sum_b += 1
    
    if sum_a == sum_b:
        for i in range(length):
            if (a[i] == 1 and b[i] == 1) or (a[i] == 0 and b[i] == 0):
                continue
            elif a[i] == 1 and b[i] == 0:
                return -1
            elif a[i] == 0 and b[i] == 1:
                return 1
            else:
                return 0
    if sum_a < sum_b:
        return 1
    
    return -1

def order_mutualistic_matrix(matrix, prods, herbs):
    rows = len(matrix)
    cols = len(matrix[0])
    
    dict_order_prods = dict()
    
    #for ordering the columns we need to implement a somewhat more complicated algorithm...
    for i in range(cols):
        colk = []
        for j in range(rows):
            colk.append(matrix[j][i])
                
        dict_order_prods[prods[i]] = colk
    
    ordered_prods = sorted(dict_order_prods.items(), compare_rows, key=itemgetter(1))
    temp_matrix = []
    producers_order = []
    for i in range(rows):
        temp_row = []
        for j in range(cols):
            temp_row.append(0)
        temp_matrix.append(temp_row)
    
    for i in range(len(ordered_prods)):
        id, col = ordered_prods[i]
        producers_order.append(id)
        for j in range(rows):
            temp_matrix[j][i] = col[j]
        
    print 'plants: ', producers_order
    
    dict_rows = dict()
    for idx in range(len(herbs)):
        dict_rows[herbs[idx]] = temp_matrix[idx]
    
    #this only line arranges the rows in descending order from top to bottom based on the number of ones they have
    ordered_rows = sorted(dict_rows.items(), compare_rows, key=itemgetter(1));
    
    matrix = []
    ordered_herbs = []
    for (idx, row) in ordered_rows:
        matrix.append(row)
        ordered_herbs.append(idx)
    
    print 'animals: ',ordered_herbs
    
    return matrix, producers_order, ordered_herbs



#### this function differs from the previous one in which here hosts are considered to be
#### represented by rows and mutualists by columns
def sort_mutualistic_matrix(matrix, prods, herbs):
    rows = len(matrix)
    cols = len(matrix[0])
    
    dict_order_muts = dict()
    
    #for ordering the columns we need to implement a somewhat more complicated algorithm...
    for i in range(cols):
        colk = []
        for j in range(rows):
            colk.append(matrix[j][i])
                
        dict_order_muts[herbs[i]] = colk
    
    ordered_muts = sorted(dict_order_muts.items(), compare_rows, key=itemgetter(1))
    #####  this creates a matrix of zeroes with the same size as the original matrix 
    temp_matrix = []
    visitors_order = []
    for i in range(rows):
        temp_row = []
        for j in range(cols):
            temp_row.append(0)
        temp_matrix.append(temp_row)
    
    #########
    
    for i in range(len(ordered_muts)):
        id, col = ordered_muts[i]
        visitors_order.append(id)
        for j in range(rows):
            temp_matrix[j][i] = col[j]
        
    #print 'visitors: ', visitors_order
    
    dict_rows = dict()
    for idx in range(len(prods)):
        dict_rows[prods[idx]] = temp_matrix[idx]
    
    #this only line arranges the rows in descending order from top to bottom based on the number of ones they have
    ordered_rows = sorted(dict_rows.items(), compare_rows, key=itemgetter(1));
  
    matrix = []
    ordered_hosts = []
    for (idx, row) in ordered_rows:
        matrix.append(row)
        ordered_hosts.append(idx)
    
    #print 'hosts: ',ordered_hosts
    
    return matrix, visitors_order, ordered_hosts

    
def calculate_nodf(matrix):
    '''
    algorithm for calculating the NODF measure for nestedness of a 0-1 matrix of interactions
    the matrix must be ordered according to number of ones in each row/column from top/left
    to bottom/right in order for the algorithm to work correctly
    (after Almeida-Neto et al.)
    '''
    n_pair_row = 0.0
    n_pair_col= 0.0
    rows = len(matrix)
    cols = len(matrix[0])
    np_row = (rows*(rows-1))/2
    np_col = (cols*(cols-1))/2
    
    paired_nested_degrees_rows = []
    #we first calculate the pairing indexes row-wise
    for i in range(rows):
        row_i = matrix[i];
        visited = False;
        mt_i = 0;
        for j in range(i+1, rows):
            row_j = matrix[j];
            mt_j = 0;
            ones_i = 0;
            n_pair_ij = 0.0;
            
            for k in range(len(row_i)):
                if not visited:
                    mt_i += row_i[k]
                mt_j += row_j[k];                    
                if row_j[k] == 1 and row_i[k] == 1:
                    ones_i += 1
            visited = True;
            
            if mt_i > mt_j and mt_j > 0:
                n_pair_ij = float(ones_i)/float(mt_j)
            
            n_pair_row += n_pair_ij;
            paired_nested_degrees_rows.append(n_pair_ij);    
     
    #we now proceed to calculate the pairing indexes column-wise
    paired_nested_degrees_cols = []
    for k in range(cols):
        colk = []
        mt_k = 0;
        for i in range(rows):
            colk.append(matrix[i][k])
            mt_k += colk[i]
            
        for l in range(k+1,cols):
            mt_l = 0
            ones_k = 0
            n_pair_kl = 0.0;
            
            for i in range(rows):
                pos = matrix[i][l]
                mt_l += pos
                if pos == 1 and colk[i] == 1:
                    ones_k += 1
            
            if mt_k > mt_l and mt_l > 0:
                n_pair_kl = float(ones_k)/float(mt_l);
            
            n_pair_col += n_pair_kl;
            paired_nested_degrees_cols.append(n_pair_kl);
    
    nodf = (n_pair_row + n_pair_col)*100/(np_row + np_col);
    
    #print 'NODF = ', nodf,' n_pair_row = ',(float((n_pair_row/len(paired_nested_degrees_rows))*100)),' n_pair_col = ', (float((n_pair_col/len(paired_nested_degrees_cols))*100))
    
    return nodf


def get_out_row(iteration, net, series_counts, offset, centroids, areas):
    out_row = dict()
    
    net.longest_path_length()
    out_row['iteration'] = iteration   
    out_row['S'] = net.order()
    out_row['L'] = net.size()
    out_row['L/S'] = net.linkage_density()
    out_row['C'] = net.connectance()
    out_row['T'], ts = net.top_predators()
    out_row['B'], bs = net.basal()
    out_row['I'], ints = net.intermediate() 
    out_row['Ca'] = net.cannibalism()
    out_row['Loop'], out_row['NCycles'] = net.fraction_in_loops()
    out_row['O'], os = net.omnivory()
    
    fractions = net.get_links_fractions_between_levels()
    
    out_row['T-B'] = fractions['tb']
    out_row['T-I'] = fractions['ti']
    out_row['I-I'] = fractions['ii']
    out_row['I-B'] = fractions['ib']
    out_row['GenSD'], out_row['VulSD'] = net.generality_vulnerability_sd()
    out_row['MxSim'] = net.maximum_similarity()
    out_row['MaxChainLength'] = net.longest_path_length() 
    tps, nops, out_row['MeanFoodChainLength'] = net.find_trophic_positions()
    chn_var, out_row['ChnSD'], out_row['ChnNo'] = net.get_path_length_feats()
    out_row['complexity'] = net.complexity()
    out_row['dynamic_complexity'] = net.dynamic_complexity()
    out_row['components'] = net.components()
    
    if net.size() > 0 and net.order() > 0:
        out_row['cc'] = nx.average_clustering(net.to_undirected())
    else:
        out_row['cc'] = 0.0
        
    out_row['compartmentalisation'] = net.degree_of_compartmentalization()
    
    out_row['mean_tp'] = numpy.mean(tps.values())
    out_row['sd_tp'] = numpy.std(tps.values())
    
    #populations stability measures
    ### we look only at the stability in terms of individuals per species since it does not make much sense
    ### to look at changes in overall numbers of individuals since it is not biomass (we cannot obtain biomass change)
    if series_counts == '' or offset == 0:
        out_row['stable'] = '';
        out_row['mean_cv'] = '';
        
        out_row['spatially_stable'] = '';
        out_row['mean_cv_centroid'] = '';
        out_row['mean_cv_area'] = '';
        out_row['mean_cv_density'] = '';
    else:
        start = iteration - offset + 1
        stop = iteration
        cvs = []
        
        if centroids != '' and areas != '' and series_counts != '':
            space = True
            cvs_centroids = []
            cvs_areas = []
            cvs_densities = []
        else:
            out_row['spatially_stable'] = '';
            out_row['mean_cv_centroid'] = '';
            out_row['mean_cv_area'] = '';
            out_row['mean_cv_density'] = '';
            
        for sp in net.nodes():
            current_sps = []
            
            if space:
                current_sps_sp = []
                current_cents = []
                current_ars = []
                current_dens = []
            
            for iterat in range(start,stop):
                current_sps.append(series_counts[iterat][sp]);
                
                if space and centroids.has_key(iterat):
                    current_cents.append(centroids[iterat][sp]);
                    current_ars.append(areas[iterat][sp]);
                    current_sps_sp.append(series_counts[iterat][sp]);
                
                
            mean = numpy.mean(current_sps);
            sd = numpy.std(current_sps);
        
            cvs.append(sd/mean);
            
            if space:
                current_dens = numpy.array(current_sps_sp)/numpy.array(current_ars)
                current_cents = numpy.array(current_cents)
                
                mean_area = numpy.mean(current_ars);
                sd_area = numpy.std(current_ars);
                cvs_areas.append(sd_area/mean_area);
                
                mean_dens = numpy.mean(current_dens);
                sd_dens = numpy.std(current_dens);
                cvs_densities.append(sd_dens/mean_dens);
                
                mean_cents = numpy.mean(current_cents, axis=0);
                sd_cents = numpy.std(current_cents, axis=0);
                cvs_centroids.append(sd_cents/mean_cents);
            
            
        mean_cv = numpy.mean(cvs);
        if mean_cv <= 0.1:
            out_row['stable'] = True;
        else:
            out_row['stable'] = False;
        
        out_row['mean_cv'] = mean_cv;
    
        if space:
            mean_cv_area = numpy.mean(cvs_areas);
            mean_cv_density = numpy.mean(cvs_densities);
            mean_cv_centroids = numpy.mean(numpy.array(cvs_centroids), axis=0);
            
            out_row['mean_cv_area'] = mean_cv_area;
            out_row['mean_cv_density'] = mean_cv_density;
            out_row['mean_cv_centroid'] = mean_cv_centroids;
            
        if mean_cv_area <= 0.1 and mean_cv_density <= 0.1 and mean_cv_centroids[0] <= 0.1 and mean_cv_centroids[1] <= 0.1:
            out_row['spatially_stable'] = True;
        else:
            out_row['spatially_stable'] = False;
        
    ####   mutualistic measures
    
    #### nestedness (using NODF algorithm)
    mutualistic_matrix, weighted_matrix, prods, muts = get_mutualistic_matrix_from_network(net);
    
    if len(mutualistic_matrix) == 0:
        out_row['nodf'] = '';
    else:
        ordered_matrix, prods, muts = sort_mutualistic_matrix(mutualistic_matrix, prods, muts);
        nodf = calculate_nodf(ordered_matrix);
        out_row['nodf'] = nodf;
    
    #### H2 measure (after Bluthgen et al. 2006)
    weighted_matrix = numpy.matrix(weighted_matrix);
    weight_sum = float(numpy.sum(weighted_matrix));
    
    if weight_sum == 0.0:
        out_row['h2'] = '';
    else:
        weighted_matrix = (weighted_matrix)/weight_sum;
        
        ### this is for applying the log function only to non-zero elements
        ### the coo scipy matrix contains the matrix of log-transformed weighted elements
        coo = coo_matrix(weighted_matrix);
        coo.data = numpy.log(coo.data);
        
        multi = numpy.multiply(weighted_matrix, coo.todense());
        
        sum_h = numpy.sum(multi);    
        out_row['h2'] = -sum_h;
    
    ####   quantitative measures
    
    # ts, bs, ints are the top, basal and intermediate species respectively
    if( (len(ts) + len(ints)) == 0 ):
        out_row['G_qi'] = '';
    else:
        out_row['G_qi'] = net.size()/float(len(ts) + len(ints));
    
    if( (len(bs) + len(ints)) == 0 ):
        out_row['V_qi'] = '';
    else:
        out_row['V_qi'] = net.size()/float(len(bs) + len(ints));
    
    if( net.size() == 0 ):
        out_row['G_q'] = '';
        out_row['V_q'] = '';
    else:
        (a,b,w) = net.edges(data=True)[0];
        if not 'is' in w:
            out_row['G_q'] = '';
            out_row['V_q'] = '';    
        else:
            total_biomass_flux = 0;
            for (v, u, atts) in net.edges(data=True):
                if 'is' in atts:
                    total_biomass_flux += atts['is'];
            
            g_q_val = 0.0;
            v_q_val = 0.0;
            
            print total_biomass_flux;
            
            for n in net.nodes():
                prey = net.predecessors(n);
                predators = net.successors(n);
                
                b_coming = 0;
                for p in prey:
                    if 'is' in net[p][n]:
                        b_coming += net[p][n]['is']; 
                
                g_q_val += (b_coming/float(total_biomass_flux)) * len(prey); 
                
                b_going = 0;
                for p in predators:
                    if 'is' in net[n][p]:
                        b_going += net[n][p]['is'];
                
                v_q_val += (b_going/float(total_biomass_flux)) * len(predators);
            
            out_row['G_q'] = g_q_val;
            out_row['V_q'] = v_q_val;
        
    return out_row
    
def get_eco_state_row(iteration, ecosystem, spatial=False):
    out_row = dict()
    
    gc = ecosystem.get_groups_counts(spatial)
    out_row['iteration'] = iteration
    out_row['total_sp'], out_row['total_count'] = gc['total']
    out_row['prod_sp'], out_row['prod_count'] = gc['prods']
    out_row['mut_prod_sp'], out_row['mut_prod_count'] = gc['mut_prods']
    out_row['herb_sp'], out_row['herb_count'] = gc['herbs']
    out_row['mut_sp'], out_row['mut_count'] = gc['muts']
    out_row['prim_pred_sp'], out_row['prim_pred_count'] = gc['prim_preds']
    out_row['sec_pred_sp'], out_row['sec_pred_count'] = gc['second_preds']
    
    indexes = ecosystem.get_shannon_biodiversity_index(gc)
    out_row['shannon_index'], out_row['shannon_eq'] = indexes['general']
    out_row['shannon_index_prods'], out_row['shannon_eq_prods'] = indexes['prods']
    out_row['shannon_index_herbs'], out_row['shannon_eq_herbs'] = indexes['herbs']
    out_row['shannon_index_interm'], out_row['shannon_eq_interm'] = indexes['ints']
    out_row['shannon_index_top'], out_row['shannon_eq_top'] = indexes['tops']
    return out_row

#def plot_populations_series(plot, i, x, u, y, v, w, z, title):
#    plot.clear()
#    plot.hold(True)
#      
#    if title == 'inmigration':
#        style1 = 'ro'
#        style2 = 'bo'
#        style3 = 'go'
#        style4 = 'yo'
#        style5 = 'co'
#        style6 = 'mo'
#    else:
#        style1 = 'r-'
#        style2 = 'b-'
#        style3 = 'g-'
#        style4 = 'y-'
#        style5 = 'c-'
#        style6 = 'm-'
#        
#    plot.plot(i, x, style1)
#    plot.plot(i, u, style6)
#    plot.plot(i, y, style2)
#    plot.plot(i, v, style3)
#    plot.plot(i, w, style4)
#    plot.plot(i, z, style5)
#    
#    if title == 'populations':
#        legend_labels = ['prods', 'mut_prods', 'herbs', 'muts', 'prim', 'sec']
#        plot.legend(legend_labels)
#    
#    for i, label in enumerate(years):
#        plt.text(x1[i], y1[i], label)
#        plt.text(x2[i], y2[i], label)
#        plt.text(x3[i], y3[i], label)
#
#    plot.set_xlabel(x_label)
#    plot.set_ylabel(y_label)
#    plot.set_title(title)
    
    #plot.set_xscale('log')
#    if title == 'populations' or title == 'reproduction' or title == 'death':
#        plot.set_yscale('log')

#    plot.hold(False)


def calculate_stability_measures(before, after, tls, mut_prod):
    
    species_before = set()
    species_after = set()
    
    producers = set()
    for sp in tls.keys():
        if tls[sp] == 0:
            producers.add(sp)
        
    values_all = []
    values_mut = []
    values_prod = []
    for i in before.keys():
        pops = before[i]
        all = 0
        prod = 0
        mut = 0
        for sp in pops:
            if pops[sp] > 0:
                species_before.add(sp)
            if sp in producers:
                all += pops[sp]
                if sp in mut_prod:
                    mut += pops[sp]
                else:
                    prod += pops[sp]
        
        values_all.append(all)
        values_mut.append(mut)
        values_prod.append(prod)
        
    mean_all = float(sum(values_all))/len(values_all)
    mean_mut = float(sum(values_mut))/len(values_mut)
    mean_prod = float(sum(values_prod))/len(values_prod)
    
    variance_all = 0.0
    variance_mut = 0.0
    variance_prod = 0.0
    
    for i in range(len(values_all)):
        variance_all += (values_all[i]-mean_all)**2
        variance_mut += (values_mut[i]-mean_mut)**2
        variance_prod += (values_prod[i]-mean_prod)**2
    
    variance_all = variance_all/len(values_all)
    variance_mut = variance_mut/len(values_mut)
    variance_prod = variance_prod/len(values_prod)
    
    temp_var_before = math.sqrt(variance_all)/mean_all
    temp_var_prod_before = math.sqrt(variance_prod)/mean_prod
    #temp_var_mut_before = math.sqrt(variance_mut)/mean_mut
    
    #after
    values_all = []
    values_mut = []
    values_prod = []
    for i in after.keys():
        pops = after[i]
        all = 0
        prod = 0
        mut = 0
        for sp in pops:
            if pops[sp] > 0:
                species_after.add(sp)
            if sp in producers:
                all += pops[sp]
                if sp in mut_prod:
                    mut += pops[sp]
                else:
                    prod += pops[sp]
        
        values_all.append(all)
        values_mut.append(mut)
        values_prod.append(prod)
        
    mean_all = float(sum(values_all))/len(values_all)
    mean_mut = float(sum(values_mut))/len(values_mut)
    mean_prod = float(sum(values_prod))/len(values_prod)
    
    variance_all = 0.0
    variance_mut = 0.0
    variance_prod = 0.0
    
    for i in range(len(values_all)):
        variance_all += (values_all[i]-mean_all)**2
        variance_mut += (values_mut[i]-mean_mut)**2
        variance_prod += (values_prod[i]-mean_prod)**2
    
    variance_all = variance_all/len(values_all)
    variance_mut = variance_mut/len(values_mut)
    variance_prod = variance_prod/len(values_prod)
    
    temp_var_after = math.sqrt(variance_all)/mean_all
    temp_var_prod_after = math.sqrt(variance_prod)/mean_prod
    #temp_var_mut_after = math.sqrt(variance_mut)/mean_mut
    
    
    extinctions = species_before - species_after
    print 'extinctions =', extinctions
    for sp in extinctions:
        print sp, ':', tls[sp]
    
    print temp_var_before, temp_var_after #, temp_var_prod_before, temp_var_prod_after #, temp_var_mut_before, temp_var_mut_after


def write_adjacency_matrix(iteration, offset, series_counts, net, output_dir):  
    init_iter = iteration - offset
    average_counts = dict()

    for i in range(init_iter+1, iteration+1):
        pops_iter = series_counts[i]
        for sp in pops_iter.keys():
            if not average_counts.has_key(sp):
                average_counts[sp] = pops_iter[sp]
            else:
                average_counts[sp] += pops_iter[sp]

    for sp in average_counts.keys():
        average_counts[sp] = math.floor(average_counts[sp]/offset) 
    

    header_names = ['predator/prey']
    for n in sorted(net.nodes()):
        header_names.append(n)
    
    file_ad = open(output_dir+'/adjacency_'+str(iteration)+'_is1.csv', 'w')
    out_adj = csv.DictWriter(file_ad, header_names)

    file_ad_is2 = open(output_dir+'/adjacency_'+str(iteration)+'_is2.csv', 'w')
    out_adj_is2 = csv.DictWriter(file_ad_is2, header_names)

    file_ad_is3 = open(output_dir+'/adjacency_'+str(iteration)+'_is3.csv', 'w')
    out_adj_is3 = csv.DictWriter(file_ad_is3, header_names)
  
    out_adj.writeheader()
    out_adj_is2.writeheader()
    out_adj_is3.writeheader()
    
    out_adj_row = dict()
    out_adj2_row = dict()
    out_adj3_row = dict()

    for predator in sorted(net.nodes()):
        out_adj_row['predator/prey'] = predator
        out_adj2_row['predator/prey'] = predator
        out_adj3_row['predator/prey'] = predator
        for prey in sorted(net.nodes()):
            if (prey,predator) in net.edges():
                if not 'is' in net[prey][predator]:
                    strength = 0.0
                else:
                    strength = net[prey][predator]['is']
                
                out_adj_row[prey] = strength

                if float(average_counts[prey]) == 0.0:
                    out_adj2_row[prey] = 0.0
                else:
                    out_adj2_row[prey] = float(strength/average_counts[prey])
                            
                fraction = float(average_counts[prey]*average_counts[predator]) 
                if fraction == 0.0:
                    weight = 0.0
                else:
                    weight = float(strength)/fraction
            
                out_adj3_row[prey] = weight
                net[prey][predator]['weight'] = weight 
            else:
                out_adj_row[prey] = 0
                out_adj2_row[prey] = 0
                out_adj3_row[prey] = 0
    
        out_adj.writerow(out_adj_row)
        out_adj_is2.writerow(out_adj2_row)
        out_adj_is3.writerow(out_adj3_row)
    
    mutualistic_matrix, weighted_matrix, prods, muts = get_mutualistic_matrix_from_network(net);
    weighted_matrix = numpy.matrix(weighted_matrix);
    numpy.savetxt(output_dir+'/mut_network_'+ str(iteration) +'_is1.csv', weighted_matrix, delimiter=",", fmt='%d')
    

def write_spatial_analysis(ecosystem, iteration, output_dir='../output'):
    header_names = ['species', 'mass_centre', 'area_dist', 'density', 'morans_i', 'gearys_c']
    file_spatial = open(output_dir+'/output_space_'+ str(iteration) +'.csv', 'w')
    out_spatial = csv.DictWriter(file_spatial, header_names)
    
    out_spatial.writeheader()
    
    out_row = dict()
    spatial_metrics = ecosystem.calculate_spatial_autocorrelation()
    for sp in sorted(ecosystem.net.nodes()):
        out_row['species'] = sp
        out_row['mass_centre'] = spatial_metrics[sp]['mass_centre']
     
        out_row['area_dist'] = spatial_metrics[sp]['area']
        out_row['density'] = spatial_metrics[sp]['density']
        out_row['morans_i'] = spatial_metrics[sp]['morans_i_fn'] 
        out_row['gearys_c'] = spatial_metrics[sp]['gearys_c_fn']
        
    
        out_spatial.writerow(out_row)
    
    file_spatial.close()

def write_spatial_state(ecosystem, iteration, output_dir='../output'):
# method added to record the spatial state of the system 
# called on same iteration as above method
# records: output_inhabitant_ITERATION#.csv and output_visitor_ITERATION#.csv
# which contain ROWSxCOLUMNS matrices with -1 -> no species present, otherwise entry=speciesID number  
    filename_inhabitant = output_dir+'/output_inhabitant_'+ str(iteration) +'.csv'
    filename_visitor = output_dir+'/output_visitor_'+ str(iteration) +'.csv'

    (inhabitant_matrix, visitor_matrix) = ecosystem.calculate_spatial_state()
    numpy.savetxt(filename_inhabitant, inhabitant_matrix, delimiter=",")
    numpy.savetxt(filename_visitor, visitor_matrix, delimiter=",")

def NetStats(writer, net, iteration, offset, series, calculate_is, out_dir, centroids=None, areas=None):
    series_counts = copy(series)
    
    if not(centroids is None):
        centroids_tmp = copy(centroids)
    if not(areas is None):
        areas_tmp = copy(areas)
    
    if calculate_is:
        write_adjacency_matrix(iteration, offset, series_counts, net, out_dir)
    
    if not(centroids is None) and not(areas is None):
        dict_stats = get_out_row(iteration, net, series_counts, offset, centroids_tmp, areas_tmp)
    else:
        dict_stats = get_out_row(iteration, net, series_counts, offset, '', '')
        
    writer.writerow(dict_stats)

def EcosystemStats(writer, ecosystem, i, series_counts, centroids=None, areas=None):
    if centroids is None:
        dict_stats = get_eco_state_row(i, ecosystem)
    else:
        dict_stats = get_eco_state_row(i, ecosystem, True)
    
    series_counts[i] = ecosystem.populations
    
    if not(centroids is None):
        centroids[i] = copy(ecosystem.centroids)
        
    if not(areas is None):
        areas[i] = copy(ecosystem.areas)
    
    writer.writerow(dict_stats)    


class ThreadNetStats(threading.Thread):
    def __init__(self, lock, writer, net, iteration, offset, series, waiting_thread, calculate_is, centroids=None, areas=None):
        threading.Thread.__init__(self)
        self.net = net
        self.lock = lock
        self.writer = writer
        self.iteration = iteration
        self.offset = offset
        self.series_counts = series
        self.waiting_thread = waiting_thread
        self.calculate_is = calculate_is
        
        self.centroids = centroids
        self.areas = areas
        
    def run(self):
        self.waiting_thread.join()
        self.series_counts = copy(self.series_counts)
        self.centroids = copy(self.centroids)
        self.areas = copy(self.areas)
        
        if self.calculate_is:
            self.write_adjacency_matrix()
        
        dict_stats = get_out_row(self.iteration, self.net, self.series_counts, self.offset, self.centroids, self.areas)
        
        self.lock.acquire()
        self.writer.writerow(dict_stats)
        self.lock.release()
    
    def write_adjacency_matrix(self):
        init_iter = self.iteration-self.offset
        average_counts = dict()
        for i in range(init_iter+1, self.iteration+1):
            pops_iter = self.series_counts[i]
            for sp in pops_iter.keys():
                if not average_counts.has_key(sp):
                    average_counts[sp] = pops_iter[sp]
                else:
                    average_counts[sp] += pops_iter[sp]
        
        for sp in average_counts.keys():
            average_counts[sp] = math.floor(average_counts[sp]/self.offset) 

        header_names = ['predator/prey']
        for n in sorted(self.net.nodes()):
            header_names.append(n)
            
        file_ad = open('../output/adjacency_'+str(self.iteration)+'_is1.csv', 'w')
        out_adj = csv.DictWriter(file_ad, header_names)
        
        file_ad_is2 = open('../output/adjacency_'+str(self.iteration)+'_is2.csv', 'w')
        out_adj_is2 = csv.DictWriter(file_ad_is2, header_names)
        
        file_ad_is3 = open('../output/adjacency_'+str(self.iteration)+'_is3.csv', 'w')
        out_adj_is3 = csv.DictWriter(file_ad_is3, header_names)
        
        ###this is for python < 2.7
        headers_dict = dict()
        for n in header_names:
            headers_dict[n] = n
            
        out_adj.writerow(headers_dict)
        out_adj_is2.writerow(headers_dict)
        out_adj_is3.writerow(headers_dict)
        
        out_adj_row = dict()
        out_adj2_row = dict()
        out_adj3_row = dict()
        
        for predator in sorted(self.net.nodes()):
            out_adj_row['predator/prey'] = predator
            out_adj2_row['predator/prey'] = predator
            out_adj3_row['predator/prey'] = predator
            for prey in sorted(self.net.nodes()):
                if (prey,predator) in self.net.edges():
                    out_adj_row[prey] = self.net[prey][predator]['is']
                    if float(average_counts[prey]) == 0.0:
                        out_adj2_row[prey] = 0.0
                    else:
                        out_adj2_row[prey] = float(self.net[prey][predator]['is'])/average_counts[prey]
                    
                    fraction = float(average_counts[prey]*average_counts[predator]) 
                    if fraction == 0.0:
                        weight = 0.0
                    else:
                        weight = float(self.net[prey][predator]['is'])/fraction
                    out_adj3_row[prey] = weight
                    self.net[prey][predator]['weight'] = weight 
                else:
                    out_adj_row[prey] = 0
                    out_adj2_row[prey] = 0
                    out_adj3_row[prey] = 0
            
            out_adj.writerow(out_adj_row)
            out_adj_is2.writerow(out_adj2_row)
            out_adj_is3.writerow(out_adj3_row)
        
        file_ad.close()
        file_ad_is2.close()
        file_ad_is3.close()

class ThreadEcosystemStats(threading.Thread):
    def __init__(self, lock, writer, ecosystem, iter, series, centroids=None, areas=None):
        threading.Thread.__init__(self)
        self.eco = ecosystem
        self.lock = lock
        self.writer = writer
        self.iteration = iter
        self.series_counts = series
        self.centroids = centroids
        self.areas = areas
    
    def run(self):
        self.lock.acquire()
        
        if self.centroids is None:
            dict_stats = get_eco_state_row(self.iteration, self.eco)
        else:
            dict_stats = get_eco_state_row(self.iteration, self.eco, True)
            
        self.series_counts[self.iteration] = copy(self.eco.populations)
        
        if not(self.centroids is None):
            self.centroids[self.iteration] = copy(self.eco.centroids)
        
        if not(self.areas is None):
            self.areas[self.iteration] = copy(self.eco.areas)
        
        self.writer.writerow(dict_stats)
        self.lock.release()

