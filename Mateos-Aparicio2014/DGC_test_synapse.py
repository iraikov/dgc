
import sys, os
import os.path
import click
import itertools
import numpy as np
from mpi4py import MPI # Must come before importing NEURON
from neuron import h
from neuroh5.io import read_tree_selection

syn_Excitatory = 0
syn_Inhibitory = 1

swc_soma = 1
swc_axon = 2
swc_apical = 4


def new_cell (template_name, local_id=0, gid=0, dataset_path="", neurotree_dict={}):
    h('objref cell, vx, vy, vz, vradius, vlayer, vsection, secnodes, vsrc, vdst')
    h.vx       = neurotree_dict['x']
    h.vy       = neurotree_dict['y']
    h.vz       = neurotree_dict['z']
    h.vradius  = neurotree_dict['radius']
    h.vlayer   = neurotree_dict['layer']
    h.vsection = neurotree_dict['section']
    h.secnodes = neurotree_dict['section_topology']['nodes']
    h.vsrc     = neurotree_dict['section_topology']['src']
    h.vdst     = neurotree_dict['section_topology']['dst']
    hstmt      = 'cell = new %s(%d, %d, "%s", vlayer, vsrc, vdst, secnodes, vx, vy, vz, vradius)' % (template_name, local_id, gid, dataset_path)
    h(hstmt)
    return h.cell


def hoc_results_to_python(hoc_results):
    results_dict = {}
    for i in xrange(0, int(hoc_results.count())):
        vect   = hoc_results.o(i)
        gid    = int(vect.x[0])
        pyvect = vect.to_python()
        results_dict[gid] = pyvect[1:]
    hoc_results.remove_all()
    return results_dict


def write_results(results, filepath, header):
    f = open(filepath,'w')
    f.write(header+'\n')
    for item in results:
        for (gid, vect) in item.iteritems():
            f.write (str(gid)+"\t")
            f.write (("\t".join(['{:0.3f}'.format(i) for i in vect])) + "\n")
    f.close()

    
def get_node_attribute (name, content, sec, x=None):
    if name in content:
        if x is None:
            return content[name]
        elif sec.n3d() == 0:
            return content[name][0]
        else:
            for i in xrange(sec.n3d()):
                if sec.arc3d(i)/sec.L >= x:
                    return content[name][i]
    else:
        return None


def synapse_relcounts(layer_density_dicts, seglist, seed, neurotree_dict=None):
    """Computes per-segment relative counts of synapse placement"""
    relcounts_dict  = {}
    relcount_total  = 0
    layers_dict     = {}
    relcount_total  = 0
    for (syn_type, layer_density_dict) in layer_density_dicts.iteritems():
        rans = {}
        for (layer,density_dict) in layer_density_dict.iteritems():
            ran = h.Random(seed)
            ran.normal(density_dict['mean'], density_dict['variance'])
            rans[layer] = ran
        relcounts = []
        layers    = []
        for seg in seglist:
            L    = seg.sec.L
            nseg = seg.sec.nseg
            if neurotree_dict is not None:
                layer = get_node_attribute('layer', neurotree_dict, seg.sec, seg.x)
            else:
                layer = -1
            layers.append(layer)
            
            ran=None
            if layer > -1:
                if layer in rans:
                    ran = rans[layer]
            else:
                ran = rans['default']
                
            if ran is not None:
                l         = L/nseg
                rc        = ran.repick()*l
                relcount_total += rc
                relcounts.append(rc)
            else:
                relcounts.append(0)
                
        relcounts_dict[syn_type] = relcounts
        layers_dict[syn_type]   = layers
    return (relcounts_dict, relcount_total, layers_dict)
    
           
def distribute_uniform_synapses(seed, sec_layer_density_dicts, sec_lists, sec_swc_types, neurotree_dicts):
    """Computes uniformly-spaced synapse locations"""

    syn_ids    = []
    syn_locs   = []
    syn_secs   = []
    syn_layers = []
    syn_types  = []
    swc_types  = []
    syn_index  = 0

    for (layer_density_dicts, sec_list, swc_type, neurotree_dict) in itertools.izip(sec_layer_density_dicts,
                                                                                    sec_lists,
                                                                                    sec_swc_types,
                                                                                    neurotree_dicts):
        
        seg_list = []
        sec_dict = {}
        sec_index = 0
        L_total   = 0
        for sec in sec_list:
            L_total += sec.L
            sec_dict[sec] = sec_index
            sec_index    += 1
            for seg in sec:
                if seg.x < 1.0 and seg.x > 0.0:
                    seg_list.append(seg)
            
    
        relcounts_dict, total, layers_dict = synapse_relcounts(layer_density_dicts, seg_list, seed, neurotree_dict=neurotree_dict)

        print 'total = ', total, ' L_total = ', L_total
        sample_size = total
        cumcount  = 0
        for (syn_type, _) in layer_density_dicts.iteritems():
            relcounts = relcounts_dict[syn_type]
            layers    = layers_dict[syn_type]
            for i in xrange(0,len(seg_list)):
                seg = seg_list[i]
                seg_start = seg.x - (0.5/seg.sec.nseg)
                seg_end   = seg.x + (0.5/seg.sec.nseg)
                seg_range = seg_end - seg_start
                rel_count = relcounts[i]
                int_rel_count = round(rel_count)
                layer = layers[i]
                syn_count = 0
                while syn_count < int_rel_count:
                    syn_loc = seg_start + seg_range * ((syn_count + 1) / rel_count)
                    syn_locs.append(syn_loc)
                    syn_ids.append(syn_index)
                    syn_secs.append(sec_dict[seg.sec])
                    syn_layers.append(layer)
                    syn_types.append(syn_type)
                    swc_types.append(swc_type)
                    syn_index += 1
                    syn_count += 1
                cumcount += syn_count

    syn_dict = {'syn_ids': np.asarray(syn_ids, dtype='uint32'),
                'syn_locs': np.asarray(syn_locs, dtype='float32'),
                'syn_secs': np.asarray(syn_secs, dtype='uint32'),
                'syn_layers': np.asarray(syn_layers, dtype='int8'),
                'syn_types': np.asarray(syn_types, dtype='uint8'),
                'swc_types': np.asarray(swc_types, dtype='uint8')}

    return syn_dict

def print_syn_summary (gid,syn_dict):
    print 'gid %d: ' % gid
    print '\t total %d synapses' % len(syn_dict['syn_ids'])
    print '\t %d excitatory synapses' % np.size(np.where(syn_dict['syn_types'] == syn_Excitatory))
    print '\t %d inhibitory synapses' % np.size(np.where(syn_dict['syn_types'] == syn_Inhibitory))
    print '\t %d apical excitatory synapses' % np.size(np.where((syn_dict['syn_types'] == syn_Excitatory) & (syn_dict['swc_types'] == swc_apical)))
    print '\t %d apical inhibitory synapses' % np.size(np.where((syn_dict['syn_types'] == syn_Inhibitory) & (syn_dict['swc_types'] == swc_apical)))
    print '\t %d soma excitatory synapses' % np.size(np.where((syn_dict['syn_types'] == syn_Excitatory) & (syn_dict['swc_types'] == swc_soma)))
    print '\t %d soma inhibitory synapses' % np.size(np.where((syn_dict['syn_types'] == syn_Inhibitory) & (syn_dict['swc_types'] == swc_soma)))
    print '\t %d ais excitatory synapses' % np.size(np.where((syn_dict['syn_types'] == syn_Excitatory) & (syn_dict['swc_types'] == swc_axon)))
    print '\t %d ais inhibitory synapses' % np.size(np.where((syn_dict['syn_types'] == syn_Inhibitory) & (syn_dict['swc_types'] == swc_axon)))



@click.command()
@click.option("--template-path", required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True))
@click.option("--forest-path", required=True, type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("--results-path", required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True))
@click.option('--selection', default=[], callback=lambda _,__,x: map(int, x.split(',')) if x else [])
@click.option("--selection-file", default=None, type=click.Path(exists=True, file_okay=True, dir_okay=False))
def main(template_path, forest_path, results_path, selection, selection_file):

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    h.load_file(template_path+"/"+"DGC_Tests_from_file_passive_na8st.hoc")

    if selection_file is not None:
        f = open(selection_file)
        lines = f.readlines()

        while lines:
            for l in lines:
                selection.append(int(l))
            lines = f.readlines()

        f.close()

    myselection = [selection[i] for i in xrange(rank,len(selection),size)]

    print 'rank %d: myselection = ' % rank, myselection
    
    pop_name = "GC"
    (trees, _) = read_tree_selection (MPI._addressof(comm), forest_path, pop_name, myselection)
    
    dend_exc_density_dict = {n: {'mean': 1.77, 'variance': 0.13} for n in [2,3,4]} # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1360197/
    dend_exc_density_dict = {n: {'mean': 2.26, 'variance': 0.07} for n in [2,3,4]} # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2538941/
    
    dend_inh_density_dict = { 2: {'mean': 1.0, 'variance': 0.01}, # https://www.ncbi.nlm.nih.gov/pubmed/20034063
                              4: {'mean': 3.5, 'variance': 0.01} }

    soma_inh_density_dict = { 'default': {'mean': 8.0, 'variance': 0.2} } # https://www.ncbi.nlm.nih.gov/pubmed/20034063
                            
    ais_inh_density_dict = { 'default': {'mean': 8.0, 'variance': 0.2} } # https://www.ncbi.nlm.nih.gov/pubmed/20034063
                            

    dend_layer_density_dict = {syn_Excitatory: dend_exc_density_dict, syn_Inhibitory: dend_inh_density_dict}
    soma_layer_density_dict = {syn_Inhibitory: soma_inh_density_dict}
    ais_layer_density_dict  = {syn_Inhibitory: ais_inh_density_dict}
    sec_layer_density_dicts = [dend_layer_density_dict, soma_layer_density_dict, ais_layer_density_dict]

    sec_swc_types = [swc_apical, swc_soma, swc_axon]
    
    for (gid, tree) in trees.iteritems():
        cell = new_cell ("DGC", neurotree_dict=tree)
        sec_lists = [cell.alldendrites, cell.soma, cell.ais]
        syn_dict = distribute_uniform_synapses(gid, sec_layer_density_dicts, sec_lists, sec_swc_types, neurotree_dicts=[tree, None, None])
        print_syn_summary(gid,syn_dict)
     

if __name__ == '__main__':
    main(args=sys.argv[(sys.argv.index("DGC_test_synapse.py")+1):])
