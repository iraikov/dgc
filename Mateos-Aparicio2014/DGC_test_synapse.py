
import sys, os
import os.path
import click
import itertools
import numpy as np
from mpi4py import MPI # Must come before importing NEURON
from neuron import h
from neuroh5.io import read_tree_selection

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


## Compute per-segment relative counts of synapse placement
def compute_synapse_relcounts(density_dict, neurotree_dict, seglist):
    relcounts       = []
    relcount_total  = 0
    layers          = []
    for seg in seglist:
        L    = seg.sec.L
        l    = 0.0
        lsum = 0.0
        relcount_dict = {}
        relcount_total = 0
        ## length of the current segment
        l = 2 * (seg.x*L - lsum)
        ## length from the beginning of the section to the end of the current segment
        lsum = lsum + l
        layer = get_node_attribute('layer', neurotree_dict, seg.sec, seg.x)
        if layer not in density_dict:
            break
        layers.append(layer)
        rp = density_dict[layer]*l
        relcount_total += rp
        relcounts.append(rp)
    return (relcounts, relcount_total, layers)
    
           
def compute_synapse_locations(syn_type, seed, density_dict, neurotree_dict, sec_list):
    """Populates a morphology with synapse locations"""

    seg_list = []
    sec_dict = {}
    sec_index = 0
    for sec in sec_list:
        sec_dict[sec] = sec_index
        sec_index    += 1
        for seg in sec:
            if seg.x < 1.0 and seg.x > 0.0:
                seg_list.append(seg)
            
    relcounts, total, layers = compute_synapse_relcounts(density_dict, neurotree_dict, seg_list)

    ran = h.Random(seed)
    ran.uniform(0, total)

    sample_size = int(total)
    sample = h.Vector(sample_size)
    sample.setrand(ran)
    sample.sort()

    syn_ids  = []
    syn_locs = []
    syn_secs = []
    syn_layers = []
    syn_types  = []

    syn_index   = 0
    cumrelcount = 0
    for i in xrange(0,len(seg_list)-1):
        seg = seg_list[i]
        seg_start = seg.x - (0.5/seg.sec.nseg)
        seg_end   = seg.x + (0.5/seg.sec.nseg)
        seg_range = seg_end - seg_start
        rel_count = relcounts[i]
        int_rel_count = int(rel_count)
        cumrelcount += rel_count
        layer = layers[i]
        syn_count = 0
        while sample[syn_index] < cumrelcount: 
            syn_loc = seg_start + seg_range * ((syn_count + 1) / rel_count)
            syn_locs.append(syn_loc)
            syn_ids.append(syn_index)
            syn_secs.append(sec_dict[seg.sec])
            syn_layers.append(layer)
            syn_types.append(syn_type)
            syn_index += 1
            syn_count += 1
            if syn_count == int_rel_count:
                break
            if syn_index == sample_size:
                break
        if syn_index == sample_size:
            break

    return({'syn_ids': syn_ids,'syn_locs': syn_locs,'syn_secs': syn_secs,'syn_layers': syn_layers,'syn_types': syn_types})
    
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

    density_dict = {'apical': {'layer': {1: 3.36, 2: 2.28, 3: 2.02}, 'default': 2.55} }  # units: synapses / um length

    syn_Exc = 0
    
    for (gid, tree) in trees.iteritems():
        cell = new_cell ("DGC", neurotree_dict=tree)
        syn_dict = compute_synapse_locations(syn_Exc, 13, density_dict['apical']['layer'], tree, cell.alldendrites)
        
     

if __name__ == '__main__':
    main(args=sys.argv[(sys.argv.index("DGC_test_synapse.py")+1):])
