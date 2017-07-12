
import sys, os
import os.path
import click
import itertools
import numpy as np
from mpi4py import MPI # Must come before importing NEURON
from neuron import h
from neuroh5.io import read_tree_selection

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

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

@click.command()
@click.option("--template-path", required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True))
@click.option("--forest-path", required=True, type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("--results-path", required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True))
@click.option('--selection', default=[], callback=lambda _,__,x: map(int, x.split(',')) if x else [])
@click.option("--selection-file", default=None, type=click.Path(exists=True, file_okay=True, dir_okay=False))
def main(template_path, forest_path, results_path, selection, selection_file):

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
    
    pop_name = "GC"
    (trees, _) = read_tree_selection (MPI._addressof(comm), forest_path, pop_name, myselection)
    
    for (gid, tree) in trees.iteritems():
        cell = new_cell ("DGC", neurotree_dict=tree)
        h.passive_test(h.cell, results_path, gid)
        h.single_ap_test(h.cell, results_path, gid)
        h.threshold_test(h.cell, results_path, gid)
        h.ap_rate_test(h.cell, results_path, gid, 0.3)
    


if __name__ == '__main__':
    main(args=sys.argv[(sys.argv.index("DGC_test_passive_na8st_run.py")+1):])
