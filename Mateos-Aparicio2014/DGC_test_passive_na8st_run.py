
import sys, os
import os.path
import click
import itertools
import numpy as np
from mpi4py import MPI # Must come before importing NEURON
from neuron import h
from neuroh5.io import read_tree_selection

def new_cell (template_name, local_id=0, gid=0, dataset_path="", neurotree_dict={}):
    h('objref cell, vx, vy, vz, vradius, vlayer, vsection, secnodes, vsrc, vdst, swc_types')
    h.vx       = neurotree_dict['x']
    h.vy       = neurotree_dict['y']
    h.vz       = neurotree_dict['z']
    h.vradius  = neurotree_dict['radius']
    h.vlayer   = neurotree_dict['layer']
    h.vsection = neurotree_dict['section']
    h.secnodes = neurotree_dict['section_topology']['nodes']
    h.vsrc     = neurotree_dict['section_topology']['src']
    h.vdst     = neurotree_dict['section_topology']['dst']
    h.swc_types  = neurotree_dict['swc_type']
    hstmt      = 'cell = new %s(%d, %d, "%s", vlayer, vsrc, vdst, secnodes, vx, vy, vz, vradius, swc_types)' % (template_name, local_id, gid, dataset_path)
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

    print "rank = ", rank
    print "size = ", size
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
    (trees, _) = read_tree_selection (comm, forest_path, pop_name, myselection)

    h('objref results_passive, results_single_ap, results_threshold, results_ap_rate')
    h.results_passive   = h.List()
    h.results_single_ap = h.List()
    h.results_threshold = h.List()
    h.results_ap_rate   = h.List()
    
    for (gid, tree) in trees.iteritems():
        cell = new_cell ("DGC", neurotree_dict=tree)
        h.passive_test(h.cell, h.results_passive, gid)
        h.single_ap_test(h.cell, h.results_single_ap, gid)
        h.threshold_test(h.cell, h.results_threshold, gid)
        h.ap_rate_test(h.cell, h.results_ap_rate, gid, 0.3)

    results_dict_passive   = hoc_results_to_python(h.results_passive)
    results_dict_single_ap = hoc_results_to_python(h.results_single_ap)
    results_dict_threshold = hoc_results_to_python(h.results_threshold)
    results_dict_ap_rate   = hoc_results_to_python(h.results_ap_rate)

    all_results_passive = comm.gather(sendobj=results_dict_passive, root=0)
    if rank==0:
        write_results(all_results_passive,
                    results_path+"/"+"DGC_tests_results_passive.dat",
                    "# Gid, input resistance, vmin,  vtau0, tau0, dendritic surface area")

    all_results_single_ap = comm.gather(sendobj=results_dict_single_ap, root=0)
    if rank==0:
        write_results(all_results_single_ap,
                    results_path+"/"+"DGC_tests_results_single_ap.dat",
                    "# Gid, vmax, vmin, vdend0, vdend1, vdend2, vdend3, vdend4")

    all_results_threshold = comm.gather(sendobj=results_dict_threshold, root=0)
    if rank==0:
        write_results(all_results_threshold,
                    results_path+"/"+"DGC_tests_results_threshold.dat",
                    "# Gid, vmax, vmin, v threshold, t threshold, amplitude, ahp")

    all_results_ap_rate = comm.gather(sendobj=results_dict_ap_rate, root=0)
    if rank==0:
        write_results(all_results_ap_rate,
                    results_path+"/"+"DGC_tests_results_ap_rate.dat",
                    "# Gid, FR mean, ISI mean, ISI variance, ISI stdev, ISI adapt 1, ISI adapt 2, ISI adapt 3, ISI adapt 4")


if __name__ == '__main__':
    main(args=sys.argv[(sys.argv.index("DGC_test_passive_na8st_run.py")+1):])
