
import sys, os, math
import os.path
import click
import itertools
from collections import defaultdict
import numpy as np
from mpi4py import MPI # Must come before importing NEURON
from neuron import h
from neuroh5.io import read_tree_selection

syn_type_excitatory = 0
syn_type_inhibitory = 1

swc_type_soma = 1
swc_type_axon = 2
swc_type_basal  = 3
swc_type_apical = 4

swc_type_dict = {
      'soma':   1,
      'axon':   2,
      'ais':    2,
      'basal':  3,
      'apical': 4,
      'trunk':  5,
      'tuft':   6
    }

v_init = -75

def get_node_attribute (name, content, sec, secnodes, x=None):
    if content.has_key(name):
        if x is None:
            return content[name]
        elif sec.n3d() == 0:
            return content[name][0]
        else:
            for i in xrange(sec.n3d()):
                if sec.arc3d(i)/sec.L >= x:
                    return content[name][secnodes[i]]
    else:
        return None

def make_neurotree_cell (template_class, local_id=0, gid=0, dataset_path="", neurotree_dict={}):
    vx       = neurotree_dict['x']
    vy       = neurotree_dict['y']
    vz       = neurotree_dict['z']
    vradius  = neurotree_dict['radius']
    vlayer   = neurotree_dict['layer']
    vsection = neurotree_dict['section']
    secnodes = neurotree_dict['section_topology']['nodes']
    vsrc     = neurotree_dict['section_topology']['src']
    vdst     = neurotree_dict['section_topology']['dst']
    swc_type = neurotree_dict['swc_type']
    cell     = template_class (local_id, gid, dataset_path, vlayer, vsrc, vdst, secnodes, vx, vy, vz, vradius, swc_type)
    return cell

def synapse_seg_density(layer_density_dicts, sec_index_dict, seglist, seed, neurotree_dict=None):
    """Computes per-segment density of synapse placement. """
    segdensity_dict  = {}
    layers_dict     = {}
    if neurotree_dict is not None:
        secnodes_dict = neurotree_dict['section_topology']['nodes']
    else:
        secnodes_dict = None
    for (syn_type, layer_density_dict) in layer_density_dicts.iteritems():
        rans = {}
        for (layer,density_dict) in layer_density_dict.iteritems():
            ran = h.Random(seed)
            ran.normal(density_dict['mean'], density_dict['variance'])
            rans[layer] = ran
        segdensity = []
        layers     = []
        for seg in seglist:
            L    = seg.sec.L
            nseg = seg.sec.nseg
            if neurotree_dict is not None:
                secindex = sec_index_dict[seg.sec]
                secnodes = secnodes_dict[secindex]
                layer = get_node_attribute('layer', neurotree_dict, seg.sec, secnodes, seg.x)
            else:
                layer = -1
            layers.append(layer)
            
            ran=None

            if layer > -1:
                if rans.has_key(layer):
                    ran = rans[layer]
                elif rans.has_key('default'):
                    ran = rans['default']
                else:
                    ran = None
            elif rans.has_key('default'):
                ran = rans['default']
            else:
                ran = None
            if ran is not None:
                dens      = ran.repick()
                segdensity.append(dens)
            else:
                segdensity.append(0)
                
        segdensity_dict[syn_type] = segdensity
        layers_dict[syn_type]     = layers
    return (segdensity_dict, layers_dict)


def synapse_seg_counts(layer_density_dicts, sec_index_dict, seglist, seed, neurotree_dict=None):
    """Computes per-segment relative counts of synapse placement. """
    segcounts_dict  = {}
    layers_dict     = {}
    segcount_total  = 0
    if neurotree_dict is not None:
        secnodes_dict = neurotree_dict['section_topology']['nodes']
    else:
        secnodes_dict = None
    for (syn_type, layer_density_dict) in layer_density_dicts.iteritems():
        rans = {}
        for (layer,density_dict) in layer_density_dict.iteritems():
            ran = h.Random(seed)
            ran.normal(density_dict['mean'], density_dict['variance'])
            rans[layer] = ran
        segcounts = []
        layers    = []
        for seg in seglist:
            L    = seg.sec.L
            nseg = seg.sec.nseg
            if neurotree_dict is not None:
                secindex = sec_index_dict[seg.sec]
                secnodes = secnodes_dict[secindex]
                layer = get_node_attribute('layer', neurotree_dict, seg.sec, secnodes, seg.x)
            else:
                layer = -1
            layers.append(layer)
            
            ran=None

            if layer > -1:
                if rans.has_key(layer):
                    ran = rans[layer]
                elif rans.has_key('default'):
                    ran = rans['default']
                else:
                    ran = None
            elif rans.has_key('default'):
                ran = rans['default']
            else:
                ran = None
            if ran is not None:
                l         = L/nseg
                dens      = ran.repick()
                rc        = dens*l
                segcount_total += rc
                segcounts.append(rc)
            else:
                segcounts.append(0)
                
        segcounts_dict[syn_type] = segcounts
        layers_dict[syn_type]    = layers
    return (segcounts_dict, segcount_total, layers_dict)
    
           
def distribute_uniform_synapses(seed, swc_type_dict, sec_layer_density_dict, neurotree_dict, sec_dict, secidx_dict):
    """Computes uniformly-spaced synapse locations. """

    syn_ids    = []
    syn_locs   = []
    syn_secs   = []
    syn_layers = []
    syn_types  = []
    swc_types  = []
    syn_index  = 0

    for (sec_name, layer_density_dict) in sec_layer_density_dict.iteritems():

        sec_index_dict = secidx_dict[sec_name]
        swc_type = swc_type_dict[sec_name]
        seg_list = []
        sec_obj_index_dict = {}
        L_total   = 0
        (seclst, maxdist) = sec_dict[sec_name]
        secidxlst         = secidx_dict[sec_name]
        for (sec, secindex) in itertools.izip(seclst, secidxlst):
            sec_obj_index_dict[sec] = int(secindex)
            if maxdist is None:
                for seg in sec:
                    if seg.x < 1.0 and seg.x > 0.0:
                        seg_list.append(seg)
            else:
                for seg in sec:
                    if seg.x < 1.0 and seg.x > 0.0 and ((L_total + sec.L * seg.x) <= maxdist):
                        seg_list.append(seg)
            L_total += sec.L
        segcounts_dict, total, layers_dict = synapse_seg_counts(layer_density_dict, sec_obj_index_dict, seg_list, seed, neurotree_dict=neurotree_dict)

        sample_size = total
        cumcount  = 0
        for (syn_type, _) in layer_density_dict.iteritems():
            segcounts = segcounts_dict[syn_type]
            layers    = layers_dict[syn_type]
            for i in xrange(0,len(seg_list)):
                seg = seg_list[i]
                seg_start = seg.x - (0.5/seg.sec.nseg)
                seg_end   = seg.x + (0.5/seg.sec.nseg)
                seg_range = seg_end - seg_start
                seg_count = segcounts[i]
                int_seg_count = math.floor(seg_count)
                layer = layers[i]
                syn_count = 0
                while syn_count < int_seg_count:
                    syn_loc = seg_start + seg_range * ((syn_count + 1) / math.ceil(seg_count))
                    assert((syn_loc <= 1) & (syn_loc >= 0))
                    if syn_loc < 1.0:
                        syn_locs.append(syn_loc)
                        syn_ids.append(syn_index)
                        syn_secs.append(sec_obj_index_dict[seg.sec])
                        syn_layers.append(layer)
                        syn_types.append(syn_type)
                        swc_types.append(swc_type)
                        syn_index += 1
                        syn_count += 1
                cumcount += syn_count

    assert(len(syn_ids) > 0)
    syn_dict = {'syn_ids': np.asarray(syn_ids, dtype='uint32'),
                'syn_locs': np.asarray(syn_locs, dtype='float32'),
                'syn_secs': np.asarray(syn_secs, dtype='uint32'),
                'syn_layers': np.asarray(syn_layers, dtype='int8'),
                'syn_types': np.asarray(syn_types, dtype='uint8'),
                'swc_types': np.asarray(swc_types, dtype='uint8')}

    return syn_dict


def syn_in_seg(seg, syns_dict):
    if seg.sec not in syns_dict:
        return False
    if any(seg.sec(x) == seg for x in syns_dict[seg.sec]): return True
    return False

def make_syn_mech(mech_name, seg):
    if mech_name == 'AMPA':
        syn = h.Exp2Syn(seg)
    elif mech_name == 'GABA_A':
        syn = h.Exp2Syn(seg)
    elif mech_name == 'GABA_B':
        syn = h.Exp2Syn(seg)
    else:
        raise ValueError("Unrecognized synaptic mechanism name %s", mech_name)
    return syn

def add_shared_synapse(mech_name, seg, syns_dict):
    """Returns the existing synapse in segment if any, otherwise creates it."""
    if not syn_in_seg(seg, syns_dict):
        syn = make_syn_mech(mech_name, seg)
        syns_dict[seg.sec][syn.get_segment().x] = syn
        return syn
    else:
        for x, syn in syns_dict[seg.sec].iteritems():
            if seg.sec(x) == seg:
               return syn

def add_unique_synapse(mech_name, seg, syns_dict):
    """Creates a synapse in the given segment."""
    syn = make_syn_mech(mech_name, seg)
    return syn
    
def mksyns(gid,cell,syn_ids,syn_types,swc_types,syn_locs,syn_sections,syn_kinetic_params,add_synapse=add_shared_synapse,spines=False):

    syns_dict_dend = defaultdict(lambda: defaultdict(lambda: {}))
    syns_dict_axon = defaultdict(lambda: defaultdict(lambda: {}))
    syns_dict_soma = defaultdict(lambda: defaultdict(lambda: {}))
    py_sections = [sec for sec in cell.sections]
    
    syn_obj_dict = {}

    for i in xrange(0, syn_ids.size):

      syn_id      = syn_ids[i]
      if not (syn_id < syn_types.size):
          print 'mksyns syn_ids for gid %i: ' % gid, syn_ids
          raise ValueError('mksyns: cell %i received invalid syn_id %d' % (gid, syn_id))
      
      syn_type    = syn_types[syn_id]
      swc_type    = swc_types[syn_id]
      syn_loc     = syn_locs[syn_id]
      syn_section = syn_sections[syn_id]
      
      sref = None
      sec = py_sections[syn_section]
      if swc_type == swc_type_apical:
        syns_dict = syns_dict_dend
        if syn_type == syn_type_excitatory: 
            if spines and h.ismembrane('spines',sec=sec):
                sec(syn_loc).count_spines += 1
      elif swc_type == swc_type_basal:
        syns_dict = syns_dict_dend
        if syn_type == syn_type_excitatory: 
            if spines and h.ismembrane('spines',sec=sec):
                sec(syn_loc).count_spines += 1
      elif swc_type == swc_type_axon:
        syns_dict = syns_dict_axon
      elif swc_type == swc_type_soma:
        syns_dict = syns_dict_soma
      else: 
        raise RuntimeError ("Unsupported synapse SWC type %d" % swc_type)
      syn_mech_dict = {}
      for (syn_mech, params) in syn_kinetic_params.iteritems():
        syn = add_synapse(syn_mech, sec(syn_loc), syns_dict)
        syn.tau1 = params['t_rise']
        syn.tau2 = params['t_decay']
        syn.e    = params['e_rev']
        cell.syns.append(syn)
        cell.syntypes.o(syn_type).append(syn)
        syn_mech_dict[syn_mech] = syn
      syn_obj_dict[syn_id] = syn_mech_dict
        
    if spines:
        cell.correct_for_spines()

    return syn_obj_dict


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

def print_syn_summary (gid,syn_dict):
    print 'gid %d: ' % gid
    print '\t total %d synapses' % len(syn_dict['syn_ids'])
    print '\t %d excitatory synapses' % np.size(np.where(syn_dict['syn_types'] == syn_type_excitatory))
    print '\t %d inhibitory synapses' % np.size(np.where(syn_dict['syn_types'] == syn_type_inhibitory))
    print '\t %d apical excitatory synapses' % np.size(np.where((syn_dict['syn_types'] == syn_type_excitatory) & (syn_dict['swc_types'] == swc_type_apical)))
    print '\t %d apical inhibitory synapses' % np.size(np.where((syn_dict['syn_types'] == syn_type_inhibitory) & (syn_dict['swc_types'] == swc_type_apical)))
    print '\t %d soma excitatory synapses' % np.size(np.where((syn_dict['syn_types'] == syn_type_excitatory) & (syn_dict['swc_types'] == swc_type_soma)))
    print '\t %d soma inhibitory synapses' % np.size(np.where((syn_dict['syn_types'] == syn_type_inhibitory) & (syn_dict['swc_types'] == swc_type_soma)))
    print '\t %d axon excitatory synapses' % np.size(np.where((syn_dict['syn_types'] == syn_type_excitatory) & (syn_dict['swc_types'] == swc_type_axon)))
    print '\t %d axon inhibitory synapses' % np.size(np.where((syn_dict['syn_types'] == syn_type_inhibitory) & (syn_dict['swc_types'] == swc_type_axon)))

    
def synapse_vclamp_test (label, syntype, cell, w, v_holding, v_init, indexes, section):

    vv = h.Vector()
    vv.append(0,0,0,0,0,0)

    se = h.SEClamp(cell.sections[section](0.5))

    h('objref synlst')
    h.synlst = h.List()
    for i in indexes:
        h.synlst.append(cell.syns.o(i))
    if syntype == syn_type_excitatory:
        v = cell.syntest_exc(h.synlst,se,w,v_holding,v_init)
    else:
        v = cell.syntest_inh(h.synlst,se,w,v_holding,v_init)
        
        
    vv = vv.add(v)
    
    amp     = vv.x[0]
    t_10_90 = vv.x[1]
    t_20_80 = vv.x[2]
    t_all   = vv.x[3]
    t_50    = vv.x[4]
    t_decay = vv.x[5]

    print("%s synapse: \n" % label)
    print("  Amplitude %f\n" % amp)
    print("  10-90 Rise Time %f\n" % t_10_90)
    print("  20-80 Rise Time %f\n" % t_20_80)
    print("  Decay Time Constant %f\n" % t_decay)

        
def synapse_iclamp_test (label, syntype, cell, w, v_init, indexes):

    h('objref synlst')
    h.synlst = h.List()
    for i in indexes:
        h.synlst.append(cell.syns.o(i))
    if syntype == syn_type_excitatory:
        v = cell.syn_iclamp_exc(h.synlst,w,v_init)
    else:
        v = cell.syn_iclamp_inh(h.synlst,w,v_init)
        
    v_amp   = v.x[0]
    v_peak  = v.x[1]
    v_pre   = v.x[2]
    t_peak  = v.x[3]
    t_pre   = v.x[4]

    print("%s synapse: \n" % label)
    print("  V Amplitude %f\n" % v_amp)
    print("  V Peak %f\n" % v_peak)
    print("  V Pre %f\n" % v_pre)

    return v_amp

def syn_group_test(syn_label, gid, cell, syn_selection, syn_ids, syn_connection_params, v_init, id_base=0):
        
    mean_v_amp = 0.
    
    for i in syn_selection:
        print 'Synapse id: ', i
        mean_v_amp += synapse_iclamp_test(syn_label, 0, cell, syn_connection_params['weight']/1000.0, 
                                          v_init, [i-id_base])

    mean_v_amp = mean_v_amp / syn_selection.size
    
    print("gid %d: %s synapse: mean V amplitude is %f\n" % (gid, syn_label, mean_v_amp))
    return mean_v_amp

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
    (trees_iter, _) = read_tree_selection (comm, forest_path, pop_name, myselection)
    
    dend_exc_density_dict = {n: {'mean': 1.77, 'variance': 0.13} for n in [2,3,4]} # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1360197/
    dend_exc_density_dict = {n: {'mean': 2.26, 'variance': 0.07} for n in [2,3,4]} # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2538941/
    
    dend_inh_density_dict = { 2: {'mean': 1.0, 'variance': 0.01}, # https://www.ncbi.nlm.nih.gov/pubmed/20034063
                              4: {'mean': 3.5, 'variance': 0.01} }

    soma_inh_density_dict = { 'default': {'mean': 8.0, 'variance': 0.2} } # https://www.ncbi.nlm.nih.gov/pubmed/20034063
                            
    ais_inh_density_dict = { 'default': {'mean': 8.0, 'variance': 0.2} } # https://www.ncbi.nlm.nih.gov/pubmed/20034063
                            

    dend_layer_density_dict = {syn_type_excitatory: dend_exc_density_dict, syn_type_inhibitory: dend_inh_density_dict}
    soma_layer_density_dict = {syn_type_inhibitory: soma_inh_density_dict}
    ais_layer_density_dict  = {syn_type_inhibitory: ais_inh_density_dict}
    sec_layer_density_dict  = {'apical': dend_layer_density_dict, 'soma': soma_layer_density_dict, 'axon': ais_layer_density_dict}

    template_name="DGC"
    np.random.seed(int(17e6))

    cell_count = 0
    LPP_mean_v_amp = 0.
    MPP_mean_v_amp = 0.
    MC_mean_v_amp  = 0.

    for (gid, tree) in trees_iter:
        template_class = eval('h.%s' % template_name)
        cell = make_neurotree_cell (template_class, neurotree_dict=tree, gid=gid)
        print cell

        cell_sec_dict = {'apical': (cell.apical, None), 'basal': (cell.basal, None), 'soma': (cell.soma, None), 'axon': (cell.axon, 50.0)}
        cell_secidx_dict = {'apical': cell.apicalidx, 'basal': cell.basalidx, 'soma': cell.somaidx, 'axon': cell.axonidx}
        syn_dict = distribute_uniform_synapses(gid, swc_type_dict, sec_layer_density_dict, tree, cell_sec_dict, cell_secidx_dict)
        print_syn_summary(gid,syn_dict)

        syn_ids      = syn_dict['syn_ids']
        syn_types    = syn_dict['syn_types']
        swc_types    = syn_dict['swc_types']
        syn_locs     = syn_dict['syn_locs']
        syn_layers   = syn_dict['syn_layers']
        syn_sections = syn_dict['syn_secs']

        exc_inds     = np.where(syn_types == syn_type_excitatory)
        exc_syn_ids  = syn_ids[exc_inds]

        exc_syn_id_base = exc_syn_ids[0]
        
        print 'size exc_syn_ids = ', exc_syn_ids.size
        
        exc_syn_kinetic_params = { 'AMPA': { 't_rise': 0.5, 't_decay': 5.5, 'e_rev': 0. } } 
        exc_syn_obj_dict = mksyns(gid,cell,exc_syn_ids,syn_types,swc_types,syn_locs,syn_sections,exc_syn_kinetic_params,add_synapse=add_shared_synapse,spines=True)
        ampa_syn_obj_dict = { k : v for k, v in exc_syn_obj_dict.iteritems() if 'AMPA' in v }

        
        syn_sample_size = 100
        
        LPP_inds     = np.where((syn_types == syn_type_excitatory) & (syn_layers == 4))
        LPP_syn_ids  = syn_ids[LPP_inds]

        LPP_syn_connection_params = { 'weight': 0.31 }

        LPP_syn_selection=LPP_syn_ids[np.random.randint(0, LPP_syn_ids.size, size=syn_sample_size)]
        
        MPP_inds     = np.where((syn_types == syn_type_excitatory) & (syn_layers == 3))
        MPP_syn_ids  = syn_ids[MPP_inds]

        MPP_syn_connection_params = { 'weight': 0.25 }

        MPP_syn_selection=MPP_syn_ids[np.random.randint(0, MPP_syn_ids.size, size=syn_sample_size)]


        MC_inds     = np.where((syn_types == syn_type_excitatory) & (syn_layers == 2))
        MC_syn_ids  = syn_ids[MC_inds]

        MC_syn_connection_params = { 'weight': 0.19 }

        MC_syn_selection=MC_syn_ids[np.random.randint(0, MC_syn_ids.size, size=syn_sample_size)]

        LPP_mean_v_amp += syn_group_test('LPP', gid, cell, LPP_syn_selection, LPP_syn_ids, LPP_syn_connection_params, v_init, id_base=exc_syn_id_base)
        MPP_mean_v_amp += syn_group_test('MPP', gid, cell, MPP_syn_selection, MPP_syn_ids, MPP_syn_connection_params, v_init, id_base=exc_syn_id_base)
        MC_mean_v_amp += syn_group_test('MC', gid, cell, MC_syn_selection, MC_syn_ids, MC_syn_connection_params, v_init, id_base=exc_syn_id_base)
        cell_count += 1

    LPP_mean_v_amp /= cell_count
    MPP_mean_v_amp /= cell_count
    MC_mean_v_amp /= cell_count
    print("%s synapse: mean V amplitude is %f\n" % ("LPP", LPP_mean_v_amp))
    print("%s synapse: mean V amplitude is %f\n" % ("MPP", MPP_mean_v_amp))
    print("%s synapse: mean V amplitude is %f\n" % ("MC", MC_mean_v_amp))
        
if __name__ == '__main__':
    main(args=sys.argv[(sys.argv.index("DGC_test_synapse.py")+1):])
