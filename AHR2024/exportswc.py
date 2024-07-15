## Instantiates a hoc cell and exports its 3d points to SWC format

import sys
from neuron import h, gui
from collections import defaultdict

h.load_file("nrngui.hoc")
h.load_file("import3d.hoc")

h('objref nil')

layer_dict = {
    'Hilus': 0,
    'GCL':   1,
    'IML':   2,
    'MML':   3,
    'OML':   4,
}

def get_layer(distance, sectype):
    layer = layer_dict['Hilus']
    if sectype == 1:
        layer = layer_dict['GCL']
    elif sectype == 2:
        layer = layer_dict['Hilus']
    elif sectype == 7:
        layer = layer_dict['GCL']
    elif sectype == 8:
        layer = layer_dict['GCL']
    elif sectype == 4: 
        if distance < 50:
            layer = layer_dict['GCL']
        elif distance < 150:
            layer = layer_dict['IML']
        elif distance < 250:
            layer = layer_dict['MML']
        else:
            layer = layer_dict['OML']
    return layer


def export_swc(cell, sections=[("soma",1),("apical",4),("ais",7),("hillock",8),("axon",2)]):
    swc_point_idx = 0
    swc_points = []
    swc_point_sec_dict = defaultdict(list)
    sec_dict = {}
    seen = set([])
    for section, sectype in sections:
        if hasattr(cell, section):
            seclist = list(getattr(cell, section))
            for secidx, sec in enumerate(seclist):
                if hasattr(sec, 'sec'):
                    sec = sec.sec
                if sec in seen:
                    continue
                seen.add(sec)
                n3d = sec.n3d()
                if n3d == 2:
                    x1 = sec.x3d(0)
                    y1 = sec.y3d(0)
                    z1 = sec.z3d(0)
                    d1 = sec.diam3d(0)
                    x2 = sec.x3d(1)
                    y2 = sec.y3d(1)
                    z2 = sec.z3d(1)
                    d2 = sec.diam3d(1)
                    mx = (x2 + x1) / 2.
                    my = (y2 + y1) / 2.
                    mz = (z2 + z1) / 2.
                    dd = d1 - (d1 - d2)/2.
                    sec.pt3dinsert(1, mx, my, mz, dd)
                    n3d = sec.n3d()
                L = sec.L
                for i in range(n3d):
                    x = sec.x3d(i)
                    y = sec.y3d(i)
                    z = sec.z3d(i)
                    d = sec.diam3d(i)
                    ll = sec.arc3d(i)
                    rad = d / 2.
                    loc = ll / L
                    first = True if i == 0 else False
                    swc_point = (swc_point_idx, section, sectype, x, y, z, rad, loc, sec, first)
                    swc_points.append(swc_point)
                    swc_point_sec_dict[sec.name()].append(swc_point)
                    swc_point_idx += 1
    soma_sec = list(cell.soma)[0]
    for swc_point in swc_points:
        (swc_point_idx, section, sectype, x, y, z, rad, loc, sec, first) = swc_point
        parent_idx = -1
        distance_to_soma = h.distance(soma_sec(0.5), sec(loc))
        if not first:
            parent_idx = swc_point_idx-1
        else:
            parent_seg = sec.parentseg()
            if parent_seg is not None:
                parent_x = parent_seg.x
                parent_sec = parent_seg.sec
                parent_points = swc_point_sec_dict[parent_sec.name()]
                for parent_point in parent_points:
                    (parent_point_idx, _, _, _, _, _, _, parent_point_loc, _, _) = parent_point
                    if parent_point_loc >= parent_x:
                        parent_idx = parent_point_idx
                        break
        layer = get_layer(distance_to_soma, sectype)
        print("%d %i %.04f %.04f %.04f %.04f %d %d" % (swc_point_idx, sectype, x, y, z, rad, parent_idx, layer))
    
                    
h.load_file("DGC.hoc")
cell = h.DGC()
#h.topology()
export_swc(cell)
