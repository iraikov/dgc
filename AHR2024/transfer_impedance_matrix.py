import os, sys, logging, click, yaml, pprint
import numpy as np
from neuron import h, gui
import seaborn as sns
import matplotlib.pyplot as plt
from neuron_utils import calcZ,  load_template, run_vclamp

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

                         
def plot_lower_tri_heatmap(df, output="transfer_impedance_matrix.png"):
    mask = np.zeros_like(df, dtype=bool)
    mask[np.triu_indices_from(mask)] = True

    # Want diagonal elements as well
    mask[np.diag_indices_from(mask)] = False

    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=(15, 9))

    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(220, 10, as_cmap=True)

    # Draw the heatmap with the mask and correct aspect ratio
    sns_plot = sns.heatmap(df, mask=mask, cmap=cmap, center=0,
                           square=True, linewidths=.5, cbar_kws={"shrink": .5})
    # save to file
    fig = sns_plot.get_figure()
    fig.savefig(output)
    

# Load the NEURON libraries
h.load_file("stdrun.hoc")
h.load_file("rn.hoc")


config_dict = None
config_path = "config/DG_GC.yaml"
with open(config_path) as f:
    config_dict = yaml.load(f, Loader=yaml.FullLoader)

template_dict = config_dict.get("Template", None)
template_name = None
template_file = None
template = None

template_name = "DGC"
template_file = "DGC.hoc"
template = load_template(template_name, template_file)

toplevel_param_key = "best"
toplevel_param_dict = config_dict.get(toplevel_param_key, None)
if toplevel_param_dict is None:
    raise RuntimeError(f"Unable to read {toplevel_param_key} configuration")

param_key = "p0"    
param_dict = toplevel_param_dict.get(param_key, None)
if param_dict is None:
    param_dict = toplevel_param_dict[int(param_key)]
    
logger.info(f"{pprint.pformat(param_dict)}")
celsius = param_dict["celsius"]

v_init = -80

cell = template(param_dict)

h.v_init = v_init
h.init()
h.finitialize(h.v_init)
h.finitialize(h.v_init)

h.topology()
h.psection(sec=list(cell.soma)[0])

soma_sec = list(cell.soma)[0]
apical_list = list(cell.apical)
prox_list = list(cell.prox)
apical_index = 5
apical_sec = apical_list[apical_index]

logger.info(f"distance to apical section {apical_index}: {h.distance(soma_sec(0.5), apical_sec(0.5))}")
h.psection(sec=list(cell.apical)[apical_index])
logger.info(f"distance to end of proximal section: {h.distance(soma_sec(0.5), prox_list[-1](1.0))}")
h.psection(sec=list(cell.prox)[-1])
h.psection(sec=list(cell.apical)[-1])
h.psection(sec=list(cell.ais)[0])
h.psection(sec=list(cell.axon)[0])




transfer_impedance_list = [soma_sec] + apical_list
N = len(transfer_impedance_list)
transfer_impedance_matrix = np.zeros((N, N))
transfer_impedance_freq = 20

for i in range(N):
    imp = h.Impedance()
    imp.loc(0.5, sec=transfer_impedance_list[i])
    imp.compute(transfer_impedance_freq)
    print(f"computing impedance in section {transfer_impedance_list[i]}")
    for j in range(N):
        if i == j:
            continue
        z = imp.transfer(0.5, sec=transfer_impedance_list[j])
        print(f"transfer impedance between section {transfer_impedance_list[i]} and {transfer_impedance_list[j]} is {z}")
        transfer_impedance_matrix[j][i] = z


        
plot_lower_tri_heatmap(transfer_impedance_matrix)

