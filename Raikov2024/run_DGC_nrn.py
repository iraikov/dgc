import os, sys, logging, click, yaml, pprint
import numpy as np
from neuron import h
from ephys_utils import detect_spikes, measure_passive, measure_rn_from_vclamp
from neuron_utils import calcZ,  load_template, run_vclamp
from scipy import optimize

import matplotlib.pyplot as plt

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def list_find(f, lst):
    """

    :param f:
    :param lst:
    :return:
    """
    i = 0
    for x in lst:
        if f(x):
            return i
        else:
            i = i + 1
    return None


def run_iclamp(
    cell,
    amp,
    t0,
    t1,
    dt=0.025,
    record_dt=0.01,
    t_stop=1000.0,
    v_init=-75.0,
    celsius=25,
    apical_index=119,
    use_coreneuron=False,
):

    # Create the recording vectors for time and voltage
    vec_t = h.Vector()
    vec_soma_v = h.Vector()
    vec_ais_v = h.Vector()
    vec_apical_v = h.Vector()
    vec_apical_ica = h.Vector()
    vec_apical_ik = h.Vector()
    vec_apical_ina = h.Vector()
    vec_soma_ik = h.Vector()
    vec_soma_ina = h.Vector()
    vec_apical_cai = h.Vector()
    vec_apical_ki = h.Vector()
    vec_soma_ki = h.Vector()
    vec_ais_ki = h.Vector()
    vec_soma_nai = h.Vector()
    vec_soma_cai = h.Vector()

    soma = list(cell.soma)
    ais = list(cell.ais)
    apical = list(cell.apical)
    
    vec_t.record(h._ref_t, record_dt)  # Time
    vec_soma_v.record(soma[0](0.5)._ref_v, record_dt)  # Voltage
    vec_ais_v.record(ais[0](0.5)._ref_v, record_dt)  # Voltage
    vec_apical_v.record(apical[apical_index](0.5)._ref_v, record_dt)  # Voltage
    #vec_apical_ica.record(apical[apical_index](0.5)._ref_ica, record_dt)
    vec_apical_ik.record(apical[apical_index](0.5)._ref_ik, record_dt)
    vec_apical_ina.record(apical[apical_index](0.5)._ref_ina, record_dt)
    vec_soma_ik.record(soma[0](0.5)._ref_ik, record_dt)
    vec_soma_ina.record(soma[0](0.5)._ref_ina, record_dt)
    #vec_apical_cai.record(apical[apical_index](0.5)._ref_cai, record_dt)
    vec_apical_ki.record(apical[apical_index](0.5)._ref_ki, record_dt)
    vec_soma_ki.record(soma[0](0.5)._ref_ki, record_dt)
    vec_ais_ki.record(ais[0](0.5)._ref_ki, record_dt)
    vec_soma_nai.record(soma[0](0.5)._ref_nai, record_dt)
    vec_soma_cai.record(soma[0](0.5)._ref_cai, record_dt)

    # Put an IClamp at the soma
    stim = h.IClamp(0.5, sec=soma[0])
    stim.delay = t0  # Stimulus stat
    stim.dur = t1 - t0  # Stimulus length
    stim.amp = amp  # strength of current injection

    # Run the Simulation
    h.dt = dt
    h.celsius = celsius
    h.v_init = v_init
    h.init()
    h.finitialize(h.v_init)
    logger.info(f"soma ek = {soma[0].ek} ena = {soma[0].ena} eca = {soma[0].eca}")
    logger.info(f"dend ek = {apical[0].ek} ena = {apical[0].ena}")
    logger.info(f"soma nao = {soma[0].nao} ko = {soma[0].ko}")

    h.tstop = t_stop
    if use_coreneuron:
        from neuron import coreneuron

        coreneuron.enable = True
    h.run()

    result_dict = {
        "t": np.array(vec_t),
        "soma_v": np.array(vec_soma_v),
        "ais_v": np.array(vec_ais_v),
        "dend_v": np.array(vec_apical_v),
        "soma_ik": np.array(vec_soma_ik),
        "soma_ki": np.array(vec_soma_ki),
        "ais_ki": np.array(vec_ais_ki),
        "soma_ina": np.array(vec_soma_ina),
        #"dend_ica": np.array(vec_apical_ica),
        #"dend_cai": np.array(vec_apical_cai),
        "dend_ki": np.array(vec_apical_ki),
        "dend_ik": np.array(vec_apical_ik),
        "dend_ina": np.array(vec_apical_ina),
        "soma_nai": np.array(vec_soma_nai),
        "soma_cai": np.array(vec_soma_cai),
    }

    return result_dict


@click.command()
@click.option("--apical-index", type=int, default=119, help="apical section number to record")
@click.option(
    "--config-path",
    "-c",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    help="path to configuration file",
)
@click.option("--model-variant", "-m", default="default", type=str)
@click.option("--dt", type=float, default=0.025, help="default simulation time step")
@click.option("--cvode/--no-cvode", default=False, help="use adaptive time step solver")
@click.option("--coreneuron/--no-coreneuron", default=False, help="use CoreNEURON")
@click.option("--passive-features", is_flag=True, help="compute passive features")
@click.option("--target-namespace", "-n", type=str, help="namespace with target information")
@click.option(
    "--stim-amp", type=float, default=0.1, help="amplitude of stimulus current"
)
@click.option("--stim-start", type=float, default=500.0, help="start time of stimulus")
@click.option("--stim-stop", type=float, default=1000.0, help="stop time of stimulus")
@click.option(
    "--t-stop", "-t", type=float, default=2000.0, help="stop time of simulation"
)
@click.option(
    "--toplevel-param-key",
    type=str,
    default="best",
    help="key for toplevel section with parameters",
)
@click.option(
    "--param-key", "-p", type=str, required=True, help="key for section with parameters"
)
@click.option(
    "--v-rest",
    "v_init_config",
    flag_value="rest",
    default=True,
    help="initialize cell to resting potential",
)
@click.option(
    "--v-hold",
    "v_init_config",
    flag_value="hold",
    help="initialize cell to holding potential",
)
@click.option(
    "--v-clamp-hold",
    type=float,
    help="voltage clamp holding potential",
)
@click.option(
    "--v-clamp-rest",
    type=float,
    default=-120,
    help="voltage clamp rest potential",
)
def main(
    apical_index,
    config_path,
    model_variant,
    dt,
    cvode,
    coreneuron,
    passive_features,
    target_namespace,
    stim_amp,
    stim_start,
    stim_stop,
    t_stop,
    toplevel_param_key,
    param_key,
    v_init_config,
    v_clamp_hold,
    v_clamp_rest,
):

    # Load the NEURON libraries
    h.load_file("stdrun.hoc")
    h.load_file("rn.hoc")

    # Enable variable time step solver
    h.cvode.use_fast_imem(1)
    h.cvode.cache_efficient(1)
    h.cvode.active(1 if cvode else 0)
    h.secondorder = 2
    h.dt = dt

    config_dict = None
    with open(config_path) as f:
        config_dict = yaml.load(f, Loader=yaml.FullLoader)

    template_dict = config_dict.get("Template", None)
    template_name = None
    template_file = None
    template = None
    if template_dict is None:
        template_name = "DGC"
        template_file = "DGC.hoc"
    else:
        if model_variant in template_dict:
            template_name = template_dict[model_variant]["name"]
            template_file = template_dict[model_variant].get("file", None)
        else:
            raise RuntimeError(f"Unknown model variant {model_variant}")
    template = load_template(template_name, template_file)

    
    target_config = config_dict["Targets"]
    if target_namespace:
        target_config = config_dict["Target namespaces"][target_namespace]
    Rin_config = target_config["Rin"]
    tau0_config = target_config["tau0"]

    toplevel_param_dict = config_dict.get(toplevel_param_key, None)
    if toplevel_param_dict is None:
        raise RuntimeError(f"Unable to read {toplevel_param_key} configuration")

    param_dict = toplevel_param_dict.get(param_key, None)
    if param_dict is None:
        param_dict = toplevel_param_dict[int(param_key)]
    logger.info(f"{pprint.pformat(param_dict)}")
    celsius = param_dict["celsius"]

    if v_init_config == "rest":
        v_init = config_dict["Targets"]["V_rest"]["val"]
    elif v_init_config == "hold":
        v_init = config_dict["Targets"]["V_hold"]["val"]
    else:
        raise RuntimeError(f"Unknown v_init configuration {v_init}")

    cell = template(param_dict)

    h.v_init = v_init
    h.init()
    h.finitialize(h.v_init)
    h.finitialize(h.v_init)

    h.psection(sec=list(cell.soma)[0])
    soma_sec = list(cell.soma)[0]
    apical_list = list(cell.apical)
    apical_sec = apical_list[apical_index]
    logger.info(f"distance to apical section {apical_index}: {h.distance(soma_sec(0.5), apical_sec(0.5))}")
    h.psection(sec=list(cell.apical)[apical_index])
    h.psection(sec=list(cell.ais)[0])
    h.psection(sec=list(cell.axon)[0])

    if passive_features:
        if "I" in Rin_config:
            Rin_amp = Rin_config["I"][0] * Rin_config.get("I_factor", 1.0)
            iclamp_results = run_iclamp(
                cell,
                Rin_amp,
                stim_start,
                stim_stop,
                t_stop=10000.0,
                v_init=v_init,
                dt=dt,
                celsius=celsius,
                use_coreneuron=coreneuron,
            )
            t = iclamp_results["t"]
            v = iclamp_results["soma_v"]
            passive_features = measure_passive(t, v, stim_start, stim_stop, Rin_amp)
            Rinp = passive_features["Rinp"]
            tau = passive_features["tau"]
            
        elif "V" in Rin_config:
            Rin_amp = np.asarray(Rin_config["V"]) * Rin_config.get("V_factor", 1.0)
            Rin_ts = Rin_config["t"]
            vclamp_results = run_vclamp(
                cell,
                Rin_amp,
                Rin_ts,
                v_init=v_init,
                dt=dt,
                use_coreneuron=coreneuron,
                celsius=celsius,
            )
            vclamp_results["t0"] = Rin_ts[0]
            vclamp_results["t1"] = Rin_ts[1]
            vclamp_results["t2"] = Rin_ts[2]
            Rinp = measure_rn_from_vclamp(**vclamp_results)
            tau_amp = tau0_config["I"][0] * tau0_config.get("I_factor", 1.0)
            iclamp_results = run_iclamp(
                cell,
                tau_amp,
                stim_start,
                stim_stop,
                celsius=celsius,
                t_stop=10000.0,
                v_init=v_init,
                dt=dt,
                use_coreneuron=coreneuron,
            )
            t = iclamp_results["t"]
            v = iclamp_results["soma_v"]
            passive_features = measure_passive(t, v, stim_start, stim_stop, tau_amp)
            tau = passive_features["tau"]
            
        else:
            raise RuntimeError("Unknown configuration for Rinp")
                
        logger.info(f"soma input resistance: {Rinp} time constant: {tau}")
        # logger.info(f'analytical input and transfer impedance per segment: {pprint.pformat(calcZ(cell))}')

    nrn_type = config_dict["Celltype"]

    if v_clamp_hold is not None:
        V_amp = np.asarray([v_clamp_rest, v_clamp_hold, v_clamp_rest])
        V_ts  = np.asarray([250, 1000, 1250])
        vclamp_results = run_vclamp(
            cell,
            V_amp,
            V_ts,
            sec=apical_sec,
            ion_current_names=["k", "na"],
            v_init=v_init,
            dt=dt,
            use_coreneuron=coreneuron,
            celsius=celsius,
        )

        vec_t = vclamp_results["t"]
        vec_v = vclamp_results["v"]
        vec_ik = vclamp_results["ik"]
        vec_ina = vclamp_results["ina"]

        index_range = np.argwhere(np.logical_and(vec_t >= V_ts[1] - 1,
                                                 vec_t <= V_ts[1],))
        
        nrows = 3
        ncols = 1
        fig, axs = plt.subplots(nrows, ncols)
#        axs[0].plot(vec_t[index_range], vec_v[index_range], linewidth=3, color="r", label="v")
#        axs[1].plot(vec_t[index_range], vec_ik[index_range], linewidth=3, color="r", label="dend_ik")
#        axs[2].plot(vec_t[index_range], vec_ina[index_range], linewidth=3, color="b", label="dend_ina")
        axs[0].plot(vec_t, vec_v, linewidth=3, color="r", label="v")
        axs[1].plot(vec_t, vec_ik, linewidth=3, color="r", label="dend_ik")
        axs[2].plot(vec_t, vec_ina, linewidth=3, color="b", label="dend_ina")
        axs[-1].set_xlabel("Time (ms)")
        axs[0].set_ylabel("V (mV)")
        axs[1].set_ylabel("I (nA)")
        
        for i in range(nrows):
            for j in range(ncols):
                axs[i].legend()
        plt.savefig(f"{nrn_type}_vclamp.svg")
        plt.show()
        
        
    else:

        iclamp_results = run_iclamp(
            cell,
            stim_amp,
            stim_start,
            stim_stop,
            apical_index=apical_index,
            t_stop=t_stop,
            v_init=v_init,
            celsius=celsius,
            dt=dt,
            use_coreneuron=coreneuron,
        )
        vec_t = iclamp_results["t"]
        vec_soma_v = iclamp_results["soma_v"]
        vec_ais_v = iclamp_results["ais_v"]
        vec_apical_v = iclamp_results["dend_v"]
        vec_apical_ik = iclamp_results["dend_ik"]
        vec_soma_ik = iclamp_results["soma_ik"]
        vec_soma_ina = iclamp_results["soma_ina"]
        vec_apical_ina = iclamp_results["dend_ina"]
        #vec_apical_ica = iclamp_results["dend_ica"]
        #vec_apical_cai = iclamp_results["dend_cai"]
        vec_apical_ki = iclamp_results["dend_ki"]
        vec_soma_ki = iclamp_results["soma_ki"]
        vec_ais_ki = iclamp_results["ais_ki"]
        vec_soma_nai = iclamp_results["soma_nai"]
        vec_soma_cai = iclamp_results["soma_cai"]

        logger.info(f"spikes: {detect_spikes(vec_t, vec_soma_v, stim_start, stim_stop+5.0)}")
        logger.info(f"dend ki min/max: {np.min(vec_apical_ki)} / {np.max(vec_apical_ki)}")
        
        nrows = 6
        ncols = 2
        fig, axs = plt.subplots(nrows, ncols)
        axs[0, 0].plot(vec_t, vec_soma_v, linewidth=3, color="r", label="soma_v")
        axs[0, 0].plot(vec_t, vec_ais_v, linewidth=3, color="b", label="ais_v")
        axs[1, 0].plot(vec_t, vec_apical_v, linewidth=3, color="r", label="dend_v")
        axs[2, 0].plot(vec_t, vec_soma_ina, linewidth=3, color="b", label="soma_ina")
        axs[3, 0].plot(vec_t, vec_soma_ik, linewidth=3, color="b", label="soma_ik")
        axs[4, 0].plot(vec_t, vec_apical_ik, linewidth=3, color="b", label="dend_ik")
        axs[5, 0].plot(vec_t, vec_apical_ina, linewidth=3, color="b", label="dend_ina")
        axs[-1, 0].set_xlabel("Time (ms)")
        axs[0, 0].set_ylabel("V (mV)")
        
        axs[0, 1].plot(vec_t, vec_soma_ki, linewidth=3, color="g", label="soma_ki")
        axs[1, 1].plot(vec_t, vec_ais_ki, linewidth=3, color="g", label="ais_ki")
        axs[2, 1].plot(vec_t, vec_apical_ki, linewidth=3, color="g", label="dend_ki")
        axs[3, 1].plot(vec_t, vec_soma_nai, linewidth=3, color="g", label="soma_nai")
        axs[4, 1].plot(vec_t, vec_soma_cai, linewidth=3, color="g", label="soma_cai")

        for i in range(nrows):
            for j in range(ncols):
                axs[i, j].legend()
        plt.savefig(f"{nrn_type}_iclamp.svg")
        plt.show()


if __name__ == "__main__":
    main(
        args=sys.argv[
            (
                list_find(
                    lambda x: os.path.basename(x) == os.path.basename(__file__),
                    sys.argv,
                )
                + 1
            ) :
        ]
    )
