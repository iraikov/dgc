import os, sys, logging, click, yaml, pprint
import numpy as np
from neuron import h
from ephys_utils import detect_spikes, measure_passive, measure_rn_from_vclamp
from neuron_utils import calcZ, ic_constant_f, load_template, run_vclamp
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
    v_init=-65.0,
    celsius=36,
    use_coreneuron=False,
):

    # Create the recording vectors for time and voltage
    vec_t = h.Vector()
    vec_soma_v = h.Vector()
    vec_dend_v = h.Vector()
    vec_dend_ica = h.Vector()
    vec_dend_ik = h.Vector()
    vec_soma_ina = h.Vector()
    vec_soma_ik = h.Vector()
    vec_dend_cai = h.Vector()
    vec_dend_g_KAHP = h.Vector()
    vec_dend_g_KCa = h.Vector()
    vec_soma_g_KM = None

    vec_t.record(h._ref_t, record_dt)  # Time
    vec_soma_v.record(cell.soma(0.5)._ref_v, record_dt)  # Voltage
    vec_dend_v.record(cell.dend(0.5)._ref_v, record_dt)  # Voltage
    vec_dend_ica.record(cell.dend(0.5)._ref_ica, record_dt)
    vec_dend_ik.record(cell.dend(0.5)._ref_ik, record_dt)
    vec_soma_ik.record(cell.soma(0.5)._ref_ik, record_dt)
    vec_soma_ina.record(cell.soma(0.5)._ref_ina, record_dt)
    vec_dend_cai.record(cell.dend(0.5)._ref_cai, record_dt)
    vec_dend_g_KAHP.record(cell.dend(0.5)._ref_g_KAHP_PR, record_dt)
    vec_dend_g_KCa.record(cell.dend(0.5)._ref_g_KCa_PR, record_dt)
    if h.ismembrane("KM", sec=cell.soma):
        vec_soma_g_KM = h.Vector()
        vec_soma_g_KM.record(cell.soma(0.5)._ref_g_KM, record_dt)

    # Put an IClamp at the soma
    stim = h.IClamp(0.5, sec=cell.soma)
    stim.delay = t0  # Stimulus stat
    stim.dur = t1 - t0  # Stimulus length
    stim.amp = amp  # strength of current injection

    # Run the Simulation
    h.dt = dt
    h.celsius = celsius
    h.v_init = v_init
    h.init()
    h.finitialize(h.v_init)
    logger.info(f"soma ek = {cell.soma.ek} ena = {cell.soma.ena}")
    logger.info(f"dend ek = {cell.dend.ek} eca = {cell.dend.eca}")
    logger.info(f"soma nao = {cell.soma.nao} ko = {cell.soma.ko}")
    logger.info(
        f"dend cao = {cell.dend.cao} cai = {cell.dend.cai} irest_CaconcPR = {cell.dend.irest_Ca_conc_PR}"
    )

    h.tstop = t_stop
    if use_coreneuron:
        from neuron import coreneuron

        coreneuron.enable = True
    h.run()

    result_dict = {
        "t": np.array(vec_t),
        "soma_v": np.array(vec_soma_v),
        "dend_v": np.array(vec_dend_v),
        "soma_ik": np.array(vec_soma_ik),
        "soma_ina": np.array(vec_soma_ina),
        "dend_ica": np.array(vec_dend_ica),
        "dend_cai": np.array(vec_dend_cai),
        "dend_ik": np.array(vec_dend_ik),
        "dend_g_KAHP": np.array(vec_dend_g_KAHP),
        "dend_g_KCa": np.array(vec_dend_g_KCa),
    }
    if vec_soma_g_KM is not None:
        result_dict["soma_g_KM"] = vec_soma_g_KM

    return result_dict


@click.command()
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
@click.option("--ic-constant", is_flag=True, help="compute ic constant")
@click.option("--passive-features", is_flag=True, help="compute passive features")
@click.option("--target-namespace", "-t", type=str, help="namespace with target information")
@click.option(
    "--stim-amp", type=float, default=0.08, help="amplitude of stimulus current"
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
def main(
    config_path,
    model_variant,
    dt,
    cvode,
    coreneuron,
    ic_constant,
    passive_features,
    target_namespace,
    stim_amp,
    stim_start,
    stim_stop,
    t_stop,
    toplevel_param_key,
    param_key,
    v_init_config,
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
        template_name = "PR_nrn"
        template_file = "PR_nrn.hoc"
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
    cell.init_ic(h.v_init)
    ic_constant_0 = cell.soma.ic_constant
    ic_constant_val = ic_constant_0

    # Obtain value for ic_constant such that RMP = v_init
    if ic_constant:
        try:
            x0, res = optimize.brentq(
                ic_constant_f,
                -0.1,
                0.1,
                args=(template, param_dict, ic_constant_0, h.v_init),
                xtol=1e-6,
                maxiter=200,
                disp=False,
                full_output=True,
            )
        except ValueError:
            x0 = 0.
        else:
            if not res.converged:
                x0 = 0.
        ic_constant_val = x0 + ic_constant_0
    else:
        if v_init_config == "rest":
            ic_constant_val = param_dict["ic_constant_rest"]
        elif v_init_config == "hold":
            ic_constant_val = param_dict["ic_constant_hold"]
        else:
            raise RuntimeError(f"Unknown v_init configuration {v_init}")

    cell.soma.ic_constant = ic_constant_val
    h.finitialize(h.v_init)
    h.finitialize(h.v_init)

    h.psection(sec=cell.soma)
    h.psection(sec=cell.dend)

    initial_v_error = ic_constant_f(
        0.0, template, param_dict, cell.soma.ic_constant, v_hold=v_init, use_cvode=cvode
    )
    logger.info(f"initial_v_error: {initial_v_error}")
    initial_v_constr = (
        1 if np.isclose(0.0, initial_v_error, rtol=1e-4, atol=1e-4) else 0
    )

    logger.info(f"ic_constant0: {ic_constant_0}")
    logger.info(f"ic_constant: {ic_constant_val}")
    logger.info(f"mean initial vm constraint: {initial_v_constr}")

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

    iclamp_results = run_iclamp(
        cell,
        stim_amp,
        stim_start,
        stim_stop,
        t_stop=t_stop,
        v_init=v_init,
        dt=dt,
        use_coreneuron=coreneuron,
    )
    vec_t = iclamp_results["t"]
    vec_soma_v = iclamp_results["soma_v"]
    vec_dend_v = iclamp_results["dend_v"]
    vec_dend_ik = iclamp_results["dend_ik"]
    vec_soma_ik = iclamp_results["soma_ik"]
    vec_soma_ina = iclamp_results["soma_ina"]
    vec_dend_ik = iclamp_results["dend_ik"]
    vec_dend_ica = iclamp_results["dend_ica"]
    vec_dend_cai = iclamp_results["dend_cai"]
    vec_dend_g_KAHP = iclamp_results["dend_g_KAHP"]
    vec_dend_g_KCa = iclamp_results["dend_g_KCa"]
    vec_soma_g_KM = iclamp_results.get("soma_g_KM", None)
    logger.info(f"spikes: {detect_spikes(vec_t, vec_soma_v, stim_start, stim_stop)}")

    nrn_type = config_dict["Celltype"]

    fig, axs = plt.subplots(5, 2)
    axs[0, 0].plot(vec_t, vec_dend_v, linewidth=3, color="r", label="dend_v")
    axs[0, 0].set_xlabel("Time (ms)")
    axs[0, 0].set_ylabel("V (mV)")
    axs[1, 0].plot(vec_t, vec_dend_ica, linewidth=3, color="r", label="dend_ica")
    axs[2, 0].plot(vec_t, vec_dend_cai, linewidth=3, color="g", label="dend_cai")
    axs[3, 0].plot(vec_t, vec_dend_ik, linewidth=3, color="b", label="dend_ik")
    axs[4, 0].plot(vec_t, vec_dend_g_KAHP, linewidth=3, color="b", label="dend_g_KAHP")
    axs[-1, 0].set_xlabel("Time (ms)")
    axs[0, 1].plot(vec_t, vec_soma_v, linewidth=3, color="r", label="soma_v")
    axs[1, 1].plot(vec_t, vec_soma_ina, linewidth=3, color="r", label="soma_ina")
    axs[2, 1].plot(vec_t, vec_soma_ik, linewidth=3, color="b", label="soma_ik")
    axs[3, 1].plot(vec_t, vec_dend_g_KCa, linewidth=3, color="b", label="dend_g_KCa")
    if vec_soma_g_KM is not None:
        axs[4, 1].plot(vec_t, vec_soma_g_KM, linewidth=3, color="b", label="soma_g_KM")
    axs[-1, 1].set_xlabel("Time (ms)")
    axs[0, 0].legend()
    axs[1, 0].legend()
    axs[2, 0].legend()
    axs[3, 0].legend()
    axs[4, 0].legend()
    axs[0, 1].legend()
    axs[1, 1].legend()
    axs[2, 1].legend()
    axs[3, 1].legend()
    axs[4, 1].legend()
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
