import os, sys, math, logging, yaml, click
from functools import partial
from neuron import h
from neuron_utils import ic_constant_f, run_iclamp, load_template
from mpi4py import MPI
import numpy as np
from numpy.random import default_rng
from dmosopt import dmosopt
import ephys_utils as ephys
from protocol import ExperimentalProtocol
from scipy import optimize, signal

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def range_distance(x, lb, ub):
    # Returns 0. if x is within the range [lb, ub], otherwise returns the smaller of the distance between x and lb, ub
    return 0.0 if (x >= lb) and (x <= ub) else min(abs(x - lb), abs(x - ub))


def lb_distance(x, lb, ub):
    # Returns 0. if x >= lb, otherwise returns the the distance between x and lb
    return 0.0 if (x >= lb) else lb - x


def init_cell(template_name, pp, v_hold=-60, celsius=36.0, ic_constant_val=None):

    h.cvode.use_fast_imem(1)
    h.cvode.cache_efficient(1)
    h.secondorder = 2

    h.celsius = celsius

    # Create the cell
    template = getattr(h, template_name)
    cell = template(pp)

    # Initialize cell
    h.v_init = v_hold
    h.init()

    if ic_constant_val is None:
        cell.init_ic(h.v_init)
        ic_constant_0 = cell.soma.ic_constant

        # Obtain value for ic_constant such that RMP = v_hold
        x0 = 0.0
        ic_constant0 = ic_constant_0
        try:
            x0, res = optimize.brentq(
                ic_constant_f,
                -0.5,
                0.5,
                args=(template, pp, ic_constant_0, h.v_init),
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
        
        ic_constant_val = ic_constant_0 + x0

    cell.soma.ic_constant = ic_constant_val
    h.finitialize(h.v_init)
    h.finitialize(h.v_init)

    return cell


def make_obj_fun(protocol_config_dict, feature_dtypes, template_name, target_namespace, worker):

    exp_protocol = ExperimentalProtocol(protocol_config_dict, target_namespace=target_namespace)

    return partial(obj_fun, exp_protocol, feature_dtypes, template_name)


# This is the function which is going to be minimized
def obj_fun(exp_protocol, feature_dtypes, template_name, pp):

    template = getattr(h, template_name)

    cell = init_cell(template_name, pp, v_hold=exp_protocol.v_hold)
    ic_constant_hold = cell.soma.ic_constant

    # Check whether the initial voltage constraint was satisfied
    initial_v_error_hold = float(
        ic_constant_f(0.0, template, pp, ic_constant_hold, v_hold=exp_protocol.v_hold)
    )

    initial_v_constr = 1 if abs(initial_v_error_hold) < 1.0 else -1
    logger.info(f'ic_constant check: {initial_v_error_hold} constr: {initial_v_constr}')

    cell = init_cell(
        template_name, pp, v_hold=exp_protocol.v_hold, ic_constant_val=ic_constant_hold
    )

    rn, tau = np.nan, np.nan
    iclamp_results = None
    vclamp_results = None
    # Measure input resistance and membrane time constant
    if initial_v_constr > 0:
        # Run single current injection to measure subthreshold features
        try:
            if exp_protocol.rn_exp_type == "iclamp":
                iclamp_results = exp_protocol.run_iclamp(cell, target="Rin", tstop=3000.0)
            elif exp_protocol.rn_exp_type == "vclamp":
                vclamp_results = exp_protocol.run_vclamp(cell, target="Rin")
                iclamp_results = exp_protocol.run_iclamp(cell, target="tau0", tstop=3000.0)
        except:
            pass
        else:
            passive_results = ephys.measure_passive(**iclamp_results)
            rn = passive_results["Rinp"]
            tau = passive_results["tau"]
            if vclamp_results is not None:
                rn = ephys.measure_rn_from_vclamp(**vclamp_results)

    target_rn = exp_protocol.target_rn
    target_tau = exp_protocol.target_tau
    rn_obj_value = range_distance(rn, target_rn[0], target_rn[1]) ** 2
    tau_obj_value = range_distance(tau, target_tau[0], target_tau[1]) ** 2

    tau_constr = 1 if ((tau > 0.0) and (tau < 1000.0)) else -1
    rn_constr = 1 if ((rn > 0.0) and (rn < 1000.0)) else -1

    # Run iclamp experiments
    iclamp_results = []
    cell = init_cell(
        template_name, pp, v_hold=exp_protocol.v_hold, ic_constant_val=ic_constant_hold
    )
    iclamp_results = exp_protocol.run_iclamp_steps(cell)


    # Measure spike features
    (
        pre_spk_cnt,
        spk_cnt,
        spk_infos,
        thresholds,
        mean_spike_amplitudes,
    ) = ephys.measure_spike_features(
        iclamp_results, exp_protocol.exp_i_inj_t0_f_I, exp_protocol.exp_i_inj_t1_f_I+2.0
    )

    pre_spk_count_constr = -1 if np.sum(pre_spk_cnt) > 0 else 1

    ISI_values = ephys.measure_ISI(exp_protocol.exp_i_inj_amp_f_I, spk_infos)

    ISI_adaptation_dists = list(
        map(
            lambda ratio, target_range: range_distance(
                ratio * 100.0, target_range[0] * 100.0, target_range[1] * 100.0
            ),
            ISI_values["ratio"],
            zip(
                exp_protocol.exp_i_lb_spk_adaptation,
                exp_protocol.exp_i_ub_spk_adaptation,
            ),
        )
    )
    ISI_adaptation_obj_value = np.mean([dist**2 for dist in ISI_adaptation_dists])
    ISI_adaptation_constr = -1 if np.isnan(ISI_adaptation_obj_value) else 1

    first_ISI_constr = (
        1 if np.all(ISI_values["first"] > exp_protocol.exp_first_ISI_lower) else -1
    )

    fI_values = ephys.measure_fI(
        spk_cnt,
        exp_protocol.exp_i_inj_t0_f_I,
        exp_protocol.exp_i_inj_t1_f_I,
        exp_protocol.exp_i_inj_amp_f_I,
    )

    fI_mean_target_rate_diff = np.mean(
        [
            (target_rate - rate) ** 2
            for rate, target_rate in zip(
                fI_values["frequency"], exp_protocol.exp_i_mean_rate_f_I
            )
        ]
    )

    fI_range_dists = list(
        map(
            lambda rate, target_range: range_distance(
                rate, target_range[0], target_range[1]
            ),
            fI_values["frequency"],
            zip(exp_protocol.exp_i_lb_rate_f_I, exp_protocol.exp_i_ub_rate_f_I),
        )
    )
    fI_obj_value = np.mean([dist**2 for dist in fI_range_dists])

    # Compute objectives
    mean_spike_amplitude_range_dists = list(
        map(
            lambda amp, target_amp: None
            if np.isnan(target_amp[0])
            else lb_distance(amp, target_amp[0], target_amp[1]),
            mean_spike_amplitudes,
            zip(exp_protocol.exp_i_lb_spk_amp, exp_protocol.exp_i_ub_spk_amp),
        )
    )
    mean_spike_amplitude_obj_value = np.mean(
        [
            dist**2
            for dist in filter(
                lambda x: False if x is None else True, mean_spike_amplitude_range_dists
            )
        ]
    )
    spike_amplitude_constr = -1 if np.isnan(mean_spike_amplitude_obj_value) else 1

    # Check for fI monotonicity
    fI_rate_diff = np.diff(fI_values["frequency"][:-1])
    monotonic_fI_constr = 1 if np.all(fI_rate_diff > 0) else -1

    # Obtain ic_constant for v_rest target
    cell = init_cell(template_name, pp, v_hold=exp_protocol.v_rest)
    ic_constant_rest = cell.soma.ic_constant

    # Pass to dmosopt
    feature_list = [
            (
                ic_constant_hold,
                ic_constant_rest,
                initial_v_error_hold, 
                rn,
                tau,
                fI_values,
                fI_mean_target_rate_diff,
                ISI_values,
                thresholds,
                mean_spike_amplitudes,
            )
        ]

    feature_values = np.asarray(
        [
            (
                ic_constant_hold,
                ic_constant_rest,
                initial_v_error_hold,
                rn,
                tau,
                fI_values,
                fI_mean_target_rate_diff,
                ISI_values,
                thresholds,
                mean_spike_amplitudes,
            )
        ],
        dtype=np.dtype(feature_dtypes),
    )
    obj_values = np.asarray(
        [
            rn_obj_value,
            tau_obj_value,
            fI_obj_value,
            mean_spike_amplitude_obj_value,
            ISI_adaptation_obj_value,
        ],
        dtype=np.float32,
    )
    constr_values = np.asarray(
        [
            monotonic_fI_constr,
            rn_constr,
            tau_constr,
            spike_amplitude_constr,
            first_ISI_constr,
            ISI_adaptation_constr,
            pre_spk_count_constr,
            initial_v_constr,
        ],
        dtype=np.float32,
    )

    return obj_values, feature_values, constr_values


script_name = os.path.basename(__file__)


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


@click.command()
@click.option(
    "--config-path",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@click.option(
    "--target-namespace", "-t",
    required=False,
    type=str,
)
@click.option(
    "--results-path",
    "-p",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
)
@click.option("--model-variant", "-m", default="default", type=str)
@click.option("--num-epochs", default=10, type=int)
@click.option("--num-initial", default=800, type=int)
@click.option("--population-size", default=400, type=int)
@click.option("--num-generations", default=400, type=int)
@click.option("--seed", default=None, type=int)
@click.option("--initial-method", default="slh", type=str)
@click.option("--optimizer", default="nsga2", type=str)
@click.option("--resample-fraction", default=1.0, type=float)
@click.option("--surrogate-method", default="gpr", type=str)
@click.option("--save-surrogate-eval", is_flag=True)
@click.option(
    "--sensitivity/--no-sensitivity",
    default=False,
    help="enable parameter sensitivity settings",
)
@click.option(
    "--feasibility/--no-feasibility",
    default=False,
    help="enable parameter feasibility prediction",
)
@click.option("--interactive", is_flag=True)
@click.option("--verbose", "-v", is_flag=True)
def main(
    config_path,
    target_namespace,
    results_path,
    model_variant,
    num_epochs,
    num_initial,
    population_size,
    num_generations,
    seed,
    initial_method,
    optimizer,
    resample_fraction,
    surrogate_method,
    save_surrogate_eval,
    sensitivity,
    feasibility,
    interactive,
    verbose,
):

    comm = MPI.COMM_WORLD
    protocol_config_dict = None
    if comm.rank == 0:
        with open(config_path) as f:
            protocol_config_dict = yaml.load(f, Loader=yaml.FullLoader)

    local_random = None
    if seed is not None:
        local_random = default_rng(seed=seed)

    protocol_config_dict = comm.bcast(protocol_config_dict)

    celltype = protocol_config_dict["Celltype"]
    template_dict = protocol_config_dict.get("Template", None)
    template_name = None
    template_file = None
    template = None
    if template_dict is None:
        template_name = "GC_nrn"
        template_file = "GC_nrn.hoc"
    else:
        if model_variant in template_dict:
            template_name = template_dict[model_variant]["name"]
            template_file = template_dict[model_variant].get("file", None)
        else:
            raise RuntimeError(f"Unknown model variant {model_variant}")
    template = load_template(template_name, template_file)

    N_exp = len(protocol_config_dict["Targets"]["f_I"]["I"])
    if target_namespace is not None:
        N_exp = len(protocol_config_dict["Target namespaces"][target_namespace]["f_I"]["I"])
    feature_dtypes = [
        (
            "ic_constant_hold",
            np.float32,
        ),
        (
            "ic_constant_rest",
            np.float32,
        ),
        (
            'initial_v_error_hold',
            np.float32,
        ), 
        (
            "rn",
            np.float32,
        ),
        (
            "tau",
            np.float32,
        ),
        ("fI", ephys.fi_value_dtype, N_exp),
        ("mean_fI_diff", np.float32),
        ("ISI", ephys.isi_value_dtype, N_exp),
        ("threshold", np.float32, N_exp),
        ("spike_amplitude", np.float32, N_exp),
    ]

    problem_parameters = protocol_config_dict["Parameters"]
    variant_parameters_dict = protocol_config_dict.get("Variant Parameters", {})
    if model_variant in variant_parameters_dict:
        variant_parameters = variant_parameters_dict[model_variant]
        for k in variant_parameters:
            problem_parameters[k] = variant_parameters[k]

    space = protocol_config_dict["Space"]
    variant_space_dict = protocol_config_dict.get("Variant Space", {})
    if model_variant in variant_space_dict:
        variant_space = variant_space_dict[model_variant]
        for k in variant_space:
            space[k] = variant_space[k]

    space_sensitivity = None
    if sensitivity:
        space_sensitivity = protocol_config_dict.get("Space Sensitivity", None)
        variant_sensitivity_dict = protocol_config_dict.get(
            "Variant Space Sensitivity", {}
        )

        if model_variant in variant_sensitivity_dict:
            variant_sensitivity = variant_sensitivity_dict[model_variant]
            for k in variant_sensitivity:
                space_sensitivity[k] = variant_sensitivity[k]

    objective_names = [
        "rn_error",
        "tau_error",
        "fI_error",
        "spike_amplitude_error",
        "ISI_adaptation_error",
    ]
    constraint_names = [
        "monotonic_fI",
        "rn_constr",
        "tau_constr",
        "spike_amplitude_constr",
        "first_ISI_constr",
        "ISI_adaptation_constr",
        "pre_spk_count",
        "initial_v_constr",
    ]

    exp_protocol = ExperimentalProtocol(protocol_config_dict,
                                        target_namespace=target_namespace)

    N_spk_amp = min(
        len(exp_protocol.exp_i_lb_spk_amp), len(exp_protocol.exp_i_inj_amp_f_I)
    )
    N_spk_adpt = min(
        len(exp_protocol.exp_i_lb_spk_adaptation), len(exp_protocol.exp_i_inj_amp_f_I)
    )
    obj_targets = {
        "rn": (np.asarray(exp_protocol.target_rn, dtype=np.float32), np.float32, 2),
        "tau": (np.asarray(exp_protocol.target_tau, dtype=np.float32), np.float32, 2),
        "ISI_adaptation": (
            np.row_stack(
                (
                    exp_protocol.exp_i_inj_amp_f_I[:N_spk_adpt],
                    exp_protocol.exp_i_lb_spk_adaptation[:N_spk_adpt],
                    exp_protocol.exp_i_ub_spk_adaptation[:N_spk_adpt],
                )
            ),
            np.float32,
            (3, N_exp),
        ),
        "fI": (
            np.row_stack(
                (
                    exp_protocol.exp_i_inj_amp_f_I,
                    exp_protocol.exp_i_lb_rate_f_I,
                    exp_protocol.exp_i_ub_rate_f_I,
                )
            ),
            np.float32,
            (3, N_exp),
        ),
        "spike_amplitude": (
            np.row_stack(
                (
                    exp_protocol.exp_i_inj_amp_f_I[:N_spk_amp],
                    exp_protocol.exp_i_lb_spk_amp[:N_spk_amp],
                    exp_protocol.exp_i_ub_spk_amp[:N_spk_amp],
                )
            ),
            np.float32,
            (3, N_exp),
        ),
    }
    problem_metadata = np.array(
        [tuple((obj_targets[k][0] for k in sorted(obj_targets)))],
        dtype=[
            (f"{k}_target", obj_targets[k][1], obj_targets[k][2])
            for k in sorted(obj_targets)
        ],
    )

    # Create an optimizer
    dmosopt_params = {
        "opt_id": f"dmosopt_{celltype}_neuron",
        "obj_fun_init_name": "make_obj_fun",
        "obj_fun_init_module": "dmosopt_GC_nrn",
        "obj_fun_init_args": {
            "protocol_config_dict": protocol_config_dict,
            "feature_dtypes": feature_dtypes,
            "template_name": template_name,
            "target_namespace": target_namespace,
        },
        "problem_parameters": problem_parameters,
        "space": space,
        "objective_names": objective_names,
        "constraint_names": constraint_names,
        "feature_dtypes": feature_dtypes,
        "optimizer": optimizer,
        "n_initial": num_initial,
        "n_epochs": num_epochs,
        "population_size": population_size,
        "num_generations": num_generations,
        "termination_conditions": True,
        "resample_fraction": resample_fraction,
        "initial_maxiter": 10,
        "initial_method": initial_method,
        "surrogate_method": surrogate_method,
        "surrogate_options": {"lengthscale_bounds": (1e-5, 100.0),
                              "min_elbo_pct_change": 5.0 },
        "feasibility_model": feasibility,
        "file_path": f"{results_path}/dmosopt_{celltype}.h5",
        "save": True,
        "save_surrogate_eval": save_surrogate_eval,
        "metadata": problem_metadata,
        "local_random": local_random,
    }

    if space_sensitivity is not None:
        dmosopt_params["di_crossover"] = space_sensitivity
        dmosopt_params["di_mutation"] = space_sensitivity

    best = dmosopt.run(dmosopt_params, verbose=True)

    if interactive and best is not None:

        import matplotlib.pyplot as plt

        bestx, besty = best

        bestx_dict = dict(bestx)

        param_dict = {k: bestx_dict[k][0] for k in space}
        for p in problem_parameters:
            param_dict[p] = problem_parameters[p]

        cell = template(param_dict)

        vec_t, vec_v = run_iclamp(cell, 0.15, 500.0, 1000.0)

        plt.plot(vec_t, vec_v, linewidth=3, color="r")
        plt.xlabel("Time (ms)")
        plt.ylabel("V (mV)")
        plt.savefig(f"{celltype}_nrn.svg")


if __name__ == "__main__":
    main(
        args=sys.argv[
            (list_find(lambda x: os.path.basename(x) == script_name, sys.argv) + 1) :
        ]
    )
