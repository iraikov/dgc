import os, sys, math, logging, yaml, click
from functools import partial
from neuron import h
from neuron_utils import ic_constant_f, run_iclamp, load_template
import numpy as np
import matplotlib.pyplot as plt
import ephys_utils as ephys
from protocol import ExperimentalProtocol
from dmosopt_PR_nrn import make_obj_fun


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)



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
    "--params-config",
    "-c",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
)
@click.option("--params-id", "-i", required=True, type=int)
@click.option("--model-variant", "-m", default="default", type=str)
@click.option("--verbose", "-v", is_flag=True)
def main(
    config_path,
    target_namespace,
    params_config,
    params_id,
    model_variant,
    verbose,
):

    comm = MPI.COMM_WORLD
    protocol_config_dict = None
    if comm.rank == 0:
        with open(config_path) as f:
            protocol_config_dict = yaml.load(f, Loader=yaml.FullLoader)

    protocol_config_dict = comm.bcast(protocol_config_dict)

    celltype = protocol_config_dict["Celltype"]
    template_dict = protocol_config_dict.get("Template", None)
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

    
    param_dict = None
    with open(params_config) as f:
        param_dict = yaml.load(f, Loader=yaml.FullLoader)

    params = param_dict[param_id]
    _, feature_vals, constr_vals = obj_fun(params)


    # Run iclamp experiments
    iclamp_results = []
    cell = init_cell(
        template_name, pp, v_hold=exp_protocol.v_hold, ic_constant_val=ic_constant_hold
    )
    iclamp_results = exp_protocol.run_iclamp_steps(cell)
    
    


if __name__ == "__main__":
    main(
        args=sys.argv[
            (list_find(lambda x: os.path.basename(x) == script_name, sys.argv) + 1) :
        ]
    )
