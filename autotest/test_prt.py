from pathlib import Path
import flopy
import pytest
from pathlib import Path
from pytest_cases import parametrize_with_cases
from simulation import TestSimulation
from test_prt_cases import PrtCases


def get_basic_flow_sim(name, ws, exe):
    """
    Get a minimal flow model simulation with two stress periods, each
    with 10 time steps and no multiplier, and constant head boundary
    which changes between periods, making for a transient simulation.
    """

    ws = Path(ws)
    ws.mkdir(exist_ok=True, parents=True)

    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=exe)
    pd = [(1.0, 10, 1.0), (1.0, 10, 1.0)]
    tdis = flopy.mf6.ModflowTdis(sim, nper=len(pd), perioddata=pd)
    ims = flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    dis = flopy.mf6.ModflowGwfdis(gwf, nrow=10, ncol=10)
    ic = flopy.mf6.ModflowGwfic(gwf)
    npf = flopy.mf6.ModflowGwfnpf(
        gwf, save_specific_discharge=True, save_saturation=True
    )
    spd = {
        0: [[(0, 0, 0), 1.0, 1.0], [(0, 9, 9), 0.0, 0.0]],
        1: [[(0, 0, 0), 0.0, 0.0], [(0, 9, 9), 1.0, 2.0]],
    }
    chd = flopy.mf6.ModflowGwfchd(
        gwf, pname="CHD-1", stress_period_data=spd, auxiliary=["concentration"]
    )
    budget_file = f"{name}.bud"
    head_file = f"{name}.hds"
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=budget_file,
        head_filerecord=head_file,
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    return sim
    

def get_basic_tracking_sim(name, ws, exe, release_timing: dict = None):
    """
    Get a minimal particle tracking simulation with a single particle
    and two stress periods each with 10 time steps and no multiplier.
    If release_timing is not provided the particle is released at the
    beginning of the first stress period. The particle is released at
    the time specified by the release_timing dict otherwise, which is
    of the form expected by perioddata in the ModflowPrtprp package.
    """

    ws = Path(ws)
    ws.mkdir(exist_ok=True, parents=True)

    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=exe)
    pd = [(1.0, 10, 1.0), (1.0, 10, 1.0)]
    tdis = flopy.mf6.ModflowTdis(sim, nper=len(pd), perioddata=pd)
    prt = flopy.mf6.ModflowPrt(sim, modelname=name)
    dis = flopy.mf6.ModflowGwfdis(prt, nrow=10, ncol=10)
    mip = flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=.2)
    releasepts = [
        # particle id, k, i, j, localx, localy, localz
        (0, 0, 0, 0, 0.5, 0.5, 0.5)
    ]
    perioddata = {0: ["FIRST"]} if release_timing is None else release_timing
    prp = flopy.mf6.ModflowPrtprp(
        prt, pname="prp1", filename=f"{name}_1.prp",
        nreleasepts=len(releasepts), packagedata=releasepts,
        perioddata=perioddata,
    )
    prt_budget_file = f"{name}.cbb"
    oc = flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        budget_filerecord=[prt_budget_file],
        saverecord=[("BUDGET", "ALL")],
    )
    gwf_budget_file = f"{name}.bud"
    gwf_head_file = f"{name}.hds"
    pd = [
        ("GWFHEAD", gwf_head_file),
        ("GWFBUDGET", gwf_budget_file),
    ]
    fmi = flopy.mf6.ModflowPrtfmi(prt, packagedata=pd)
    ems = flopy.mf6.ModflowEms(
        sim, pname="ems",
        filename=f"{name}.ems",
    )
    sim.register_solution_package(ems, [prt.name])

    return sim


def get_basic_combined_sim(name, ws, exe, release_timing: dict = None):
    ws = Path(ws)
    ws.mkdir(exist_ok=True, parents=True)

    gwf_name = f"{name}_gwf"
    prt_name = f"{name}_prt"

    # create simulation and GWF model
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=exe)
    pd = [(1.0, 10, 1.0), (1.0, 10, 1.0)]
    tdis = flopy.mf6.ModflowTdis(sim, nper=len(pd), perioddata=pd)
    ims = flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwf_name, save_flows=True)
    dis = flopy.mf6.ModflowGwfdis(gwf, nrow=10, ncol=10)
    ic = flopy.mf6.ModflowGwfic(gwf)
    npf = flopy.mf6.ModflowGwfnpf(
        gwf, save_specific_discharge=True, save_saturation=True
    )
    spd = {
        0: [[(0, 0, 0), 1.0, 1.0], [(0, 9, 9), 0.0, 0.0]],
        1: [[(0, 0, 0), 0.0, 0.0], [(0, 9, 9), 1.0, 2.0]],
    }
    chd = flopy.mf6.ModflowGwfchd(
        gwf, pname="CHD-1", stress_period_data=spd, auxiliary=["concentration"]
    )
    budget_file = f"{name}.bud"
    head_file = f"{name}.hds"
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=budget_file,
        head_filerecord=head_file,
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # create PRT model
    prt = flopy.mf6.ModflowPrt(sim, modelname=prt_name)
    dis = flopy.mf6.ModflowGwfdis(prt, nrow=10, ncol=10)
    mip = flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=.2)
    releasepts = [
        # particle id, k, i, j, localx, localy, localz
        (0, 0, 0, 0, 0.5, 0.5, 0.5)
    ]
    prp = flopy.mf6.ModflowPrtprp(
        prt, pname="prp1", filename=f"{prt_name}_1.prp",
        nreleasepts=len(releasepts), packagedata=releasepts,
        perioddata=release_timing,
    )
    prt_budget_file = f"{name}.cbb"
    oc = flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        budget_filerecord=[prt_budget_file],
        saverecord=[("BUDGET", "ALL")],
    )
    gwf_budget_file = f"{name}.bud"
    gwf_head_file = f"{name}.hds"
    pd = [
        ("GWFHEAD", gwf_head_file),
        ("GWFBUDGET", gwf_budget_file),
    ]
    fmi = flopy.mf6.ModflowPrtfmi(prt, packagedata=pd)
    exg = flopy.mf6.ModflowGwfprt(
        sim, exgtype="GWF6-PRT6",
        exgmnamea=gwf_name, exgmnameb=prt_name,
        filename=f"{name}.gwfprt",
    )
    ems = flopy.mf6.ModflowEms(
        sim, pname="ems",
        filename=f"{prt_name}.ems",
    )
    sim.register_solution_package(ems, [prt.name])

    return sim


def test_fmi_basic(function_tmpdir, targets):
    """
    Tests ability to run GWF model first, then PRT model
    in separate simulations with a flow model interface.
    """

    name = "prtfmi0"
    mf6 = targets.mf6
    ws = function_tmpdir

    sim = get_basic_flow_sim(name, ws, mf6)
    sim.write_simulation()
    sim.run_simulation()
    assert (ws / f"{name}.bud").is_file()
    assert (ws / f"{name}.hds").is_file()

    sim = get_basic_tracking_sim(name, ws, mf6)
    sim.write_simulation()
    sim.run_simulation()
    assert (ws / f"{name}.cbb").is_file()


invalid_cases = {
    "frac_array": {0: ["FRACTION", 0.5, 0.2]},
    "frac_frac_array": {0: ["FRACTION", 0.5, 0.2], 1: ["FRACTION", 0.5, 0.1]},
}
@pytest.mark.parametrize("release_timing", invalid_cases.values(), ids=invalid_cases.keys())
def test_invalid_releasesetting(function_tmpdir, targets, release_timing):
    """
    Tests invalid fractional timestep release settings. All cases should fail.
    """

    name = "prtrel"
    mf6 = targets.mf6
    ws = function_tmpdir

    sim = get_basic_combined_sim(name, ws, mf6, release_timing)
    sim.write_simulation()
    success, buff = sim.run_simulation(report=True)
    assert not success
    assert any("FRACTION must be a single value if no other release setting is provided" in l for l in buff)


# dict insertion order is preserved since py3.6 so
# cases can be unzipped with .keys() and .values() 
valid_cases = {
    "none": None,
    "first": {0: ["FIRST"]},
    "first_first": {0: ["FIRST"], 1: ["FIRST"]},
    # "first_frac": {0: [("FIRST"), ("FRACTION", 0.5)]},  # todo debug flopy hanging
    "all": {0: ["ALL"]},
    "all_all": {0: ["ALL"], 1: ["ALL"]},
    # "all_frac": {0: [("ALL"), ("FRACTION", 0.5)]},  # todo debug flopy hanging
    "freq": {0: ["FREQUENCY", 2]},
    "freq_freq": {0: ["FREQUENCY", 2], 1: ["FREQUENCY", 5]},
    "freq_frac": {0: [("FREQUENCY", 2), ("FRACTION", 0.5)]},
    "steps": {0: ["STEPS", 1, 2, 3]},
    "steps_steps": {0: ["STEPS", 1, 2, 3], 1: ["STEPS", 4, 5, 6, 7]},
    "steps_frac": {0: [("STEPS", 1, 2, 3), ("FRACTION", 0.5, 0.4, 0.2)]},
    "frac": {0: ["FRACTION", 0.5]},
    "frac_frac": {0: ["FRACTION", 0.5], 1: ["FRACTION", 0.5]},
}
# define total particle mass (# of particles currently
# since mass is hardcoded 1 for now) expected for each
# release timing test case
valid_cmass = {
    "none": 1.,  # if no period data provided,  default is FIRST in first stress period
    "first": 2.,  # if only first stress period is provided, settings apply to all sp's
    "first_first": 2.,  # should be totally equivalent to 'first'
    # "first_frac": 2.,  # mass should be equal to 'first' even though timing is different
    "all": 20.,
    "all_all": 20.,  # should be equivalent to 'all' since 1st sp settings apply to all sp's
    # "all_frac": 20.,  # mass should be equal to 'all' and 'all_all'
    "freq": 10.,
    "freq_freq": 7.,  # 1st sp uses freq 2, 2nd sp uses freq 5, so (10/2) + (10/5) = 7
    "freq_frac": 10.,
    "steps": 6.,
    "steps_steps": 7.,  # different steps specified in 1st and 2nd sp's (1 more in 2nd)
    "steps_frac": 6.,  # 3 time steps selected in both sp's, mass equal to 'steps'
    "frac": 2.,  # mass should be equal to 'first' and 'first_first'
    "frac_frac": 2.,  # mass should be equal to 'first' and 'first_first'
}
@pytest.mark.parametrize("release_timing", valid_cases.values(), ids=valid_cases.keys())
def test_valid_releasesetting(request, function_tmpdir, targets, release_timing):
    """
    Tests fractional timestep release configuration, where particles 
    are released some fraction of the way through selected timesteps
    """

    name = "prtrel"
    mf6 = targets.mf6
    ws = function_tmpdir

    sim = get_basic_combined_sim(name, ws, mf6, release_timing)
    sim.write_simulation()
    success, buff = sim.run_simulation(report=True)
    assert success
    assert (ws / f"{name}_prt.lst").is_file()
    assert (ws / f"{name}.bud").is_file()  # GWF budget binary output
    assert (ws / f"{name}.hds").is_file()  # GWF head binary output
    assert (ws / f"{name}.cbb").is_file()  # PRT budget binary output
    # assert (ws / f"{name}.prk").is_file()  # particle track binary output
    # assert (ws / f"{name}.prh").is_file()  # particle track ascii headers/dtypes for .prk binary output

    # get expected particle mass for this case
    case_name = request.node.name.rpartition("[")[2].rpartition("]")[0]
    expected_mass = valid_cmass.get(case_name, None)
    if expected_mass is None:
        print(f"No expected mass for case {case_name}, skipping check")
        return

    # parse particle mass from list file
    lines = (ws / f"{name}_prt.lst").open().readlines()
    line = next(reversed([l for l in lines if "TOTAL OUT" in l]), None)
    assert line
    mass = float(line.split("TOTAL OUT")[1].replace("=", "").strip())
    assert mass == expected_mass


@parametrize_with_cases("case", cases=PrtCases)
def test_prt_models(case, targets):
    # cases are tuples (context, simulation, optional comparison simulation, evaluation function)
    ctx, sim, cmp, evl = case
    sim.write_simulation()
    if cmp:
        cmp.write_simulation()
    
    test = TestSimulation(
        name=ctx.name,
        exe_dict=targets,
        exfunc=lambda s: evl(ctx, s),  # hack the context into the evaluation function for now
        idxsim=0,
        mf6_regression=True,
        require_failure=ctx.xfail,
        make_comparison=False,
    )

    test.set_model(sim.simulation_data.mfpath.get_sim_path(), testModel=False)
    test.run()
    test.compare()