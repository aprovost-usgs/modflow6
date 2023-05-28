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


# dict insertion order is preserved since py3.6 so
# cases can be unzipped with .keys() and .values() 
cases = {
    "none": None,
    "first": {0: ["FIRST"]},
    "first_first": {0: ["FIRST"], 1: ["FIRST"]},
    # "fraction": {0: ["FRACTION", 0.5], 1: ["FRACTION", 0.5]},
    # "first_frac": {0: [("FIRST"), ("FRACTION", 0.5)]},  # todo debug flopy, hangs and doesn't write prp file
    "all": {0: ["ALL"]},
    "all_all": {0: ["ALL"], 1: ["ALL"]},
    # "all_frac": {0: ["ALL", ("FRACTION", 0.5)]},  # todo debug flopy, hangs and doesn't write prp file
    "freq": {0: ["FREQUENCY", 2], 1: ["FREQUENCY", 5]},
    # "freq_frac": {0: [("FREQUENCY", 2), ("FRACTION", 0.5)]},  # todo debug flopy, KeyError 'fractionrecord'
    "steps": {0: ["STEPS", 1, 2, 3], 1: ["STEPS", 4, 5, 6]},
    # "steps_frac": {0: [("STEPS", 1, 2, 3), ("FRACTION", 0.5)]},  # todo debug flopy, KeyError 'fractionrecord'
}
# define total particle mass (# of particles currently
# since mass is hardcoded 1 for now) expected for each
# release timing test case
cmass = {
    "none": 1.,
    "first": 2.,
    "first_first": 2.,
    "all": 20.,
    "all_all": 20.,
    "freq": 7.,
    "steps": 6.
}
@pytest.mark.parametrize("release_timing", cases.values(), ids=cases.keys())
def test_releasesetting(request, function_tmpdir, targets, release_timing):
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
    assert (ws / f"{name}.bud").is_file()
    assert (ws / f"{name}.hds").is_file()
    assert (ws / f"{name}.cbb").is_file()

    # get expected particle mass for this case
    case_name = request.node.name.rpartition("[")[2].rpartition("]")[0]
    expected_mass = cmass.get(case_name, None)
    if expected_mass is None:
        print(f"No expected mass for case {case_name}, skipping check")
        return

    # parse particle mass from list file
    lines = (ws / f"{name}_prt.lst").open().readlines()
    line = next(reversed([l for l in lines if "TOTAL OUT" in l]), None)
    assert line
    mass = float(line.split("TOTAL OUT")[1].replace("=", "").strip())
    assert mass == expected_mass

    # check release timing by logs in list file
    # lines = (ws / f"{name}_1.prp").open().readlines()
    # assert any(l for l in lines if "PARTICLE RELEASE SCHEDULED AT TIME STEP FRACTION" in l)

    # check release timing by pathline output
    # lines = (ws / f"{name}.pathline").open().readlines()
    # assert any(l for l in lines if "0.500000" in l)


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