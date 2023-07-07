"""
Tests ability to run a GWF model then a PRT model
in separate simulations via flow model interface.

The grid is a 10x10 square with a single layer.
Particles are released from the top left cell.
"""


import os
from pathlib import Path

import flopy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
from flopy.utils.binaryfile import HeadFile

from prt_track_utils import check_track_data, get_track_dtype

# model names
name = "prtfmi0"
gwfname = f"{name}"
prtname = f"{name}_prt"

# output file names
gwf_budget_file = f"{name}.bud"
gwf_head_file = f"{name}.hds"
prp_track_file = f"{name}.prp.trk"
prp_track_csv_file = f"{name}.prp.trk.csv"
prt_budget_file = f"{name}.cbb"
prt_track_file = f"{name}.trk"
prt_track_csv_file = f"{name}.trk.csv"

# problem info
nlay=1
nrow=10
ncol=10
top=1.0
botm=[0.0]
nper=1
perlen=1.0
nstp=1
tsmult=1.0
porosity=0.1
releasepts = [
    # index, k, i, j, x, y, z
    # (0-based indexing converted to 1-based for mf6 by flopy)
    (i, 0, 0, 0, float(f"0.{i}"), float(f"9.{i}"), 0.5)
    for i in range(10)
]


def build_gwf_sim(ws, mf6):
    # create simulation
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        exe_name=mf6,
        version="mf6",
        sim_ws=ws,
    )

    # create tdis package
    pd = (perlen, nstp, tsmult)
    flopy.mf6.modflow.mftdis.ModflowTdis(
        sim,
        pname="tdis",
        time_units="DAYS",
        nper=nper,
        perioddata=[pd],
    )

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, save_flows=True)

    # create gwf discretization
    flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        gwf,
        pname="dis",
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
    )

    # create gwf initial conditions package
    flopy.mf6.modflow.mfgwfic.ModflowGwfic(gwf, pname="ic")

    # create gwf node property flow package
    flopy.mf6.modflow.mfgwfnpf.ModflowGwfnpf(
        gwf,
        pname="npf",
        save_saturation=True,
        save_specific_discharge=True,
    )

    # create gwf chd package
    spd = {
        0: [[(0, 0, 0), 1.0, 1.0], [(0, 9, 9), 0.0, 0.0]],
        1: [[(0, 0, 0), 0.0, 0.0], [(0, 9, 9), 1.0, 2.0]],
    }
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        pname="CHD-1",
        stress_period_data=spd,
        auxiliary=["concentration"],
    )

    # create gwf output control package
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=gwf_budget_file,
        head_filerecord=gwf_head_file,
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # create iterative model solution for gwf model
    ims = flopy.mf6.ModflowIms(sim)

    return sim


def build_prt_sim(ws, mf6):
    # create simulation
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        exe_name=mf6,
        version="mf6",
        sim_ws=ws,
    )

    # create tdis package
    pd = (perlen, nstp, tsmult)
    flopy.mf6.modflow.mftdis.ModflowTdis(
        sim,
        pname="tdis",
        time_units="DAYS",
        nper=nper,
        perioddata=[pd],
    )

    # create prt model
    prt = flopy.mf6.ModflowPrt(sim, modelname=prtname)

    # create prt discretization
    flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        prt,
        pname="dis",
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
    )

    # create mip package
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=porosity)

    # create prp package
    
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prtname}_1.prp",
        nreleasepts=len(releasepts),
        packagedata=releasepts,
        track_filerecord=[prp_track_file],
        trackcsv_filerecord=[prp_track_csv_file],
        perioddata={0: ["FIRST"]},
    )

    # create output control package
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        budget_filerecord=[prt_budget_file],
        track_filerecord=[prt_track_file],
        trackcsv_filerecord=[prt_track_csv_file],
        saverecord=[("BUDGET", "ALL")],
    )

    # create the flow model interface
    pd = [
        ("GWFHEAD", gwf_head_file),
        ("GWFBUDGET", gwf_budget_file),
    ]
    flopy.mf6.ModflowPrtfmi(prt, packagedata=pd)

    # add explicit model solution
    ems = flopy.mf6.ModflowEms(
        sim,
        pname="ems",
        filename=f"{prtname}.ems",
    )
    sim.register_solution_package(ems, [prt.name])

    return sim


def check_budget_data(lst: os.PathLike, cbb: os.PathLike):
    # load PRT model's list file
    mflist = flopy.utils.mflistfile.ListBudget(
        lst, budgetkey="MASS BUDGET FOR ENTIRE MODEL"
    )
    names = mflist.get_record_names()
    entries = mflist.entries

    # check timesteps
    inc = mflist.get_incremental()
    v = inc["totim"][-1]
    assert v == perlen * nper, f"Last time should be {perlen}.  Found {v}"

    # entries should be a subset of names
    assert all(e in names for e in entries)

    # todo what other record names should we expect?
    expected_entries = [
        "PRP_IN",
        "PRP_OUT",
    ]
    assert all(en in names for en in expected_entries)

    # load and check cell budget file
    # todo: reinstate below after mass budget is saved to budget file
    # mfbud = flopy.utils.binaryfile.CellBudgetFile(cbb)
    # assert mfbud.nlay == nlay
    # assert mfbud.nrow == nrow
    # assert mfbud.ncol == ncol
    # assert len(mfbud.times) == 1
    # assert mfbud.times[0] == perlen


def test_mf6model(function_tmpdir, targets):
    ws = function_tmpdir
    gwfsim = build_gwf_sim(ws, targets.mf6)
    prtsim = build_prt_sim(ws, targets.mf6)

    for sim in [gwfsim, prtsim]:
        sim.write_simulation()
        success, _ = sim.run_simulation()
        assert success

    # check output files exist
    assert (ws / gwf_budget_file).is_file()
    assert (ws / gwf_head_file).is_file()
    assert (ws / prt_budget_file).is_file()
    assert (ws / prt_track_file).is_file()
    assert (ws / prt_track_csv_file).is_file()
    assert (ws / prp_track_file).is_file()
    assert (ws / prp_track_csv_file).is_file()

    # check cell budget file
    check_budget_data(ws / f"{name}_prt.lst", ws / prt_budget_file)

    # check track data
    check_track_data(
        track_bin=ws / prt_track_file,
        track_hdr=ws / Path(prt_track_file.replace(".trk", ".trk.hdr")),
        track_csv=ws / prt_track_csv_file,
    )

    # extract head, budget, and specific discharge results from GWF model
    gwf = gwfsim.get_model(gwfname)
    prt = prtsim.get_model(prtname)
    hds = HeadFile(ws / gwf_head_file).get_data()
    bud = gwf.output.budget()
    spdis = bud.get_data(text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    # debugging: plot pathlines in map view
    mg = gwf.modelgrid
    pmv = flopy.plot.PlotMapView(modelgrid=mg)
    pmv.plot_grid()
    pmv.plot_array(hds[0], alpha=0.1)
    pmv.plot_vector(qx, qy, normalize=True, color="white")
    csvdata = pd.read_csv(ws / prt_track_csv_file)
    plines = csvdata.groupby(['iprp', 'irpt', 'trelease'])
    for ipl, (pl_name, pl) in enumerate(plines):
        plt.plot(pl['x'], pl['y'], '.-', lw=.01, color=cm.viridis(ipl / len(plines)))
    # plt.show()
