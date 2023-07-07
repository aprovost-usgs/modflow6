"""Test gwf+prt models in same simulation via exchange"""


import os
from pathlib import Path

import flopy
from flopy.utils.binaryfile import write_budget, write_head
from flopy.utils.gridutil import uniform_flow_field
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import numpy as np
import pytest
from flopy.utils.binaryfile import HeadFile
from framework import TestFramework
from simulation import TestSimulation

from prt_track_utils import check_track_data, get_track_dtype

# model names
name="prtexg1"
gwfname = f"{name}"
prtname = f"{name}_prt"

# output file names
gwf_budget_file = f"{gwfname}.bud"
gwf_head_file = f"{gwfname}.hds"
prp_track_file = f"{prtname}.prp.trk"
prp_track_csv_file = f"{prtname}.prp.trk.csv"
prt_budget_file = f"{prtname}.cbb"
prt_track_file = f"{prtname}.trk"
prt_track_csv_file = f"{prtname}.trk.csv"

# model info
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

# expected particle track solution locations
exp_locs = np.array(
    [
        (0.00000000e00, 9.0, 0.5),
        (9.97537086e-21, 8.0, 0.5),
        (1.45058368e-20, 7.0, 0.5),
        (1.95476429e-20, 6.0, 0.5),
        (2.54566727e-20, 5.0, 0.5),
        (3.31945641e-20, 4.0, 0.5),
        (4.50762642e-20, 3.0, 0.5),
        (6.79747993e-20, 2.0, 0.5),
        (1.35948624e-19, 1.0, 0.5),
        (1.00000000e00, 0.0, 0.5),
        (2.00000000e00, 0.0, 0.5),
        (3.00000000e00, 0.0, 0.5),
        (4.00000000e00, 0.0, 0.5),
        (5.00000000e00, 0.0, 0.5),
        (6.00000000e00, 0.0, 0.5),
        (7.00000000e00, 0.0, 0.5),
        (8.00000000e00, 0.0, 0.5),
        (9.00000000e00, 0.0, 0.5),
        (1.00000000e-01, 9.1, 0.5),
        (1.11111111e-01, 9.0, 0.5),
        (1.84020109e-01, 8.0, 0.5),
        (2.67595632e-01, 7.0, 0.5),
        (3.60604074e-01, 6.0, 0.5),
        (4.69610579e-01, 5.0, 0.5),
        (6.12354908e-01, 4.0, 0.5),
        (8.31541922e-01, 3.0, 0.5),
        (1.00000000e00, 2.49994832, 0.5),
        (1.25597089e00, 2.0, 0.5),
        (2.00000000e00, 1.2558761, 0.5),
        (2.49970272e00, 1.0, 0.5),
        (3.00000000e00, 0.83148486, 0.5),
        (4.00000000e00, 0.61228533, 0.5),
        (5.00000000e00, 0.46947412, 0.5),
        (6.00000000e00, 0.36041286, 0.5),
        (7.00000000e00, 0.26740321, 0.5),
        (8.00000000e00, 0.18388997, 0.5),
        (9.00000000e00, 0.1110318, 0.5),
        (2.00000000e-01, 9.2, 0.5),
        (2.50000000e-01, 9.0, 0.5),
        (4.14045245e-01, 8.0, 0.5),
        (6.02090173e-01, 7.0, 0.5),
        (8.11359167e-01, 6.0, 0.5),
        (1.00000000e00, 5.18731402, 0.5),
        (1.06083307e00, 5.0, 0.5),
        (1.39371125e00, 4.0, 0.5),
        (1.88967622e00, 3.0, 0.5),
        (2.00000000e00, 2.8357984, 0.5),
        (2.83548743e00, 2.0, 0.5),
        (3.00000000e00, 1.88940228, 0.5),
        (4.00000000e00, 1.39328875, 0.5),
        (5.00000000e00, 1.06028978, 0.5),
        (5.18544956e00, 1.0, 0.5),
        (6.00000000e00, 0.81077304, 0.5),
        (7.00000000e00, 0.60154156, 0.5),
        (8.00000000e00, 0.41367289, 0.5),
        (9.00000000e00, 0.24977352, 0.5),
        (3.00000000e-01, 9.3, 0.5),
        (4.28571429e-01, 9.0, 0.5),
        (7.09791849e-01, 8.0, 0.5),
        (1.00000000e00, 7.07079686, 0.5),
        (1.04000660e00, 7.0, 0.5),
        (1.44194733e00, 6.0, 0.5),
        (1.87178010e00, 5.0, 0.5),
        (2.00000000e00, 4.7211973, 0.5),
        (2.44542502e00, 4.0, 0.5),
        (3.00000000e00, 3.28698638, 0.5),
        (3.28692505e00, 3.0, 0.5),
        (4.00000000e00, 2.44511796, 0.5),
        (4.72039333e00, 2.0, 0.5),
        (5.00000000e00, 1.87135558, 0.5),
        (6.00000000e00, 1.44147854, 0.5),
        (7.00000000e00, 1.03952535, 0.5),
        (7.06996554e00, 1.0, 0.5),
        (8.00000000e00, 0.70953937, 0.5),
        (9.00000000e00, 0.42841615, 0.5),
        (4.00000000e-01, 9.4, 0.5),
        (6.66666667e-01, 9.0, 0.5),
        (1.00000000e00, 8.15867673, 0.5),
        (1.15868934e00, 8.0, 0.5),
        (1.75342419e00, 7.0, 0.5),
        (2.00000000e00, 6.5110226, 0.5),
        (2.37704205e00, 6.0, 0.5),
        (3.00000000e00, 5.07937977, 0.5),
        (3.06875198e00, 5.0, 0.5),
        (3.94003343e00, 4.0, 0.5),
        (4.00000000e00, 3.94002526, 0.5),
        (5.00000000e00, 3.06865624, 0.5),
        (5.07927134e00, 3.0, 0.5),
        (6.00000000e00, 2.37700418, 0.5),
        (6.51110124e00, 2.0, 0.5),
        (7.00000000e00, 1.75351945, 0.5),
        (8.00000000e00, 1.1588233, 0.5),
        (8.15880985e00, 1.0, 0.5),
        (9.00000000e00, 0.66671565, 0.5),
        (5.00000000e-01, 9.5, 0.5),
        (1.00000000e00, 9.0, 0.5),
        (2.00000000e00, 8.00007947, 0.5),
        (2.00012544e00, 8.0, 0.5),
        (3.00000000e00, 7.00015821, 0.5),
        (3.00020767e00, 7.0, 0.5),
        (4.00000000e00, 6.00013079, 0.5),
        (4.00014941e00, 6.0, 0.5),
        (5.00000000e00, 5.0001141, 0.5),
        (5.00011410e00, 5.0, 0.5),
        (6.00000000e00, 4.00012377, 0.5),
        (6.00010833e00, 4.0, 0.5),
        (7.00000000e00, 3.00010831, 0.5),
        (7.00008250e00, 3.0, 0.5),
        (8.00000000e00, 2.00008248, 0.5),
        (8.00005226e00, 2.0, 0.5),
        (9.00000000e00, 1.00005227, 0.5),
        (9.00002071e00, 1.0, 0.5),
        (6.00000000e-01, 9.6, 0.5),
        (1.00000000e00, 9.33333333, 0.5),
        (1.84132327e00, 9.0, 0.5),
        (2.00000000e00, 8.84133588, 0.5),
        (3.00000000e00, 8.24660575, 0.5),
        (3.48903677e00, 8.0, 0.5),
        (4.00000000e00, 7.62300176, 0.5),
        (4.92068497e00, 7.0, 0.5),
        (5.00000000e00, 6.9313041, 0.5),
        (6.00000000e00, 6.06003039, 0.5),
        (6.06003857e00, 6.0, 0.5),
        (6.93139984e00, 5.0, 0.5),
        (7.00000000e00, 4.92079341, 0.5),
        (7.62303964e00, 4.0, 0.5),
        (8.00000000e00, 3.48895816, 0.5),
        (8.24651050e00, 3.0, 0.5),
        (8.84120193e00, 2.0, 0.5),
        (9.00000000e00, 1.84121538, 0.5),
        (9.33329435e00, 1.0, 0.5),
        (7.00000000e-01, 9.7, 0.5),
        (1.00000000e00, 9.57142857, 0.5),
        (2.00000000e00, 9.29020815, 0.5),
        (2.92920314e00, 9.0, 0.5),
        (3.00000000e00, 8.9599934, 0.5),
        (4.00000000e00, 8.55805267, 0.5),
        (5.00000000e00, 8.1282199, 0.5),
        (5.27880270e00, 8.0, 0.5),
        (6.00000000e00, 7.55457498, 0.5),
        (6.71301362e00, 7.0, 0.5),
        (7.00000000e00, 6.71307495, 0.5),
        (7.55488204e00, 6.0, 0.5),
        (8.00000000e00, 5.27960667, 0.5),
        (8.12864442e00, 5.0, 0.5),
        (8.55852146e00, 4.0, 0.5),
        (8.96047465e00, 3.0, 0.5),
        (9.00000000e00, 2.93003446, 0.5),
        (9.29046063e00, 2.0, 0.5),
        (9.57158385e00, 1.0, 0.5),
        (8.00000000e-01, 9.8, 0.5),
        (1.00000000e00, 9.75, 0.5),
        (2.00000000e00, 9.58595475, 0.5),
        (3.00000000e00, 9.39790983, 0.5),
        (4.00000000e00, 9.18864083, 0.5),
        (4.81268598e00, 9.0, 0.5),
        (5.00000000e00, 8.93916693, 0.5),
        (6.00000000e00, 8.60628875, 0.5),
        (7.00000000e00, 8.11032378, 0.5),
        (7.16420160e00, 8.0, 0.5),
        (8.00000000e00, 7.16451257, 0.5),
        (8.11059772e00, 7.0, 0.5),
        (8.60671125e00, 6.0, 0.5),
        (8.93971022e00, 5.0, 0.5),
        (9.00000000e00, 4.81455044, 0.5),
        (9.18922696e00, 4.0, 0.5),
        (9.39845844e00, 3.0, 0.5),
        (9.58632711e00, 2.0, 0.5),
        (9.75022648e00, 1.0, 0.5),
        (9.00000000e-01, 9.9, 0.5),
        (1.00000000e00, 9.88888889, 0.5),
        (2.00000000e00, 9.81597989, 0.5),
        (3.00000000e00, 9.73240437, 0.5),
        (4.00000000e00, 9.63939593, 0.5),
        (5.00000000e00, 9.53038942, 0.5),
        (6.00000000e00, 9.38764509, 0.5),
        (7.00000000e00, 9.16845808, 0.5),
        (7.50005168e00, 9.0, 0.5),
        (8.00000000e00, 8.74402911, 0.5),
        (8.74412390e00, 8.0, 0.5),
        (9.00000000e00, 7.50029728, 0.5),
        (9.16851514e00, 7.0, 0.5),
        (9.38771467e00, 6.0, 0.5),
        (9.53052588e00, 5.0, 0.5),
        (9.63958714e00, 4.0, 0.5),
        (9.73259679e00, 3.0, 0.5),
        (9.81611003e00, 2.0, 0.5),
        (9.88896820e00, 1.0, 0.5),
    ],
    dtype=[("x", "<f8"), ("y", "<f8"), ("z", "<f8")],
).view(np.recarray)
exp_locs.sort(order=["x", "y", "z"])

# test cases
ex = [name]


def build_sim(ws, mf6):
    workspace = ws
    
    # create simulation
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        exe_name=mf6,
        version="mf6",
        sim_ws=workspace,
    )

    # create tdis package
    perioddata = (perlen, nstp, tsmult)
    flopy.mf6.modflow.mftdis.ModflowTdis(
        sim,
        pname="tdis",
        time_units="DAYS",
        nper=nper,
        perioddata=[perioddata],
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
    flopy.mf6.ModflowPrtfmi(prt, packagedata=[
        ("GWFHEAD", gwf_head_file),
        ("GWFBUDGET", gwf_budget_file),
    ])

    # create exchange
    flopy.mf6.ModflowGwfprt(
        sim,
        exgtype="GWF6-PRT6",
        exgmnamea=gwfname,
        exgmnameb=prtname,
        filename=f"{gwfname}.gwfprt",
    )

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


def eval_results(sim):
    print(f"Evaluating results for sim {sim.name}")
    simpath = Path(sim.simpath)

    # check budget data
    check_budget_data(
        simpath / f"{sim.name}_prt.lst",
        simpath / f"{sim.name}_prt.cbb",
    )

    # check particle track data
    prt_track_file = simpath / f"{sim.name}_prt.trk"
    prt_track_hdr_file = simpath / f"{sim.name}_prt.trk.hdr"
    prt_track_csv_file = simpath / f"{sim.name}_prt.trk.csv"
    assert prt_track_file.exists()
    assert prt_track_hdr_file.exists()
    assert prt_track_csv_file.exists()
    check_track_data(
        track_bin=prt_track_file,
        track_hdr=prt_track_hdr_file,
        track_csv=prt_track_csv_file,
    )


@pytest.mark.parametrize("name", ex)
def test_mf6model(name, function_tmpdir, targets):
    ws = function_tmpdir
    sim = build_sim(str(ws), targets.mf6)
    sim.write_simulation()

    test = TestFramework()
    test.run(
        TestSimulation(
            name=name, exe_dict=targets, exfunc=eval_results, 
            idxsim=0, make_comparison=False,
        ),
        str(ws),
    )

    # extract head, budget, and specific discharge results from GWF model
    gwf = sim.get_model(name)
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
    ax = plt.gca()

    # plot color-coded pathlines
    for ipl, (pl_name, pl) in enumerate(plines):
        data = csvdata[(csvdata["iprp"] == pl_name[0]) & (csvdata["irpt"] == pl_name[1]) & (csvdata["trelease"] == pl_name[2])]
        data.plot(
            kind="line",
            x='x',
            y='y',
            ax=ax,
            legend=False,
            color=cm.plasma(ipl / len(plines)))

    # plt.show()

    # cols = csvdata.columns
    cols = ["x", "y", "z", "t", "ilay"]
    locs = (
        pd.DataFrame(csvdata, columns=cols)
        .drop_duplicates()
        .to_records(index=False)
    )
    locs.sort()
    # assert np.allclose(locs, exp_locs)  # why does this fail with promotion error?
    for col in cols:
        if col in exp_locs.dtype.names:
            assert np.allclose(locs[col], exp_locs[col])
        
        if col == "ilay":
            assert np.all(locs[col] == 1)
