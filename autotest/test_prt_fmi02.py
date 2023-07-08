"""
Tests ability to run a GWF model then a PRT model
in separate simulations via flow model interface.

The grid is a 10x10 square with a single layer,
the same flow system shown on the FloPy readme,
except for 2 inactive cells in the bottom left
and top right corners.

Particles are released from the bottom left and
top right cells of the model grid.

One motivation for this test case is to check 
cell numbers reported in pathline data - they
should have been converted from reduced node
numbers to user node numbers 
"""


import os
from pathlib import Path

import flopy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd
from flopy.utils.binaryfile import HeadFile

from prt_track_utils import check_track_data

# model names
name = "prtfmi02"
gwfname = f"{name}"
prtname = f"{name}_prt"

# output file names
gwf_budget_file = f"{name}.bud"
gwf_head_file = f"{name}.hds"
prp_track_file_a = f"{name}_a.prp.trk"
prp_track_csv_file_a = f"{name}_a.prp.trk.csv"
prp_track_file_b = f"{name}_b.prp.trk"
prp_track_csv_file_b = f"{name}_b.prp.trk.csv"
prt_budget_file = f"{name}.cbb"
prt_track_file = f"{name}.trk"
prt_track_csv_file = f"{name}.trk.csv"

# problem info
nlay = 1
nrow = 10
ncol = 10
top = 1.0
botm = [0.0]
nper = 1
perlen = 1.0
nstp = 1
tsmult = 1.0
porosity = 0.1
releasepts_a = [
    # index, k, i, j, x, y, z
    # (0-based indexing converted to 1-based for mf6 by flopy)
    (i, 0, 0, 0, float(f"0.{i + 1}"), float(f"9.{i + 1}"), 0.5)
    for i in range(4)
]
releasepts_b = [
    # index, k, i, j, x, y, z
    # (0-based indexing converted to 1-based for mf6 by flopy)
    (i, 0, 0, 0, float(f"0.{i + 5}"), float(f"9.{i + 5}"), 0.5)
    for i in range(5)
]
idomain = np.ones((nlay, nrow, ncol), dtype=int)
idomain[0, 0, 9] = 0
idomain[0, 9, 0] = 0
# idomain = idomain.ravel()


# expected particle track solution locations
exp_locs = np.array(
    [

    ],
    dtype=[
        ("x", "<f8"),
        ("y", "<f8"),
        ("z", "<f8"),
        ("t", "<f8"),
        ("ilay", "<i8"),
        ("icell", "<i8"),
    ],
).view(np.recarray)
exp_locs.sort()


def build_gwf_sim(ws, mf6):
    # create simulation
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        exe_name=mf6,
        version="mf6",
        sim_ws=ws,
    )

    # create tdis package
    flopy.mf6.modflow.mftdis.ModflowTdis(
        sim,
        pname="tdis",
        time_units="DAYS",
        nper=nper,
        perioddata=[(perlen, nstp, tsmult)],
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
        idomain=idomain,
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
    flopy.mf6.modflow.mftdis.ModflowTdis(
        sim,
        pname="tdis",
        time_units="DAYS",
        nper=nper,
        perioddata=[(perlen, nstp, tsmult)],
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
        idomain=idomain,
    )

    # create mip package
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=porosity)

    # create prp packages
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp",
        filename=f"{prtname}_a.prp",
        nreleasepts=len(releasepts_a),
        packagedata=releasepts_a,
        track_filerecord=[prp_track_file_a],
        trackcsv_filerecord=[prp_track_csv_file_a],
        perioddata={0: ["FIRST"]},
    )
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp_b",
        filename=f"{prtname}_b.prp",
        nreleasepts=len(releasepts_b),
        packagedata=releasepts_b,
        track_filerecord=[prp_track_file_b],
        trackcsv_filerecord=[prp_track_csv_file_b],
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


def test_prt_fmi02(function_tmpdir, targets):
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
    assert (ws / prp_track_file_a).is_file()
    assert (ws / prp_track_csv_file_a).is_file()
    assert (ws / prp_track_file_b).is_file()
    assert (ws / prp_track_csv_file_b).is_file()

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

    # plot pathlines in map view
    mg = gwf.modelgrid
    pmv = flopy.plot.PlotMapView(modelgrid=mg)
    pmv.plot_grid()
    pmv.plot_array(hds[0], alpha=0.1)
    pmv.plot_vector(qx, qy, normalize=True, color="white")
    csvdata = pd.read_csv(ws / prt_track_csv_file)
    plines = csvdata.groupby(["iprp", "irpt", "trelease"])
    ax = plt.gca()

    # plot color-coded pathlines
    for ipl, (pl_name, pl) in enumerate(plines):
        data = csvdata[
            (csvdata["iprp"] == pl_name[0])
            & (csvdata["irpt"] == pl_name[1])
            & (csvdata["trelease"] == pl_name[2])
        ]
        data.plot(
            kind="line",
            x="x",
            y="y",
            ax=ax,
            legend=False,
            color=cm.plasma(ipl / len(plines)),
        )

    # cols = csvdata.columns
    cols = ["x", "y", "z", "t", "ilay", "icell"]
    locs = (
        pd.DataFrame(csvdata, columns=cols)
        .drop_duplicates()
        .to_records(index=False)
    )
    locs.sort()
    # assert np.allclose(locs, exp_locs)  # why does this fail with promotion error?
    # for col in cols:
    #     if col in exp_locs.dtype.names:
    #         assert np.allclose(locs[col], exp_locs[col])

    #     if col == "ilay":
    #         assert np.all(locs[col] == 1)

    # save the plot to file
    plt.savefig(ws / f"test_{name}.png")

    # debugging
    # plt.show()
    # from pprint import pprint
    # pprint(locs)

    for x, y, z, t, ilay, icell in list(locs):
        # debugging
        # print(x, y, z, t, ilay, icell)

        k, i, j = mg.intersect(x, y, z)
        nn = mg.get_node([k, i, j]) + 1
        neighbors = mg.neighbors(nn)
        # todo: figure out why cell numbers are off when idomain is used
        # assert np.isclose(nn + 1, icell, atol=1) or any(nn == n for n in neighbors)
        assert ilay == (k + 1) == 1
