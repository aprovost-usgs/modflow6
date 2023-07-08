"""
Tests ability to run a GWF model then a PRT model
in separate simulations via flow model interface.

The grid is a 10x10 square with a single layer,
the same flow system shown on the FloPy readme.

Particles are released from the top left cell.
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
name = "prtfmi01"
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
releasepts = [
    # index, k, i, j, x, y, z
    # (0-based indexing converted to 1-based for mf6 by flopy)
    (i, 0, 0, 0, float(f"0.{i + 1}"), float(f"9.{i + 1}"), 0.5)
    for i in range(9)
]

# expected particle track solution locations
exp_locs = np.array(
    [
        (0.1, 9.1, 0.5, 0.0, 1, 1),
        (0.11111111, 9.0, 0.5, 0.06345987, 1, 1),
        (0.18402011, 8.0, 0.5, 0.83043083, 1, 11),
        (0.2, 9.2, 0.5, 0.0, 1, 1),
        (0.25, 9.0, 0.5, 0.13440196, 1, 1),
        (0.26759563, 7.0, 0.5, 2.02638995, 1, 21),
        (0.3, 9.3, 0.5, 0.0, 1, 1),
        (0.36060407, 6.0, 0.5, 3.70426544, 1, 31),
        (0.4, 9.4, 0.5, 0.0, 1, 1),
        (0.41404525, 8.0, 0.5, 0.90137293, 1, 11),
        (0.42857143, 9.0, 0.5, 0.21482947, 1, 1),
        (0.46961058, 5.0, 0.5, 5.92893723, 1, 41),
        (0.5, 9.5, 0.5, 0.0, 1, 1),
        (0.6, 9.6, 0.5, 0.0, 1, 1),
        (0.60209017, 7.0, 0.5, 2.09733204, 1, 21),
        (0.61235491, 4.0, 0.5, 8.82784228, 1, 51),
        (0.66666667, 9.0, 0.5, 0.30767622, 1, 1),
        (0.7, 9.7, 0.5, 0.0, 1, 1),
        (0.70979185, 8.0, 0.5, 0.98180044, 1, 11),
        (0.8, 9.8, 0.5, 0.0, 1, 1),
        (0.81135917, 6.0, 0.5, 3.77520753, 1, 31),
        (0.83154192, 3.0, 0.5, 12.68120183, 1, 61),
        (0.9, 9.9, 0.5, 0.0, 1, 1),
        (1.0, 2.49994832, 0.5, 15.14939194, 1, 71),
        (1.0, 5.18731402, 0.5, 5.53595905, 1, 41),
        (1.0, 7.07079686, 0.5, 2.07667244, 1, 21),
        (1.0, 8.15867673, 0.5, 0.92407092, 1, 11),
        (1.0, 9.0, 0.5, 0.41749062, 1, 1),
        (1.0, 9.0, 0.5, 0.41749062, 1, 2),
        (1.0, 9.33333333, 0.5, 0.30767622, 1, 1),
        (1.0, 9.57142857, 0.5, 0.21482947, 1, 1),
        (1.0, 9.75, 0.5, 0.13440196, 1, 1),
        (1.0, 9.88888889, 0.5, 0.06345987, 1, 1),
        (1.0400066, 7.0, 0.5, 2.20361733, 1, 22),
        (1.06083307, 5.0, 0.5, 6.03726831, 1, 42),
        (1.15868934, 8.0, 0.5, 1.16531307, 1, 12),
        (1.25597089, 2.0, 0.5, 18.21418588, 1, 72),
        (1.39371125, 4.0, 0.5, 9.10642069, 1, 52),
        (1.44194733, 6.0, 0.5, 4.19757573, 1, 32),
        (1.75342419, 7.0, 0.5, 2.81982124, 1, 22),
        (1.84132327, 9.0, 0.5, 0.92407092, 1, 2),
        (1.8717801, 5.0, 0.5, 6.65181663, 1, 42),
        (1.88967622, 3.0, 0.5, 13.08389058, 1, 62),
        (2.0, 1.2558761, 0.5, 24.57630394, 1, 82),
        (2.0, 2.8357984, 0.5, 13.86103603, 1, 72),
        (2.0, 4.7211973, 0.5, 7.43761494, 1, 52),
        (2.0, 6.5110226, 0.5, 3.74450869, 1, 32),
        (2.0, 8.00007947, 0.5, 1.93770701, 1, 12),
        (2.0, 8.84133588, 0.5, 1.1652939, 1, 12),
        (2.0, 9.29020815, 0.5, 0.98180044, 1, 2),
        (2.0, 9.58595475, 0.5, 0.90137293, 1, 2),
        (2.0, 9.81597989, 0.5, 0.83043083, 1, 2),
        (2.00012544, 8.0, 0.5, 1.93789772, 1, 13),
        (2.37704205, 6.0, 0.5, 5.06055662, 1, 33),
        (2.44542502, 4.0, 0.5, 9.90819793, 1, 53),
        (2.49970272, 1.0, 0.5, 27.64079914, 1, 83),
        (2.83548743, 2.0, 0.5, 18.71985648, 1, 73),
        (2.92920314, 9.0, 0.5, 2.07667244, 1, 3),
        (3.0, 0.83148486, 0.5, 30.11083496, 1, 93),
        (3.0, 1.88940228, 0.5, 19.49876924, 1, 83),
        (3.0, 3.28698638, 0.5, 12.75938985, 1, 63),
        (3.0, 5.07937977, 0.5, 7.65478139, 1, 43),
        (3.0, 7.00015821, 0.5, 4.33732448, 1, 23),
        (3.0, 8.24660575, 0.5, 2.81980207, 1, 13),
        (3.0, 8.9599934, 0.5, 2.20361733, 1, 13),
        (3.0, 9.39790983, 0.5, 2.09733204, 1, 3),
        (3.0, 9.73240437, 0.5, 2.02638995, 1, 3),
        (3.00020767, 7.0, 0.5, 4.33782286, 1, 24),
        (3.06875198, 5.0, 0.5, 7.92230807, 1, 44),
        (3.28692505, 3.0, 0.5, 14.07202724, 1, 64),
        (3.48903677, 8.0, 0.5, 3.74460723, 1, 14),
        (3.94003343, 4.0, 0.5, 11.52558563, 1, 54),
        (4.0, 0.61228533, 0.5, 33.96453147, 1, 94),
        (4.0, 1.39328875, 0.5, 23.47729073, 1, 84),
        (4.0, 2.44511796, 0.5, 16.92425856, 1, 74),
        (4.0, 3.94002526, 0.5, 11.75758164, 1, 64),
        (4.0, 6.00013079, 0.5, 7.48733924, 1, 34),
        (4.0, 7.62300176, 0.5, 5.06050603, 1, 24),
        (4.0, 8.55805267, 0.5, 4.19757573, 1, 14),
        (4.0, 9.18864083, 0.5, 3.77520753, 1, 4),
        (4.0, 9.63939593, 0.5, 3.70426544, 1, 4),
        (4.00014941, 6.0, 0.5, 7.48780988, 1, 35),
        (4.72039333, 2.0, 0.5, 19.39266633, 1, 75),
        (4.81268598, 9.0, 0.5, 5.53595905, 1, 5),
        (4.92068497, 7.0, 0.5, 7.65492574, 1, 25),
        (5.0, 0.46947412, 0.5, 36.8634256, 1, 95),
        (5.0, 1.06028978, 0.5, 26.54694563, 1, 85),
        (5.0, 1.87135558, 0.5, 20.18090394, 1, 85),
        (5.0, 3.06865624, 0.5, 15.36106403, 1, 65),
        (5.0, 5.0001141, 0.5, 11.08573482, 1, 45),
        (5.0, 6.9313041, 0.5, 7.92223486, 1, 35),
        (5.0, 8.1282199, 0.5, 6.65181663, 1, 15),
        (5.0, 8.93916693, 0.5, 6.03726831, 1, 15),
        (5.0, 9.53038942, 0.5, 5.92893723, 1, 5),
        (5.0001141, 5.0, 0.5, 11.0861454, 1, 46),
        (5.07927134, 3.0, 0.5, 15.62822007, 1, 66),
        (5.18544956, 1.0, 0.5, 27.04342109, 1, 86),
        (5.2788027, 8.0, 0.5, 7.43761494, 1, 16),
        (6.0, 0.36041286, 0.5, 39.08799504, 1, 96),
        (6.0, 0.81077304, 0.5, 28.80857582, 1, 96),
        (6.0, 1.44147854, 0.5, 22.63520794, 1, 86),
        (6.0, 2.37700418, 0.5, 18.22255708, 1, 76),
        (6.0, 4.00012377, 0.5, 14.68403545, 1, 56),
        (6.0, 6.06003039, 0.5, 11.52551242, 1, 36),
        (6.0, 7.55457498, 0.5, 9.90819793, 1, 26),
        (6.0, 8.60628875, 0.5, 9.10642069, 1, 16),
        (6.0, 9.38764509, 0.5, 8.82784228, 1, 6),
        (6.00010833, 4.0, 0.5, 14.68442526, 1, 57),
        (6.06003857, 6.0, 0.5, 11.75775673, 1, 37),
        (6.51110124, 2.0, 0.5, 19.53853999, 1, 77),
        (6.71301362, 7.0, 0.5, 12.75938985, 1, 27),
        (6.93139984, 5.0, 0.5, 15.36123911, 1, 47),
        (7.0, 0.26740321, 0.5, 40.76579383, 1, 97),
        (7.0, 0.60154156, 0.5, 30.48637461, 1, 97),
        (7.0, 1.03952535, 0.5, 24.62883003, 1, 87),
        (7.0, 1.75351945, 0.5, 20.46286516, 1, 87),
        (7.0, 3.00010831, 0.5, 17.83379884, 1, 67),
        (7.0, 4.92079341, 0.5, 15.62817756, 1, 57),
        (7.0, 6.71307495, 0.5, 14.07202724, 1, 37),
        (7.0, 8.11032378, 0.5, 13.08389058, 1, 17),
        (7.0, 9.16845808, 0.5, 12.68120183, 1, 7),
        (7.0000825, 3.0, 0.5, 17.83405871, 1, 68),
        (7.06996554, 1.0, 0.5, 24.754262, 1, 88),
        (7.1642016, 8.0, 0.5, 13.86103603, 1, 18),
        (7.50005168, 9.0, 0.5, 15.14939194, 1, 8),
        (7.55488204, 6.0, 0.5, 16.92425856, 1, 38),
        (7.62303964, 4.0, 0.5, 18.22270956, 1, 58),
        (8.0, 0.18388997, 0.5, 41.96175836, 1, 98),
        (8.0, 0.41367289, 0.5, 31.68233914, 1, 98),
        (8.0, 0.70953937, 0.5, 25.85031011, 1, 98),
        (8.0, 1.1588233, 0.5, 22.11719375, 1, 88),
        (8.0, 2.00008248, 0.5, 20.23326379, 1, 78),
        (8.0, 3.48895816, 0.5, 19.5385433, 1, 68),
        (8.0, 5.27960667, 0.5, 19.39266633, 1, 48),
        (8.0, 7.16451257, 0.5, 18.71985648, 1, 28),
        (8.0, 8.74402911, 0.5, 18.21418588, 1, 18),
        (8.00005226, 2.0, 0.5, 20.2333892, 1, 79),
        (8.11059772, 7.0, 0.5, 19.49876924, 1, 29),
        (8.12864442, 5.0, 0.5, 20.18090394, 1, 49),
        (8.15880985, 1.0, 0.5, 22.35864633, 1, 89),
        (8.2465105, 3.0, 0.5, 20.46298621, 1, 69),
        (8.55852146, 4.0, 0.5, 22.63520794, 1, 59),
        (8.60671125, 6.0, 0.5, 23.47729073, 1, 39),
        (8.7441239, 8.0, 0.5, 24.57630394, 1, 19),
        (8.84120193, 2.0, 0.5, 22.1173148, 1, 79),
        (8.93971022, 5.0, 0.5, 26.54694563, 1, 49),
        (8.96047465, 3.0, 0.5, 24.62883003, 1, 69),
        (9.0, 0.1110318, 0.5, 42.72875485, 1, 99),
        (9.0, 0.24977352, 0.5, 32.44933563, 1, 99),
        (9.0, 0.42841615, 0.5, 26.61730659, 1, 99),
        (9.0, 0.66671565, 0.5, 22.97494177, 1, 99),
        (9.0, 1.00005227, 0.5, 21.75363331, 1, 89),
        (9.0, 1.84121538, 0.5, 22.35872903, 1, 89),
        (9.0, 2.93003446, 0.5, 24.754262, 1, 79),
        (9.0, 4.81455044, 0.5, 27.04342109, 1, 59),
        (9.0, 7.50029728, 0.5, 27.64079914, 1, 29),
        (9.00002071, 1.0, 0.5, 21.75366479, 1, 90),
        (9.16851514, 7.0, 0.5, 30.11083496, 1, 30),
        (9.18922696, 4.0, 0.5, 28.80857582, 1, 60),
        (9.29046063, 2.0, 0.5, 25.85031011, 1, 80),
        (9.33329435, 1.0, 0.5, 22.97504726, 1, 90),
        (9.38771467, 6.0, 0.5, 33.96453147, 1, 40),
        (9.39845844, 3.0, 0.5, 30.48637461, 1, 70),
        (9.53052588, 5.0, 0.5, 36.8634256, 1, 50),
        (9.57158385, 1.0, 0.5, 26.61730659, 1, 90),
        (9.58632711, 2.0, 0.5, 31.68233914, 1, 80),
        (9.63958714, 4.0, 0.5, 39.08799504, 1, 60),
        (9.73259679, 3.0, 0.5, 40.76579383, 1, 70),
        (9.75022648, 1.0, 0.5, 32.44933563, 1, 90),
        (9.81611003, 2.0, 0.5, 41.96175836, 1, 80),
        (9.8889682, 1.0, 0.5, 42.72875485, 1, 90),
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


def test_prt_fmi01(function_tmpdir, targets):
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
    for col in cols:
        if col in exp_locs.dtype.names:
            assert np.allclose(locs[col], exp_locs[col])

        if col == "ilay":
            assert np.all(locs[col] == 1)

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
        nn = mg.get_node([k, i, j])
        neighbors = mg.neighbors(nn)
        assert np.isclose(nn + 1, icell, atol=1) or any(nn == n for n in neighbors)
        assert ilay == (k + 1) == 1
