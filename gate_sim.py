"""Gate 10 (opengate) simulation: a point source in a water cylinder, inside a toy
cylindrical PET scanner.

Produces REAL event kinematics -- true trajectories, real scatter angles and energies
-- for the 3D event viewer and for true / scattered / random statistics. Every
dimension of the scanner, the phantom, the source and the digitizer is a command-line
argument; the defaults are a 300 mm bore with 40 x 5 crystals of 20 mm.

Two modes, because a run that feeds the 3D view and a run that measures rates want
opposite things:

    python gate_sim.py --mode visu  --activity 1 --duration 0.002 --outdir runs/x/visu
    python gate_sim.py --mode stats --activity 1 --duration 1.0   --outdir runs/x/stats

"visu" adds a PhaseSpaceActor on the world with steps_to_store="all" -- every step of
every track, with its event and track number attached. That is what the viewer draws
and it is the expensive actor: 3 MB at 2000 decays, well over a gigabyte at a million.
"stats" leaves it out; counting coincidences needs no trajectories, since scatter comes
from the interaction history inside the phantom.

Every run writes into its own --outdir and refuses to overwrite one that already has
output in it, so runs accumulate instead of replacing each other. run_visu.sh puts a
matching pair under one timestamped folder and points the viewer at it.

The energy window is applied by Gate itself, during the simulation. To try a different
one, run the simulation again -- it takes a minute.
"""
import json
import math
import os, sys, time, argparse
import numpy as np
import uproot
import opengate as gate
from opengate.geometry.utility import get_circular_repetition

mm, cm, ns, ps, sec, Bq, keV, MBq = (gate.g4_units.mm, gate.g4_units.cm, gate.g4_units.ns,
                                     gate.g4_units.ps, gate.g4_units.s, gate.g4_units.Bq,
                                     gate.g4_units.keV, 1e6 * gate.g4_units.Bq)

# ---------------------------------------------------------------- the scanner
# Defaults only: every one of these is a command-line argument (see main()), and
# configure() below recomputes everything derived from them. The defaults are a
# coarse stand-in for a whole-body scanner: 700 mm bore, 240 mm axial field of
# view. The crystals are deliberately large -- a real WB scanner has ~4 mm ones --
# because that keeps the count rate high enough for a 2000-decay run to show
# something, and nothing here depends on resolving individual crystals.
BORE_MM = 700.0              # face-to-face diameter of the detector ring
CRYSTAL = (20.0, 20.0, 20.0)  # radial depth, tangential, axial [mm]
MATERIAL = "LSO"             # crystal material, one of CRYSTAL_MATERIALS below
# What fills the bore. Air is what a real scanner has; vacuum is ~5 % faster and
# makes the run photon-transport only, at the cost of the ~0.5 % of photons that
# air would have scattered or absorbed between phantom and detector.
WORLD_MATERIALS = {"air": "G4_AIR", "vacuum": "G4_Galactic"}
WORLD = "air"
N_TANG, N_AXIAL = 100, 12    # crystals around the ring, and along the axis
AXIAL_GAP_MM = 0.0           # edge-to-edge gap between neighbouring axial rings

# The two scintillators worth comparing. Neither is in Geant4's NIST database, so
# both are defined by composition and density. BGO is the denser and much
# higher-Z of the two (bismuth, Z = 83, against lutetium's 71), so it stops more
# 511 keV photons and stops a larger share of them photoelectrically -- more
# full-energy deposits, and more of them inside the energy window. What it gives up
# is light yield, hence the poorer energy and timing resolution of a real BGO
# scanner: if you switch, raise --eres and --ctr to match.
CRYSTAL_MATERIALS = {
    #        elements          atoms per molecule   density [g/cm3]
    "LSO": (["Lu", "Si", "O"], [2, 1, 5], 7.40),      # Lu2SiO5
    "BGO": (["Bi", "Ge", "O"], [4, 3, 12], 7.13),     # Bi4Ge3O12
}

# ---------------------------------------------------------------- the phantom
PHANTOM_D_MM = 200.0         # water cylinder diameter [mm]
# None means "as long as the scanner", which is the usual thing to want: the
# phantom then fills the axial field of view exactly. A number overrides it.
PHANTOM_L_MM = None

# ----------------------------------------------------------------- the source
# One point source at (transaxial, 0, axial); both 0 puts it at the isocentre.
SOURCE_T_MM, SOURCE_A_MM = 0.0, 0.0

# --------------------------------------------------------------- the digitizer
CTR_PS = 200.0               # coincidence timing resolution (FWHM) [ps]
E_RES = 0.10                 # energy resolution at 511 keV (fraction, 0.10 = 10 %)
E_WINDOW = (425.0, 650.0)    # energy acceptance window [keV]
# None means "as wide as the physics needs and no wider" -- see configure().
COINC_WINDOW_NS = None

# Stop tracking a photon in the WATER once it is too degraded to ever pass the
# energy window. None means "work it out from the window and the energy
# resolution" -- see configure(); 0 turns it off.
KILL_BELOW_KEV = None
KILL_SIGMA = 3.0

# How many sigma of timing blur to leave beyond the longest real delay when the
# coincidence window is chosen automatically. 3 sigma keeps ~99.7 % of the true
# coincidences that sit right at the geometric limit.
WINDOW_SIGMA = 3.0

# Derived; configure() keeps them in step with the above.
RING_RADIUS = PHANTOM_R = PHANTOM_HL = RING_PITCH = 0.0
N_DET = N_RINGS = 0
CRYSTAL_CENTRE_R = ADJACENT_CRYSTAL_MM = 0.0
AXIAL_FOV_MM = LONGEST_LOR_MM = KILL_BELOW = 0.0
COINC_WINDOW = 0.0           # the resolved window [ns], auto or as given
SOURCE_POS = (0.0, 0.0, 0.0)
C_MM_PER_NS = 299.792458


def configure(**kw):
    """Set any of the constants above and recompute what depends on them.

    Everything the simulation, the analysis and the viewer need about the scanner
    comes from here, so a re-parameterised ring carries through to all three without
    an edit anywhere else.
    """
    g = globals()
    for k, v in kw.items():
        if k not in g:
            raise KeyError(f"unknown setting {k!r}")
        g[k] = v

    g["RING_RADIUS"] = BORE_MM / 2.0                 # inner radius: crystal faces
    g["N_DET"], g["N_RINGS"] = N_TANG, N_AXIAL
    g["RING_PITCH"] = CRYSTAL[2] + AXIAL_GAP_MM      # centre-to-centre, axially
    # the axial field of view: outer face to outer face of the end crystal rings
    g["AXIAL_FOV_MM"] = (N_AXIAL - 1) * g["RING_PITCH"] + CRYSTAL[2]
    g["PHANTOM_R"] = PHANTOM_D_MM / 2.0
    g["PHANTOM_HL"] = (AXIAL_FOV_MM if PHANTOM_L_MM is None
                       else PHANTOM_L_MM) / 2.0
    # centre-to-centre distance between two neighbouring crystals: the scale of "the
    # same detector" for the coincidence sorter. Derived, so changing the ring size
    # or the crystal count carries through instead of leaving a stale constant.
    g["CRYSTAL_CENTRE_R"] = RING_RADIUS + CRYSTAL[0] / 2.0
    g["ADJACENT_CRYSTAL_MM"] = (2.0 * g["CRYSTAL_CENTRE_R"]
                                * math.sin(math.pi / N_TANG))

    # The longest line of response the scanner can record: straight across the ring
    # AND from one end of the axial field of view to the other.
    g["LONGEST_LOR_MM"] = math.hypot(2.0 * g["CRYSTAL_CENTRE_R"],
                                     (N_AXIAL - 1) * g["RING_PITCH"])
    # A real coincidence can be delayed by at most the time light takes to cross
    # that LOR -- an annihilation at one end, so one photon arrives immediately and
    # the other travels the whole length. Anything beyond that is a random, apart
    # from what the timing blur adds, so the window is that limit plus a few sigma.
    # Widening it further buys no trues and linearly more randoms.
    g["COINC_WINDOW"] = (COINC_WINDOW_NS if COINC_WINDOW_NS is not None else
                         LONGEST_LOR_MM / C_MM_PER_NS
                         + WINDOW_SIGMA * (CTR_PS * 1e-3) / 2.355)

    # The deposit E that the energy blurring could still push up to the lower
    # threshold, at KILL_SIGMA sigma. With the digitizer's InverseSquare law the
    # resolution is R(E) = R_ref*sqrt(E_ref/E), so sigma(E) = R_ref*sqrt(E_ref*E)/2.355
    # and E + k*sigma(E) = Elow is a quadratic in sqrt(E). A photon carrying less
    # than that can never produce an accepted single on its own, so tracking it
    # further through the water is wasted work.
    k = KILL_SIGMA * E_RES * math.sqrt(511.0) / 2.355
    g["KILL_BELOW"] = (KILL_BELOW_KEV if KILL_BELOW_KEV is not None else
                       ((-k + math.sqrt(k * k + 4.0 * E_WINDOW[0])) / 2.0) ** 2)
    g["SOURCE_POS"] = (float(SOURCE_T_MM), 0.0, float(SOURCE_A_MM))
    return g


configure()


def read_run_info(outdir):
    try:
        with open(os.path.join(outdir, "run_info.json")) as fh:
            return json.load(fh)
    except OSError:
        return {}


def scatter_counts(outdir):
    """How many times each annihilation photon scattered in the water.

    {(EventID, TrackID): n}, from the slim phantom_scatter.npz when it is there and
    the raw phantom_hits.root otherwise. One row per Compton or Rayleigh vertex, so
    counting rows counts scatters. Rayleigh is included: it is elastic, so the energy
    is unchanged, but the photon has been deflected all the same and the line of
    response is just as wrong.

    This is the object-scatter history, the only correct basis for the true/scattered
    split -- UnscatteredPrimaryFlag is NOT: at 511 keV most interactions in LSO are
    Compton, so that flag is 0 for most perfectly good trues.

    It counts every vertex in the water, including the rare one a photon makes on its
    way back through after scattering in a crystal. Measured against the stricter
    "before it first left the water" definition on a 4016-photon run: 3 photons
    differ, 0.07 %.
    """
    import numpy as _np
    slim = os.path.join(outdir, "phantom_scatter.npz")
    if os.path.exists(slim):
        z = _np.load(slim, allow_pickle=False)
        ev, tr = z["EventID"], z["TrackID"]
    else:
        f = uproot.open(os.path.join(outdir, "phantom_hits.root"))
        d = f[f.keys()[0]].arrays(library="np")
        proc = _np.asarray(d["ProcessDefinedStep"], dtype=object)
        name = _np.asarray(d["ParticleName"], dtype=object)
        m = _np.array([(p in ("compt", "Rayl")) and (n == "gamma")
                       for p, n in zip(proc, name)])
        ev, tr = d["EventID"][m], d["TrackID"][m]
    import collections as _c
    return _c.Counter(zip(ev.tolist(), tr.tolist()))


def slim_phantom_hits(outdir):
    """Reduce phantom_hits.root to the scatter vertices and drop the raw file.

    The phantom hits collection records every step inside the water -- 285 MB for a
    million decays -- but the only thing anything downstream asks of it is which
    photons scattered. Keep that, throw the rest away.
    """
    import numpy as _np
    raw = os.path.join(outdir, "phantom_hits.root")
    if not os.path.exists(raw):
        return
    f = uproot.open(raw)
    d = f[f.keys()[0]].arrays(library="np")
    proc = _np.asarray(d["ProcessDefinedStep"], dtype=object)
    name = _np.asarray(d["ParticleName"], dtype=object)
    m = _np.array([(p in ("compt", "Rayl")) and (n == "gamma")
                   for p, n in zip(proc, name)])
    _np.savez_compressed(os.path.join(outdir, "phantom_scatter.npz"),
                         EventID=d["EventID"][m].astype(_np.int32),
                         TrackID=d["TrackID"][m].astype(_np.int16))
    before = os.path.getsize(raw)
    os.remove(raw)
    after = os.path.getsize(os.path.join(outdir, "phantom_scatter.npz"))
    print(f"  slimmed: {int(m.sum())} scatter vertices of {m.size} phantom steps, "
          f"{before/1e6:.0f} MB -> {after/1e6:.1f} MB")


def _n(x):
    """1001215 -> '1 001 215'. Thin-spaced groups read far better in a column."""
    return f"{x:,}".replace(",", " ")


def format_stats(st, info=None):
    """The counting result as a table.

    Two numbers per row, because they answer different questions: the count is what
    this run produced, and the count per emission is a property of the scanner and
    the phantom -- the absolute sensitivity, comparable across runs of any length.
    The share of coincidences is the third, and the one that decides how much of the
    data is background.
    """
    info = info or {}
    n_dec = st["n_decays"] or 1
    n_cc = st["n_coincidences"] or 1
    W = 66
    rows = [("emissions", st["n_decays"], None),
            ("singles in window", st["n_singles"], None),
            ("coincidences", st["n_coincidences"], st["n_coincidences"]),
            ("  true", st["n_true"], st["n_true"]),
            ("  scattered", st["n_scattered"], st["n_scattered"]),
            ("    single scatter", st.get("n_single_scatter"), None),
            ("    higher order", st.get("n_multi_scatter"), None),
            ("  random", st["n_random"], st["n_random"])]
    out = ["", "=" * W]
    act, dur = info.get("activity_MBq"), info.get("duration_s")
    head = f" {_n(st['n_decays'])} emissions"
    if act and dur:
        head += f"   ({act:g} MBq x {dur:g} s)"
    out.append(head)
    geo = info.get("crystal_mm", list(CRYSTAL))
    out.append(f" {info.get('bore_mm', BORE_MM):.0f} mm bore, "
               f"{info.get('n_tangential', N_TANG)} x {info.get('n_axial', N_AXIAL)}"
               f" crystals of {geo[0]:.0f}x{geo[1]:.0f}x{geo[2]:.0f} mm"
               f" ({info.get('material', MATERIAL)}),"
               f" {info.get('axial_fov_mm', AXIAL_FOV_MM):.0f} mm axial FOV")
    out.append(f" {st['energy_window_keV'][0]:.0f}-{st['energy_window_keV'][1]:.0f}"
               f" keV window, {st['coincidence_window_ns']:.3f} ns coincidence window")
    out.append("=" * W)
    out.append(f"{'':<20}{'count':>12}{'per emission':>16}{'of coinc.':>12}")
    out.append("-" * W)
    for label, n, share in rows:
        if n is None:                     # a run from before this was recorded
            continue
        # the two scatter-order rows are a breakdown of the row above, so they get
        # their share of the coincidences too -- it is what they are usually quoted as
        share = n if share is None and label.startswith("    ") else share
        per = f"{n / n_dec:.6f}"
        pct = "" if share is None else f"{100.0 * share / n_cc:8.2f} %"
        out.append(f"  {label:<18}{_n(n):>12}{per:>16}{pct:>12}")
    out.append("=" * W)
    return "\n".join(out)


def analyse_stats_run(outdir):
    """Coincidence rates from a long acquisition, with no track data at all.

    The step-level PhaseSpaceActor that feeds the 3D view is the expensive actor --
    a few MB at 2000 decays, over a GB at a million -- and none of it is needed for
    counting. Scatter is taken from the phantom interaction history instead.
    """
    import collections as _c
    scat = scatter_counts(outdir)
    sw = uproot.open(os.path.join(outdir, "singles_window.root"))
    tree = sw[sw.keys()[0]]
    cc = sort_coincidences(tree, scat)
    k = _c.Counter(c["kind"] for c in cc)
    # scattered coincidences split by order: exactly one deflection over the two
    # photons, or more than one
    n_single = sum(1 for c in cc if c["kind"] == "scattered" and c["n_scatter"] == 1)
    n_multi = k["scattered"] - n_single
    info = read_run_info(outdir)
    n_dec = read_n_events(outdir)
    tot = max(len(cc), 1)
    out = {"n_decays": n_dec,
           "activity_MBq": info.get("activity_MBq"),
           "duration_s": info.get("duration_s"),
           "material": MATERIAL, "world": WORLD,
           "energy_window_keV": list(E_WINDOW),
           "coincidence_window_ns": COINC_WINDOW,
           "n_singles": int(tree.num_entries),
           "n_coincidences": len(cc), "n_true": k["true"],
           "n_scattered": k["scattered"],
           "n_single_scatter": n_single, "n_multi_scatter": n_multi,
           "n_random": k["random"],
           "frac_true": round(k["true"] / tot, 5),
           "frac_scattered": round(k["scattered"] / tot, 5),
           "frac_single_scatter": round(n_single / tot, 5),
           "frac_multi_scatter": round(n_multi / tot, 5),
           "frac_random": round(k["random"] / tot, 5),
           # per emission: the scanner's absolute sensitivity to each class,
           # comparable across runs of any length
           "per_emission": {
               "singles": round(int(tree.num_entries) / max(n_dec, 1), 8),
               "coincidences": round(len(cc) / max(n_dec, 1), 8),
               "true": round(k["true"] / max(n_dec, 1), 8),
               "scattered": round(k["scattered"] / max(n_dec, 1), 8),
               "single_scatter": round(n_single / max(n_dec, 1), 8),
               "multi_scatter": round(n_multi / max(n_dec, 1), 8),
               "random": round(k["random"] / max(n_dec, 1), 8)}}
    with open(os.path.join(outdir, "stats.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print(format_stats(out, info))
    print(f"written to {os.path.join(outdir, 'stats.json')}")
    return out


def read_n_events(outdir, fallback=None):
    """How many annihilations Geant4 actually simulated, from the statistics actor.

    activity x duration is only the MEAN: the decays are sampled, so the run almost
    never contains exactly that many. Falls back to counting events in the track file
    when stats.txt is missing (an outdir written before this actor existed).
    """
    path = os.path.join(outdir, "stats.txt")
    try:
        with open(path) as fh:
            txt = fh.read()
    except OSError:
        return fallback
    try:                                  # opengate 10 writes json despite the name
        return int(json.loads(txt)["events"]["value"])
    except Exception:
        pass
    for line in txt.splitlines():         # older "# NumberOfEvents = N" format
        if "NumberOfEvents" in line:
            return int(float(line.split("=")[-1].strip()))
    return fallback


def sort_coincidences(windowed_singles_tree, scatter_of):
    """Coincidences from opengate's own sorter, on the singles Gate accepted.

    The tree is singles_window.root -- the output of the DigitizerEnergyWindowsActor,
    so the energy window has already been applied by Gate itself. scatter_of maps
    (EventID, TrackID) to the number of times that photon scattered in the water.

    "TakeAllGoods" keeps every good pair inside a multiple instead of discarding the
    whole multiple, which matters here: a random that steals a photon from a true
    coincidence is exactly the case worth showing, and RemoveMultiples deletes it.
    min_transaxial_distance is the sorter's way of saying "not the same detector":
    adjacent crystal centres are 2*R*sin(pi/N) apart, so a shade more than that rules
    out a pair that cannot be a line of response, and nothing else. Computed from the
    ring rather than hardcoded, so it follows the scanner.

    Each pair gets "n_scatter", the total number of water scatters over its two
    photons. 1 is a single-scatter coincidence -- one photon, deflected once, the
    case a single-scatter simulation is built to estimate -- and anything more is
    higher order.
    """
    from opengate.actors.coincidences import coincidences_sorter
    cc = coincidences_sorter(windowed_singles_tree, COINC_WINDOW,
                             policy="TakeAllGoods",
                             min_transaxial_distance=ADJACENT_CRYSTAL_MM * 1.05,
                             transaxial_plane="xy",
                             max_axial_distance=1e9, return_type="pd")
    out = []
    for r in cc.itertuples():
        e1, t1 = int(r.EventID1), int(r.TrackID1)
        e2, t2 = int(r.EventID2), int(r.TrackID2)
        n_sc = scatter_of.get((e1, t1), 0) + scatter_of.get((e2, t2), 0)
        if e1 != e2:                                         # two different decays
            kind = "random"
        elif n_sc:
            kind = "scattered"
        else:
            kind = "true"
        out.append({"e1": e1, "t1": t1, "e2": e2, "t2": t2, "kind": kind,
                    "n_scatter": n_sc,
                    "dt_ns": round(float(r.GlobalTime1 - r.GlobalTime2), 6)})
    return out


def build(activity_MBq, duration_s, seed, outdir, mode="stats"):
    """mode "visu" adds the step-level actor that feeds the 3D view; "stats" leaves
    it out, which is the whole difference between a 2000-decay run and a million."""
    ctr = CTR_PS
    visu = mode == "visu"
    sim = gate.Simulation()
    sim.random_seed = seed
    sim.number_of_threads = 1
    sim.progress_bar = False
    sim.output_dir = outdir
    sim.g4_verbose = False

    # 1.2 m is plenty for the scanner in this file; the max() only matters if you
    # enlarge the ring past about a 0.8 m bore, where a fixed world would start
    # cutting through the detectors.
    sim.world.size = [max(1200.0,
                          3.0 * (RING_RADIUS + CRYSTAL[0]),
                          1.5 * (N_RINGS * RING_PITCH + CRYSTAL[2]),
                          3.0 * PHANTOM_HL) * mm] * 3
    sim.world.material = WORLD_MATERIALS[WORLD]

    phantom = sim.add_volume("Tubs", "phantom")
    phantom.rmin, phantom.rmax, phantom.dz = 0, PHANTOM_R * mm, PHANTOM_HL * mm
    phantom.material = "G4_WATER"
    phantom.color = [0.3, 0.6, 1.0, 0.25]

    crystal = sim.add_volume("Box", "crystal")
    crystal.size = [CRYSTAL[0] * mm, CRYSTAL[1] * mm, CRYSTAL[2] * mm]
    gcm3 = gate.g4_units.g_cm3
    elements, atoms, rho = CRYSTAL_MATERIALS[MATERIAL]
    sim.volume_manager.material_database.add_material_nb_atoms(
        MATERIAL, elements, atoms, rho * gcm3)
    crystal.material = MATERIAL
    tr, rot = get_circular_repetition(
        N_DET, [(RING_RADIUS + 0.5 * CRYSTAL[0]) * mm, 0, 0],
        angular_step_deg="auto_full_circle", axis=(0, 0, 1))
    translations, rotations = [], []
    zs = (np.arange(N_RINGS) - 0.5 * (N_RINGS - 1)) * RING_PITCH
    for z in zs:
        for t, r in zip(tr, rot):
            translations.append([t[0], t[1], z * mm])
            rotations.append(r)
    crystal.translation = translations
    crystal.rotation = rotations
    crystal.color = [1, 0.6, 0.2, 1]

    sim.physics_manager.physics_list_name = "G4EmStandardPhysics_option4"
    # A 1 mm RANGE cut, which Geant4 turns into a per-material energy threshold.
    # It is why this is a photon-transport simulation in all but name:
    #   * in water 1 mm is a ~351 keV electron, and the most energetic Compton
    #     electron a 511 keV photon can make is 340.7 keV, so no Compton electron
    #     is ever created in the phantom -- the energy is deposited on the spot.
    #     Photoelectric absorption does clear the threshold, but it is ~0.02 % of
    #     511 keV interactions in water: a 2000-decay vacuum-bore run produced
    #     exactly ONE electron track, a 510.5 keV photoelectron.
    #   * in LSO the threshold is far higher (7.4x the density), so nothing is
    #     created there at all -- zero electron steps in the crystals, measured --
    #     which is what a digitizer summing deposits per crystal wants.
    #   * air is the exception: 1 mm of range there is a few keV, so it is where
    #     charged tracks actually appear (88 of them in that same run). --world
    #     vacuum removes them and with them the last of the charged transport.
    sim.physics_manager.global_production_cuts.all = 1.0 * mm

    # Stop following a photon in the WATER once it is too degraded to pass the
    # energy window even with the blurring in its favour (see configure()). This is
    # Geant4's G4UserLimits/G4UserSpecialCuts, set per volume.
    #
    # The phantom ONLY, and that restriction is the whole point. G4UserSpecialCuts
    # does not merely stop the track: it deposits the remaining kinetic energy
    # locally. In the water that is harmless -- nothing scores energy there. Inside
    # a crystal it would be a bug that flatters the scanner: a photon that Compton-
    # scatters and would have escaped instead dumps its remainder into the same
    # crystal, turning escapes into full-energy events and inflating the photopeak.
    # Regions propagate to daughter volumes, which is also why the cut is not put on
    # the world.
    if visu:
        # the 3D view is 2000 decays, so the time saved is nothing, and truncated
        # trajectories would be a real loss -- the low-energy tracks are exactly
        # what it draws as "detected but rejected" and "never detected"
        pass
    elif KILL_BELOW > 0:
        # Gammas only, and that is not a shorthand for "and charged particles are
        # cut too". opengate attaches G4UserSpecialCuts -- which is what enforces
        # min_ekine -- only to the particles named here; the default "all_charged"
        # gets a G4StepLimiter instead, which enforces max_step_size and nothing
        # else. Since no max_step_size is set anywhere in this file, dropping
        # "all_charged" changes nothing: checked at 1e6 decays, with and without,
        # the run gives bit-identical counts. (If you ever do set a max step size
        # for charged particles, add "all_charged" back here.)
        sim.physics_manager.user_limits_particles = ["gamma"]
        sim.physics_manager.set_min_ekine(phantom.name, KILL_BELOW * keV)

    src = sim.add_source("GenericSource", "src")
    src.particle = "back_to_back"
    src.energy.mono = 511 * keV
    src.position.type = "point"
    src.position.translation = [c * mm for c in SOURCE_POS]
    src.direction.type = "iso"
    src.activity = activity_MBq * MBq

    # Gate's own event counter. The requested activity x duration is only the MEAN
    # number of decays; what actually gets simulated is Poisson around it, so this is
    # the honest answer to "how many annihilations were simulated".
    stats = sim.add_actor("SimulationStatisticsActor", "Stats")
    stats.track_types_flag = True
    stats.output_filename = "stats.txt"

    # ---- hits in the crystals -------------------------------------------------
    hc = sim.add_actor("DigitizerHitsCollectionActor", "Hits")
    hc.attached_to = crystal.name
    hc.authorize_repeated_volumes = True
    hc.attributes = ["EventID", "TrackID", "GlobalTime", "TotalEnergyDeposit",
                     "PostPosition", "PreStepUniqueVolumeID", "UnscatteredPrimaryFlag",
                     "EventPosition", "PreKineticEnergy"]
    hc.output_filename = None          # nothing downstream reads raw hits

    ad = sim.add_actor("DigitizerAdderActor", "Singles0")
    ad.attached_to = crystal.name
    ad.authorize_repeated_volumes = True
    ad.input_digi_collection = "Hits"
    ad.policy = "EnergyWeightedCentroidPosition"
    ad.output_filename = None

    eb = sim.add_actor("DigitizerBlurringActor", "SinglesE")
    eb.attached_to = crystal.name
    eb.authorize_repeated_volumes = True
    eb.input_digi_collection = "Singles0"
    eb.blur_attribute = "TotalEnergyDeposit"
    eb.blur_method = "InverseSquare"
    eb.blur_resolution = E_RES
    eb.blur_reference_value = 511 * keV
    eb.output_filename = None

    tb = sim.add_actor("DigitizerBlurringActor", "SinglesRaw")
    tb.attached_to = crystal.name
    tb.authorize_repeated_volumes = True
    tb.input_digi_collection = "SinglesE"
    tb.blur_attribute = "GlobalTime"
    tb.blur_method = "Gaussian"
    # per-single time resolution: the coincidence resolution is the quadratic sum of two
    tb.blur_fwhm = max(ctr, 1e-6) / np.sqrt(2.0) * ps
    # every single, before the energy window: the 3D view needs it to draw the
    # photons the scanner detected and then REJECTED. Off for a stats run.
    tb.output_filename = "singles_all.root" if visu else None

    ew = sim.add_actor("DigitizerEnergyWindowsActor", "EnergyWindow")
    ew.attached_to = crystal.name
    ew.authorize_repeated_volumes = True
    ew.input_digi_collection = "SinglesRaw"
    ew.channels = [{"name": "Singles", "min": E_WINDOW[0] * keV,
                    "max": E_WINDOW[1] * keV}]
    ew.output_filename = "singles_window.root"

    # ---- hits inside the water phantom (to reconstruct the real photon paths) --
    ph = sim.add_actor("DigitizerHitsCollectionActor", "PhantomHits")
    ph.attached_to = phantom.name
    ph.attributes = ["EventID", "TrackID", "ParentID", "PrePosition", "PostPosition",
                     "ProcessDefinedStep", "PreKineticEnergy", "PostKineticEnergy",
                     "ParticleName", "GlobalTime"]
    # the object-scatter history, the only correct basis for true vs scattered.
    # A visu run gets the same information, with identity attached, from the
    # Tracks actor below, so it does not need this one.
    ph.output_filename = "phantom_hits.root" if not visu else None

    if visu:
        # THE expensive actor: every step of every track, with its EventID and
        # TrackID attached. This is what the 3D view draws, and it is also why the
        # visualisation run is kept to a couple of thousand decays -- at a million
        # it would write well over a gigabyte, and none of it is needed to count
        # coincidences. "stats" mode leaves it out.
        tracks = sim.add_actor("PhaseSpaceActor", "Tracks")
        tracks.attached_to = "world"
        tracks.steps_to_store = "all"
        tracks.attributes = ["EventID", "TrackID", "ParentID", "ParticleName",
                             "PrePosition", "PostPosition", "KineticEnergy",
                             "PostKineticEnergy", "TotalEnergyDeposit",
                             "ProcessDefinedStep", "PreStepUniqueVolumeID",
                             # times, so the flight from the decay to the crystal can
                             # be read off rather than divided out of a path length
                             "PreGlobalTime", "GlobalTime",
                             "TimeFromBeginOfEvent"]
        tracks.output_filename = "tracks.root"

    sim.run_timing_intervals = [[0, duration_s * sec]]
    return sim


def write_tracks_json(outdir):
    """Everything the 3D viewer needs, as one plain JSON file.

    One entry per photon: its true EventID/TrackID, its exact path, the energy it had
    when it first left the water, whether it interacted in the water on the way out,
    and what the scanner made of it.  Nothing here is inferred from geometry -- it all
    comes straight out of the hit files, joined on (EventID, TrackID).

    Two definitions worth stating, because the obvious ones are subtly wrong:

    * "left the water" means the FIRST exit.  A photon can leave the phantom, Compton
      scatter in a crystal, and fly back through the water on its way out the far
      side; taking the last phantom step would then report the post-detector energy.
    * "scattered in the water" likewise counts only interactions before that first
      exit -- the ones that actually corrupt the line of response.  A Rayleigh scatter
      counts: it is elastic, so the exit energy stays 511 keV, but the photon has been
      deflected all the same.
    """
    import collections
    import uproot

    f = uproot.open(os.path.join(outdir, "singles_all.root"))
    sg = f[f.keys()[0]].arrays(library="np")
    # (event, track) -> the (time, energy) of every single the digitizer reported for
    # that photon. The time is already blurred by the CTR, so it is what the scanner
    # measures rather than what happened.
    per_photon = collections.defaultdict(list)
    for ev, tr, t, e, v in zip(sg["EventID"], sg["TrackID"], sg["GlobalTime"],
                               sg["TotalEnergyDeposit"],
                               sg["PreStepUniqueVolumeID"]):
        per_photon[(int(ev), int(tr))].append((float(t), float(1000.0 * e), str(v)))

    # ... and which of those singles Gate's energy window actually kept. Taking this
    # from the actor's own output rather than re-deriving it means the viewer and the
    # coincidence sorter can never disagree about what "accepted" means.
    w = uproot.open(os.path.join(outdir, "singles_window.root"))
    wg = w[w.keys()[0]].arrays(library="np")
    accepted = {(int(ev), int(tr), str(v))
                for ev, tr, v in zip(wg["EventID"], wg["TrackID"],
                                     wg["PreStepUniqueVolumeID"])}

    f = uproot.open(os.path.join(outdir, "tracks.root"))
    st = f[f.keys()[0]].arrays(library="np")
    name = np.asarray(st["ParticleName"], dtype=object)
    proc = np.asarray(st["ProcessDefinedStep"], dtype=object)
    vol = np.asarray(st["PreStepUniqueVolumeID"], dtype=object)
    pre = np.stack([st[f"PrePosition_{c}"] for c in "XYZ"], -1)
    post = np.stack([st[f"PostPosition_{c}"] for c in "XYZ"], -1)

    # Only the two annihilation photons. ParentID == 0 is what makes them primaries:
    # Geant4 also tracks the recoil electrons, and occasionally a tertiary gamma (a
    # few keV of bremsstrahlung from a Compton electron, say), and those are neither
    # p1 nor p2 of anything. Event 634 of this run has one.
    steps = collections.defaultdict(list)
    for i in range(len(name)):
        if name[i] == "gamma" and int(st["ParentID"][i]) == 0:
            steps[(int(st["EventID"][i]), int(st["TrackID"][i]))].append(i)

    # when each annihilation happened: both photons start at the same instant, so the
    # earliest step time of the event is the decay time
    t_emit = {}
    for (ev, tr), idx in steps.items():
        t0 = float(st["PreGlobalTime"][idx[0]])
        t_emit[ev] = min(t_emit.get(ev, t0), t0)

    tracks_out = []
    for (ev, tr), idx in sorted(steps.items()):
        pts, n_scatter, exit_keV, left = [pre[idx[0]]], 0, None, False
        # a photon can miss the water entirely -- the source may sit outside the
        # cylinder -- and then "the energy it left the water with" is simply the
        # energy it started with, not "absorbed in the water"
        entered, t_hit = False, None
        # the true (unblurred) deposit per crystal, keyed by the same volume id the
        # singles carry, so a recorded single can be matched back to the physics
        # that made it -- which matters because some photons deposit in TWO
        # crystals and only one of them becomes the recorded single.
        in_crystal = {}
        for i in idx:
            if np.linalg.norm(pre[i] - pts[-1]) > 1e-6:
                pts.append(pre[i])
            pts.append(post[i])
            in_phantom = str(vol[i]).startswith("phantom")
            if str(vol[i]).startswith("crystal"):
                c = in_crystal.setdefault(str(vol[i]), {"dep": 0.0})
                c["dep"] += float(1000.0 * st["TotalEnergyDeposit"][i])
                # when the photon reached the first crystal that took energy off
                # it: the START of that step, which is when it got there
                if t_hit is None and st["TotalEnergyDeposit"][i] > 0:
                    t_hit = float(st["PreGlobalTime"][i])
            if not left:
                if in_phantom:
                    entered = True
                    if str(proc[i]) in ("compt", "Rayl"):
                        n_scatter += 1
                    exit_keV = float(1000.0 * st["PostKineticEnergy"][i])
                elif entered:               # first step outside: it is out
                    left = True
        if entered and not left:            # absorbed before it ever got out
            exit_keV = None
        elif not entered:                   # never in the water: nothing taken off
            exit_keV = float(1000.0 * st["KineticEnergy"][idx[0]])

        # Every single the digitizer made from this photon, with the physics behind
        # each one attached. This -- not a single pre-selected number -- is what lets
        # the energy window be changed afterwards: pick_single() below re-runs the
        # choice for any window, and so can a viewer that only has tracks.json.
        # What the scanner made of this photon. A photon can leave a single in more
        # than one crystal -- Compton in the first, absorbed in the second -- so the
        # one that counts is the one Gate's energy window kept; failing that, the
        # biggest of the ones it rejected.
        all_s = per_photon.get((ev, tr), [])
        kept = [x for x in all_s if (ev, tr, x[2]) in accepted]
        if kept:
            klass, pick = "accepted", min(kept, key=lambda x: x[0])
        elif all_s:
            klass, pick = "rejected", max(all_s, key=lambda x: x[1])
        else:
            klass, pick = "undetected", None
        keV = None if pick is None else round(pick[1], 3)
        t_meas = None if pick is None else round(pick[0] - t_emit[ev], 6)
        c = in_crystal.get(pick[2], {}) if pick is not None else {}
        dep_keV = None if c.get("dep") is None else round(c["dep"], 3)
        tracks_out.append({
            # water: did the phantom deflect it at all; n_scatter: how many times.
            # Only deflections BEFORE it first left the water are counted -- those
            # are the ones that corrupt the line of response.
            "event": ev, "track": tr, "water": n_scatter > 0,
            "n_scatter": n_scatter, "exit_keV": exit_keV,
            "klass": klass, "keV": keV,
            # keV. exit_keV above is what the photon carried out of the water --
            # which, bar the odd scatter in the air gap, is also what reached the
            # crystal, so there is no separate "arrived" number. dep is what it
            # actually left in the crystal the recorded single came from, less than
            # exit when it Compton-scatters there and the scattered photon escapes;
            # keV is that same deposit after the digitizer's energy blurring, which
            # is the only one of the three a real scanner ever sees.
            "dep_keV": dep_keV,
            # ns. tof is the true flight from the decay to the first crystal that took
            # energy; t_meas_ns is the blurred timestamp the scanner actually records,
            # both measured from the decay so the numbers are small and comparable.
            "tof_ns": None if t_hit is None else round(t_hit - t_emit[ev], 6),
            "t_meas_ns": t_meas,
            "pts": [[round(float(c), 4) for c in p] for p in pts]})

    # ---- what the scanner actually records ------------------------------------
    scat = {(t["event"], t["track"]): t["n_scatter"] for t in tracks_out}
    n_decays = read_n_events(outdir, fallback=len({t["event"] for t in tracks_out}))
    coincidences = sort_coincidences(w[w.keys()[0]], scat)
    k = collections.Counter(c["kind"] for c in coincidences)

    n = {x: sum(1 for t in tracks_out if t["klass"] == x)
         for x in ("accepted", "rejected", "undetected")}
    stats = {"n_decays": n_decays,
             "n_photon_tracks": len(tracks_out),
             "n_singles": int(w[w.keys()[0]].num_entries),
             "n_accepted_photons": n["accepted"],
             "n_scattered_in_water": int(sum(t["water"] for t in tracks_out)),
             "n_single_scatter": sum(1 for c in coincidences
                                     if c["kind"] == "scattered"
                                     and c["n_scatter"] == 1),
             "n_multi_scatter": sum(1 for c in coincidences
                                    if c["kind"] == "scattered"
                                    and c["n_scatter"] > 1),
             "n_absorbed_in_water": sum(1 for t in tracks_out
                                        if t["exit_keV"] is None),
             "n_coincidences": len(coincidences), "n_true": k["true"],
             "n_scattered": k["scattered"], "n_random": k["random"]}
    doc = {"energy_window_keV": list(E_WINDOW),
           "coincidence_window_ns": COINC_WINDOW,
           "coincidences": coincidences,
           "stats": stats,
           "source_mm": list(SOURCE_POS),
           "phantom_radius_mm": PHANTOM_R,
           "phantom_halflength_mm": PHANTOM_HL,
           # enough to redraw the scanner from this file alone, so the viewer needs
           # no geometry export and nothing from Geant4
           "scanner": {"ring_radius_mm": RING_RADIUS,
                       "crystal_mm": list(CRYSTAL),
                       "material": MATERIAL,
                       "n_det": N_DET, "n_rings": N_RINGS,
                       "ring_pitch_mm": RING_PITCH,
                       "axial_gap_mm": AXIAL_GAP_MM},
           "digitizer": {"energy_resolution": E_RES, "ctr_ps": CTR_PS},
           "tracks": tracks_out}
    with open(os.path.join(outdir, "tracks.json"), "w") as fh:
        json.dump(doc, fh)
    print(f"tracks.json: {len(tracks_out)} photon tracks from {n_decays} simulated "
          f"annihilations -- {n['accepted']} accepted, {n['rejected']} detected but "
          f"rejected, {n['undetected']} never detected; "
          f"{stats['n_scattered_in_water']} scattered in the water on the way out, "
          f"{stats['n_absorbed_in_water']} absorbed before leaving it")
    print(f"  {E_WINDOW[0]:.0f}-{E_WINDOW[1]:.0f} keV: {stats['n_singles']} singles "
          f"-> {len(coincidences)} coincidences = {k['true']} true + "
          f"{k['scattered']} scattered + {k['random']} random")


def main():
    ap = argparse.ArgumentParser(
        description="Gate 10 PET simulation: a point source in a water cylinder, "
                    "inside a toy cylindrical scanner.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    g = ap.add_argument_group("run")
    g.add_argument("--mode", choices=["visu", "stats"], default="stats",
                   help="visu adds the step-level actor that records every "
                        "trajectory and writes tracks.json for the 3D viewer -- use "
                        "a couple of thousand decays. stats leaves it out, which is "
                        "what makes a million decays affordable")
    g.add_argument("--activity", type=float, default=1.0, help="MBq")
    g.add_argument("--duration", type=float, default=0.002, help="s")
    g.add_argument("--seed", type=int, default=12345)
    g.add_argument("--outdir", default=None,
                   help="where this run writes. Default: "
                        "runs/<timestamp>/<mode>. An existing non-empty folder is "
                        "refused, so runs accumulate rather than overwrite")
    g.add_argument("--force", action="store_true",
                   help="write into --outdir even if it already holds a run")

    g = ap.add_argument_group("scanner")
    g.add_argument("--bore", type=float, default=BORE_MM,
                   help="face-to-face diameter of the detector ring [mm]")
    g.add_argument("--crystal", type=float, nargs=3, default=list(CRYSTAL),
                   metavar=("DEPTH", "TANG", "AXIAL"),
                   help="crystal size [mm]: radial depth, tangential, axial")
    g.add_argument("--n-tangential", type=int, default=N_TANG,
                   help="crystals around the ring")
    g.add_argument("--n-axial", type=int, default=N_AXIAL,
                   help="crystal rings along the axis")
    g.add_argument("--axial-gap", type=float, default=AXIAL_GAP_MM,
                   help="edge-to-edge gap between axial crystal rings [mm]")
    g.add_argument("--world", choices=sorted(WORLD_MATERIALS), default=WORLD,
                   help="what fills the bore. vacuum is ~5 %% faster and leaves no "
                        "charged particle anywhere in the run, but drops the ~0.5 %% "
                        "of photons air would have scattered or absorbed on the way "
                        "to the detector")
    g.add_argument("--material", choices=sorted(CRYSTAL_MATERIALS), default=MATERIAL,
                   help="crystal material. BGO stops more 511 keV photons than LSO "
                        "and more of them photoelectrically, but a real BGO scanner "
                        "has poorer energy and timing resolution -- raise --eres and "
                        "--ctr to match")

    g = ap.add_argument_group("phantom and source")
    g.add_argument("--phantom-diameter", type=float, default=PHANTOM_D_MM,
                   help="water cylinder diameter [mm]")
    g.add_argument("--phantom-length", type=float, default=None,
                   help="water cylinder length [mm]. Default: the scanner's axial "
                        "field of view, so the phantom fills it exactly")
    g.add_argument("--source-transaxial", type=float, default=SOURCE_T_MM,
                   help="point source offset from the axis, along x [mm]")
    g.add_argument("--source-axial", type=float, default=SOURCE_A_MM,
                   help="point source offset along the axis, z [mm]")

    g = ap.add_argument_group("digitizer")
    g.add_argument("--elow", type=float, default=E_WINDOW[0],
                   help="lower energy acceptance threshold [keV]")
    g.add_argument("--ehigh", type=float, default=E_WINDOW[1],
                   help="upper energy acceptance threshold [keV]")
    g.add_argument("--eres", type=float, default=E_RES,
                   help="energy resolution at 511 keV, as a fraction (0.10 = 10 %%)")
    g.add_argument("--ctr", type=float, default=CTR_PS,
                   help="coincidence timing resolution, FWHM [ps]")
    g.add_argument("--kill-below", type=float, default=None, metavar="KEV",
                   help="stop tracking a photon in the water once it drops below "
                        "this energy [keV]. Default: the deposit that the energy "
                        "blurring could still push up to --elow, at 3 sigma, so "
                        "nothing detectable is lost. 0 disables it. Never applied "
                        "inside the crystals, and not applied at all in visu mode")
    g.add_argument("--coinc-window", type=float, default=None,
                   help="coincidence time window [ns]. Default: the time light "
                        "takes to cross the longest LOR the scanner can record, "
                        "plus 3 sigma of timing blur -- as wide as the physics "
                        "needs and no wider")

    a = ap.parse_args()
    configure(BORE_MM=a.bore, CRYSTAL=tuple(a.crystal), N_TANG=a.n_tangential,
              N_AXIAL=a.n_axial, AXIAL_GAP_MM=a.axial_gap, MATERIAL=a.material,
              WORLD=a.world,
              PHANTOM_D_MM=a.phantom_diameter, PHANTOM_L_MM=a.phantom_length,
              SOURCE_T_MM=a.source_transaxial, SOURCE_A_MM=a.source_axial,
              E_WINDOW=(a.elow, a.ehigh), E_RES=a.eres, CTR_PS=a.ctr,
              COINC_WINDOW_NS=a.coinc_window, KILL_BELOW_KEV=a.kill_below)
    print(f"scanner: {BORE_MM:.0f} mm bore, {N_TANG} x {N_AXIAL} crystals of "
          f"{CRYSTAL[0]:.0f}x{CRYSTAL[1]:.0f}x{CRYSTAL[2]:.0f} mm {MATERIAL}, "
          f"{AXIAL_FOV_MM:.0f} mm axial FOV, {WORLD}-filled bore")
    print(f"phantom: {PHANTOM_D_MM:.0f} mm diameter x {2*PHANTOM_HL:.0f} mm long"
          + ("  (= the axial FOV)" if a.phantom_length is None else ""))
    if a.mode == "stats" and KILL_BELOW > 0:
        print(f"tracking cut: photons below {KILL_BELOW:.1f} keV are stopped in the "
              f"water ({KILL_SIGMA:.0f} sigma below the {E_WINDOW[0]:.0f} keV "
              f"threshold; --kill-below 0 disables)")
    print(f"coincidence window: {COINC_WINDOW:.3f} ns"
          + (f"  (longest LOR {LONGEST_LOR_MM:.0f} mm -> "
             f"{LONGEST_LOR_MM / C_MM_PER_NS:.3f} ns, plus "
             f"{WINDOW_SIGMA:.0f} sigma of {CTR_PS:.0f} ps blur)"
             if a.coinc_window is None else "  (given)"))

    outdir = a.outdir or os.path.join(
        "runs", time.strftime("%Y%m%d-%H%M%S"), a.mode)
    if os.path.isdir(outdir) and os.listdir(outdir) and not a.force:
        raise SystemExit(
            f"{outdir} already holds a run. Pick another --outdir, or pass --force "
            f"to overwrite it.")
    os.makedirs(outdir, exist_ok=True)

    sim = build(a.activity, a.duration, a.seed, outdir, mode=a.mode)
    with open(os.path.join(outdir, "run_info.json"), "w") as fh:
        json.dump({"mode": a.mode,
                   "activity_MBq": a.activity, "duration_s": a.duration,
                   # nominal: activity x duration. The number actually simulated
                   # is Poisson around it and is added below, once Geant4 has run.
                   "n_decays": a.activity * 1e6 * a.duration,
                   "seed": a.seed,
                   "source_mm": list(SOURCE_POS),
                   "bore_mm": BORE_MM, "crystal_mm": list(CRYSTAL),
                   "material": MATERIAL, "world": WORLD,
                   "n_tangential": N_TANG, "n_axial": N_AXIAL,
                   "axial_gap_mm": AXIAL_GAP_MM,
                   "phantom_diameter_mm": PHANTOM_D_MM,
                   "phantom_length_mm": 2 * PHANTOM_HL,
                   "axial_fov_mm": AXIAL_FOV_MM,
                   "longest_lor_mm": LONGEST_LOR_MM,
                   "ctr_ps": CTR_PS, "energy_resolution": E_RES,
                   "kill_below_keV": KILL_BELOW if a.mode == "stats" else 0.0,
                   "coinc_window_ns": COINC_WINDOW,
                   "energy_window_keV": list(E_WINDOW)}, fh, indent=2)
    sim.run()
    n_sim = read_n_events(outdir)
    info = read_run_info(outdir)
    info["n_decays_simulated"] = n_sim
    with open(os.path.join(outdir, "run_info.json"), "w") as fh:
        json.dump(info, fh, indent=2)
    print("DONE", a.activity, "MBq x", a.duration, "s  ->",
          f"{a.activity * 1e6 * a.duration:.3e} decays nominal, {n_sim} simulated")
    print(f"  output in {outdir}")

    if a.mode == "visu":
        write_tracks_json(outdir)
    else:
        slim_phantom_hits(outdir)
        analyse_stats_run(outdir)


if __name__ == "__main__":
    main()
