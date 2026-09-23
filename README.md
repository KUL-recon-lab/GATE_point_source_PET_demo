# Gate 10 — true, scattered and random coincidences

A Gate 10 (Geant4) simulation of a point source in a water cylinder, inside a
cylindrical PET scanner, and a 3D viewer for the events it produces. Real
trajectories, real scatter angles and energies — the point is to *see* what each class
of coincidence is, rather than to be told.

The default scanner is a coarse stand-in for a whole-body system: **700 mm bore,
240 mm axial field of view**, 100 × 12 crystals. The crystals are deliberately large
(20 mm, against ~4 mm in a real one) because that keeps the count rate high enough for
a 2000-decay run to show something, and nothing here depends on resolving individual
crystals.

| file | what it is |
|---|---|
| `gate_sim.py` | the simulation: geometry, digitizer chain, coincidence sorting |
| `event_viewer.py` | the 3D event viewer (plotly) |
| `run_visu.sh` | one command: both acquisitions, then the viewer |
| `runs/` | one folder per run, and the page built from it (git-ignored, ~15 MB each) |

Nothing under `runs/` is committed: a run takes about a minute, so the repo holds the
code that produces the events, not the events.

## Running it

```bash
mamba env create -f environment.yml      # once -- pulls Geant4, a few hundred MB
mamba activate pet_gate

bash run_visu.sh                               # both acquisitions + viewer (~1 min)
python event_viewer.py                         # re-open the newest run
python event_viewer.py runs/20260922-093658    # ... or a particular one
```

`opengate` is PyPI-only and its binary wheels are tagged `manylinux_2_34` and
`macosx_15_0_arm64`, so it needs glibc ≥ 2.34 (Ubuntu 22.04 / RHEL 9 or newer) or
macOS ≥ 15 on Apple silicon. The page it produces needs none of that: `event_view.html`
is a single self-contained file that opens in any browser, so it travels on its own.

Every run writes into its own folder under `runs/`, named by timestamp, and
`gate_sim.py` **refuses to write into a folder that already holds a run** (pass
`--force` to override). Runs accumulate rather than replacing each other, so an
earlier result stays reproducible and comparable after you have moved on.

## Two acquisitions, because they want opposite things

`run_visu.sh` calls the simulation twice with the same geometry, into one run folder:

```bash
RUN=runs/$(date +%Y%m%d-%H%M%S)
python gate_sim.py --mode visu  --activity 1 --duration 0.002 --outdir $RUN/visu
python gate_sim.py --mode stats --activity 1 --duration 1.0   --outdir $RUN/stats
python event_viewer.py $RUN
```

**`--mode visu`** adds a `PhaseSpaceActor` on the world with `steps_to_store="all"` —
every step of every track, with its event and track number attached. That is what the
3D view draws, and it is the expensive actor: 3 MB for 2000 decays, so well over a
gigabyte at a million. Hence 2000 decays. It writes `tracks.json`.

**`--mode stats`** leaves it out. Counting coincidences needs no trajectories at all:
scatter comes from the interaction history inside the phantom, a cheap second hits
collection that is reduced to its scatter vertices as soon as the run ends
(285 MB → 2 MB). A million decays takes about 35 s. It writes `stats.json`, and prints

```
==================================================================
 1 001 215 emissions   (1 MBq x 1 s)
 700 mm bore, 100 x 12 crystals of 20x20x20 mm (LSO), 240 mm axial FOV
 425-650 keV window, 2.766 ns coincidence window
==================================================================
                           count    per emission   of coinc.
------------------------------------------------------------------
  emissions            1 001 215        1.000000
  singles in window      168 736        0.168531
  coincidences            20 497        0.020472    100.00 %
    true                  14 758        0.014740     72.00 %
    scattered              5 666        0.005659     27.64 %
      single scatter       4 663        0.004657     22.75 %
      higher order         1 003        0.001002      4.89 %
    random                    73        0.000073      0.36 %
==================================================================
```

The two number columns answer different questions. The **count** is what this run
produced. The **count per emission** is a property of the scanner and the phantom
rather than of the run — the absolute sensitivity — so it is what to compare across
runs of different length or activity. The **share of coincidences** is what decides
how much of the data is background.

**Single scatter** means exactly one deflection across the two photons of the pair —
the case a single-scatter simulation is built to estimate. **Higher order** is
everything else: two or more, whether both in one photon or one in each. Rayleigh
counts alongside Compton, since it is elastic but deflects the photon just the same,
and the line of response is equally wrong.

The split is worth watching against the lower energy threshold, because each extra
scatter costs more energy and so is rejected sooner:

| `--elow` | scattered | single | higher order | higher : single |
|---|---|---|---|---|
| 470 keV | 15.9 % | 14.6 % | **1.3 %** | 0.09 |
| **425 keV** | 27.0 % | 22.6 % | **4.4 %** | 0.20 |
| 350 keV | 38.8 % | 27.8 % | **11.0 %** | 0.39 |
| 250 keV | 49.4 % | 26.9 % | **22.5 %** | 0.84 |

Single scatter roughly doubles across that range while higher order grows seventeen-
fold. A high threshold does not reject scatter evenly — it rejects the multiply
scattered first, which is also why single-scatter simulation is a workable scatter
estimate at a clinical threshold and much less so at a low one.

The viewer draws the short run's events and quotes the long run's **rates**, and says
which is which — a handful of coincidences is no basis for a fraction.

## The 3D view

```bash
python event_viewer.py RUN                  # writes RUN/event_view.html and opens it
python event_viewer.py RUN --no-open        # ... or just writes it
python event_viewer.py RUN --kind true      # only one class
python event_viewer.py RUN --events 470 15  # two named events, one figure each
python event_viewer.py RUN --max-items 100  # the slider is capped
python event_viewer.py RUN --inline-js      # bundle plotly.js in (~5 MB, offline)
```

`RUN` can be the run folder, its `visu/` subfolder, or a `tracks.json` directly;
omitted, it takes the newest run under `runs/`.

`event_view.html` is one self-contained file: a section per coincidence class, each
with a **slider** that steps through the events of that class one at a time. Drag to
rotate, scroll to zoom, hover a track for its energies. A class with nothing in it is
left out and named in the header instead — at 1 MBq a 2000-decay run usually has no
randoms at all.

The classes are the ones the stats table counts, scatter included:

| section | what is in it |
|---|---|
| True coincidence | neither photon deflected in the water |
| Single-scatter coincidence | exactly one deflection across the pair |
| Higher-order scatter coincidence | two or more |
| Random coincidence | the two photons came from different decays |
| Not a coincidence | the event the scanner never recorded as a pair |

`event_view.html` loads plotly.js from the CDN, which is what keeps it under 1 MB;
`--inline-js` bundles it instead (~5 MB) for a copy that works with no network.

The report in the left column carries, for both photons of the pair:

* the energy it had **when it left the water** (`exit`) — 511 keV means it came out
  untouched, less means it scattered on the way, and `(nx scat)` says how many times;
* `dep / rec` — what the crystal actually absorbed, and what the digitizer chain
  reported after energy blurring. `dep` is below `exit` when the photon Compton-scatters
  in the crystal and the scattered photon escapes; `rec` is the only one of the three a
  real scanner ever sees. (There is no separate "arrived at the crystal" number: it
  equalled `exit` for 965 of the 966 photons in a reference run, the exception being one
  that scattered in the air gap.)
* `flight` — emission to first energy deposit, per photon;
* `nominal dt` — the difference of those two flight times, and the position it puts the
  annihilation at along the LOR (`c·Δt/2`);
* `measured dt` — the same after the timing blur, i.e. what a real scanner sees.

With the source at the isocentre the nominal Δt of a true coincidence is 0 by
construction, so the measured spread is pure timing resolution — which makes the
default geometry a clean demonstration of exactly that.

For a random pair there is no nominal Δt (the two photons come from different decays,
so the viewer says so instead), and the two photons that were *not* recorded are drawn
thin and labelled `partner`, which shows why the pair was formed.

Colour says what the water did to that photon: **green** = straight through,
**pink** = deflected once, **purple** = deflected twice or more. Weight says what the
scanner did with it: **thick and opaque** = accepted into the energy window,
**thinner and faded** = detected but rejected, **thin and pale** = never detected.

Because the hue is per photon and the class is per pair, a higher-order section often
shows two *pink* tracks rather than a purple one — 1 + 1 rather than 2 + 0. The report
spells the sum out so there is nothing to infer:

```
ev0207 p1  exit  490.8 keV (1x scat)   flight 1.234 ns
           dep  490.8 / rec  498.7 keV
ev0207 p2  exit  441.3 keV (1x scat)   flight 1.249 ns
           dep  441.3 / rec  430.2 keV
nominal  dt = -0.016 ns ->   -2.3 mm on LOR
measured dt = -0.142 ns ->  -21.3 mm on LOR
water scatters: 1 + 1 = 2  (higher order)
```

The scanner is rebuilt from the constants in `tracks.json`, so a re-parameterised ring
carries through to the viewer with no edit.

## Everything the scanner is, from the command line

No constant has to be edited. `python gate_sim.py --help` lists them all.

```bash
python gate_sim.py --mode visu --outdir runs/wide/visu \
    --bore 830 --crystal 25 4 4 --n-tangential 400 --n-axial 60 --axial-gap 1.0 \
    --phantom-diameter 300 --phantom-length 400 \
    --source-transaxial 60 --source-axial -40 \
    --elow 470 --ehigh 620 --eres 0.12 --ctr 400
python event_viewer.py runs/wide
```

| group | arguments | default |
|---|---|---|
| scanner | `--bore` face-to-face diameter | 700 mm |
| | `--crystal DEPTH TANG AXIAL` | 20 20 20 mm |
| | `--n-tangential` around the ring | 100 |
| | `--n-axial` rings along the axis | 12 |
| | `--axial-gap` between those rings | 0 mm |
| | `--material` LSO or BGO | LSO |
| | `--world` air or vacuum in the bore | air |
| phantom | `--phantom-diameter` | 200 mm |
| | `--phantom-length` | **the axial FOV** (240 mm) |
| source | `--source-transaxial`, `--source-axial` | 0, 0 (isocentre) |
| digitizer | `--elow`, `--ehigh` | 425, 650 keV |
| | `--eres` at 511 keV, as a fraction | 0.10 |
| | `--ctr` coincidence timing resolution FWHM | 200 ps |
| | `--coinc-window` | **from the geometry** (2.77 ns) |
| | `--kill-below` stop tracking in the water below | **from the window** (370 keV) |
| run | `--mode`, `--activity`, `--duration`, `--seed`, `--outdir`, `--force` | |

Two defaults are computed rather than fixed, because the sensible value depends on the
scanner:

* **the phantom is as long as the axial field of view**, so it fills it exactly. Give
  `--phantom-length` to override.
* **the coincidence window is as wide as the physics needs and no wider.** A real
  coincidence can be delayed by at most the time light takes to cross the longest line
  of response the scanner can record — straight across the ring *and* from one end of
  the axial FOV to the other, √((2R)² + Z²) = 753 mm here, so 2.511 ns — plus what the
  timing blur adds, taken as 3σ of the CTR. That gives 2.766 ns by default, and the
  script prints the arithmetic every run. Give `--coinc-window` to override.

The world box, the "same detector" distance the sorter uses, and everything the viewer
draws are derived from these too, so nothing is left stale by changing one. `run_visu.sh`
passes extras to both acquisitions through `GEOM`:

```bash
GEOM="--bore 830 --n-axial 40 --ctr 400" bash run_visu.sh
```

### LSO or BGO

`--material` switches the scintillator. Both are defined by composition and density
(neither is in Geant4's NIST database): LSO is Lu₂SiO₅ at 7.40 g/cm³, BGO is
Bi₄Ge₃O₁₂ at 7.13 g/cm³. Bismuth's Z of 83 against lutetium's 71 is what matters —
BGO stops more 511 keV photons, and stops a larger share of them photoelectrically,
so more deposits land inside the energy window. Same geometry, same seed, 10⁵ decays:

| per emission | LSO | BGO | |
|---|---|---|---|
| singles in window | 0.1701 | 0.1909 | ×1.122 |
| coincidences | 0.02064 | 0.02625 | ×1.272 |
| true | 0.01491 | 0.01993 | ×1.337 |
| scatter fraction | 27.2 % | 23.5 % | |

Worth noticing: the coincidence gain is the **square** of the single-photon gain
(1.122² = 1.260 against the 1.272 measured, 1 % apart), because a coincidence needs
both photons detected. That the simulation reproduces it is a free check on the
digitizer chain.

What this comparison does *not* include is the reason most scanners use LSO anyway:
BGO's light yield is far lower, so a real BGO system has worse energy and timing
resolution. `--eres` and `--ctr` are left independent of `--material` rather than
guessed at — raise them when you switch, or the comparison flatters BGO.

### Not tracking photons that cannot be detected

A photon that has lost enough energy in the water can no longer produce a single
inside the window, however lucky the blurring is. Following it any further is wasted
work, and Geant4 can be told to stop: `G4UserLimits::SetUserMinEkine`, enforced by
`G4UserSpecialCuts`, reached in opengate as

```python
sim.physics_manager.user_limits_particles = ["all_charged", "gamma"]
sim.physics_manager.set_min_ekine(phantom.name, 369.6 * keV)
```

`--kill-below` sets the threshold; the default is computed, not guessed. With the
digitizer's `InverseSquare` law the resolution is R(E) = R₀·√(E₀/E), so
σ(E) = R₀·√(E₀·E)/2.355, and the threshold is the deposit that 3σ of blurring could
still lift to `--elow` — solving E + 3σ(E) = E_low, which is 369.6 keV for the
defaults. Measured against the digitizer rather than assumed: σ fitted from the true
deposit and the blurred readout of a real run agrees with that formula to ~3 % over
250–520 keV.

**It is worth about 30 %.** A million decays with the default geometry:

| | runtime | singles/emission | coinc/emission |
|---|---|---|---|
| no cut | 65.2 s | 0.168452 | 0.020294 |
| cut at 369.6 keV | **45.2 s** | 0.167476 | 0.020277 |
| | | −1.7 σ | −0.1 σ |

Every class agrees within Poisson error — the largest discrepancy across singles,
coincidences, trues, both scatter orders and randoms is 1.7 σ, and all but one are
under 1 σ.

**Two things about it are load-bearing, and both were checked by breaking them.**

*The cut is applied to the phantom only.* `G4UserSpecialCuts` does not merely stop the
track — it deposits the remaining kinetic energy **locally**. In water that is
harmless, since nothing scores energy there. Inside a crystal it is a bug that
flatters the scanner: a photon that Compton-scatters and would have escaped instead
dumps its remainder into the same crystal, turning escapes into full-energy events.
Adding the same cut to the crystals, same seed:

| | singles/emission | coinc/emission |
|---|---|---|
| phantom only | 0.168218 | 0.020267 |
| phantom + crystals | 0.199561 | 0.028405 |
| | **+19 %** | **+40 %** |

Geant4 regions propagate to daughter volumes, which is also why the cut is not put on
the world.

*The 3σ margin is not decoration.* Cutting at `--elow` itself, with no margin, already
loses counts — and precisely the scattered ones, which is the worst thing to lose
silently:

| `--kill-below` | singles/emission | vs the default |
|---|---|---|
| **369.6 keV (3σ, default)** | 0.168218 | — |
| 425 keV (= `--elow`) | 0.164832 | −2.0 % |
| 475 keV | 0.145704 | −13.4 % |

`--kill-below 0` turns it off. It is never applied in `--mode visu`: that run is 2000
decays, so the time saved is nothing, and truncated trajectories would cost the 3D
view exactly the low-energy tracks it draws as "detected but rejected" and "never
detected".

**The cut is on gammas only, and that is not shorthand.** opengate attaches
`G4UserSpecialCuts` — which is what enforces `min_ekine` — only to the particles named
in `user_limits_particles`; its default `"all_charged"` gets a `G4StepLimiter`
instead, which enforces `max_step_size` and nothing else. So charged particles are
untouched by the energy cut either way. Since no max step size is set anywhere here,
`["gamma"]` and `["all_charged", "gamma"]` give bit-identical counts at 10⁶ decays.

**And there is almost nothing charged to track.** The 1 mm production cut is a *range*
cut, which Geant4 turns into a per-material energy threshold, and that makes this a
photon-transport simulation in all but name. In water 1 mm is a ~351 keV electron
while the most energetic Compton electron a 511 keV photon can produce is 340.7 keV —
so no Compton electron is ever created in the phantom; the energy is deposited on the
spot. In LSO the threshold is far higher again, so nothing is created there at all,
which is what a digitizer summing deposits per crystal wants. Counted on a 2000-decay
run:

| | gamma steps | electron steps |
|---|---|---|
| water | 11 815 | 13 |
| crystals | 2 064 | **0** |

All 88 charged tracks in that run are born in the **air**, where 1 mm of range is a
threshold of a few keV — which is what `--world vacuum` below is for. Photoelectric
absorption in water does clear the threshold, but it is ~0.02 % of 511 keV
interactions there: a vacuum-bore run of the same size produced exactly **one**
electron track, a 510.5 keV photoelectron.

Raising the global production cut to 1 m would also remove the air electrons and is
worth a further ~9 % (48.3 s → 44.1 s at 10⁶ decays), at the price of a ~1.2 σ shift
in the counts that is probably but not certainly noise. It is left at 1 mm because
that is a physics knob rather than a free optimisation.

### Air or vacuum in the bore

`--world vacuum` (against the default `air`) fills the world with `G4_Galactic`. It
is worth ~5 %, and it leaves the run with essentially no charged transport anywhere.
What it costs is the air itself: 10⁶ decays, same seed,

| per emission | air | vacuum | change | |
|---|---|---|---|---|
| singles | 0.168213 | 0.169030 | +0.49 % | +1.4 σ |
| coincidences | 0.020430 | 0.020706 | +1.35 % | +1.4 σ |
| true | 0.014608 | 0.014748 | +0.95 % | +0.8 σ |
| scattered | 0.005760 | 0.005895 | +2.34 % | +1.3 σ |

and 44.8 s → 42.7 s. The singles gain is the air that is no longer there: at 511 keV
the attenuation of 25 cm of air is 0.26 %, against the +0.49 % measured, and the
coincidence gain is close to its square, as it must be when a coincidence needs both
photons. Both are around 1.4 σ, so the effect is real in direction and plausible in
size rather than sharply resolved at this statistics.

The default is **air**, because that is what a real scanner has and 5 % is not much
to pay for it. `--world vacuum` is there for when it is: the bias it introduces is an
overestimate of absolute sensitivity by about a percent, far smaller than the
idealisations already in the geometry, so for relative comparisons it costs nothing.

The **energy window is applied by Gate itself**, during the simulation, by the
`DigitizerEnergyWindowsActor`. To try a different one, run the simulation again with
`--elow` / `--ehigh`: a visualisation run takes a second and a full stats run under a
minute, which is cheaper than any machinery for changing it afterwards.

## What the coincidence sorter does

The pairs come from Gate's own `coincidences_sorter`, run on the singles that passed
the energy window, with:

* **`policy="TakeAllGoods"`** — when a multiple occurs, keep every good pair in it
  rather than discarding the lot. A random that steals a photon from a true
  coincidence is exactly the case worth showing, and `RemoveMultiples` deletes it.
* **`min_transaxial_distance`** — the sorter's way of saying "not the same detector".
  There are no blocks in this scanner, so every crystal is its own detector and
  adjacent centres are 2·R·sin(π/N) apart; a shade more than that (×1.05) rules out
  pairs that cannot be a line of response, and nothing else. It is a guard, not a
  filter: in the default geometry no recorded pair comes anywhere near it.

Both are computed from the ring, so they follow `--bore`, `--n-tangential` and
`--crystal` without an edit.

Classification is the standard one: **random** — the two singles belong to different
annihilations; **scattered** — same annihilation, but at least one photon Compton- or
Rayleigh-scattered inside the phantom; **true** — same annihilation, neither did.

## What the runs show

**Randoms are made by the activity, not by the physics.** A random is two photons from
different decays that happen to arrive within the coincidence window; nothing extra is
simulated to produce them. At a fixed number of decays the rate goes as N²/T:

| activity | decays | coincidences | true | scattered | random |
|---|---|---|---|---|---|
| 1 MBq | 10⁶ | 20 497 | 72.0 % | 27.6 % | **0.4 %** |
| 5 MBq | 2·10⁵ | 4 095 | 69.8 % | 28.7 % | **1.5 %** |
| 20 MBq | 2·10⁵ | 4 358 | 65.8 % | 27.1 % | **7.0 %** |
| 80 MBq | 2·10⁵ | 5 174 | 54.5 % | 22.5 % | **23.0 %** |

The scatter fraction barely moves across that range while the random fraction moves by
a factor of fifteen — the two backgrounds have entirely different origins, and only one
of them is a property of the patient. Reproduce any row with, for example,

```bash
python gate_sim.py --mode stats --activity 80 --duration 0.0025 --outdir runs/a80/stats
```

**The coincidence window is already as narrow as it can be.** Re-sorting the 1 MBq run
at three widths:

| window | coincidences | true | scattered | random |
|---|---|---|---|---|
| **2.766 ns (the default)** | 20 497 | 14 758 | 5 666 | **73** |
| 4 ns | 20 531 | 14 758 | 5 666 | **107** |
| 6 ns | 20 597 | 14 758 | 5 666 | **173** |

Widening it buys **not one true and not one scattered** coincidence, and randoms
linearly — they are uniform in Δt, because nothing correlates them. That is the whole
argument for sizing the window from the geometry, and it is also exactly what
delayed-window randoms estimation exploits.

In this particular run the largest │Δt│ of a real coincidence is only 0.386 ns, far
inside the 2.511 ns geometric limit, because a point source at the isocentre sits
equidistant from both detectors. The limit is what a source anywhere in the field of
view could produce, which is what the window has to allow for.

**And by TOF.** Taking `c·Δt/2` as a position along the line of response puts a random
nowhere near the LOR at all, and often outside the phantom entirely. The viewer says so
in as many words when you step onto one: *different decays -> that position is
meaningless*.

**Why a given random happened.** Stepping onto a random draws four tracks, not two: the
recorded pair, and **each decay's other photon** — because a random can only form when
those two were *not* recorded, so they are the explanation. Typically both partners
were scattered in the water and then rejected or lost, which is why a heavier patient
raises scatter and randoms together. When a partner was itself accepted, that decay was
recorded as its own coincidence too and the random is a **multiple** — the case
`RemoveMultiples` would have thrown away.

**Widening the energy window buys the wrong counts.** The 2000-decay visualisation run,
re-simulated at three windows:

| window (keV) | singles | coincidences | true | scattered |
|---|---|---|---|---|
| 470–650 | 298 | 33 | 31 | 2 |
| **425–650 (default)** | **347** | **42** | **32** | **10** |
| 350–700 | 432 | 54 | 33 | 21 |

One extra true for eight extra scattered: that is what the lower threshold buys.

## Three things that are easy to get wrong

1. **`UnscatteredPrimaryFlag` is not an object-scatter flag.** At 511 keV roughly two
   thirds of the interactions in LSO are Compton, so the flag is 0 for most perfectly
   good true coincidences. Classification here uses the interaction history inside the
   phantom (a second hits collection attached to the water cylinder).
2. **`coincidences_sorter` orders the two singles in time**, so `GlobalTime1 <=
   GlobalTime2` always and the raw `c·(t1−t2)/2` is never positive. Histogramming it
   gives a folded distribution and an apparent timing resolution a factor of two too
   good. The displacement has to be signed by geometry.
3. **The step-level actor is what makes a run expensive**, not the number of decays. A
   `PhaseSpaceActor` on the world with `steps_to_store="all"` writes every step of
   every track: 3 MB for 2000 decays, and over a gigabyte at a million. Nothing about
   counting coincidences needs it, which is the whole reason `--mode visu` and
   `--mode stats` are separate runs.
