"""An interactive 3D event viewer for a Gate 10 run: plotly and nothing else.

    python event_viewer.py runs/20260922-101500     -> writes event_view.html and opens it
    python event_viewer.py                          -> the most recent run under runs/
    python event_viewer.py RUN --kind true          -> only one class of coincidence
    python event_viewer.py RUN --events 470 15      -> only those two events

A run folder holds the two acquisitions gate_sim.py produces:

    <run>/visu/tracks.json     the events this page draws (a short run)
    <run>/stats/stats.json     the true / scattered / random rates (a long run)

Either path can also be given directly. Everything drawn comes out of tracks.json --
the scanner is redrawn from the constants the simulation recorded there, so this needs
no geometry export and nothing from Geant4.

Shown for each event:

  * the ring, the water cylinder and the source, rotatable;
  * both photons of the annihilation, coloured by what the scanner made of them;
  * a report of what the photon carried out of the phantom, what a crystal recorded,
    and the two time differences -- the true one and the one the scanner measured.

The events drawn come from a short run of a couple of thousand decays. The true /
scattered / random RATES in the header come from the long run, because a handful of
coincidences is no basis for a fraction; the header says which is which.
"""
import argparse
import glob
import json
import math
import os
import webbrowser

import numpy as np
import plotly.graph_objects as go

HERE = os.path.dirname(os.path.abspath(__file__))
C_MM_PER_NS = 299.792458

# The html is committed, so it must not change just because two machines have
# different plotly versions: include_plotlyjs="cdn" writes whatever plotly.js
# the local package pins, and the file then shows up dirty in git after every
# run. Pin it here instead. Any 4.x renders these figures; --inline-js still
# bundles the locally installed one, for working offline.
PLOTLY_JS = "https://cdn.plot.ly/plotly-4.0.0.min.js"

# hue = what the water did, weight = whether the scanner kept it
# hue by how many times the water deflected the photon, so the scatter order of a
# pair is visible in the tracks and not only in the section heading
HUE = {0: "#00C48A",        # straight through
       1: "#E040A0",        # deflected once
       2: "#8B3FD6"}        # deflected twice or more
FATE_STYLE = {                                      # opacity, width
    "accepted":   (1.00, 7),
    "rejected":   (0.55, 4),
    # an undetected photon is often the ONLY thing in the frame (a "not a
    # coincidence" event where both flew straight out), so it has to stay
    # visible on its own while still reading as subordinate
    "undetected": (0.55, 3),
}
# The classes the page is split into, in display order. Scatter is separated by
# order because the two behave differently: a single scatter is what a
# single-scatter simulation models and is the bulk of it at a clinical threshold,
# while higher order is the part that estimate does not describe.
KIND_LABEL = {"true": "True coincidence",
              "single_scatter": "Single-scatter coincidence",
              "multi_scatter": "Higher-order scatter coincidence",
              # only for a run written before the scatter order was recorded
              "scattered": "Scattered coincidence",
              "random": "Random coincidence",
              "lost": "Not a coincidence"}
KINDS = list(KIND_LABEL)


def class_of(c):
    """The viewer's class for one coincidence, finer than the sorter's.

    The sorter says true / scattered / random; here "scattered" is split by the
    pair's total number of water deflections, which it records as n_scatter.
    """
    if c["kind"] != "scattered":
        return c["kind"]
    n = c.get("n_scatter")
    if n is None:
        return "scattered"
    return "single_scatter" if n == 1 else "multi_scatter"


# --------------------------------------------------------------------- data
def find_run(path=None):
    """(tracks.json, stats.json-or-None) for a run folder, or the newest run."""
    if path is None:
        runs = sorted(glob.glob(os.path.join(HERE, "runs", "*")))
        if not runs:
            raise SystemExit(
                "no run to show. Make one with\n    bash run_visu.sh\n"
                "or point me at a folder: python event_viewer.py <run>")
        path = runs[-1]
    if os.path.isfile(path):                         # a tracks.json given directly
        tracks = path
        run = os.path.dirname(os.path.dirname(os.path.abspath(path)))
    elif os.path.isfile(os.path.join(path, "tracks.json")):     # the visu folder
        tracks = os.path.join(path, "tracks.json")
        run = os.path.dirname(os.path.abspath(path))
    else:                                                        # the run folder
        tracks = os.path.join(path, "visu", "tracks.json")
        run = os.path.abspath(path)
        if not os.path.isfile(tracks):
            raise SystemExit(f"no visu/tracks.json in {path}")
    stats = os.path.join(run, "stats", "stats.json")
    return tracks, (stats if os.path.isfile(stats) else None), run


def load(path=None):
    """The run, ready to draw: tracks indexed by event, classes, randoms, rates."""
    tracks_json, stats_json, run = find_run(path)
    with open(tracks_json) as fh:
        doc = json.load(fh)
    doc["_run"] = run
    doc["_rates"] = None
    if stats_json:
        with open(stats_json) as fh:
            doc["_rates"] = json.load(fh)

    doc["_by_event"] = by = {}
    for t in doc["tracks"]:
        by.setdefault(t["event"], []).append(t)
    for g in by.values():
        g.sort(key=lambda t: t["track"])
    recorded = {c["e1"]: class_of(c) for c in doc["coincidences"]
                if c["e1"] == c["e2"]}
    doc["_class"] = {ev: recorded.get(ev, "lost") for ev in by}
    doc["_randoms"] = [c for c in doc["coincidences"] if c["kind"] == "random"]
    return doc


def events_of(doc, kind):
    """Event numbers of one class. For randoms the unit is a PAIR, so this returns
    the index into the random list rather than an event number."""
    if kind == "random":
        return list(range(len(doc["_randoms"])))
    return sorted(ev for ev, k in doc["_class"].items() if k == kind)


def kinds_present(doc):
    """The coincidence classes this run actually has, in display order.

    A class with nothing in it is left out rather than drawn empty -- at a low
    activity a short run often has no randoms at all.
    """
    return [k for k in KINDS if events_of(doc, k)]


def window_label(doc):
    lo, hi = doc["energy_window_keV"]
    return f"{lo:.0f}-{hi:.0f} keV window"


def summary(doc):
    """One line of run statistics, saying which run each number comes from."""
    g = doc["stats"]
    drawn = (f"drawn: {g['n_decays']} annihilations, "
             f"{g['n_coincidences']} coincidences")
    r = doc.get("_rates")
    if not r:
        return drawn + " &nbsp;|&nbsp; no long run alongside it"
    # a real thin space, not &thinsp;: this string goes into the plotly figure
    # title as well as the page header, and plotly's SVG text renders the entity
    # literally
    n = f"{r['n_decays']:,}".replace(",", "\u2009")
    return (drawn + f" &nbsp;|&nbsp; rates from a {n}-decay run: "
            f"{100 * r['frac_true']:.1f} % true, "
            f"{100 * r['frac_scattered']:.1f} % scattered, "
            f"{100 * r['frac_random']:.1f} % random")


# ----------------------------------------------------------------- geometry
def _box(centre, half, u, v, w):
    """Eight corners of a box with the given half-extents along u, v, w."""
    c = np.asarray(centre, float)
    return np.array([c + a * half[0] * u + b * half[1] * v + d * half[2] * w
                     for a in (-1, 1) for b in (-1, 1) for d in (-1, 1)])


def scanner_traces(doc, crystal_opacity=0.16):
    """The detector ring, the water cylinder and the point source."""
    sc = doc["scanner"]
    R = sc["ring_radius_mm"] + sc["crystal_mm"][0] / 2.0        # centre radius
    dr, dt, dz = (x / 2.0 for x in sc["crystal_mm"])
    n, nz, pitch = sc["n_det"], sc["n_rings"], sc["ring_pitch_mm"]
    z0 = -(nz - 1) * pitch / 2.0

    verts, faces = [], []
    quads = [(0, 1, 3, 2), (4, 5, 7, 6), (0, 1, 5, 4),
             (2, 3, 7, 6), (0, 2, 6, 4), (1, 3, 7, 5)]
    for i in range(n):
        a = 2 * math.pi * i / n
        u = np.array([math.cos(a), math.sin(a), 0.0])           # radial
        v = np.array([-math.sin(a), math.cos(a), 0.0])          # tangential
        w = np.array([0.0, 0.0, 1.0])
        for k in range(nz):
            base = len(verts)
            verts.extend(_box(R * u + (z0 + k * pitch) * w, (dr, dt, dz), u, v, w))
            for q in quads:
                faces.append((base + q[0], base + q[1], base + q[2]))
                faces.append((base + q[0], base + q[2], base + q[3]))
    verts = np.array(verts); faces = np.array(faces)
    crystals = go.Mesh3d(
        x=verts[:, 0], y=verts[:, 1], z=verts[:, 2],
        i=faces[:, 0], j=faces[:, 1], k=faces[:, 2],
        color="#D9C7A3", opacity=crystal_opacity, flatshading=True,
        name=f"{n * nz} {sc.get('material', '')} crystals".replace("  ", " "),
        showlegend=True, legendgroup="scanner", hoverinfo="skip")

    # The water cylinder as a WIREFRAME, not a surface. Plotly sorts transparent
    # surfaces per-object, so a translucent cylinder in front of a translucent ring
    # paints over it whatever the depth -- the crystals disappear and so do the
    # tracks inside. A wireframe has no such problem and reads just as well.
    Rp, HL = doc["phantom_radius_mm"], doc["phantom_halflength_mm"]
    th = np.linspace(0, 2 * np.pi, 121)
    wx, wy, wz = [], [], []
    for z in (-HL, 0.0, HL):                                   # three hoops
        wx += list(Rp * np.cos(th)) + [None]
        wy += list(Rp * np.sin(th)) + [None]
        wz += [z] * len(th) + [None]
    for a in np.linspace(0, 2 * np.pi, 13)[:-1]:               # a few uprights
        wx += [Rp * math.cos(a)] * 2 + [None]
        wy += [Rp * math.sin(a)] * 2 + [None]
        wz += [-HL, HL, None]
    water = go.Scatter3d(x=wx, y=wy, z=wz, mode="lines",
                         line=dict(color="#4C8BF5", width=2), opacity=0.55,
                         name="water phantom", legendgroup="scanner",
                         hoverinfo="skip")

    src = np.array([doc["source_mm"]])
    sources = go.Scatter3d(
        x=src[:, 0], y=src[:, 1], z=src[:, 2], mode="markers",
        marker=dict(size=6, color="#FFD166", line=dict(color="#8A6D1F", width=1)),
        name="source", legendgroup="scanner",
        hovertemplate="source (%{x:.0f}, %{y:.0f}, %{z:.0f}) mm<extra></extra>")
    return [crystals, water, sources]


# ------------------------------------------------------------------- tracks
def _role(kind, n):
    """For randoms the first two tracks are the pair the sorter recorded; the
    rest are the other photons of the same two decays, drawn to show why."""
    if kind != "random":
        return ""
    return "  recorded" if n < 2 else "  partner"


def _track_trace(t, visible=True, role=""):
    pts = np.array(t["pts"])
    op, width = FATE_STYLE[t["klass"]]
    n = t.get("n_scatter", 1 if t["water"] else 0)
    tag = ("straight out of the water" if not n else
           "scattered once in the water" if n == 1 else
           f"scattered {n} times in the water")
    dep = ("never detected" if t["keV"] is None
           else f"{t['keV']:.0f} keV recorded ({t['klass']})")
    return go.Scatter3d(
        x=pts[:, 0], y=pts[:, 1], z=pts[:, 2], mode="lines",
        line=dict(color=HUE[min(t.get("n_scatter", 1 if t["water"] else 0), 2)],
                  width=width),
        opacity=op, visible=visible,
        name=f"ev{t['event']:04d} p{t['track']}{role}",
        legendgroup=f"p{t['track']}",
        hovertemplate=(f"<b>event {t['event']}, photon {t['track']}</b><br>"
                       f"{tag}<br>left the water with "
                       f"{'-' if t['exit_keV'] is None else format(t['exit_keV'], '.0f')} keV<br>"
                       f"{dep}<extra></extra>"))


def event_items(doc, kind):
    """[(label, [tracks], report)] for one class, in order."""
    out = []
    if kind == "random":
        for c in doc["_randoms"]:
            a = next(t for t in doc["_by_event"][c["e1"]] if t["track"] == c["t1"])
            b = next(t for t in doc["_by_event"][c["e2"]] if t["track"] == c["t2"])
            partners = [t for t in doc["_by_event"][c["e1"]] + doc["_by_event"][c["e2"]]
                        if t is not a and t is not b]
            out.append((f"ev{c['e1']:04d} + ev{c['e2']:04d}",
                        [a, b] + partners, report([a, b], c["dt_ns"])))
    else:
        for ev in events_of(doc, kind):
            g = doc["_by_event"][ev]
            out.append((f"ev{ev:04d}", g, report(g)))
    return out


def report(pair, measured_dt=None):
    """The four numbers, as plain text lines."""
    lines = []
    for t in pair:
        # Kept narrow on purpose: this text sits in the left column of the
        # figure and must not run under the 3D scene.
        n = t.get("n_scatter", 1 if t["water"] else 0)
        ex = ("absorbed in water" if t["exit_keV"] is None
              else f"exit {t['exit_keV']:6.1f} keV"
                   + ("" if not n else f" ({n}x scat)"))
        fly = "" if t["tof_ns"] is None else f"   flight {t['tof_ns']:.3f} ns"
        det = ("never detected" if t["keV"] is None else
               f"dep {t['dep_keV']:6.1f} / rec {t['keV']:6.1f} keV")
        lines.append(f"ev{t['event']:04d} p{t['track']}  {ex}{fly}")
        lines.append(f"           {det}")
    a, b = pair[0], pair[1]
    same_decay = a["event"] == b["event"]
    if same_decay and a["tof_ns"] is not None and b["tof_ns"] is not None:
        dt = a["tof_ns"] - b["tof_ns"]
        lines.append(f"nominal  dt = {dt:+.3f} ns -> {C_MM_PER_NS * dt / 2:+6.1f} mm on LOR")
    if a["t_meas_ns"] is not None and b["t_meas_ns"] is not None:
        dt = measured_dt if measured_dt is not None else a["t_meas_ns"] - b["t_meas_ns"]
        lines.append(f"measured dt = {dt:+.3f} ns -> {C_MM_PER_NS * dt / 2:+6.1f} mm on LOR")
    if same_decay:
        # Say the pair's scatter order explicitly. A "higher order" pair is often
        # 1 + 1 -- each photon deflected once -- rather than one photon deflected
        # twice, and the two tracks then look identical to a single-scatter pair.
        n1 = a.get("n_scatter", 1 if a["water"] else 0)
        n2 = b.get("n_scatter", 1 if b["water"] else 0)
        if n1 + n2:
            order = "single scatter" if n1 + n2 == 1 else "higher order"
            lines.append(f"water scatters: {n1} + {n2} = {n1 + n2}  ({order})")
    else:
        lines.append("different decays -> that position is meaningless")
    return "\n".join(lines)


def _layout(doc, title, report_text=""):
    """Scene on the right, report column on the left, so the two never overlap."""
    sc = doc["scanner"]
    L = sc["ring_radius_mm"] + sc["crystal_mm"][0] + 15
    # The scanner is a short wide barrel. aspectmode="cube" would stretch z to
    # match x and y and spend half the frame on empty space -- but "data" is
    # worse: it takes its proportions from the traces that happen to be VISIBLE,
    # so every event of a class drew the ring at a different shape. Fix the box
    # to the geometry with aspectmode="manual" and it never moves again.
    Lz = max(0.5 * sc["n_rings"] * sc["ring_pitch_mm"] + sc["crystal_mm"][2],
             doc["phantom_halflength_mm"]) + 15
    ax = dict(showbackground=False, showgrid=True,
              gridcolor="#EEF1F5", zeroline=False)
    return dict(
        title=dict(text=f"{title}<br><span style='font-size:11.5px;color:#64748B'>"
                        f"{window_label(doc)} &nbsp;·&nbsp; {summary(doc)}</span>",
                   x=0.01, y=0.975, font=dict(size=17)),
        scene=dict(xaxis=dict(ax, range=[-L, L], title="x [mm]"),
                   yaxis=dict(ax, range=[-L, L], title="y [mm]"),
                   zaxis=dict(ax, range=[-Lz, Lz], title="z [mm]"),
                   aspectmode="manual",
                   aspectratio=dict(x=1.0, y=1.0, z=Lz / L),
                   domain=dict(x=[0.27, 1.0], y=[0.0, 0.92]),
                   camera=dict(eye=dict(x=1.30, y=-1.30, z=0.85))),
        margin=dict(l=8, r=8, t=62, b=8), height=720,
        paper_bgcolor="white",
        legend=dict(y=0.88, x=0.0, xanchor="left", font=dict(size=11),
                    bgcolor="rgba(255,255,255,0)"),
        annotations=[dict(text=report_text.replace("\n", "<br>").replace(" ", "&nbsp;"),
                          xref="paper", yref="paper", x=0.0, y=0.42,
                          showarrow=False, align="left",
                          xanchor="left", yanchor="top",
                          font=dict(family="monospace", size=10.5, color="#334155"))])


def figure(doc, event=None, kind=None):
    """One event, as a standalone figure."""
    if kind is None:
        kind = doc["_class"].get(event, "lost")
    items = event_items(doc, kind)
    label, tracks, rep = (next((it for it in items
                                if it[0].startswith(f"ev{event:04d}")), items[0])
                          if event is not None else items[0])
    fig = go.Figure(scanner_traces(doc)
                    + [_track_trace(t, role=_role(kind, n))
                       for n, t in enumerate(tracks)])
    fig.update_layout(**_layout(doc, f"{KIND_LABEL[kind]} — {label}", rep))
    return fig


def figure_with_slider(doc, kind, max_items=40):
    """Every event of one class, with a slider to step through them."""
    items = event_items(doc, kind)[:max_items]
    if not items:
        return None
    geo = scanner_traces(doc)
    traces, spans = list(geo), []
    for _, tracks, _ in items:
        start = len(traces)
        traces += [_track_trace(t, visible=False, role=_role(kind, n))
                   for n, t in enumerate(tracks)]
        spans.append((start, len(traces)))
    for i in range(*spans[0]):
        traces[i].visible = True

    steps = []
    for (label, _, rep), (lo, hi) in zip(items, spans):
        vis = [True] * len(geo) + [False] * (len(traces) - len(geo))
        for i in range(lo, hi):
            vis[i] = True
        steps.append(dict(method="update", label=label,
                          args=[{"visible": vis},
                                {"title.text": _layout(doc, f"{KIND_LABEL[kind]} — "
                                                            f"{label}")["title"]["text"],
                                 "annotations[0].text": rep.replace("\n", "<br>")}]))
    fig = go.Figure(traces)
    fig.update_layout(**_layout(doc, f"{KIND_LABEL[kind]} — {items[0][0]}",
                                items[0][2]))
    fig.update_layout(sliders=[dict(active=0, steps=steps, x=0.08, len=0.88,
                                    y=0, pad=dict(t=0, b=10),
                                    currentvalue=dict(prefix="event: "))])
    return fig


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run", nargs="?", default=None,
                    help="the run folder gate_sim.py wrote (or its visu/ subfolder, "
                         "or a tracks.json). Default: the newest under runs/")
    ap.add_argument("--out", default=None,
                    help="where to write the page (default: <run>/event_view.html)")
    ap.add_argument("--kind", default="all",
                    choices=["all"] + KINDS)
    ap.add_argument("--events", nargs="*", type=int,
                    help="draw only these event numbers, in one figure")
    ap.add_argument("--max-items", type=int, default=40,
                    help="cap the slider length (the lost class has thousands)")
    ap.add_argument("--no-open", action="store_true",
                    help="just write the file; by default it is opened in "
                         "your browser")
    ap.add_argument("--inline-js", action="store_true",
                    help="bundle plotly.js into the html (~5 MB, works offline); "
                         f"the default loads {PLOTLY_JS.rsplit('/', 1)[1]} from the "
                         "CDN and keeps the file under 1 MB")
    a = ap.parse_args()
    doc = load(a.run)
    out = a.out or os.path.join(doc["_run"], "event_view.html")
    js = True if a.inline_js else PLOTLY_JS

    if a.events:
        parts = [figure(doc, ev).to_html(full_html=False,
                                         include_plotlyjs=(js if n == 0 else False))
                 for n, ev in enumerate(a.events)]
        nav = ""
    else:
        kinds = kinds_present(doc) if a.kind == "all" else [a.kind]
        parts, links = [], []
        for k in kinds:
            fig = figure_with_slider(doc, k, a.max_items)
            if fig is None:
                continue
            div = f"fig_{k}"
            parts.append(fig.to_html(full_html=False, div_id=div,
                                     include_plotlyjs=(js if not parts else False)))
            links.append((k, div))
        nav = " &nbsp;·&nbsp; ".join(f'<a href="#{d}">{KIND_LABEL[k]}</a>'
                                     for k, d in links)
        missing = [KIND_LABEL[k].lower()
                   for k in ("true", "single_scatter", "multi_scatter", "random")
                   if not events_of(doc, k)]
        if missing:
            nav += ("<br><span style='font-size:12px'>none of: "
                    + ", ".join(missing) + " in this run</span>")

    html = ("<html><head><meta charset='utf-8'><title>Gate 10 event viewer</title>"
            "</head><body style='font-family:sans-serif;margin:18px'>"
            "<h2 style='margin:0 0 2px'>Gate 10 event viewer</h2>"
            f"<p style='color:#64748B;margin:0 0 4px'>{window_label(doc)} "
            f"&nbsp;·&nbsp; {summary(doc)}</p>"
            f"<p style='color:#64748B;margin:0 0 14px'>{nav}</p>"
            + "".join(parts) + "</body></html>")
    with open(out, "w") as fh:
        fh.write(html)

    st = doc["stats"]
    print(f"{out}  ({os.path.getsize(out)/1e6:.2f} MB)")
    print(f"  drawn: {st['n_coincidences']} coincidences = {st['n_true']} true + "
          f"{st.get('n_single_scatter', 0)} single-scatter + "
          f"{st.get('n_multi_scatter', 0)} higher-order + "
          f"{st['n_random']} random; "
          f"{len(events_of(doc, 'lost'))} events not a coincidence")
    r = doc.get("_rates")
    if r:
        print(f"  rates from the {r['n_decays']}-decay run: "
              f"{100*r['frac_true']:.1f} % true, "
              f"{100*r['frac_scattered']:.1f} % scattered, "
              f"{100*r['frac_random']:.1f} % random")
    else:
        print("  no stats run in this folder -- see gate_sim.py --mode stats")

    if not a.no_open:
        url = "file://" + os.path.abspath(out)
        # webbrowser returns False on a headless box instead of raising, and
        # raises if it finds no browser at all -- neither is worth a traceback
        # over, the file is already written.
        try:
            opened = webbrowser.open(url)
        except Exception:
            opened = False
        print("  opened in your browser" if opened else f"  open it yourself: {url}")


if __name__ == "__main__":
    main()
