#!/usr/bin/env bash
# One complete run: the two acquisitions and the viewer built from them.
#     mamba activate pet_gate && bash run_visu.sh        (~1 min)
#
# Two acquisitions, because they want opposite things:
#
#   visu   2000 decays WITH the step-level actor that records every trajectory.
#          That actor is the expensive one -- 3 MB here, well over a GB at a
#          million decays -- and it is the only reason this run is kept short.
#   stats  1e6 decays WITHOUT it. Nothing about counting coincidences needs
#          trajectories; scatter comes from the phantom interaction history.
#          This is where the true / scattered / random rates come from, and the
#          viewer quotes them rather than the handful of coincidences it draws.
#
# Both go into one timestamped folder under runs/, which is never reused, so
# every run is kept and can be reopened later.
set -euo pipefail
cd "$(dirname "$0")"

# e.g. GEOM="--bore 680 --n-tangential 80 --elow 470" bash run_visu.sh
GEOM=${GEOM:-}
ACTIVITY=${ACTIVITY:-1}                 # MBq, both acquisitions
RUN=${RUN:-runs/$(date +%Y%m%d-%H%M%S)}

python gate_sim.py --mode visu  --activity "$ACTIVITY" --duration 0.002 --seed 7 \
    --outdir "$RUN/visu" $GEOM

python gate_sim.py --mode stats --activity "$ACTIVITY" --duration 1.0 --seed 11 \
    --outdir "$RUN/stats" $GEOM

python event_viewer.py "$RUN"

echo
echo "run kept in $RUN"
