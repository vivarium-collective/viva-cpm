"""Model-building loop driver: climb the recruitment mechanism ladder from the
LOCKED contract (workspace/studies/recruitment-adaptive/study.yaml) and capture
a real, emergent trajectory (Phase 2, Task 7).

The trajectory is NOT scripted. `navigate.next_mechanism` reads the graded
behavior tests each iteration, takes the failing HARD test with the
most-negative margin, and installs the mechanism `mechanisms._FIXES` maps it
to. Every observation is `mechanisms.simulate_condition` running the real
Rust CPM engine (`pbg_cpm_studies.composites.chemotaxis.recruitment`); every
verdict is computed by `viva_superpowers.test_contract`.

Also captures: a two-round sufficiency audit (a one-sided `recruits_high`
draft flagged as an unverifiable, gameable discriminator -> revised to the
two-sided banded contract that mirrors study.yaml -> re-audit clean), a
GIVE_UP companion trace (adaptive_receptor withheld from the library, so the
loop climbs to hill_occupancy and terminates HONESTLY instead of faking a
pass), and perturbation-stability of the final pass under a different seed
set. Writes trajectory.json + a real .pbg/loop/recruitment-adaptive.json.
"""
from __future__ import annotations

import json
import os

from pbg_cpm_studies.model_building import mechanisms as M
from pbg_cpm_studies.model_building import navigate
from viva_superpowers import loop_state as ls, test_audit, test_contract as tc

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir))
OUT = os.path.join(ROOT, "workspace", "investigations", "recruitment-model-building")

STUDY = "recruitment-adaptive"

# The installable rungs (Task 7 drives THIS library; navigate never sees
# anything outside it).
LIBRARY = M.MECHANISMS  # ('static_lambda', 'hill_occupancy', 'adaptive_receptor')

QUESTION = (
    "Across a range of background cue levels, can a receptor mechanism keep "
    "recruiting responders to a rising source cue, where a fixed half-occupancy "
    "mechanism saturates and collapses?"
)

CITES = {
    "recruitment-baseline": "measured recruitment_index, faithful mechanism-ladder calibration (3-12 seeds)",
    "nasser2009": "Nasser & Simoni 2009 -- receptor-gated response abolished when downstream signaling is blocked",
    "barkai1997": "Barkai & Leibler 1997, Nature -- robust adaptation via a re-centering set-point",
    "tu2008": "Tu, Shimizu & Berg 2008, PNAS -- fold-change detection in bacterial chemotaxis",
}

# KEEP IT FAST knobs. steps=30 (as first tried) gave noisy, unconverged
# recruitment_index values that did not reliably clear the locked bands (see
# task-7-report.md); steps=40 (still far short of a sweep) converges cleanly
# and the whole driver (real run + GIVE_UP companion + perturbation) still
# completes in well under 6 minutes, so steps was raised rather than the bands
# weakened, per the task's own escape hatch.
DEFAULT_SEEDS = (17, 29)
DEFAULT_STEPS = 40
PERTURB_SEEDS = (7, 101)   # a different seed pair for the final stability check

# id / label / condition / locked band (low, high) OR one-sided (op, value) /
# classification (behavior_tests shape) / control / cites / provenance prose.
_TEST_DEFS = [
    dict(id="recruits_low", label="Recruits responders at low background",
         condition="low_bg", low=0.4, high=1.0, classification="primary",
         cites=["recruitment-baseline"],
         prov="recruitment_index in [0.4, 1.0] at background=0.0 -- baseline "
              "competence, not the discriminator (measured adaptive 0.82, hill 0.58; "
              "%s)" % CITES["recruitment-baseline"]),
    dict(id="receptor_gating", label="Blocking the response abolishes recruitment "
                                      "(receptor-gated mechanisms only)",
         condition="high_bg_blocked", op="<=", value=0.15,
         classification="diagnostic", control="negative", cites=["nasser2009"],
         prov="recruitment_index <= 0.15 at background=400.0, blocked=True -- "
              "receptor-mediated recruitment is abolished; the non-receptor static "
              "mechanism ignores the block (measured 0.72) and FAILS this test by "
              "design (%s)" % CITES["nasser2009"]),
    dict(id="recruits_high", label="Recruits responders at high background "
                                    "(the discriminator)",
         condition="high_bg", low=0.35, high=1.0, classification="primary",
         cites=["barkai1997", "tu2008"],
         prov="recruitment_index in [0.35, 1.0] at background=400.0 -- fold-change "
              "detection: an adaptive kd rescues recruitment (measured ~0.56, "
              "12-seed mean) where a fixed-kd hill mechanism collapses to 0.00 "
              "(%s; %s)" % (CITES["barkai1997"], CITES["tu2008"])),
]


def _expected(t):
    return tc.value(t["value"], op=t["op"]) if "op" in t else tc.band(t["low"], t["high"])


TESTS = [dict(t, expected=_expected(t)) for t in _TEST_DEFS]
BYID = {t["id"]: t for t in TESTS}


def _pass_if(t, *, oneside_high):
    if t["id"] == "recruits_high" and oneside_high:
        # the DRAFT form: a bare one-sided floor, no upper band -- unverifiable,
        # gameable tightness that test_audit.one_sided_loose_primary flags.
        return {"op": ">=", "value": t["low"], "provenance": t["prov"]}
    if "op" in t:
        return {"op": t["op"], "value": t["value"], "provenance": t["prov"]}
    return {"op": "in_range", "low": t["low"], "high": t["high"], "provenance": t["prov"]}


def _build_spec(*, oneside_high):
    bt = []
    for t in TESTS:
        entry = {"name": t["id"], "classification": t["classification"],
                 "description": t["label"],
                 "measure": {"kind": "recruitment_index", "condition": t["condition"]},
                 "pass_if": _pass_if(t, oneside_high=oneside_high), "cites": t["cites"]}
        if t.get("control"):
            entry["control"] = t["control"]
        bt.append(entry)
    return {"name": STUDY, "question": QUESTION, "behavior_tests": bt}


def draft_spec():
    """One-sided draft: `recruits_high` is `>= 0.35` with no upper band -- the
    audit flags it as an unverifiable primary threshold (discrimination drift),
    so `audit_gate(build_audit_report(draft_spec()))` returns "warn" or "fail"."""
    return _build_spec(oneside_high=True)


def locked_spec():
    """Two-sided, locked contract -- mirrors study.yaml's behavior_tests bands.
    `audit_gate(build_audit_report(locked_spec()))` returns "pass"."""
    return _build_spec(oneside_high=False)


def grade(mechanism, *, seeds=DEFAULT_SEEDS, steps=DEFAULT_STEPS):
    """Run `mechanism` against every locked test's condition (the real CPM
    engine, via `mechanisms.simulate_condition`) and grade each observation
    against its locked Expected band/value. Returns the list of v2 axes."""
    axes = []
    for t in TESTS:
        obs = M.simulate_condition(mechanism, t["condition"], seeds=seeds, steps=steps)
        axis = tc.check(t["id"], t["label"], obs, t["expected"], severity="hard",
                         detail={"observed": round(obs, 4), "condition": t["condition"],
                                 "provenance": t["prov"]})
        axes.append(axis)
    return axes


def run_loop(library, *, seeds=DEFAULT_SEEDS, steps=DEFAULT_STEPS, withhold=None, st_box=None):
    """Climb `library` from static_lambda, driven by the graded locked contract.

    Each iteration: grade the active mechanism -> navigate.next_mechanism picks
    the next rung. None -> DONE. The target is in (library - withhold) -> install
    it (`st_box[0]`, if given, gets `ls.record_iteration`d). Otherwise the target
    is unavailable -> GIVE_UP. Capped at len(mechanisms.MECHANISMS) + 1 iterations
    so a policy bug can never spin forever. Returns
    (outcome, iterations, final_mechanism).
    """
    avail = tuple(m for m in library if m != withhold)
    active = "static_lambda"
    iterations = []
    cap = len(M.MECHANISMS) + 1
    prev_margins = {}
    for step in range(cap):
        axes = grade(active, seeds=seeds, steps=steps)
        margins = {a["id"]: a.get("margin") for a in axes}
        target = navigate.next_mechanism(active, axes)
        record = {
            "iteration": step, "active": active,
            "verdicts": [{"id": a["id"], "label": a["label"], "verdict": a["verdict"],
                          "observed": a["detail"]["observed"], "margin": a.get("margin")}
                         for a in axes],
            "target": target,
        }
        iterations.append(record)
        if target is None:
            return "DONE", iterations, active
        if target not in avail:
            record["blocked"] = True
            return "GIVE_UP", iterations, active
        if st_box is not None:
            gate = "pass" if all(a["verdict"] == "within_tol" for a in axes) else "fail"
            deltas = {k: (None if margins[k] is None or prev_margins.get(k) is None
                          else round(margins[k] - prev_margins[k], 4)) for k in margins}
            st_box[0] = ls.record_iteration(st_box[0], edit=f"climb {active} -> {target}",
                                             target="mechanism", margin_deltas=deltas, gate=gate)
        prev_margins = margins
        active = target
    return "GIVE_UP", iterations, active


def _fmt_iter(it):
    tag = f"install {it['target']}" if it["target"] else "DONE"
    if it.get("blocked"):
        tag = f"BLOCKED (no available mechanism provides `{it['target']}`)"
    verdicts = " ".join(f"{v['id']}={v['verdict']}" for v in it["verdicts"])
    return f"  iter {it['iteration']} active={it['active']:<17} {verdicts}  -> {tag}"


def main():
    os.makedirs(OUT, exist_ok=True)

    # AUTHOR -> AUDIT round 1 (one-sided recruits_high) -> revise -> AUDIT round 2
    # (two-sided, locked) -> SELECT -> LOCK.
    draft, locked = draft_spec(), locked_spec()
    rep1 = test_audit.build_audit_report(draft)
    gate1 = test_audit.audit_gate(rep1)
    rep2 = test_audit.build_audit_report(locked)
    gate2 = test_audit.audit_gate(rep2)
    assert gate1 in ("warn", "fail"), f"draft spec audit expected warn/fail, got {gate1!r}"
    assert gate2 == "pass", f"locked spec audit expected pass, got {gate2!r}"

    st = ls.create(ROOT, STUDY, QUESTION, max_iterations=len(M.MECHANISMS) + 1)
    st = ls.advance(st, "AUDIT", audit={"gate": gate2, "reopened_from": gate1})
    st = ls.advance(st, "SELECT", sourcing={
        "decision": "compose", "modules": ["chemotaxis_receptor"],
        "rationale": "receptor occupancy is catalogued; adaptation is composed on top via "
                     "occupancy-space chemotaxis with an adaptive kd", "gate": "pass"})
    st = ls.lock_tests(st, locked["behavior_tests"])
    locked_hash = st.get("locked_tests_hash")

    # the real, emergent climb
    st_box = [st]
    outcome, iterations, final_mechanism = run_loop(LIBRARY, st_box=st_box)
    st = st_box[0]
    violations = ls.validate(st, locked["behavior_tests"])
    st = ls.advance(st, outcome, last_verdict={"gate": "pass" if outcome == "DONE" else "fail"})
    ls.save(ROOT, STUDY, st)

    # GIVE_UP companion: withhold adaptive_receptor -- the loop must terminate
    # honestly rather than fake a pass with a mechanism the library can't supply.
    gu_outcome, gu_iterations, gu_final = run_loop(LIBRARY, withhold="adaptive_receptor")

    # perturbation: re-grade the final mechanism at a different seed set
    perturb_axes = grade(final_mechanism, seeds=PERTURB_SEEDS, steps=DEFAULT_STEPS)
    perturb_pass = all(a["verdict"] == "within_tol" for a in perturb_axes)

    final_numbers = {v["id"]: v["observed"] for v in iterations[-1]["verdicts"]}
    per_mechanism_numbers = {it["active"]: {v["id"]: v["observed"] for v in it["verdicts"]}
                             for it in iterations}

    traj = {
        "schema": "model_build_trajectory/v2",
        "study": STUDY,
        "contract": {"question": QUESTION,
                     "observables": ["recruitment_index (fraction of responders within radius)"],
                     "success": "all 3 hard behavior tests pass for the installed mechanism"},
        "draft": {"description": "One-sided draft acceptance contract -- recruits_high is a bare "
                                  ">= 0.35 floor with no upper band, an unverifiable, gameable "
                                  "discriminator the audit flags before the loop is driven.",
                  "spec": draft},
        "select": {"decision": "compose", "modules": ["chemotaxis_receptor"],
                   "rationale": "receptor occupancy is catalogued; adaptation is composed on top via "
                                "occupancy-space chemotaxis with an adaptive kd",
                   "library": list(LIBRARY)},
        "audit": {"round1_gate": gate1, "round2_gate": gate2,
                  "revision": "recruits_high's one-sided >=0.35 threshold flagged (discrimination "
                              "drift, unverifiable tightness) -> revised to the two-sided band "
                              "[0.35, 1.0] mirroring study.yaml -> re-audit clean"},
        "lock": {"tests_hash": locked_hash, "n_tests": len(locked["behavior_tests"])},
        "tests": [{"id": t["id"], "label": t["label"], "condition": t["condition"],
                   "expected": ({"op": t["op"], "value": t["value"]} if "op" in t
                               else {"low": t["low"], "high": t["high"]}),
                   "severity": "hard", "provenance": t["prov"]} for t in TESTS],
        "iterations": iterations,
        "result": {"state": outcome, "final_mechanism": final_mechanism,
                   "n_iterations": len(iterations), "final_numbers": final_numbers,
                   "violations": violations},
        "giveup_companion": {"outcome": gu_outcome, "final_mechanism": gu_final,
                             "n_iterations": len(gu_iterations),
                             "note": "same loop with adaptive_receptor withheld from the library -- it "
                                     "climbs to hill_occupancy (clears recruits_low and receptor_gating) "
                                     "then stalls on recruits_high (fixed-kd occupancy collapses at high "
                                     "background) with no further mechanism available, and terminates "
                                     "HONESTLY (GIVE_UP) instead of faking a pass."},
        "perturbation": {"seeds": list(PERTURB_SEEDS), "mechanism": final_mechanism,
                         "pass": perturb_pass,
                         "axes": [{"id": a["id"], "verdict": a["verdict"], "observed": a["detail"]["observed"]}
                                  for a in perturb_axes]},
        "per_mechanism_numbers": per_mechanism_numbers,
    }

    path = os.path.join(OUT, "trajectory.json")
    with open(path, "w") as fh:
        json.dump(traj, fh, indent=2)

    print(f"AUDIT: draft={gate1} -> revised -> locked={gate2}   "
          f"LOCK: {len(locked['behavior_tests'])} tests ({locked_hash[:18]}...)")
    print("REAL CLIMB:")
    for it in iterations:
        print(_fmt_iter(it))
    print(f"RESULT: {outcome} at {final_mechanism} in {len(iterations)} iterations; violations={violations}")
    print(f"PERTURBATION (seeds={PERTURB_SEEDS}): {final_mechanism} still passes = {perturb_pass}")
    print("GIVE_UP COMPANION (adaptive_receptor withheld):")
    for it in gu_iterations:
        print(_fmt_iter(it))
    print(f"GIVE_UP companion result: {gu_outcome} at {gu_final} in {len(gu_iterations)} iterations")
    print(f"trajectory -> {os.path.relpath(path, ROOT)}")


if __name__ == "__main__":
    main()
