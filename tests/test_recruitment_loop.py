"""tests/test_recruitment_loop.py -- the recruitment model-building loop
driver (Phase 2, Task 7) must: (1) flag a one-sided draft contract and pass a
two-sided locked one under test_audit; (2) climb the real ladder to
adaptive_receptor and reach DONE; (3) terminate honestly as GIVE_UP when
adaptive_receptor is withheld from the library. These run the real Rust CPM
engine (mechanisms.simulate_condition) -- slow relative to unit tests, but
each full loop completes in seconds, not minutes.
"""
import importlib.util
import pathlib

_s = importlib.util.spec_from_file_location("brl", pathlib.Path(__file__).parent.parent / "scripts" / "build_recruitment_loop.py")
brl = importlib.util.module_from_spec(_s)
_s.loader.exec_module(brl)


def test_audit_two_rounds():
    import viva_superpowers.test_audit as ta
    assert ta.audit_gate(ta.build_audit_report(brl.draft_spec())) in ("warn", "fail")
    assert ta.audit_gate(ta.build_audit_report(brl.locked_spec())) == "pass"


def test_loop_reaches_adaptive():
    outcome, iters, final = brl.run_loop(brl.LIBRARY, seeds=(17, 29), steps=40)
    assert outcome == "DONE" and final == "adaptive_receptor" and len(iters) >= 2


def test_giveup_companion_open():
    outcome, iters, final = brl.run_loop(brl.LIBRARY, seeds=(17, 29), steps=40, withhold="adaptive_receptor")
    assert outcome == "GIVE_UP" and final != "adaptive_receptor"
