import time
import pytest
from pbg_cpm_studies.influenza import build, sheet


@pytest.mark.perf
def test_mcs_throughput_1mm():
    w = build.world_from_spec(sheet.build_sheet_spec(1.0))  # 500x500, 10k cells
    w.step_parallel(5, 64)          # warm up
    t0 = time.perf_counter()
    steps = 50
    w.step_parallel(steps, 64)
    rate = steps / (time.perf_counter() - t0)
    print(f"\n1mm^2 throughput: {rate:.2f} MCS/s")
    # Budget floor: the 2-week capstone is ~20160 steps. The minimum needed
    # for feasibility is >=2 MCS/s (one replica in <~2.8 h). This floor is a
    # HARDWARE-PORTABLE feasibility + gross-regression gate, deliberately not a
    # machine-specific benchmark: measured baselines vary widely by hardware
    # (~410-422 MCS/s on a dev laptop, ~124 MCS/s on GitHub CI runners), so a
    # tight "50% of the fastest machine" floor false-fails on slower CI. We set
    # 20 MCS/s -- a 10x margin over the 2.0 feasibility minimum that still trips
    # on a catastrophic engine regression, and holds on the slowest CI seen.
    # The actual measured rate is printed (and recorded in the study) for
    # tracking; it is informational, not asserted.
    assert rate >= 20.0, f"throughput {rate:.2f} MCS/s below feasibility floor"
