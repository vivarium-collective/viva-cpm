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
    # is >=2 MCS/s (one replica in <~2.8 h). Measured baseline on this
    # machine (500x500 lattice, 10k cells, 64-thread block): ~410-422 MCS/s
    # across 3 runs -- far above the 2.0 minimum. Per the perf-gate plan,
    # the floor is raised to ~50% of the measured baseline (~200 MCS/s) so
    # the test is a real regression gate rather than a no-op at 2.0.
    assert rate >= 200.0, f"throughput {rate:.2f} MCS/s below budget floor"
