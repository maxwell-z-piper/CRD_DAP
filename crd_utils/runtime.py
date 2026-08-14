"""Runtime measurement and wall-time forecasting utilities."""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
import time


@dataclass(frozen=True)
class TimingResult:
    label: str
    elapsed_seconds: float


@contextmanager
def timed_block(label: str, sink: list[TimingResult] | None = None):
    """Context manager for lightweight timing instrumentation."""
    start = time.perf_counter()
    try:
        yield
    finally:
        result = TimingResult(label=label, elapsed_seconds=time.perf_counter() - start)
        if sink is not None:
            sink.append(result)


def estimate_parallel_walltime(serial_seconds: float, n_tasks: int, n_workers: int, efficiency: float = 0.85) -> float:
    """Simple planning estimate for embarrassingly parallel independent tasks.

    This is intentionally only a forecast. Final Script-7 planning should use
    measured one-realization timings and, when possible, empirical scaling tests
    at multiple worker counts.
    """
    if serial_seconds < 0 or n_tasks < 0 or n_workers < 1:
        raise ValueError("Invalid timing/task/worker input.")
    if not 0 < efficiency <= 1:
        raise ValueError("efficiency must lie in (0, 1].")
    return float(serial_seconds * n_tasks / (n_workers * efficiency))
