"""Runtime measurement and wall-time forecasting utilities."""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
import os
from shutil import get_terminal_size
import sys
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


def format_duration(seconds):
    minutes, seconds = divmod(max(0, round(seconds)), 60)
    hours, minutes = divmod(minutes, 60)
    return f'{hours}:{minutes:02d}:{seconds:02d}' if hours else f'{minutes:02d}:{seconds:02d}'


class BinProgress:
    """Script 03 style terminal heartbeat for checkpointed initial scans."""

    def __init__(self, total, saved, workers):
        self.enabled = sys.stderr.isatty()
        self.total = total
        self.saved = saved
        self.workers = workers
        self.started = time.perf_counter()
        self.last_bin = self.started
        self.done = 0
        self.fit_seconds = 0.
        self.width = 0
        self.spinner_index = 0

    def checkpoint(self, fit_seconds):
        self.done += 1
        self.fit_seconds += fit_seconds
        self.last_bin = time.perf_counter()

    def eta(self):
        remaining = self.total - self.saved - self.done
        if not self.done:
            return None
        return self.fit_seconds / self.done * remaining / min(self.workers, max(1, remaining))

    def refresh(self):
        if not self.enabled:
            return
        now = time.perf_counter()
        done = self.saved + self.done
        fraction = done / self.total
        spinner = '|/-\\'[self.spinner_index % 4]
        self.spinner_index += 1
        elapsed = format_duration(now - self.started)
        estimate = self.eta()
        eta = 'ETA warming up' if estimate is None else f'ETA {format_duration(estimate)}'
        last = (f'last new bin {format_duration(now - self.last_bin)} ago'
                if self.done else 'last new bin: none yet')

        def line(bar_width):
            filled = int(fraction * bar_width)
            bar = '#' * filled + '-' * (bar_width - filled)
            return f'{spinner} RUNNING [{bar}] {done}/{self.total} bins ({fraction:5.1%})'

        columns = max(1, get_terminal_size().columns - 1)
        status = f'{line(24)} | workers={self.workers} | elapsed {elapsed} | {last} | {eta}'
        if len(status) > columns:
            status = f'{line(12)} | {elapsed} | {eta}'
        status = status[:columns]
        self.width = min(columns, max(self.width, len(status)))
        sys.stderr.write('\r' + status.ljust(self.width))
        sys.stderr.flush()

    def clear(self):
        if self.width:
            sys.stderr.write('\r' + ' ' * self.width + '\r')
            sys.stderr.flush()
            self.width = 0


@contextmanager
def output_lock(directory, lock_name):
    """Release the run lock even if the process is interrupted or killed."""
    with (directory / lock_name).open('a+b') as stream:
        stream.seek(0)
        try:
            if os.name == 'nt':
                import msvcrt
                msvcrt.locking(stream.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            raise RuntimeError(f'Another process holds the run lock: {directory}') from exc
        yield
