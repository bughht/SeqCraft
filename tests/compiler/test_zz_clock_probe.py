"""TEMPORARY probe: effective granularity of the CPU clock on each CI platform.  Not for merge."""
from __future__ import annotations

import platform
import time


def test_clock_granularity_probe() -> None:
    def effective(clock):
        """Smallest non-zero increment the clock actually reports."""
        deltas = []
        for _ in range(200):
            a = clock()
            while True:
                b = clock()
                if b != a:
                    deltas.append(b - a)
                    break
        return min(deltas)

    lines = [f'platform={platform.system()} {platform.machine()} python={platform.python_version()}']
    for name in ('perf_counter', 'process_time', 'thread_time'):
        clock = getattr(time, name, None)
        if clock is None:
            lines.append(f'{name}: ABSENT')
            continue
        info = time.get_clock_info(name)
        lines.append(f'{name}: reported={info.resolution * 1e9:.1f} ns  '
                     f'effective={effective(clock) * 1e9:.0f} ns')
    # Also: can process_time resolve a ~1.3 ms workload at all?
    def busy() -> int:
        return sum(i * i for i in range(60000))

    for name in ('perf_counter', 'process_time'):
        clock = getattr(time, name)
        samples = []
        for _ in range(8):
            t0 = clock()
            busy()
            samples.append((clock() - t0) * 1e3)
        lines.append(f'{name}: 8 samples of a ~ms workload (ms) = '
                     + ', '.join(f'{s:.3f}' for s in samples))
    # Only Windows is the unknown, and `fail-fast: true` would cancel it if a Linux lane went
    # red first, so only Windows reports.
    if platform.system() == 'Windows':
        raise AssertionError('CLOCK PROBE\n' + '\n'.join(lines))
