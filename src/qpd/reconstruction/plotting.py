"""Quick-look plots of a readout trace with reconstructed tunnelling times.

The point of these helpers is to put the *input* and the *output* on the same
time axis: the raw I and Q voltages that went in, and the flip times that came
out, so a reconstruction can be eyeballed rather than taken on trust.

Two things about this data make a naive plot useless, and both are handled here.

* **Volume.** A 5 s trace at 100 kHz is 500,000 samples per quadrature. Drawn
  whole it is a solid band. :func:`plot_trace_with_flips` therefore decimates to
  a bounded number of points, and takes a ``window`` so you can zoom to the few
  milliseconds around an event.
* **Modulation.** When ``n_g`` is swept, I and Q are dominated by the ramp --
  at the reference scenario a 5.5 kHz oscillation -- and individual parity
  flips are invisible underneath it. That is expected: in the swept case the
  parity lives in the *spread* of the trace, not its mean (see
  ``docs/reconstruction.md`` §7.2), so pass ``projected=True`` to plot the
  learned discrimination axis with the branch means overlaid, which is where
  the flips actually become visible. At a fixed bias the raw I/Q is already
  informative, and ``smooth_hz`` makes the telegraph plain.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

__all__ = ["plot_trace_with_flips", "plot_iq_plane"]

_STYLE = Path(__file__).resolve().parents[1] / "theory" / "qpd.mplstyle"


def _decimate(t, *arrays, max_points):
    """Uniformly thin a set of equal-length arrays for plotting."""
    n = t.size
    if max_points is None or n <= max_points:
        return (t,) + arrays
    step = int(np.ceil(n / max_points))
    return (t[::step],) + tuple(a[::step] for a in arrays)


def _lowpass(x, sample_rate, cutoff_hz):
    """Zero-phase Butterworth low-pass, for making a telegraph legible."""
    from scipy.signal import butter, filtfilt

    b, a = butter(N=4, Wn=cutoff_hz, btype="low", fs=sample_rate)
    return filtfilt(b, a, x)


def plot_trace_with_flips(
    iq,
    sample_rate,
    flip_times=None,
    truth_times=None,
    t0=0.0,
    window=None,
    projected=False,
    emission=None,
    smooth_hz=None,
    baseline=False,
    max_points=20000,
    confidence=None,
    figsize=(9.5, 4.6),
    title=None,
    axes=None,
):
    """Plot a readout trace with reconstructed flip times overlaid.

    Parameters
    ----------
    iq : complex array
        The readout trace, ``I + 1j*Q``, exactly as handed to the
        reconstruction (e.g. ``SimResult.iq``).
    sample_rate : float
        Sampling rate [Hz].
    flip_times : array, optional
        Reconstructed tunnelling times [s] -- ``rec.flip_times``. Drawn as
        dashed red verticals.
    truth_times : array, optional
        Known tunnelling times [s] -- ``SimResult.flip_times`` in simulation.
        Drawn as solid green verticals, so detected-vs-truth is visible at a
        glance. Omit for measured data, where no truth exists.
    t0 : float
        Time of the first sample [s]; the time axis and both event lists are
        interpreted on this axis.
    window : (float, float), optional
        ``(t_start, t_end)`` in seconds. Strongly recommended: a whole trace is
        far too dense to read. Defaults to the first 4 ms, or the full trace if
        it is shorter.
    projected : bool
        Plot the trace projected onto the learned discrimination axis instead
        of raw I and Q. This is the view in which flips are visible for a
        **swept** ``n_g``; requires ``emission``.
    emission : EmissionModel or StaticBlobModel, optional
        Supplies the projection axis and the two branch means, which are
        overlaid when ``projected`` is set. Pass ``rec.emission`` (swept) or
        ``rec.model`` (fixed bias).
    smooth_hz : float, optional
        Overlay a zero-phase low-pass at this cutoff. Useful at a fixed bias,
        where it makes the telegraph obvious; misleading on a swept trace,
        where it would filter away the ramp modulation that carries the signal.
    baseline : bool
        Plot each quadrature as ``x / mean(x) - 1`` (the notebook's convention)
        rather than in raw units.
    max_points : int or None
        Decimate to at most this many plotted points per curve.
    confidence : array, optional
        Per-flip confidence in ``[0, 1]`` (``rec.confidence``), matched to
        ``flip_times``; scales the opacity of each detected line so
        low-confidence flips are visibly tentative.
    figsize, title, axes
        Usual matplotlib escapes. ``axes`` must hold two axes when
        ``projected`` is False, one when it is True.

    Returns
    -------
    (fig, axes)

    Examples
    --------
    >>> rec = reconstruct_parity_flips_ramped(result.iq, 1e5)
    >>> plot_trace_with_flips(result.iq, 1e5, rec.flip_times,
    ...                       truth_times=result.flip_times,
    ...                       window=(0.020, 0.024))
    """
    iq = np.asarray(iq)
    if iq.ndim != 1:
        raise ValueError("iq must be a 1-D complex trace")
    if projected and emission is None:
        raise ValueError("projected=True requires an emission/blob model")

    n = iq.size
    dt = 1.0 / float(sample_rate)
    t = t0 + np.arange(n) * dt

    if window is None:
        window = (t0, min(t0 + 4e-3, t0 + n * dt))
    lo, hi = float(window[0]), float(window[1])
    sel = (t >= lo) & (t <= hi)
    if not np.any(sel):
        raise ValueError(f"window {window} selects no samples "
                         f"(trace spans {t[0]:.6g}..{t[-1]:.6g} s)")

    def _events(times):
        if times is None:
            return np.empty(0)
        times = np.asarray(times, dtype=float)
        return times[(times >= lo) & (times <= hi)]

    det, tru = _events(flip_times), _events(truth_times)
    conf_in = None
    if confidence is not None and flip_times is not None:
        c = np.asarray(confidence, dtype=float)
        ft = np.asarray(flip_times, dtype=float)
        if c.size == ft.size:
            conf_in = c[(ft >= lo) & (ft <= hi)]

    def _mark(ax):
        for tt in tru:
            ax.axvline(tt * 1e3, color="tab:green", lw=1.0, alpha=0.85, zorder=1)
        for j, tt in enumerate(det):
            a = 0.9 if conf_in is None else float(np.clip(conf_in[j], 0.15, 1.0))
            ax.axvline(tt * 1e3, color="tab:red", lw=1.0, ls="--", alpha=a,
                       zorder=2)

    with plt.style.context(_STYLE):
        if projected:
            x = emission.project(iq)
            tt, xs = _decimate(t[sel], x[sel], max_points=max_points)
            if axes is None:
                fig, ax = plt.subplots(figsize=(figsize[0], figsize[1] * 0.62))
            else:
                ax = np.atleast_1d(axes)[0]
                fig = ax.figure
            ax.plot(tt * 1e3, xs, lw=0.4, color="0.6", label="projected trace")
            mu = getattr(emission, "branch_means", None)
            if callable(mu):
                a_, b_ = mu()
                a_ = np.asarray(a_)
                b_ = np.asarray(b_)
                if a_.size == n:
                    _, ma, mb = _decimate(t[sel], a_[sel], b_[sel],
                                          max_points=max_points)
                else:  # static model: two scalars
                    ma = np.full(tt.size, float(a_))
                    mb = np.full(tt.size, float(b_))
                ax.plot(tt * 1e3, ma, lw=1.0, color="tab:blue",
                        label=r"$\mu_A(t)$")
                ax.plot(tt * 1e3, mb, lw=1.0, color="tab:orange",
                        label=r"$\mu_B(t)$")
            elif hasattr(emission, "mu_a"):
                ax.axhline(emission.mu_a, lw=1.0, color="tab:blue",
                           label=r"$\mu_A$")
                ax.axhline(emission.mu_b, lw=1.0, color="tab:orange",
                           label=r"$\mu_B$")
            ax.set_ylabel("projection [a.u.]")
            ax.set_xlabel("Time [ms]")
            ax.legend(fontsize=6, ncol=3, loc="upper left")
            _mark(ax)
            out_axes = np.array([ax])
        else:
            i_raw, q_raw = iq.real, iq.imag
            if baseline:
                i_raw = i_raw / np.mean(i_raw) - 1.0
                q_raw = q_raw / np.mean(q_raw) - 1.0
            if axes is None:
                fig, out_axes = plt.subplots(2, 1, figsize=figsize, sharex=True)
            else:
                out_axes = np.atleast_1d(axes)
                fig = out_axes[0].figure
            unit = "baselined [a.u.]" if baseline else "[a.u.]"
            for ax, raw, name, colour in ((out_axes[0], i_raw, "I", "tab:blue"),
                                          (out_axes[1], q_raw, "Q", "tab:orange")):
                tt, xs = _decimate(t[sel], raw[sel], max_points=max_points)
                ax.plot(tt * 1e3, xs, lw=0.4,
                        color="0.6" if smooth_hz else colour,
                        label=f"{name}(t)")
                if smooth_hz:
                    sm = _lowpass(raw, sample_rate, smooth_hz)
                    _, sms = _decimate(t[sel], sm[sel], max_points=max_points)
                    ax.plot(tt * 1e3, sms, lw=1.1, color=colour,
                            label=f"{name} < {smooth_hz / 1e3:g} kHz")
                ax.set_ylabel(f"{name} {unit}")
                ax.legend(fontsize=6, loc="upper left")
                _mark(ax)
            out_axes[-1].set_xlabel("Time [ms]")

        head = np.atleast_1d(out_axes)[0]
        if title is None:
            bits = []
            if tru.size:
                bits.append(f"{tru.size} true")
            if det.size:
                bits.append(f"{det.size} detected")
            bits.append("green = truth, red dashed = reconstructed"
                        if tru.size and det.size else "")
            title = "  |  ".join(b for b in bits if b)
        if title:
            head.set_title(title, fontsize=8)
        fig.tight_layout()
    return fig, out_axes


def plot_iq_plane(iq, sample_rate=None, branch=None, model=None,
                  max_points=20000, figsize=(4.2, 4.0), title=None, ax=None):
    """Scatter the trace in the complex plane, coloured by decoded branch.

    The companion view to :func:`plot_trace_with_flips`: it shows *why* the
    problem is one-dimensional. The cloud is an isotropic noise blob smeared
    along a single line, and at a fixed bias it is visibly two blobs.

    ``branch`` is the decoded label per sample (``rec.branch``); ``model`` may
    be a ``StaticBlobModel``, whose two fitted centres are then marked.
    """
    iq = np.asarray(iq)
    idx = np.arange(iq.size)
    if max_points is not None and iq.size > max_points:
        idx = idx[:: int(np.ceil(iq.size / max_points))]
    z = iq[idx]

    with plt.style.context(_STYLE):
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.figure
        if branch is not None:
            b = np.asarray(branch)[idx]
            for val, colour, name in ((0, "tab:blue", "branch A"),
                                      (1, "tab:orange", "branch B")):
                m = b == val
                if np.any(m):
                    ax.plot(z[m].real, z[m].imag, ".", ms=1.0, alpha=0.25,
                            color=colour, label=name)
            ax.legend(fontsize=6, markerscale=6)
        else:
            ax.plot(z.real, z.imag, ".", ms=1.0, alpha=0.25, color="0.4")
        if model is not None and hasattr(model, "mu_a"):
            d, o = model.direction, model.origin
            for mu, colour in ((model.mu_a, "tab:blue"),
                               (model.mu_b, "tab:orange")):
                c = o + d * mu
                ax.plot([c.real], [c.imag], "x", ms=9, mew=2, color=colour,
                        zorder=5)
        ax.set_xlabel("I [a.u.]")
        ax.set_ylabel("Q [a.u.]")
        ax.set_aspect("equal", adjustable="datalim")
        if title:
            ax.set_title(title, fontsize=8)
        fig.tight_layout()
    return fig, ax
