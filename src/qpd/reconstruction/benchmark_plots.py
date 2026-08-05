"""Diagnostic plots for the surrogate-replay benchmark.

Three questions a detector paper has to answer about its reconstruction, and
one plot each:

:func:`plot_efficiency_vs_rate`
    *How fast can the background tunnel before flips start being lost?* The
    background is Poisson, so the rate alone sets how often two tunnels land
    unresolvably close -- and an unresolved pair is invisible, not merely
    mistimed, because two toggles return the parity to where it started.

:func:`plot_burst_efficiency`
    *How big must a burst be before it is found at all?* This is burst-level,
    not flip-level: the object is the energy deposition, and the question is
    whether it is distinguishable from a background coincidence.

:func:`plot_burst_multiplicity`
    *Once found, how many of its quasiparticles are counted?* The answer
    saturates -- a large burst crowds its tunnels into a few milliseconds and
    the pairs that cancel are never seen -- so this is the plot that says how
    far a quoted multiplicity can be trusted.

All three take the objects the sweeps in :mod:`.benchmark` return, and all
three return ``(fig, ax)`` so the caller can annotate or save.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

__all__ = ["plot_efficiency_vs_rate", "plot_burst_efficiency",
           "plot_burst_multiplicity"]

_STYLE = Path(__file__).resolve().parents[1] / "theory" / "qpd.mplstyle"


def _injected_rate(report) -> float:
    """The rate a swept report actually injected."""
    if report.trials:
        return float(np.mean([t.rate_injected for t in report.trials]))
    return float(report.fidelity.rate_hz)


def plot_efficiency_vs_rate(
    reports,
    show_purity: bool = True,
    figsize=(6.4, 4.2),
    title=None,
    ax=None,
):
    """Flip detection efficiency against background tunnelling rate.

    Parameters
    ----------
    reports : list of BenchmarkReport
        Output of :func:`~qpd.reconstruction.benchmark.sweep_rate`.
    show_purity : bool
        Also draw purity, which behaves quite differently: crowding destroys
        recall long before it costs precision, because the flips that are lost
        are lost silently rather than replaced by spurious ones.

    Returns
    -------
    (fig, ax)
    """
    reports = list(reports)
    if not reports:
        raise ValueError("no reports to plot")
    rate = np.array([_injected_rate(r) for r in reports], float)
    eff = np.array([r.efficiency[0] for r in reports], float)
    eff_e = np.array([r.efficiency[1] for r in reports], float)
    pur = np.array([r.purity[0] for r in reports], float)
    pur_e = np.array([r.purity[1] for r in reports], float)
    order = np.argsort(rate)

    with plt.style.context(_STYLE):
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.figure
        ax.errorbar(rate[order], eff[order], yerr=eff_e[order], marker="o",
                    ms=4, capsize=2.5, lw=1.4, label="efficiency (recall)")
        if show_purity:
            ax.errorbar(rate[order], pur[order], yerr=pur_e[order], marker="s",
                        ms=3.5, capsize=2.5, lw=1.2, ls="--",
                        label="purity (precision)")
        ax.set_xscale("log")
        ax.set_xlabel("background tunnelling rate [Hz per state]")
        ax.set_ylabel("detection performance")
        ax.set_ylim(-0.02, 1.04)
        ax.axhline(1.0, color="0.6", lw=0.7, zorder=0)

        f = reports[0].fidelity
        # The measured trace's own rate: the one point on this curve that is
        # not a hypothetical.
        ax.axvline(f.rate_hz, color="0.4", lw=0.9, ls=":", zorder=0)
        # Annotate at the top: the curve starts flat at 1.0 on the left, so the
        # bottom-left corner belongs to the legend.
        ax.annotate(f"measured\n{f.rate_hz:.3g} Hz", xy=(f.rate_hz, 0.90),
                    xytext=(4, 0), textcoords="offset points",
                    fontsize="x-small", color="0.35", va="top")
        ax.legend(loc="lower left", frameon=False)
        ax.set_title(title or
                     f"Flip detection vs background rate "
                     f"(contrast {f.contrast_median:.2f}, "
                     f"{reports[0].fidelity.sample_rate / 1e3:.0f} kSa/s)")

        # The matching tolerance is not necessarily constant along this curve
        # -- `sweep_rate` shrinks it with the dwell -- and a reader comparing
        # points has to know that.
        tols = np.array([r.tol for r in reports], float)
        if tols.size and not np.allclose(tols, tols[0]):
            ax.annotate(
                f"match tolerance shrinks with dwell: "
                f"{tols.max() * 1e6:.0f} to {tols.min() * 1e6:.1f} $\\mu$s",
                xy=(0.5, -0.20), xycoords="axes fraction", ha="center",
                fontsize="x-small", color="0.35")
        fig.tight_layout()
    return fig, ax


def plot_burst_efficiency(
    points,
    figsize=(6.4, 4.2),
    title=None,
    ax=None,
):
    """Burst-level detection efficiency against burst quasiparticle number.

    Parameters
    ----------
    points : list of BurstSizePoint
        Output of :func:`~qpd.reconstruction.benchmark.sweep_burst_size`.

    Returns
    -------
    (fig, ax)
    """
    points = list(points)
    if not points:
        raise ValueError("no sweep points to plot")
    n_qp = np.array([p.n_qp_expected for p in points], float)
    eff = np.array([p.efficiency for p in points], float)
    err = np.array([p.efficiency_err for p in points], float)
    order = np.argsort(n_qp)

    with plt.style.context(_STYLE):
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.figure
        ax.errorbar(n_qp[order], eff[order], yerr=err[order], marker="o",
                    ms=4, capsize=2.5, lw=1.4, color="C0")
        ax.axhline(1.0, color="0.6", lw=0.7, zorder=0)
        ax.axhline(0.5, color="0.85", lw=0.7, ls=":", zorder=0)

        # Where the curve crosses 0.5 is the useful summary number: the burst
        # size this reconstruction is half-sensitive to.
        e, q = eff[order], n_qp[order]
        cross = np.flatnonzero((e[:-1] < 0.5) & (e[1:] >= 0.5))
        if cross.size:
            i = int(cross[0])
            span = e[i + 1] - e[i]
            n50 = q[i] + (0.5 - e[i]) * (q[i + 1] - q[i]) / span if span else q[i]
            ax.axvline(n50, color="C3", lw=0.9, ls="--", zorder=0)
            ax.annotate(f"50% at {n50:.1f} qp", xy=(n50, 0.5),
                        xytext=(5, -12), textcoords="offset points",
                        fontsize="x-small", color="C3")

        ax.set_xscale("log")
        ax.set_xlabel("burst quasiparticle number (Poisson mean)")
        ax.set_ylabel("burst detection efficiency")
        ax.set_ylim(-0.03, 1.05)
        bg = points[0].background_rate_hz
        ax.set_title(title or
                     f"Burst detection vs multiplicity "
                     f"(background {bg:.3g} Hz, measured from the data)")
        fig.tight_layout()
    return fig, ax


def plot_burst_multiplicity(
    points,
    figsize=(6.4, 5.4),
    title=None,
    axes=None,
):
    """True against reconstructed burst multiplicity, with the residual.

    Upper panel: one point per detected burst, true quasiparticle count against
    the count the reconstruction recovered, over the identity line. Lower
    panel: the mean residual, which is the bias a quoted multiplicity carries.

    Only *detected* bursts appear -- an undetected burst has no measured
    multiplicity. That makes the low-multiplicity end a selected sample: a
    2-quasiparticle burst is generally found only when it fluctuated upward,
    so the points there sit above where an unbiased sample would. The
    saturation at the high end is the real effect and is not a selection.

    Parameters
    ----------
    points : list of BurstSizePoint
        Output of :func:`~qpd.reconstruction.benchmark.sweep_burst_size`.

    Returns
    -------
    (fig, (ax_main, ax_resid))
    """
    points = list(points)
    tr = np.concatenate([p.n_qp_true for p in points if p.n_qp_true.size]
                        or [np.empty(0)])
    de = np.concatenate([p.n_qp_detected for p in points
                         if p.n_qp_detected.size] or [np.empty(0)])
    if tr.size == 0:
        raise ValueError("no detected bursts to plot")

    with plt.style.context(_STYLE):
        if axes is None:
            fig, (ax, axr) = plt.subplots(
                2, 1, figsize=figsize, sharex=True,
                gridspec_kw={"height_ratios": [2.4, 1.0], "hspace": 0.08})
        else:
            ax, axr = axes
            fig = ax.figure

        hi = float(max(tr.max(), de.max())) * 1.05
        ax.plot([0, hi], [0, hi], color="0.55", lw=0.9, ls="--",
                label="perfect recovery", zorder=1)
        ax.plot(tr, de, "o", ms=3.0, alpha=0.35, mew=0, color="C0",
                label="detected burst", zorder=2)

        # Binned mean: the trend the scatter is hiding.
        edges = np.unique(np.round(
            np.geomspace(max(tr.min(), 1.0), max(tr.max(), 2.0), 12)))
        if edges.size >= 2:
            idx = np.clip(np.digitize(tr, edges) - 1, 0, edges.size - 2)
            bx, by, be = [], [], []
            for b in range(edges.size - 1):
                sel = idx == b
                if np.count_nonzero(sel) >= 3:
                    bx.append(float(np.mean(tr[sel])))
                    by.append(float(np.mean(de[sel])))
                    be.append(float(np.std(de[sel], ddof=1)
                                    / np.sqrt(np.count_nonzero(sel))))
            if bx:
                ax.errorbar(bx, by, yerr=be, marker="o", ms=4.5, lw=1.5,
                            capsize=2.5, color="C3", label="binned mean",
                            zorder=3)
                axr.errorbar(bx, np.array(by) - np.array(bx), yerr=be,
                             marker="o", ms=4.5, lw=1.5, capsize=2.5,
                             color="C3", zorder=3)

        ax.set_ylabel("reconstructed quasiparticles")
        ax.legend(loc="upper left", frameon=False)
        ax.set_title(title or "Burst multiplicity: truth vs reconstruction")

        axr.axhline(0.0, color="0.55", lw=0.9, ls="--", zorder=1)
        axr.plot(tr, de - tr, "o", ms=2.5, alpha=0.25, mew=0, color="C0",
                 zorder=2)
        axr.set_xlabel("true quasiparticles in burst")
        axr.set_ylabel("bias")
        fig.tight_layout()
    return fig, (ax, axr)
