"""Multivariable Mendelian randomization (summary-data IVW).

Estimates the direct effect of several exposures jointly, holding the others
fixed. Used here to separate SHBG from bioavailable testosterone, whose genetic
instruments overlap, which is where the two prior single-exposure MR studies
disagree.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats


@dataclass
class MVMRResult:
    names: list[str]
    beta: list[float]
    se: list[float]
    n_snps: int

    def as_dict(self) -> dict:
        out = {"n_snps": self.n_snps, "exposures": {}}
        for i, nm in enumerate(self.names):
            b, s = self.beta[i], self.se[i]
            z = b / s
            out["exposures"][nm] = {
                "or": float(np.exp(b)),
                "or_ci_low": float(np.exp(b - 1.96 * s)),
                "or_ci_high": float(np.exp(b + 1.96 * s)),
                "p": float(2 * stats.norm.sf(abs(z))),
            }
        return out


def mvmr_ivw(Bx, by, se_out, names) -> MVMRResult:
    """Multivariable IVW with multiplicative random-effects standard errors.

    Bx: (n_snp x k) matrix of SNP->exposure effects, harmonised to a common allele.
    by: (n_snp,) SNP->outcome effects; se_out: (n_snp,) their standard errors.
    """
    Bx = np.asarray(Bx, dtype=float)
    by = np.asarray(by, dtype=float)
    w = 1.0 / np.asarray(se_out, dtype=float) ** 2
    W = np.diag(w)
    xtwx = Bx.T @ W @ Bx
    beta = np.linalg.solve(xtwx, Bx.T @ W @ by)
    resid = by - Bx @ beta
    dof = max(len(by) - Bx.shape[1], 1)
    phi = max(1.0, float((resid * w) @ resid) / dof)   # random-effects scaling
    cov = np.linalg.inv(xtwx) * phi
    se = np.sqrt(np.diag(cov))
    return MVMRResult(list(names), [float(b) for b in beta], [float(s) for s in se], len(by))
