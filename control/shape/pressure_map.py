"""The pressure-to-shape map for one finger design and material.

Loads a table exported by `export_map.py` and interpolates it with monotone
(PCHIP) splines, so no wiggles appear between FEM points. Outside the table
the map does not extrapolate: pressure is clipped to the table's top level and
to the design's pressure cap, whichever is lower.

Units: kPa (gauge), degrees, 1/m, ml.
"""

import json
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

P_ATM_KPA = 101.325


class ShapeMap:
    def __init__(self, table, p_cap_kPa=None):
        self.table = table
        p = np.asarray(table["pressure_kPa"], float)
        kappa = np.asarray(table["curvature_per_m"], float)
        if np.any(np.diff(p) <= 0) or np.any(np.diff(kappa) <= 0):
            raise ValueError("pressure and curvature must both increase along the table")
        self.p = p
        self._kappa = PchipInterpolator(p, kappa)
        self.p_max = float(p[-1] if p_cap_kPa is None else min(p[-1], p_cap_kPa))
        self.kappa_max = float(self._kappa(self.p_max))
        # inverse by dense tabulation of the monotone kappa(P), so the round trip is exact
        self._p_fine = np.linspace(0.0, p[-1], 4001)
        self._k_fine = self._kappa(self._p_fine)
        self._bend = PchipInterpolator(p, np.asarray(table["bend_deg"], float))
        vol = table.get("cavity_volume_mm3")
        self.has_volume = vol is not None
        if self.has_volume:
            self._vol = PchipInterpolator(p, np.asarray(vol, float) / 1000.0)   # ml
            self._dvol = self._vol.derivative()

    @classmethod
    def load(cls, path, p_cap_kPa=None):
        return cls(json.loads(Path(path).read_text()), p_cap_kPa)

    def _clip(self, P):
        return np.clip(P, 0.0, self.p_max)

    def kappa(self, P):
        """Curvature (1/m) at pressure P (kPa)."""
        return self._kappa(self._clip(P))

    def bend_deg(self, P):
        return self._bend(self._clip(P))

    def pressure_for(self, kappa):
        """P*(kappa*): the pressure (kPa) that gives curvature kappa* (1/m).

        Clipped to [0, p_max], so a target above the reachable curvature gives
        the highest allowed pressure.
        """
        k = np.clip(kappa, 0.0, self.kappa_max)
        return float(self._clip(np.interp(k, self._k_fine, self._p_fine)))

    def volume_ml(self, P):
        return self._vol(self._clip(P))

    def compliance_ml_per_kPa(self, P, v_extra_ml=0.0):
        """C(P) = V/P_atm + (1 + P/P_atm) dV/dP, plus V_extra/P_atm of tube gas.

        With flow Q in standard ml/s (volume at atmosphere), C(P) dP/dt = Q.
        """
        if not self.has_volume:
            raise ValueError("this map has no cavity volume, so it has no compliance")
        P = self._clip(P)
        return (self._vol(P) + v_extra_ml) / P_ATM_KPA + (1 + P / P_ATM_KPA) * self._dvol(P)
