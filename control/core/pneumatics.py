"""Pump and bleed flow laws, shared by the lumped plant and the controller's prior.

Flows are in standard ml/s (volume at atmosphere), pressures are gauge kPa.
"""

import math

RHO_AIR = 1.2          # kg/m^3


def orifice_k(diameter_mm, cd=0.6):
    """Bleed constant k in Q = k * P^0.5 (ml/s, P in kPa) for a round orifice.

    Q = Cd * A * sqrt(2 dP / rho). 0.3 mm gives about 1.7 ml/s at 1 kPa.
    """
    area = math.pi * (diameter_mm * 1e-3) ** 2 / 4
    return cd * area * math.sqrt(2000.0 / RHO_AIR) * 1e6


def bleed_flow(k, n, P_kPa):
    """Q_bleed = k * P^n for P >= 0."""
    return k * max(P_kPa, 0.0) ** n


def duty_to_drive(u, u_min):
    """Pump drive h(u) in [0, 1]: 0 below the start duty u_min, linear above."""
    if u <= u_min:
        return 0.0
    return min(1.0, (u - u_min) / (1.0 - u_min))


def drive_to_duty(h, u_min):
    """Inverse of duty_to_drive for h in (0, 1]."""
    return u_min + max(0.0, h) * (1.0 - u_min)


def pump_flow(h, P_kPa, q0, p_stall):
    """One-way pump: Q = Q0 * h * (1 - P/P_stall), never negative (no backflow)."""
    return max(0.0, q0 * h * (1.0 - P_kPa / p_stall))
