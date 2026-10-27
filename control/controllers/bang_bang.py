"""bang_bang: pump at u_on below P* - h, off above P* + h, hold in between.

The first closed loop to try. It tolerates a noisy or biased sensor because it
only needs the sign of the error, and the band h absorbs sensor noise: choose
h at least about the sensor's noise. The pressure ripples around P* by about
h plus the overshoot during one control tick and the sensor latency.
"""


class BangBang:
    def __init__(self, h_kPa=0.15, u_on=0.75, **_unused):
        self.h, self.u_on = h_kPa, u_on
        self.on = False

    def reset(self):
        self.on = False

    def step(self, t, P_hat, dP_hat, P_target):
        if P_target <= 0.05:
            self.on = False
        elif P_hat < P_target - self.h:
            self.on = True
        elif P_hat > P_target + self.h:
            self.on = False
        return self.u_on if self.on else 0.0
