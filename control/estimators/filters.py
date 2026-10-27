"""Small streaming filters. The EMA and alpha-beta `update(t, x)` return (x_hat, dx_hat)."""

from collections import deque


class Median3:
    """Median of the last three values; passes the first two through unchanged."""

    def __init__(self):
        self.buf = deque(maxlen=3)

    def update(self, x):
        self.buf.append(x)
        return sorted(self.buf)[1] if len(self.buf) == 3 else x

    def reset(self):
        self.buf.clear()


class EMA:
    """Exponential moving average, with a smoothed finite-difference slope."""

    def __init__(self, alpha=0.4):
        self.alpha = alpha
        self.reset()

    def reset(self):
        self.x = self.v = self.t = None

    def update(self, t, z):
        if self.x is None:
            self.x, self.v, self.t = z, 0.0, t
            return self.x, self.v
        dt = max(t - self.t, 1e-6)
        x_new = self.x + self.alpha * (z - self.x)
        self.v += self.alpha * ((x_new - self.x) / dt - self.v)
        self.x, self.t = x_new, t
        return self.x, self.v


class AlphaBeta:
    """Alpha-beta filter: position and velocity, tolerant of uneven sample times."""

    def __init__(self, alpha=0.5, beta=0.15):
        self.alpha, self.beta = alpha, beta
        self.reset()

    def reset(self):
        self.x = self.v = self.t = None

    def update(self, t, z):
        if self.x is None:
            self.x, self.v, self.t = z, 0.0, t
            return self.x, self.v
        dt = max(t - self.t, 1e-6)
        pred = self.x + self.v * dt
        r = z - pred
        self.x = pred + self.alpha * r
        self.v += self.beta * r / dt
        self.t = t
        return self.x, self.v


def make_filter(kind, cfg):
    if kind == "alpha_beta":
        return AlphaBeta(cfg.alpha, cfg.beta)
    if kind == "ema":
        return EMA(cfg.ema_alpha)
    raise ValueError(f"unknown filter {kind!r}; use 'alpha_beta' or 'ema'")
