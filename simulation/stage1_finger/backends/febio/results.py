"""Read FEBio text output: node_data log files and the main .log.

A node_data record written to its own file (file="...") looks like
(FECore/DataRecord.cpp, v4.13):

    *Step  = 3
    *Time  = 0.3
    *Data  = tip
    1,12.5,0.0,19.0,...

One record per output point (the writer sets OUTPUT_MUST_POINTS). The first column is the 1-based position
in the node set, not the node id. Rows follow node-set order, which is the
order the writer used. The binary .xplt is not read.
"""

import re
from pathlib import Path

import numpy as np

_HEADER = re.compile(r"^\*(Step|Time|Data)\s*=\s*(.*)$")
# FEBio prints a spaced banner at the end of the log. The exact text is not in
# the docs; these patterns match the "N O R M A L   T E R M I N A T I O N" banner.
_NORMAL = re.compile(r"N\s*O\s*R\s*M\s*A\s*L\s+T\s*E\s*R\s*M\s*I\s*N\s*A\s*T\s*I\s*O\s*N")
_ERROR = re.compile(r"E\s*R\s*R\s*O\s*R\s+T\s*E\s*R\s*M\s*I\s*N\s*A\s*T\s*I\s*O\s*N")


def read_node_data(path):
    """Records of one node_data file: list of {step, time, values (n, k)}."""
    records, cur, rows = [], None, []

    def flush():
        if cur is not None:
            cur["values"] = np.array(rows, float) if rows else np.zeros((0, 0))
            records.append(cur)

    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            m = _HEADER.match(line)
            if m:
                key, val = m.groups()
                if key == "Step":
                    flush()
                    cur, rows = {"step": int(val)}, []
                elif cur is not None:
                    cur[key.lower()] = float(val) if key == "Time" else val
                continue
            parts = [p for p in re.split(r"[,\s]+", line) if p]
            rows.append([float(p) for p in parts[1:]])   # drop the item number
    flush()
    return records


def log_status(log_path):
    """'normal', 'error', or 'unknown' from the main FEBio log."""
    path = Path(log_path)
    if not path.exists():
        return "missing"
    text = path.read_text(encoding="utf-8", errors="replace")
    if _NORMAL.search(text):
        return "normal"
    if _ERROR.search(text):
        return "error"
    return "unknown"


def read_results(run_dir, log_files, n_nodes):
    """Per-step tip data, clamp reactions, and full displacement fields."""
    run_dir = Path(run_dir)
    tip = read_node_data(run_dir / log_files["tip"])
    clamp = read_node_data(run_dir / log_files["clamp"])
    nodes = read_node_data(run_dir / log_files["all_nodes"])
    times = np.array([r["time"] for r in nodes])
    u = np.stack([r["values"] for r in nodes]) if nodes else np.zeros((0, n_nodes, 3))
    if u.shape[1:] != (n_nodes, 3):
        raise ValueError(f"all-node log has shape {u.shape}, expected (steps, {n_nodes}, 3)")
    by_time = lambda recs: {round(r["time"], 9): r["values"] for r in recs}
    tip_t, clamp_t = by_time(tip), by_time(clamp)
    steps = []
    for t in times:
        key = round(t, 9)
        tv, cv = tip_t.get(key), clamp_t.get(key)
        steps.append({
            "time": float(t),
            "tip_reaction_N": tv[:, 6:9].sum(0).tolist() if tv is not None else None,
            "clamp_reaction_N": cv.sum(0).tolist() if cv is not None else None,
        })
    return times, u, steps


_STATS = {
    "time_steps": r"Number of time steps completed \.+ : (\d+)",
    "equilibrium_iterations": r"Total number of equilibrium iterations \.+ : (\d+)",
    "rhs_evaluations": r"Total number of right hand evaluations \.+ : (\d+)",
    "reformations": r"Total number of stiffness reformations \.+ : (\d+)",
    "linear_solver_s": r"time in linear solver \.+ : \S+ \(([\d.]+) sec\)",
    "solve_s": r"Solve time \.+ : \S+ \(([\d.]+) sec\)",
    "equations": r"Nr of equations \.+ : (\d+)",
}


def log_stats(log_path):
    """Iteration and timing totals from the end of an FEBio log (empty if missing)."""
    path = Path(log_path)
    if not path.exists():
        return {}
    text = path.read_text(encoding="utf-8", errors="replace")
    out = {}
    for key, pat in _STATS.items():
        m = re.search(pat, text)
        if m:
            out[key] = float(m.group(1)) if "." in m.group(1) else int(m.group(1))
    out["retries"] = len(re.findall(r"AUTO STEPPER: retry step", text))
    out["negative_jacobian_events"] = len(re.findall(r"negative jacobians? detected", text))
    times = re.findall(r"------- converged at time : ([\d.eE+-]+)", text)
    out["last_converged_time"] = float(times[-1]) if times else 0.0
    return out
