#!/usr/bin/env python3
# demos/booster_reset_demo.py
# Monte Carlo booster attitude simulation with reset shim
# Now includes automatic MP4 summary video creation
# Requires: pip install numpy matplotlib
# Quick demo (smaller run)
# python demos/booster_reset_demo.py --mode fast --vector-pdf --out results/booster_fast
# Full run (big Monte Carlo + MP4)
# python demos/booster_reset_demo.py --mode full --vector-pdf --out results/booster_full
# it creates folder logs/ plots/ and videos/ under the specified output folder
# If none of the runs ever reach the threshold (default = 1°), then every t_to_thr_* = np.nan.
# When matplotlib gets all-NaN or zero-height data, it just draws an empty set of axes.
# Make recovery easier so some runs record finite times:
# python demos\booster_reset_demo.py --mode full --thr 5

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'python'))
import numpy as np, matplotlib.pyplot as plt, argparse, os, csv, time
from pathlib import Path
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.animation import FuncAnimation

from so3_reset import estimate_lambda_and_R

import traceback, sys, time
try:
    print("[DEBUG] booster_reset_demo starting up...")
except Exception as e:
    with open("debug_startup_log.txt", "w") as f:
        f.write("Startup failed:\n" + traceback.format_exc())
    sys.exit(1)

# ---------- Helper Functions ----------
def quat_mul(a, b):
    w1, x1, y1, z1 = a; w2, x2, y2, z2 = b
    return np.array([
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2
    ])

def quat_conj(q): return np.array([q[0], -q[1], -q[2], -q[3]])
def quat_normalize(q): return q / (np.linalg.norm(q) + 1e-15)

def quat_from_axis_angle(n, th):
    n = np.asarray(n, float); n /= (np.linalg.norm(n) + 1e-12)
    return np.array([np.cos(th/2), *(np.sin(th/2)*n)])

def quat_to_axis_angle(q):
    q = quat_normalize(q); w = np.clip(q[0], -1.0, 1.0)
    th = 2*np.arccos(w)
    if th < 1e-12: return np.array([1.0, 0.0, 0.0]), 0.0
    s = np.sqrt(max(1 - w*w, 1e-12))
    n = q[1:] / s
    return n, th

def quat_error(q_ref, q): return quat_mul(q_ref, quat_conj(q))
def compose_seq(seq):
    q = np.array([1.0, 0.0, 0.0, 0.0])
    for n, th in seq:
        q = quat_mul(quat_from_axis_angle(n, th), q)
        q = quat_normalize(q)
    return q

# ---------- Plant & Controllers ----------
def omega_dot(omega, torque, I):
    return np.linalg.inv(I) @ (torque - np.cross(omega, I @ omega))

def quat_dot(q, omega):
    w, x, y, z = q; ox, oy, oz = omega
    O = np.array([[0, -ox, -oy, -oz],
                  [ox, 0, oz, -oy],
                  [oy, -oz, 0, ox],
                  [oz, oy, -ox, 0]])
    return 0.5 * (O @ np.array([w, x, y, z]))

class PIDAttitude:
    def __init__(self, Kp=35.0, Kd=6.0, I=np.diag([800., 900., 700.]), torque_limit=250.0):
        self.Kp, self.Kd, self.I, self.torque_limit = Kp, Kd, I, torque_limit

    def torque(self, q, omega, q_ref=np.array([1.0, 0.0, 0.0, 0.0])):
        q_err = quat_error(q_ref, q)
        n, th = quat_to_axis_angle(q_err)
        tau = -self.Kp * th * n - self.Kd * omega
        nrm = np.linalg.norm(tau)
        return tau if nrm <= self.torque_limit else tau * (self.torque_limit / (nrm + 1e-12))

class ResetShim:
    def __init__(self, dt, window_sec=0.2, r_thresh=0.1, lam_min=0.5, lam_max=3.0,
                 T_reset=0.18, omega_ref=1.0, min_impulse_dt=0.005):
        self.dt = dt; self.Nwin = int(window_sec / dt)
        self.window = []
        self.r_thresh, self.lam_min, self.lam_max = r_thresh, lam_min, lam_max
        self.T_reset, self.omega_ref, self.min_impulse_dt = T_reset, omega_ref, min_impulse_dt
        self.queue = []; self.active = False
        self.last_report = (1.0, 1.0, 0.0)

    def update_window(self, omega_b):
        th = np.linalg.norm(omega_b) * self.dt
        n = np.array([1.0, 0.0, 0.0]) if th == 0 else (omega_b / (np.linalg.norm(omega_b) + 1e-12))
        self.window.append((n, th))
        if len(self.window) > self.Nwin: self.window.pop(0)

    def plan_if_applicable(self, qbar=1.0, qbar_max=5.0):
        if len(self.window) < max(4, self.Nwin // 4):
            self.active = False; return
        lam, R, th = estimate_lambda_and_R(self.window)
        self.last_report = (lam, R, th)
        if qbar < qbar_max and R < self.r_thresh and th > np.deg2rad(0.2):
            lam = float(np.clip(lam, self.lam_min, self.lam_max))
            seq = [(n, lam * thk) for (n, thk) in self.window] * 2
            sched = []
            for n, thk in seq:
                dur = max(abs(thk) / self.omega_ref, self.min_impulse_dt)
                steps = max(1, int(dur / self.dt))
                tau = -np.sign(thk) * n
                for _ in range(steps): sched.append(tau)
            if sched:
                desired = max(1, int(self.T_reset / self.dt))
                idx = np.linspace(0, len(sched) - 1, desired).astype(int)
                self.queue = [sched[i] for i in idx]
                self.active = True
        else:
            self.active = False; self.queue = []

    def torque_ff(self, scale=25.0):
        if self.active and self.queue:
            tau = scale * np.array(self.queue.pop(0))
            if not self.queue: self.active = False
            return tau
        return np.zeros(3)

def disturbance_profile(t, p):
    tau = np.zeros(3)
    if p["gust_t0"] < t < p["gust_t1"]:
        tau += p["gust_amp"] * np.array([0.0, 0.01 * np.sin(p["gust_freq"] * t), 0.0])
    if p["bias_t0"] < t < p["bias_t1"]:
        tau += np.array([p["bias_px"], p["bias_py"], 0.0])
    if p["yaw_t0"] < t < p["yaw_t1"]:
        tau += np.array([0.0, 0.0, p["yaw_bump"]])
    return tau

def simulate(params, use_reset=True):
    I = np.diag([800., 900., 700.])
    dt = params["dt"]; T = params["T_total"]; steps = int(T / dt)
    pid = PIDAttitude(Kp=params["Kp"], Kd=params["Kd"], I=I, torque_limit=params["torque_limit"])
    shim = ResetShim(dt=dt, window_sec=0.2, r_thresh=0.1, T_reset=0.18, omega_ref=1.0)
    q = np.array([1.0, 0.0, 0.0, 0.0]); w = np.zeros(3); q_ref = np.array([1.0, 0.0, 0.0, 0.0])
    tlog = np.zeros(steps); err = np.zeros(steps); Rlog = np.full(steps, np.nan); lamlog = np.full(steps, np.nan); active = np.zeros(steps, int)

    for i in range(steps):
        t = i * dt; tlog[i] = t
        if use_reset:
            shim.update_window(w); shim.plan_if_applicable(qbar=1.0, qbar_max=5.0)
            lam, R, th = shim.last_report; Rlog[i] = R; lamlog[i] = lam; active[i] = int(shim.active)
        tau_ctrl = pid.torque(q, w, q_ref)
        tau_ff = shim.torque_ff(scale=params["ff_scale"]) if use_reset else np.zeros(3)
        tau_tot = tau_ctrl + tau_ff + disturbance_profile(t, params)
        w = w + omega_dot(w, tau_tot, I) * dt
        q = quat_normalize(q + quat_dot(q, w) * dt)
        _, th_err = quat_to_axis_angle(quat_error(q_ref, q)); err[i] = np.degrees(th_err)

    thr = params.get("thr_deg", 1.0); win = int(0.5 / dt); t_to = np.nan
    for k in range(steps - win):
        if np.all(err[k:k+win] <= thr): t_to = tlog[k]; break
    return {"t": tlog, "err": err, "R": Rlog, "lam": lamlog, "active": active,
            "t_to_thr": t_to, "residual": float(err[-1])}

# ---------- Monte Carlo ----------
def run_monte_carlo(args):
    """
    Monte Carlo booster attitude simulation with batching, auto-resume,
    and progress diagnostics.
    """
    print(f"[DEBUG] Starting Monte Carlo mode={args.mode}, runs={args.runs}, dt={args.dt}, T={args.t}")
    out = Path(args.out)
    (out / "plots").mkdir(parents=True, exist_ok=True)
    (out / "logs").mkdir(parents=True, exist_ok=True)
    print("[DEBUG] Output directories created:", out)

    # ---- Base parameters ----
    base = {"dt": args.dt, "T_total": args.t, "Kp": 35.0, "Kd": 6.0,
            "torque_limit": 250.0, "ff_scale": 25.0,
            "gust_t0": 1.0, "gust_t1": 3.0, "gust_amp": 20.0, "gust_freq": 4.0,
            "bias_t0": 4.0, "bias_t1": 4.3, "bias_px": 10.0, "bias_py": -15.0,
            "yaw_t0": 7.0, "yaw_t1": 7.2, "yaw_bump": 8.0, "thr_deg": args.thr}
    rng = np.random.default_rng(args.seed)

    # ---- Example run ----
    print("[DEBUG] Running single example simulation...")
    ex = simulate(base, use_reset=True)
    print("[DEBUG] Example simulation done, creating plot...")
    plt.figure(figsize=(8, 4))
    plt.plot(ex["t"], ex["err"])
    plt.xlabel("Time [s]"); plt.ylabel("Attitude error [deg]")
    plt.title("Example run (with reset)")
    plt.tight_layout()
    ex_fig = out / "plots/fig_attitude_example.png"
    plt.savefig(ex_fig, dpi=150); plt.close()

    # ---- Batch & resume configuration ----
    batch_size = 50
    total_runs = args.runs
    batches = (total_runs + batch_size - 1) // batch_size
    summary = out / "logs/summary.csv"

    # Determine resume point
    completed_runs = 0
    if summary.exists():
        try:
            with open(summary, newline="") as f:
                completed_runs = sum(1 for _ in csv.DictReader(f))
            print(f"[DEBUG] Resuming from existing summary.csv ({completed_runs} runs found).")
        except Exception as e:
            print(f"[WARN] Could not read summary.csv ({e}), starting fresh.")
            completed_runs = 0
    else:
        print("[DEBUG] No existing summary.csv, starting fresh.")

    start_batch = completed_runs // batch_size
    header_written = completed_runs > 0

    print(f"[DEBUG] Monte Carlo configured: {batches} batches × {batch_size} runs = {total_runs} total")
    print(f"[DEBUG] Resuming from batch {start_batch + 1}/{batches} (run {completed_runs})")
    t0 = time.time()

    for batch_idx in range(start_batch, batches):
        start = batch_idx * batch_size
        end = min(total_runs, start + batch_size)
        print(f"[DEBUG] --- Batch {batch_idx + 1}/{batches}: runs {start}–{end - 1} ---")

        rows = []
        for run in range(start, end):
            if run % 10 == 0:
                print(f"[DEBUG] Run {run}/{total_runs} in progress...")
            pset = dict(base)
            pset["gust_amp"] *= rng.uniform(0.8, 1.3)
            pset["bias_px"]  *= rng.uniform(0.8, 1.3)
            pset["bias_py"]  *= rng.uniform(0.8, 1.3)
            pset["yaw_bump"] *= rng.uniform(0.8, 1.3)
            s = rng.uniform(-0.2, 0.2)
            for k in ["gust_t0", "gust_t1", "bias_t0", "bias_t1", "yaw_t0", "yaw_t1"]:
                pset[k] = max(0.0, pset[k] + s)

            b = simulate(pset, use_reset=False)
            r = simulate(pset, use_reset=True)
            R_active = np.nanmedian(r["R"][r["active"] == 1]) if np.any(r["active"] == 1) else np.nan
            rows.append({
                "run": run, "R_active": R_active,
                "t_to_thr_base": b["t_to_thr"], "t_to_thr_reset": r["t_to_thr"],
                "residual_base": b["residual"], "residual_reset": r["residual"]
            })
            if run % 25 == 0:
                print(f"[DEBUG] Completed run {run}/{total_runs}  Residual={r['residual']:.2f}  t_to_thr={r['t_to_thr']}")

        mode = "a" if header_written else "w"
        with open(summary, mode, newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            if not header_written:
                w.writeheader()
                header_written = True
            w.writerows(rows)
        print(f"[DEBUG] Saved batch {batch_idx + 1}/{batches} → {summary}")

    print(f"[DEBUG] All batches complete in {time.time() - t0:.1f}s; compiling results...")

    # ---- Read all results back ----
    with open(summary, newline="") as f:
        reader = csv.DictReader(f)
        rows = [r for r in reader]

    # ---- Helper functions ----
    def binR(R):
        if np.isnan(R): 
            return "unknown (no reset)"
        if R < 0.05: 
            return "low(<0.05)"
        if R < 0.20: 
            return "mid(0.05-0.20)"
        return "high(>0.20)"

    bins = {
        "low(<0.05)": [],
        "mid(0.05-0.20)": [],
        "high(>0.20)": [],
        "unknown (no reset)": []   # clearer label
    }

    for r in rows:
        try:
            R_val = float(r["R_active"]) if r["R_active"] not in ["", "nan"] else np.nan
        except Exception:
            R_val = np.nan
        bins[binR(R_val)].append({
            k: (float(v) if isinstance(v, (int, float, str)) and str(v).replace('.', '', 1).isdigit() else v)
            for k, v in r.items()
        })

    def stats(B, key):
        vals = [x[key] for x in B if not np.isnan(x[key])]
        if not vals: return np.nan, np.nan, np.nan
        a = np.array(vals)
        return float(np.median(a)), float(np.percentile(a, 25)), float(np.percentile(a, 75))

    # ---- Plotting section ----
    labels = list(bins.keys()); x = np.arange(len(labels)); width = 0.35

    base_med, reset_med = [], []
    for lab in labels:
        m, _, _ = stats(bins[lab], "residual_base"); base_med.append(m)
        m, _, _ = stats(bins[lab], "residual_reset"); reset_med.append(m)
    print("[DEBUG] Generating residual vs R plot...")
    plt.figure(figsize=(9, 4))
    plt.bar(x - width / 2, base_med, width, label="baseline")
    plt.bar(x + width / 2, reset_med, width, label="reset")
    plt.xticks(x, labels, rotation=20)
    plt.ylabel("Residual [deg]")
    plt.title("Residual vs R bin (median)")
    plt.legend()
    fig_resid = out / "plots/fig_residual_vs_R.png"
    plt.tight_layout(); plt.savefig(fig_resid, dpi=150); plt.close()

    base_med, reset_med = [], []
    for lab in labels:
        m, _, _ = stats(bins[lab], "t_to_thr_base"); base_med.append(m)
        m, _, _ = stats(bins[lab], "t_to_thr_reset"); reset_med.append(m)
    print("[DEBUG] Generating time-to-threshold plot...")
    plt.figure(figsize=(9, 4))
    plt.bar(x - width / 2, base_med, width, label="baseline")
    plt.bar(x + width / 2, reset_med, width, label="reset")
    plt.xticks(x, labels, rotation=20)
    plt.ylabel(f"Time to ≤{base['thr_deg']}° [s]")
    plt.title("Recovery time vs R bin (median)")
    plt.legend()
    fig_time = out / "plots/fig_time_to_1deg_vs_R.png"
    plt.tight_layout(); plt.savefig(fig_time, dpi=150); plt.close()

    # ---- MP4 animation summary ----
    print("[DEBUG] Starting MP4 animation block...")
    try:
        print("🎥 Generating summary video...")
        out_vid = out / "videos"; out_vid.mkdir(exist_ok=True)
        mp4_path = out_vid / f"booster_reset_summary_{int(time.time())}.mp4"
        imgs = [plt.imread(fig_resid), plt.imread(fig_time)]
        fig, ax = plt.subplots(figsize=(9, 5))
        img = ax.imshow(imgs[0]); ax.axis("off")

        def update(i):
            img.set_data(imgs[i % len(imgs)])
            return [img]

        ani = FuncAnimation(fig, update, frames=20, interval=600, blit=True)
        ani.save(str(mp4_path), fps=1)
        print(f"✅ Video saved: {mp4_path}")
    except Exception as e:
        import traceback
        print("⚠️ Could not generate MP4:", repr(e))
        traceback.print_exc()

    print("✅ Monte Carlo complete.")

# ---- Export simplified results for validator ----
    import pandas as pd
    summary_for_validation = []
    with open(summary, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                t_recover = float(row["t_to_thr_reset"]) if row["t_to_thr_reset"] not in ["", "nan"] else np.nan
                summary_for_validation.append({
                    "R": float(row["R_active"]),
                    "residual": float(row["residual_reset"]),
                    "t_recover": t_recover,
                    "domain": "booster"
                })
            except Exception:
                continue

    if len(summary_for_validation) > 0:
        out_path = out / "booster_results.csv"
        pd.DataFrame(summary_for_validation).to_csv(out_path, index=False)
        print(f"✅ Booster results exported for validator → {out_path}")
    else:
        print("⚠️ No valid booster runs to export.")

# ---------- CLI ----------
def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["fast", "full"], default="fast")
    ap.add_argument("--runs", type=int, default=None)
    ap.add_argument("--dt", type=float, default=None)
    ap.add_argument("--t", type=float, default=None)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", type=str, default=".")
    ap.add_argument("--thr", type=float, default=1.0, help="Recovery threshold (deg)")
    ap.add_argument("--vector-pdf", action="store_true", help="Also save vector PDF plots")
    args = ap.parse_args()
    if args.mode == "fast":
        args.runs = args.runs or 50; args.dt = args.dt or 0.01; args.t = args.t or 12.0
    else:
        args.runs = args.runs or 500; args.dt = args.dt or 0.005; args.t = args.t or 12.0
    return args

if __name__ == "__main__":
    run_monte_carlo(parse_args())
