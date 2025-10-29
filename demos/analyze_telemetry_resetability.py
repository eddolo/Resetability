#!/usr/bin/env python3
# demos/analyze_telemetry_resetability.py
# Analyze rocket/satellite telemetry for resetability (R) opportunities + MP4 animation

import pandas as pd, numpy as np, matplotlib.pyplot as plt, argparse, os, time
from pathlib import Path
import sys
from matplotlib.animation import FuncAnimation

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'python'))
from so3_reset import estimate_lambda_and_R, quat_mul, quat_to_axang

def quat_conj(q):
    q = np.asarray(q, dtype=float)
    return np.array([q[0], -q[1], -q[2], -q[3]])

def quat_error(q_ref, q):
    return quat_mul(q_ref, quat_conj(q))

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=str, default="data/telemetry.csv", help="Input telemetry CSV file")
    ap.add_argument("--window", type=int, default=50, help="Sliding window size")
    ap.add_argument("--fps", type=int, default=10, help="Frame rate for timestamps")
    ap.add_argument("--animate", action="store_true", help="Create MP4 animation")
    return ap.parse_args()

if __name__ == "__main__":
    args = parse_args()
    fpath = Path(args.input)
    if not fpath.exists():
        print(f"❌ Telemetry file not found: {fpath}")
        sys.exit(1)

    df = pd.read_csv(fpath)
    if not {"qw","qx","qy","qz"}.issubset(df.columns):
        print("❌ Telemetry file must contain quaternion columns: qw,qx,qy,qz")
        sys.exit(1)

    out_dir = Path("results"); out_dir.mkdir(exist_ok=True)
    window = args.window
    fps = args.fps

    Rs, thetas, ts = [], [], []
    quats = df[["qw","qx","qy","qz"]].to_numpy()
    timestamps = df["timestamp"].to_numpy() if "timestamp" in df.columns else np.arange(len(quats)) / fps

    for i in range(window, len(quats)):
        seq = []
        for j in range(i - window, i):
            q1, q2 = quats[j], quats[j + 1]
            dq = quat_error(q2, q1)
            n, th = quat_to_axang(dq)
            seq.append((n, th))
        lam, R, th_net = estimate_lambda_and_R(seq)
        Rs.append(R)
        thetas.append(np.degrees(th_net))
        ts.append(timestamps[i])

    out_csv = out_dir / "telemetry_analysis.csv"
    pd.DataFrame({"timestamp": ts, "R": Rs, "theta_net_deg": thetas}).to_csv(out_csv, index=False)
    print(f"✅ Saved telemetry resetability results to {out_csv}")

    # --- Detect reset opportunities ---
    df_det = pd.DataFrame({"timestamp": ts, "R": Rs, "theta_net_deg": thetas})
    candidates = df_det[(df_det["R"] < 0.05) & (df_det["theta_net_deg"] > 1.0)]
    if not candidates.empty:
        opp_path = out_dir / "telemetry_reset_opportunities.csv"
        candidates.to_csv(opp_path, index=False)
        print(f"🚀 {len(candidates)} reset opportunities detected → {opp_path}")
    else:
        print("⚠️ No reset opportunities detected (no low-R, high-angle segments).")

    # --- Static Plot ---
    fig, ax1 = plt.subplots(figsize=(10, 5))
    ax1.plot(ts, Rs, color="blue", label="R")
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("R", color="blue")

    ax2 = ax1.twinx()
    ax2.plot(ts, thetas, color="orange", label="θ_net [deg]")
    ax2.set_ylabel("θ_net [deg]", color="orange")
    plt.title("Resetability Analysis from Telemetry")

    if not candidates.empty:
        ax1.scatter(candidates["timestamp"], candidates["R"], color="red", s=20, label="Reset opportunity")
        ax2.scatter(candidates["timestamp"], candidates["theta_net_deg"], color="red", s=10)

    plt.tight_layout()
    plot_path = out_dir / "telemetry_resetability_plot.png"
    plt.savefig(plot_path, dpi=150)
    plt.close(fig)
    print(f"🖼️ Saved static plot to {plot_path}")

    # --- MP4 Animation ---
    if args.animate:
        fig, ax1 = plt.subplots(figsize=(10, 5))
        ax2 = ax1.twinx()
        line1, = ax1.plot([], [], color="blue", label="R")
        line2, = ax2.plot([], [], color="orange", label="θ_net [deg]")
        scatter1 = ax1.scatter([], [], color="red", s=25)
        ax1.set_xlabel("Time [s]")
        ax1.set_ylabel("R", color="blue")
        ax2.set_ylabel("θ_net [deg]", color="orange")
        ax1.set_xlim(ts[0], ts[-1])
        ax1.set_ylim(min(Rs), max(Rs)*1.1)
        ax2.set_ylim(min(thetas)*0.9, max(thetas)*1.1)
        plt.title("Resetability Analysis Animation")

        def update(frame):
            line1.set_data(ts[:frame], Rs[:frame])
            line2.set_data(ts[:frame], thetas[:frame])
            reset_points = candidates[candidates["timestamp"] <= ts[frame]]
            scatter1.set_offsets(np.c_[reset_points["timestamp"], reset_points["R"]])
            return line1, line2, scatter1

        ani = FuncAnimation(fig, update, frames=len(ts), interval=1000/fps, blit=True)
        out_vid = out_dir / f"telemetry_resetability_{int(time.time())}.mp4"
        ani.save(out_vid, fps=fps)
        print(f"🎥 MP4 animation saved: {out_vid}")
