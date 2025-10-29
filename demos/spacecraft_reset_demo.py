# demos/spacecraft_reset_demo.py
# Spacecraft attitude reset simulation (SO3 Resetability)
# MP4 + optional 3D viz + CSV logging for cross-domain validation
# Usage:
#   python demos/spacecraft_reset_demo.py
#   python demos/spacecraft_reset_demo.py --record
#   python demos/spacecraft_reset_demo.py --3d --record

import sys, os, time
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

# make sure Python can find the SO(3) reset module
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'python'))
from so3_reset import (
    compose_seq,
    estimate_lambda_and_R,
    axang_to_quat,
    quat_mul,
    quat_normalize,
    quat_to_axang,
)

# --- Optional visualization (PyBullet 3D orientation view) ---
use_3d = "--3d" in sys.argv
if use_3d:
    import pybullet as p, pybullet_data
    p.connect(p.GUI)
    p.setAdditionalSearchPath(pybullet_data.getDataPath())
    p.resetSimulation()
    p.loadURDF("plane.urdf")

    # Big red “spacecraft” box so rotation is obvious
    visual_shape = p.createVisualShape(
        p.GEOM_BOX,
        halfExtents=[0.3, 0.15, 0.05],
        rgbaColor=[0.8, 0.2, 0.2, 1.0],
    )
    collision_shape = p.createCollisionShape(p.GEOM_BOX, halfExtents=[0.3, 0.15, 0.05])
    sat = p.createMultiBody(
        baseMass=5.0,
        baseCollisionShapeIndex=collision_shape,
        baseVisualShapeIndex=visual_shape,
        basePosition=[0, 0, 0.3],
    )

    p.setGravity(0, 0, 0)
    p.resetDebugVisualizerCamera(
        cameraDistance=1.2,
        cameraYaw=40,
        cameraPitch=-25,
        cameraTargetPosition=[0, 0, 0.3],
    )

# --- Helpers ---
def omega_dot(omega, tau, I):
    return np.linalg.inv(I) @ (tau - np.cross(omega, I @ omega))

def quat_dot(q, omega):
    w, x, y, z = q
    ox, oy, oz = omega
    O = np.array([
        [0, -ox, -oy, -oz],
        [ox,   0,  oz, -oy],
        [oy, -oz,   0,  ox],
        [oz,  oy, -ox,   0],
    ])
    return 0.5 * (O @ np.array([w, x, y, z]))

def quat_conj(q):
    q = np.asarray(q, dtype=float)
    return np.array([q[0], -q[1], -q[2], -q[3]])

def integrate(q0, w0, torque_fun, I, dt, steps):
    q, w = q0.copy(), w0.copy()
    logs = []
    for i in range(steps):
        t = i * dt
        tau = torque_fun(t)
        w = w + omega_dot(w, tau, I) * dt
        q = quat_normalize(q + quat_dot(q, w) * dt)
        _, th = quat_to_axang(q)
        logs.append((t, np.degrees(th), q.copy()))
    return np.array(logs, dtype=object)

def torque_pattern(seq, Tmax=0.05, scale=1.0, repeat=1):
    pats = []
    for _ in range(repeat):
        for n, th in seq:
            dur = abs(th) / 0.1
            tau = Tmax * n * np.sign(th) * scale
            pats.append((dur, tau))
    return pats

def make_fun(pats):
    T = sum(d for d, _ in pats)
    def f(t):
        acc = 0.0
        for d, tau in pats:
            if acc <= t < acc + d:
                return tau
            acc += d
        return np.zeros(3)
    return f, T

# --- Main simulation ---
if __name__ == "__main__":
    out = Path("videos"); out.mkdir(exist_ok=True)
    video_mode = "--record" in sys.argv

    if video_mode:
        vid_path = out / f"spacecraft_reset_{int(time.time())}.mp4"
        print(f"🎥 Recording MP4 to {vid_path}")

    # Random rotation sequence
    rng = np.random.default_rng(int(time.time()))
    N = 100
    axes = rng.normal(size=(N, 3)); axes /= np.linalg.norm(axes, axis=1, keepdims=True)
    angs = np.deg2rad(rng.normal(scale=3, size=N))
    seq = list(zip(axes, angs))

    lam, R, th = estimate_lambda_and_R(seq)
    print(f"λ = {lam:.3f}   R = {R:.3f}   θ_net = {np.degrees(th):.2f}°")

    # Stronger “disturbance”, λ-scaled reset (twice)
    dist  = torque_pattern(seq, Tmax=0.05, scale=4.0,        repeat=1)
    reset = torque_pattern(seq, Tmax=0.05, scale=lam * 4.0,  repeat=2)

    f1, T1 = make_fun(dist)
    f2, T2 = make_fun(reset)
    I  = np.diag([10, 12,  8])
    q0 = np.array([1, 0, 0, 0])
    w0 = np.zeros(3)
    dt = 0.01

    # Simulate both phases
    a1 = integrate(q0,       w0, f1, I, dt, int(T1/dt))
    a2 = integrate(a1[-1, 2], w0, f2, I, dt, int(T2/dt))

    # Plot (attitude error vs time)
    t1   = a1[:, 0];       err1 = a1[:, 1]
    t2   = T1 + a2[:, 0];  err2 = a2[:, 1]
    plt.figure(figsize=(8, 4))
    plt.plot(t1, err1, label="disturbance")
    plt.plot(t2, err2, label=f"reset (λ = {lam:.2f})")
    plt.axvline(T1, color='k', ls='--')
    plt.xlabel("Time [s]")
    plt.ylabel("Attitude error [deg]")
    plt.legend(); plt.tight_layout()
    fig_path = out / f"spacecraft_reset_plot_{int(time.time())}.png"
    plt.savefig(fig_path, dpi=150)
    print(f"🖼️ Plot saved to {fig_path}")

    # === Metrics for validation ===
    # Residual error at end of reset phase:
    residual_deg = float(err2[-1])

    # Time to reach small error after reset starts (≤ 1°); fall back to full T2
    thr_deg = 1.0
    idx = np.where(err2 <= thr_deg)[0]
    t_recover = float(a2[idx[0], 0]) if idx.size > 0 else float(T2)

    # === Optional 3D playback ===
    if use_3d:
        # Build quaternion list for viz
        q_all = np.array([np.array(q, dtype=float) for q in np.concatenate((a1[:, 2], a2[:, 2]))])

        if video_mode:
            log_id = p.startStateLogging(p.STATE_LOGGING_VIDEO_MP4, str(vid_path))

        print(f"[3D] Starting simulation with {len(q_all)} frames...")
        FAST_MODE = True   # flip to False to see trails, vectors, text

        for i, q in enumerate(q_all):
            if i % 500 == 0:
                print(f"[3D] Frame {i}/{len(q_all)}")

            # Apply quaternion to satellite (PyBullet uses [x,y,z,w])
            p.resetBasePositionAndOrientation(
                sat, [0, 0, 0.3], [float(q[1]), float(q[2]), float(q[3]), float(q[0])]
            )

            if not FAST_MODE and i % 10 == 0:
                yaw = 40 + 10 * np.sin(i * 0.002)
                p.resetDebugVisualizerCamera(
                    cameraDistance=1.2, cameraYaw=yaw, cameraPitch=-25,
                    cameraTargetPosition=[0, 0, 0.3]
                )

            p.stepSimulation()

            # Throttle a touch while recording so ffmpeg keeps up
            if FAST_MODE and video_mode and (i % 10 == 0):
                time.sleep(0.02)
            elif not FAST_MODE:
                time.sleep(0.01)

        if video_mode:
            p.stopStateLogging(log_id)
            print(f"✅ Video saved: {vid_path}")

        print("✅ 3D visualization complete. Close window to exit.")
    else:
        if video_mode:
            # Fallback: simple MP4 of the 2D plot (already saved as PNG)
            import matplotlib.animation as animation
            fig, ax = plt.subplots()
            line, = ax.plot([], [], lw=2)
            ax.set_xlim(0, T1 + T2)
            ax.set_ylim(0, max(err1.max(), err2.max()) * 1.1)
            ax.set_xlabel("Time [s]"); ax.set_ylabel("Attitude error [deg]")
            ax.set_title("Spacecraft Reset Simulation")
            all_t = np.concatenate((t1, t2))
            all_e = np.concatenate((err1, err2))
            def init(): line.set_data([], []); return (line,)
            def update(frame): line.set_data(all_t[:frame], all_e[:frame]); return (line,)
            ani = animation.FuncAnimation(fig, update, frames=len(all_t), init_func=init, interval=dt*1000)
            ani.save(str(vid_path), fps=int(1/dt))
            print(f"✅ MP4 saved: {vid_path}")
        print("✅ Simulation complete.")

    # === Append results for validator ===
    results_dir = Path("results/spacecraft")
    results_dir.mkdir(parents=True, exist_ok=True)
    results_file = results_dir / "spacecraft_results.csv"

    row = {
        "R": R,
        "residual": residual_deg,
        "t_recover": t_recover,
        "domain": "spacecraft",
    }
    df = pd.DataFrame([row])
    if results_file.exists():
        df.to_csv(results_file, mode="a", header=False, index=False)
        print(f"✅ Appended spacecraft results to {results_file}")
    else:
        df.to_csv(results_file, index=False)
        print(f"✅ Created {results_file} and wrote first row")
