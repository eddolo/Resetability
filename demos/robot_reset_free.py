# demos/robot_reset_free.py
# SO(3) Reset Demo – Floating zero-gravity robot version (λ-reset with thruster bursts)
# Requires: pip install pybullet numpy pandas
# Run visually:
#   python demos/robot_reset_free.py --gui
# Record video:
#   python demos/robot_reset_free.py --gui --record

import pybullet as p, pybullet_data, numpy as np, time, sys, os
import pandas as pd
from pathlib import Path

# --- PATH FIX ---
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'python'))
from so3_reset import estimate_lambda_and_R
# --- END PATH FIX ---

def random_seq(N=50, spread_deg=5):
    ths = np.deg2rad(np.random.randn(N) * spread_deg)
    ns = np.random.randn(N, 3)
    ns /= np.linalg.norm(ns, axis=1, keepdims=True)
    return list(zip(ns, ths))

def orbit_camera(t, radius=1.5, height=0.3, speed=0.3):
    yaw = np.degrees(speed * t)
    p.resetDebugVisualizerCamera(cameraDistance=radius, cameraYaw=yaw,
                                 cameraPitch=-30, cameraTargetPosition=[0, 0, 0.2])

if __name__ == "__main__":
    use_gui = "--gui" in sys.argv
    record = "--record" in sys.argv

    if use_gui:
        p.connect(p.GUI)
        p.configureDebugVisualizer(p.COV_ENABLE_GUI, 0)
    else:
        p.connect(p.DIRECT)

    p.setAdditionalSearchPath(pybullet_data.getDataPath())
    p.resetSimulation()
    p.setGravity(0, 0, 0)  # 🛰️ ZERO GRAVITY
    p.setRealTimeSimulation(0)

    # Floating cube
    visual_shape = p.createVisualShape(p.GEOM_BOX, halfExtents=[0.25,0.15,0.05],
                                       rgbaColor=[0.2,0.7,1,1])
    collision_shape = p.createCollisionShape(p.GEOM_BOX, halfExtents=[0.25,0.15,0.05])
    cube = p.createMultiBody(baseMass=1.0, baseCollisionShapeIndex=collision_shape,
                             baseVisualShapeIndex=visual_shape, basePosition=[0,0,0.3])

    dt = 1 / 240
    initial_pos, initial_orn = p.getBasePositionAndOrientation(cube)

    if record:
        os.makedirs("videos", exist_ok=True)
        vid_path = f"videos/robot_reset_free_{int(time.time())}.mp4"
        vid_id = p.startStateLogging(p.STATE_LOGGING_VIDEO_MP4, vid_path)
        print(f"🎥 Recording to {vid_path}")

    seq = random_seq(60, spread_deg=5)
    lam, R, th = estimate_lambda_and_R(seq)
    print(f"λ = {lam:.3f}, R = {R:.3f}, θ_net = {np.degrees(th):.2f}°")

    t = 0.0
    debug_id = None

    # --- PHASE 1: Disturbance ---
    for n, th_k in seq:
        torque = 2.5 * np.array(n)
        p.applyExternalTorque(cube, -1, torque, p.WORLD_FRAME)
        if use_gui:
            orbit_camera(t)
            if debug_id: p.removeUserDebugItem(debug_id)
            debug_id = p.addUserDebugText("Phase 1 – Disturbance", [0.6, 0, 0.7],
                                          textColorRGB=[1,1,1], textSize=1.5, lifeTime=0.1)
        p.stepSimulation(); t += dt
        if use_gui: time.sleep(0.002)

    pre_reset_t = t

    # --- PHASE 2: Thruster-based λ reset ---
    reset_seq = [(n, lam * thk) for (n, thk) in seq] * 2
    for n, th_k in reset_seq:
        torque = -12.0 * np.array(n) * np.sign(th_k)  # stronger bursts
        burst_len = np.random.randint(3, 8)
        for _ in range(burst_len):
            p.applyExternalTorque(cube, -1, torque, p.WORLD_FRAME)
            if use_gui:
                orbit_camera(t)
                p.addUserDebugLine([0,0,0.3], np.add([0,0,0.3], 0.3*torque),
                                   [0,1,0], lineWidth=3, lifeTime=0.05)
            p.stepSimulation(); t += dt
            if use_gui: time.sleep(0.002)
        # Cooldown (valve closed)
        for _ in range(int(burst_len * 0.5)):
            p.stepSimulation(); t += dt
            if use_gui: time.sleep(0.002)

    post_reset_t = t
    final_pos, final_orn = p.getBasePositionAndOrientation(cube)

    diff_q = p.getDifferenceQuaternion(initial_orn, final_orn)
    _, resid_rad = p.getAxisAngleFromQuaternion(diff_q)
    resid_deg = np.degrees(resid_rad)
    print(f"Residual attitude error ≈ {resid_deg:.2f}°")

    # --- Save results ---
    results_dir = Path("results"); results_dir.mkdir(exist_ok=True)
    results_file = results_dir / "robot_results.csv"
    environment = "zeroG"

    results_df = pd.DataFrame([{
        "R_active": lam,
        "residual_error": resid_deg,
        "t_to_thr": post_reset_t - pre_reset_t,
        "environment": environment
    }])

    if results_file.exists():
        results_df.to_csv(results_file, mode="a", header=False, index=False)
        print(f"✅ Appended new results to {results_file}")
    else:
        results_df.to_csv(results_file, index=False)
        print(f"✅ Created new results file: {results_file}")

    if record:
        p.stopStateLogging(vid_id)
        print(f"✅ Video saved: {vid_path}")

    print("✅ Simulation complete – close window to exit.")
