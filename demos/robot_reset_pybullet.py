# demos/robot_reset_pybullet.py
# Visual robot reset demo (SO(3) Resetability) – gravity version
# Includes logging, torque visualization, environment tagging
# Run visually:
#   python demos/robot_reset_pybullet.py --gui
# Record video:
#   python demos/robot_reset_pybullet.py --gui --record

import pybullet as p, pybullet_data, numpy as np, time, sys, os
import pandas as pd
from pathlib import Path

# --- PATH FIX ---
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'python'))
from so3_reset import estimate_lambda_and_R
# --- END PATH FIX ---

# --- LOGGING SETUP ---
os.makedirs("logs", exist_ok=True)
imu_log_file = open("logs/imu_log.csv", "w")
imu_log_file.write("timestamp,q_x,q_y,q_z,q_w,ang_vel_x,ang_vel_y,ang_vel_z\n")
reset_log_file = open("logs/reset_log.csv", "w")
reset_log_file.write("lambda,R,theta_net\n")
# --- END LOGGING SETUP ---

def random_seq(N=50, spread_deg=5):
    ths = np.deg2rad(np.random.randn(N) * spread_deg)
    ns = np.random.randn(N, 3)
    ns /= np.linalg.norm(ns, axis=1, keepdims=True)
    return list(zip(ns, ths))

def orbit_camera(t, radius=2.0, height=0.8, speed=0.3):
    yaw = np.degrees(speed * t)
    p.resetDebugVisualizerCamera(
        cameraDistance=radius,
        cameraYaw=yaw,
        cameraPitch=-25,
        cameraTargetPosition=[0, 0, 0.4]
    )

def log_imu_data(timestamp, body_id):
    _, orn_quat = p.getBasePositionAndOrientation(body_id)
    _, ang_vel = p.getBaseVelocity(body_id)
    q_x, q_y, q_z, q_w = orn_quat
    ang_vel_x, ang_vel_y, ang_vel_z = ang_vel
    imu_log_file.write(f"{timestamp},{q_x},{q_y},{q_z},{q_w},{ang_vel_x},{ang_vel_y},{ang_vel_z}\n")

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
    p.setGravity(0, 0, -9.81)
    p.loadURDF("plane.urdf")

    # --- Robot body ---
    visual_shape = p.createVisualShape(
        p.GEOM_BOX, halfExtents=[0.25, 0.15, 0.05],
        rgbaColor=[0.3, 0.6, 1.0, 1.0]
    )
    collision_shape = p.createCollisionShape(p.GEOM_BOX, halfExtents=[0.25, 0.15, 0.05])
    cube = p.createMultiBody(
        baseMass=1.0,
        baseCollisionShapeIndex=collision_shape,
        baseVisualShapeIndex=visual_shape,
        basePosition=[0, 0, 0.5]
    )

    dt = 1 / 240
    initial_pos, initial_orn = p.getBasePositionAndOrientation(cube)

    results_dir = Path("results"); results_dir.mkdir(exist_ok=True)
    results_file = results_dir / "robot_results.csv"

    if record:
        os.makedirs("videos", exist_ok=True)
        vid_name = f"videos/robot_reset_{int(time.time())}.mp4"
        log_id = p.startStateLogging(p.STATE_LOGGING_VIDEO_MP4, vid_name)
        print(f"🎥 Recording video to {vid_name}")

    seq = random_seq(60, spread_deg=5)
    lam, R, th = estimate_lambda_and_R(seq)
    print(f"λ={lam:.3f}   R={R:.3f}   θ_net={np.degrees(th):.2f}°")
    reset_log_file.write(f"{lam},{R},{th}\n")

    t = 0.0
    debug_text = None

    # --- Phase 1: Disturbance torques ---
    for n, th_k in seq:
        torque = 2.0 * np.array(n)
        p.applyExternalTorque(cube, -1, torque, p.WORLD_FRAME)
        p.addUserDebugLine([0,0,0.5], np.add([0,0,0.5], 0.3*torque), [1,0,0], lineWidth=3, lifeTime=0.05)
        log_imu_data(t, cube)
        if use_gui:
            orbit_camera(t)
            if debug_text is not None: p.removeUserDebugItem(debug_text)
            debug_text = p.addUserDebugText("Phase 1: Disturbance", [0.8, 0, 0.7],
                                            textColorRGB=[1,1,1], textSize=1.5, lifeTime=0.1)
        p.stepSimulation(); t += dt
        if use_gui: time.sleep(0.001)

    pre_reset_t = t

    # --- Phase 2: λ-Scaled reset (x2) ---
    reset_seq = [(n, lam * thk) for (n, thk) in seq] * 2
    for n, th_k in reset_seq:
        torque = -2.5 * np.array(n) * np.sign(th_k)
        for _ in range(8):
            p.applyExternalTorque(cube, -1, torque, p.WORLD_FRAME)
            p.addUserDebugLine([0,0,0.5], np.add([0,0,0.5], 0.3*torque), [0,1,0], lineWidth=3, lifeTime=0.05)
            log_imu_data(t, cube)
            if use_gui:
                orbit_camera(t)
                if debug_text is not None: p.removeUserDebugItem(debug_text)
                debug_text = p.addUserDebugText("Phase 2: Reset", [0.8,0,0.7],
                                                textColorRGB=[1,1,1], textSize=1.5, lifeTime=0.1)
            p.stepSimulation(); t += dt
            if use_gui: time.sleep(0.001)

    post_reset_t = t
    final_pos, final_orn = p.getBasePositionAndOrientation(cube)

    if record:
        p.stopStateLogging(log_id)
        print(f"✅ Video saved: {vid_name}")

    imu_log_file.close(); reset_log_file.close()
    print("✅ Log files saved to logs/imu_log.csv and logs/reset_log.csv")

    diff_quat = p.getDifferenceQuaternion(initial_orn, final_orn)
    _, residual_angle_rad = p.getAxisAngleFromQuaternion(diff_quat)
    residual_angle_deg = np.degrees(residual_angle_rad)
    recovery_time = post_reset_t - pre_reset_t

    results_df = pd.DataFrame([{
        "R": lam,
        "residual": residual_angle_deg,
        "t_recover": recovery_time,
        "domain": "gravity"
    }])

    results_dir = Path("results")
    results_dir.mkdir(exist_ok=True)
    results_file = results_dir / "robot_results.csv"

    if results_file.exists():
        results_df.to_csv(results_file, mode="a", header=False, index=False)
        print(f"✅ Appended new results to {results_file}")
    else:
        results_df.to_csv(results_file, index=False)
        print(f"✅ Created new results file: {results_file}")
