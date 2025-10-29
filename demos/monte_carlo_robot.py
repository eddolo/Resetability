# demos/monte_carlo_robot.py
# Monte Carlo simulation for robot reset
# Runs the robot reset simulation in a headless batch mode to collect data.

import pybullet as p
import pybullet_data
import numpy as np
import time
import sys
import os

# --- PATH FIX ---
# Add the 'python' directory to the system path to find the so3_reset module
sys.path.append(os.path.abspath('python'))
# --- END PATH FIX ---

from so3_reset import estimate_lambda_and_R

def random_seq(N=50, spread_deg=5):
    """Generates a random sequence of rotations."""
    ths = np.deg2rad(np.random.randn(N) * spread_deg)
    ns = np.random.randn(N, 3)
    ns /= np.linalg.norm(ns, axis=1, keepdims=True)
    return list(zip(ns, ths))

def run_single_simulation():
    """Runs one full headless simulation and returns the key results."""
    p.resetSimulation()
    p.setGravity(0, 0, -9.81)
    p.loadURDF("plane.urdf")
    cube = p.loadURDF("r2d2.urdf", [0, 0, 0.5])
    dt = 1 / 240
    p.setTimeStep(dt)

    # Estimate reset parameters
    seq = random_seq(60, spread_deg=5)
    lam, R, th = estimate_lambda_and_R(seq)

    # --- Simulation Phases ---
    # Phase 1: small wobble
    for n, th_k in seq:
        p.applyExternalTorque(cube, -1, 0.2 * np.array(n), p.WORLD_FRAME)
        p.stepSimulation()

    # Phase 2: external shove
    for _ in range(60):
        p.applyExternalTorque(cube, -1, [5, 0, 0], p.WORLD_FRAME)
        p.stepSimulation()

    # Phase 3: reset sequence
    reset_seq = [(n, lam * thk) for (n, thk) in seq] * 2
    for n, th_k in reset_seq:
        torque = -0.4 * np.array(n) * np.sign(th_k)
        for _ in range(10):
            p.applyExternalTorque(cube, -1, torque, p.WORLD_FRAME)
            p.stepSimulation()
            
    return lam, R, th

if __name__ == "__main__":
    num_repetitions = 100
    results_filename = "logs/robot_results.csv"

    print("Starting Monte Carlo simulation for robot reset...")
    print(f"Running {num_repetitions} repetitions...")

    # Setup environment
    p.connect(p.DIRECT)  # Headless mode
    p.setAdditionalSearchPath(pybullet_data.getDataPath())

    # Create log directory and prepare results file
    os.makedirs("logs", exist_ok=True)
    with open(results_filename, "w") as f:
        f.write("run,lambda,R,theta_net_deg\n")

    start_time = time.time()
    for i in range(num_repetitions):
        lam, R, th = run_single_simulation()
        theta_net_deg = np.degrees(th)
        
        # Append results to the CSV
        with open(results_filename, "a") as f:
            f.write(f"{i+1},{lam},{R},{theta_net_deg}\n")
            
        # Print progress
        print(f"  Completed run {i+1}/{num_repetitions} -> lambda={lam:.2f}, R={R:.3f}, theta_net={theta_net_deg:.1f} deg")

    p.disconnect()
    end_time = time.time()
    
    print("-" * 30)
    print(f"✅ Monte Carlo simulation complete.")
    print(f"Total time: {end_time - start_time:.2f} seconds.")
    print(f"Aggregated results saved to {results_filename}")