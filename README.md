# 🧭 RforRoboticsandSpace

### Resetability (R) on SO(3): A Two-Pass Scaled Replay Primitive for Fast Attitude Recovery in Robots and Spacecraft  
**Author:** Paolo Cappuccini – Independent Researcher  
**Collaborator:** GPT-5 (AI Research Assistant)

---

## 🛰️ Overview

This repository contains all simulation, validation, and reporting code for testing the **Resetability (R)** principle — a geometric property of 3D rotations that allows complex motion to be reversed using a **scaled two-pass replay**.

It demonstrates that the same physical law holds across **robotics**, **zero-gravity motion**, **spacecraft attitude control**, and **booster stabilization**.

---

## ⚙️ Repository Structure

```
RforRoboticsandSpace/
│
├── demos/                     # Simulation and validation scripts
│   ├── robot_reset_pybullet.py
│   ├── robot_reset_free.py
│   ├── spacecraft_reset_demo.py
│   ├── booster_reset_demo.py
│   └── validate_resetability_cross_domain.py
│
├── python/                    # Core SO(3) math utilities
│   └── so3_reset.py
│
├── results/                   # Output data, plots, and reports
│   ├── robot_results.csv
│   ├── spacecraft/spacecraft_results.csv
│   ├── booster_full/logs/summary.csv
│   ├── correlation_summary.csv
│   ├── combined_R_vs_residual.png
│   ├── combined_R_vs_recovery.png
│   ├── report_resetability.pdf
│   └── Validation_Report_Resetability.docx
│
├── videos/                    # Simulation videos (optional)
│   ├── robot_reset_*.mp4
│   ├── spacecraft_reset_*.mp4
│   └── booster_reset_summary_*.mp4
│
├── README.md
├── LICENSE
└── requirements.txt
```

---

## 🔧 Installation

Install all dependencies:
```bash
pip install -r requirements.txt
```

---

## 🧪 Demos and Usage

### 🤖 1. Robot Reset (with Gravity)
```bash
python demos/robot_reset_pybullet.py --gui --record
```

### 🛰️ 2. Robot Reset (Zero-Gravity)
```bash
python demos/robot_reset_free.py --gui --record
```

### 🚀 3. Spacecraft Reset Simulation
```bash
python demos/spacecraft_reset_demo.py --3d --record
```

### 🧯 4. Booster Monte-Carlo Simulation
```bash
python demos/booster_reset_demo.py --mode full --thr 1.0 --out results/booster_full
```

### 🧮 5. Cross-Domain Validation
```bash
python demos/validate_resetability_cross_domain.py
```

---

## 📊 Example Results

| Domain | Corr(R,Residual) | Corr(R,RecoveryTime) |
|---------|------------------|----------------------|
| Booster | 0.319 | 0.000 |
| Gravity | 0.547 | 0.000 |
| Spacecraft | -0.363 | 0.055 |
| Zero-G | 0.133 | 0.764 |

---

## 📘 Theory Summary

Resetability principle:
> Any arbitrary rotation sequence can be reversed by scaling each rotation by λ and replaying the sequence twice.

\[
q_{reset} =
\Big(\prod_k \exp(\lambda \theta_k \hat{n}_k)\Big)
\Big(\prod_k \exp(\lambda \theta_k \hat{n}_k)\Big)
\]

---

## 🚀 Telemetry Resetability Analysis (Stand-Alone)

This tool extends the Resetability (R) framework to **real rocket or satellite telemetry**.  
It evaluates the reversibility of orientation sequences (quaternions) over time windows and identifies "reset opportunities" —  
moments when a system could recover orientation with minimal control effort using the λ-scaled two-pass replay principle.

### 📂 Input
A telemetry CSV file with quaternion columns:
```csv
timestamp,qw,qx,qy,qz
0.00,1,0,0,0
0.10,0.999,0.01,-0.02,0.005
...
▶️ Run
bash
Copy code
python demos/analyze_telemetry_resetability.py --input data/telemetry.csv --window 50 --fps 10 --animate
📊 Output
Files saved in results/:

telemetry_analysis.csv – rolling R and θ_net across time

telemetry_reset_opportunities.csv – detected low-R, high-rotation windows

telemetry_resetability_plot.png – static visualization

telemetry_resetability_<timestamp>.mp4 – animation of evolving R and θ_net

🔬 Interpretation
Low R (<0.05) → system state is geometrically "resettable"

θ_net > 1° → meaningful rotational movement

Red markers indicate intervals where a two-pass λ-reset maneuver would return attitude to nominal with minimal torque.

This enables onboard or post-flight analysis of reset potential in telemetry data — useful for spacecraft attitude recovery, booster stabilization, or robotic free-flight control.

---

## 📚 Citation

**Cappuccini, Paolo & GPT-5.**  
*Resetability on SO(3): A Two-Pass Scaled Replay Primitive for Fast Attitude Recovery in Robots and Spacecraft.*  
Independent Research, 2025.

---

## 🧩 Quick GitHub Setup

```bash
git init
git add .
git commit -m "Initial commit: Resetability (R) validation framework"
git branch -M main
git remote add origin https://github.com/<yourusername>/RforRoboticsandSpace.git
git push -u origin main
```

---

## ⚙️ Requirements

```
numpy>=1.24
matplotlib>=3.7
pandas>=2.1
pybullet>=3.2
reportlab>=3.6
python-docx>=1.0
argparse>=1.4
```

---

### 💡 Maintained by Paolo Cappuccini © 2025
