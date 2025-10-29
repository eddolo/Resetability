# create_fake_telemetry.py
import pandas as pd, numpy as np
from pathlib import Path

out = Path("data"); out.mkdir(exist_ok=True)
N = 1000
t = np.linspace(0, 100, N)
# simple oscillatory quaternion-like data
qw = np.cos(t * 0.01)
qx = np.sin(t * 0.01)
qy = np.sin(t * 0.015)
qz = np.sin(t * 0.005)
wx = 0.05 * np.sin(t * 0.1)
wy = 0.04 * np.cos(t * 0.1)
wz = 0.03 * np.sin(t * 0.08)

df = pd.DataFrame({
    "timestamp": t,
    "qw": qw, "qx": qx, "qy": qy, "qz": qz,
    "wx": wx, "wy": wy, "wz": wz
})
df.to_csv("data/telemetry.csv", index=False)
print("✅ Fake telemetry generated → data/telemetry.csv")
