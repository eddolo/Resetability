# fix_booster_summary.py
import pandas as pd
from pathlib import Path

space = Path("results/booster_full/logs/summary.csv")
if not space.exists():
    raise SystemExit("Missing summary.csv – run booster_reset_demo.py first")

df = pd.read_csv(space)
# Rename columns so validate_resetability_cross_domain.py recognizes them
df = df.rename(columns={
    "R_active": "R",
    "residual_reset": "residual",
    "t_to_thr_reset": "t_recover"
})
df["domain"] = "spacecraft"
df.to_csv(space, index=False)
print("✅ Normalized spacecraft CSV saved:", space)
