# =====================================================
# validate_resetability_cross_domain.py
# Empirical validation of Resetability (R) concept
# across gravity, zero-G, spacecraft, and booster domains.
# Author: Paolo Cappuccini
# =====================================================

import pandas as pd, matplotlib.pyplot as plt, numpy as np
from pathlib import Path
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet

root = Path(".")

# --- Canonical files (explicit) ---
robot_csv       = root / "results" / "robot_results.csv"
spacecraft_csv  = root / "results" / "spacecraft" / "spacecraft_results.csv"
booster_summary = root / "results" / "booster_full" / "logs" / "summary.csv"

frames = []

# ---------- Load ROBOT ----------
if robot_csv.exists():
    robot = pd.read_csv(robot_csv)
    # Normalize columns
    rename_map = {}
    for c in robot.columns:
        lc = c.strip().lower()
        if lc == "r_active":       rename_map[c] = "R"
        elif lc == "residual_error": rename_map[c] = "residual"
        elif lc == "t_to_thr":       rename_map[c] = "t_recover"
        elif lc == "environment":    rename_map[c] = "domain"  # if env column used
        elif lc == "domain":         rename_map[c] = "domain"
    robot = robot.rename(columns=rename_map)

    # If no domain provided, default to "gravity" (legacy)
    if "domain" not in robot.columns:
        robot["domain"] = "gravity"

    keep_cols = [c for c in ["R","residual","t_recover","domain"] if c in robot.columns]
    robot = robot[keep_cols]
    frames.append(robot)
    print(f"[INFO] Loaded robot entries: {len(robot)}")
else:
    print("[WARN] No robot_results.csv found.")

# ---------- Load SPACECRAFT ----------
if spacecraft_csv.exists():
    try:
        space = pd.read_csv(spacecraft_csv)
        # Ensure required columns exist
        needed = {"R","residual","t_recover","domain"}
        present = set(space.columns)
        # If missing, try to map
        if "R" not in present and "r" in {c.lower() for c in space.columns}:
            for c in space.columns:
                if c.lower() == "r": space = space.rename(columns={c:"R"})
        if "residual" not in present and "residual_error" in space.columns:
            space = space.rename(columns={"residual_error":"residual"})
        if "t_recover" not in present and "t_to_thr" in space.columns:
            space = space.rename(columns={"t_to_thr":"t_recover"})
        if "domain" not in space.columns:
            space["domain"] = "spacecraft"

        space = space[["R","residual","t_recover","domain"]]
        frames.append(space)
        print(f"[INFO] Loaded spacecraft entries: {len(space)}")
    except Exception as e:
        print(f"[WARN] Could not parse spacecraft CSV: {e}")
else:
    print("[WARN] No spacecraft_results.csv found at results/spacecraft/.")

# ---------- Load BOOSTER (summary.csv) ----------
if booster_summary.exists():
    try:
        booster = pd.read_csv(booster_summary, on_bad_lines="skip")
        # Normalize
        rename_map = {}
        for c in booster.columns:
            lc = c.strip().lower()
            if lc == "r_active":            rename_map[c] = "R"
            elif lc == "residual_reset":    rename_map[c] = "residual"
            elif lc == "t_to_thr_reset":    rename_map[c] = "t_recover"
        booster = booster.rename(columns=rename_map)

        # Keep only rows that have at least R and residual
        keep_cols = [c for c in ["R","residual","t_recover"] if c in booster.columns]
        booster = booster[keep_cols].copy()
        booster["domain"] = "booster"
        frames.append(booster)
        print(f"[INFO] Loaded booster entries: {len(booster)}")
    except Exception as e:
        print(f"[WARN] Booster file unreadable: {e}")
else:
    print("[WARN] No booster summary at results/booster_full/logs/summary.csv.")

# ---------- Auto-discover any extra *_results.csv (without re-adding known ones) ----------
known_files = {robot_csv.resolve(), spacecraft_csv.resolve()}
extra_files = []
for f in (root / "results").rglob("*_results.csv"):
    if f.resolve() not in known_files:
        extra_files.append(f)

for f in extra_files:
    try:
        extra = pd.read_csv(f)
        rename_map = {}
        for c in extra.columns:
            lc = c.strip().lower()
            if lc == "r_active":       rename_map[c] = "R"
            elif lc == "r":            rename_map[c] = "R"
            elif lc == "residual_error": rename_map[c] = "residual"
            elif lc == "residual":       rename_map[c] = "residual"
            elif lc == "t_to_thr":       rename_map[c] = "t_recover"
            elif lc == "t_recover":      rename_map[c] = "t_recover"
            elif lc == "environment":    rename_map[c] = "domain"
            elif lc == "domain":         rename_map[c] = "domain"
        extra = extra.rename(columns=rename_map)
        if "domain" not in extra.columns:
            extra["domain"] = f.stem.split("_")[0]
        extra = extra[[c for c in ["R","residual","t_recover","domain"] if c in extra.columns]]
        if not extra.empty:
            frames.append(extra)
            print(f"[INFO] Loaded extra results from: {f}  ({len(extra)})")
    except Exception as e:
        print(f"[WARN] Could not parse {f}: {e}")

# ---------- Combine ----------
if not frames:
    print("⚠️ No data found. Run simulations first.")
    raise SystemExit

df = pd.concat(frames, ignore_index=True)
# keep numeric
for c in ["R","residual","t_recover"]:
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")
# drop rows with no R or no residual
df = df.dropna(subset=["R","residual"], how="any")

# ---------- Correlations ----------
corr_summary = []

for dom, dsub in df.groupby("domain"):
    dsub = dsub.dropna(subset=["R","residual"])
    n = len(dsub)
    if n < 2:
        print(f"[WARN] Domain '{dom}' has <2 samples; skipping corr.")
        c_r_resid, c_r_t = 0.0, 0.0
    elif np.std(dsub["R"]) == 0 or np.std(dsub["residual"]) == 0:
        print(f"[WARN] Domain '{dom}' has constant R or residual; corr=0.")
        c_r_resid, c_r_t = 0.0, 0.0
    else:
        c_r_resid = float(np.corrcoef(dsub["R"], dsub["residual"])[0,1])
        # recovery time correlation (only if present and non-constant)
        if "t_recover" in dsub.columns and dsub["t_recover"].notna().sum() > 1 and np.std(dsub["t_recover"].dropna()) > 0:
            c_r_t = float(np.corrcoef(dsub["R"], dsub["t_recover"].fillna(0))[0,1])
        else:
            c_r_t = 0.0
            print(f"[INFO] Domain '{dom}' has no/constant t_recover → Corr(R,Recovery)=0.")

    corr_summary.append({
        "Domain": dom,
        "Corr(R,Residual)": c_r_resid,
        "Corr(R,RecoveryTime)": c_r_t
    })
    print(f"[INFO] {dom}: Corr(R,Residual)={c_r_resid:.3f}, Corr(R,Recovery)={c_r_t:.3f}")

corr_df = pd.DataFrame(corr_summary)
corr_path = root / "results" / "correlation_summary.csv"
corr_df.to_csv(corr_path, index=False)

# ---------- Plots ----------
plt.figure(figsize=(7,5))
for dom, dsub in df.groupby("domain"):
    plt.scatter(dsub["R"], dsub["residual"], alpha=0.6, label=dom)
plt.xlabel("Resetability (R)")
plt.ylabel("Residual orientation error [deg]")
plt.title("Residual vs Resetability (R)")
plt.legend()
plt.tight_layout()
plot_resid = root / "results" / "combined_R_vs_residual.png"
plt.savefig(plot_resid, dpi=150)
plt.close()

plt.figure(figsize=(7,5))
have_any_t = False
for dom, dsub in df.groupby("domain"):
    if "t_recover" in dsub.columns and dsub["t_recover"].notna().any():
        have_any_t = True
        plt.scatter(dsub["R"], dsub["t_recover"], alpha=0.6, label=dom)
plt.xlabel("Resetability (R)")
plt.ylabel("Recovery time [s]")
plt.title("Recovery Time vs Resetability (R)")
plt.legend()
plt.tight_layout()
plot_time = root / "results" / "combined_R_vs_recovery.png"
plt.savefig(plot_time, dpi=150)
plt.close()

# ---------- PDF ----------
pdf_path = root / "results" / "report_resetability.pdf"
doc = SimpleDocTemplate(str(pdf_path), pagesize=A4, rightMargin=36, leftMargin=36, topMargin=36, bottomMargin=36)
styles = getSampleStyleSheet()
story = []

story.append(Paragraph("<b>Validation of the Resetability (R) Concept</b>", styles["Title"]))
story.append(Spacer(1, 12))
story.append(Paragraph(
    "This report aggregates simulation outputs across domains (robotics: gravity & zero-G, spacecraft, booster) "
    "and evaluates how the Resetability metric (R) correlates with residual attitude error and recovery time.",
    styles["BodyText"]))
story.append(Spacer(1, 12))

table_data = [["Domain", "Corr(R,Residual)", "Corr(R,RecoveryTime)"]]
for _, r in corr_df.iterrows():
    cr = "" if pd.isna(r["Corr(R,RecoveryTime)"]) else f"{r['Corr(R,RecoveryTime)']:.3f}"
    table_data.append([r["Domain"], f"{r['Corr(R,Residual)']:.3f}", cr])

tbl = Table(table_data, hAlign='LEFT')
tbl.setStyle(TableStyle([
    ('BACKGROUND', (0,0), (-1,0), colors.grey),
    ('TEXTCOLOR', (0,0), (-1,0), colors.whitesmoke),
    ('ALIGN', (0,0), (-1,-1), 'CENTER'),
    ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
    ('GRID', (0,0), (-1,-1), 0.5, colors.black),
]))
story.append(tbl)
story.append(Spacer(1, 18))

story.append(Paragraph("<b>Residual vs R</b>", styles["Heading2"]))
story.append(Image(str(plot_resid), width=5.5*inch, height=3.5*inch))
story.append(Spacer(1, 12))

if have_any_t:
    story.append(Paragraph("<b>Recovery Time vs R</b>", styles["Heading2"]))
    story.append(Image(str(plot_time), width=5.5*inch, height=3.5*inch))
    story.append(Spacer(1, 12))

story.append(Paragraph(
    "Interpretation: gravity-bound runs often show strong R–residual correlation; zero-G/spacecraft emphasize "
    "timing sensitivity (λ-scaled sequences). Booster Monte Carlo provides population-level behavior.", styles["BodyText"]))
story.append(Spacer(1, 12))
story.append(Paragraph("Generated automatically by the R–SO(3) validation pipeline.", styles["Italic"]))

doc.build(story)

print(f"✅ PDF report generated: {pdf_path}")
print(f"📊 Correlation CSV saved: {corr_path}")
