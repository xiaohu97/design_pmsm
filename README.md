# Humanoid Robot Joint Motor / PMSM Simulation

[中文版 (Chinese Version)](README_zh.md)

**A comprehensive, open-source simulation toolkit for rapid robot joint motor analysis, finite element evaluation, and control system integration.**

This repository bridges the gap between electromagnetic design and control engineering, providing a **FOC-ready model** of an **SPMSM** (Surface Permanent Magnet Synchronous Motor) that can be easily integrated into **MATLAB/Simulink** or Python-based systems.

## 🎯 Use Cases
- **Humanoid joint / actuator pre-design**: Rapid evaluation, sizing, and electromagnetic performance estimation.
- **Control simulation**: Algorithm validation using a fast, analytical **dq model**.
- System-level design: Deriving **torque-speed** characteristics and efficiency maps for **robot joint motor** systems.

## ✨ Simulation Results
*(Example: Electromagnetic Field Density mapped via FEMM FEA)*
![Field Density](output_femm/field_density_rpm1000_iq15.png)

## 🚀 5-Minute Quick Start

**1. Setup Environment**
```powershell
conda create -n motor python=3.11 -y
conda activate motor
pip install -r requirements.txt
```

**2. Fast Analytical dq Model Simulation**
```powershell
python .\design_pmsm.py --save-csv --save-png --out .\output
# Or simply execute the quick script: .\run.ps1
```

**3. FEMM 2D Finite Element Analysis (Windows)**
- Install [FEMM 4.2](https://www.femm.info/wiki/Download).
- Register the COM server (Run CMD as Administrator):
  ```bat
  "C:\FEMM42\femm.exe" /regserver
  ```
- Run a basic FEA simulation:
  ```powershell
  python .\femm_spm_template.py --rpm-list "1000,1500" --iq-list "15,25" --out .\output_femm
  ```

---
## 📁 Repository Structure
```text
motor design/
  design_pmsm.py        ← Fast SPMSM simulation (dq model)
  femm_spm_template.py  ← Comprehensive 2D FEA template (FEMM)
  requirements.txt
  run.ps1 / run.bat
  README.md
  output/               ← dq model output
  output_femm/          ← FEMM output
```

## 🔍 Detailed FEMM Analysis Features

Use the `--analysis` parameter to select the analysis types to execute (multiple can be selected):

| Analysis Type | Description | Output Files |
|--------|------|----------|
| `basic` | Torque / Loss / Temp Rise Waveforms | `femm_waveforms_*.png/.csv`, `femm_summary.csv` |
| `field` | Flux Density Map & Lines | `field_density_*.png`, `field_lines_*.png` |
| `airgap` | Airgap Flux Density (Bn/Bt) | `airgap_B_*.png/.csv` |
| `cogging` | Cogging Torque | `cogging_torque.png/.csv` |
| `inductance` | Ld/Lq Inductance Parameters | `inductance.png/.csv` |
| `emap` | Efficiency Map | `efficiency_map.png/.csv` |
| `tncurve` | Torque-Speed / Power-Speed Curves | `torque_speed_curve.png` |
| `all` | Run all the above | All files |

### FEMM Running Examples

**Basic Waveform Analysis (Default):**
```powershell
python .\femm_spm_template.py --rpm-list "1000,1500,2000" --iq-list "15,25,35" --out .\output_femm
```

**Flux Density + Magnetic Lines + Airgap (Single Operating Point):**
```powershell
python .\femm_spm_template.py --analysis field airgap --rpm-list "1000" --iq-list "25" --out .\output_femm
```

**Cogging Torque Analysis:**
```powershell
python .\femm_spm_template.py --analysis cogging --cogging-steps 72 --out .\output_femm
```

**Inductance Analysis (Ld/Lq):**
```powershell
python .\femm_spm_template.py --analysis inductance --ind-current 1.0 --ind-steps 37 --out .\output_femm
```

**Efficiency Map + Torque-Speed Curve:**
```powershell
python .\femm_spm_template.py --analysis emap tncurve --emap-rpm "500,1000,1500,2000,2500,3000" --emap-iq "5,10,15,20,25,30,35" --emap-pts 12 --out .\output_femm
```

**Run All Analysis Types at Once:**
```powershell
python .\femm_spm_template.py --analysis all --rpm-list "1000" --iq-list "15" --points 36 --emap-rpm "500,1000,2000" --emap-iq "10,20,30" --emap-pts 6 --out .\output_femm
```

### Common Parameters
| Parameter | Default | Description |
|------|--------|------|
| `--analysis` | `basic` | Selected analysis items |
| `--rpm-list` | `1000,1500,2000` | Speed list (rpm) |
| `--iq-list` | `15,25,35` | Peak current list (A) |
| `--points` | `72` | Sampling points per electrical cycle |
| `--out` | `output_femm` | Output directory |
| `--show` | off | Show matplotlib figures interactively |
| `--show-femm` | off | Show FEMM GUI |
| `--emap-rpm` | `500,...,3000`| Speed list for Efficiency Map |
| `--emap-iq` | `5,...,35`  | Current list for Efficiency Map |
| `--emap-pts` | `12`       | Sampling points for Efficiency Map |
| `--cogging-steps`| `72`   | Steps for Cogging Torque |
| `--ind-steps`  | `37`     | Steps for Inductance Analysis |
| `--ind-current`| `1.0`    | Test current for Inductance (A) |

### Runtime Reference (Approximations)
| Analysis Type | Est. Time |
|--------|----------|
| basic (1 condition, 72 steps) | 2 ~ 5 mins |
| field + airgap | basic + 1 min |
| cogging (72 steps) | 2 ~ 5 mins |
| inductance (37 steps × 2 axes) | 3 ~ 8 mins |
| emap (42 conditions, 12 steps) | 15 ~ 40 mins |
| all (simplified parameters) | 10 ~ 20 mins |

> **Tip:** Reducing sampling points (e.g., `--points 12` or `--emap-pts 6`) significantly accelerates the simulation.

---

## 📊 Extensive Simulation Gallery

### 4.1 Fast dq Model Waveforms
![dq Model Waveforms](output/pmsm_waveforms.png)

### 4.2 Flux Density Map
Polar representation of |B| magnitude cross-section.
![Flux Density Map](output_femm/field_density_rpm1000_iq15.png)

### 4.3 Magnetic Field Lines
Equipotential (Az) lines manifesting as magnetic field line distribution.
![Field Lines Map](output_femm/field_lines_rpm1000_iq15.png)

### 4.4 Airgap Flux Density
Radial (Bn) and Tangential (Bt) components extracted along the airgap mid-line.
![Airgap Flux](output_femm/airgap_B_rpm1000_iq15.png)

### 4.5 Waveforms (Torque / Loss / Temperature Rise)
One electrical cycle displaying electromagnetic torque, various core/winding losses, and temperature rise.
![Waveforms](output_femm/femm_waveforms_rpm1000_iq15.png)

### 4.6 Cogging Torque
Torque ripple produced by rotor rotating across one slot pitch at zero current.
![Cogging Torque](output_femm/cogging_torque.png)

### 4.7 Inductance Parameters (Ld / Lq)
Variation of d-axis and q-axis inductances against rotor position (For SPMSM, Ld ≈ Lq).
![Inductance](output_femm/inductance.png)

### 4.8 Efficiency Map
Efficiency contours on the speed × torque plane.
![Efficiency Map](output_femm/efficiency_map.png)

### 4.9 Torque-Speed Curves
FEA max torque operating points alongside theoretical voltage constraint envelopes.
![Torque-Speed Curve](output_femm/torque_speed_curve.png)

---

## 💡 Technical Notes

- `design_pmsm.py` is highly suitable for rapid control law testing and parameter trend estimation.
- `femm_spm_template.py` offers closer-to-truth field-circuit co-simulation but serves as a generic template:
  - Slot geometries and winding arrangements are simplified.
  - Iron loss and magnet loss are based on semi-empirical models.
  - If used for rigorous engineering validation, please replace it with real geometry, material properties, and authentic B-H curve data.

## 🛠 Troubleshooting (FEMM)

If you encounter an `Invalid class string` or `femm.ActiveFEMM` error, the FEMM COM server is not properly registered. Run CMD as Administrator:
```bat
"C:\FEMM42\femm.exe" /regserver
```
Check if it succeeded via:
```bat
reg query HKCR\femm.ActiveFEMM
```

If PowerShell throws a parameter list parsing error, enclose the list elements in quotes:
```powershell
python .\femm_spm_template.py --rpm-list "1000,1500,2000" --iq-list "15,25,35"
```
