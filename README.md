# Snail Morphometrics – Anopolis Lafka Ori  

This project analyzes morphometric traits of snails collected across altitudinal gradients in Anopolis, Lafka Ori.  
The pipeline summarizes trait distributions, performs statistical analyses, and visualizes patterns including PCA.  

---

## 📂 Data  
- **Input:** `../data/snails.xlsx`  
  - **Sheet:** `RAW_data`  
  - **Columns include:**  
    - `Individual_ID`  
    - `Station`  
    - `ALT` (altitude)  
    - Morphometric measurements (e.g., Shell height, Shell diameter, Aperture measures, etc.)  

---

## 🔎 Workflow  

### 1. Pre-processing  
- Calculates new ratio: `sh_sd = Shell_height / Shell_diameter`.  
- Removes empty/unwanted columns.  
- Converts dataset to **long format** for easier statistics and plotting.  
- Bins altitude every **200 m** into categorical `alt_bin` (e.g., `400–600`).  

### 2. Summary statistics  
- Mean and standard error per variable × altitude bin.  
- **Output:** `../results/variable_summary.tsv`  

### 3. ANOVA & Tukey tests  
- One-way ANOVA of each morphometric variable vs. altitude bin.  
- Tukey HSD post-hoc comparisons.  
- Results written to:  
  - `../results/stats_results_anova.tsv`  
  - `../results/tukey_results.tsv`  
  - Significant results only → `../results/tukey_results_sig.tsv`  
- Letter groupings (`a`, `b`, `c` …) generated for significance annotation.  

### 4. Boxplots  
- Per trait × altitude bin:  
  - Boxplot with jittered raw points.  
  - Tukey letters placed above boxes.  
- **Saved to:** `../figures/snails_<variable>_boxplot.png`  

### 5. Principal Component Analysis (PCA)  
- Uses only **morphometric numeric variables** (excludes ID, station, altitude).  
- Standardized before PCA.  
- **Outputs:**  
  - PCA scatterplot of specimens colored by altitude bin → `../figures/snails_pca.png`  
  - Scree plot (variance explained by PCs) → `../figures/pca_scree.png`  
  - Biplot (individuals + loadings, ellipses per altitude bin) → `../figures/pca_biplot.png`  

---

## 📦 Dependencies  
The script requires the following R packages:  

- `tidyverse`  
- `readxl`  
- `multcompView`  
- `factoextra`  
- `broom`  

---

## ▶️ Run  
Make the script executable and run from terminal:  

```bash
chmod +x analysis_snails.R
./analysis_snails.R
