# Generalized Variable-Weighted Adjacency Framework for Fuzzy Spectral Clustering

This repository provides fully reproducible code for the paper:

**A Generalized Variable-Weighted Adjacency Framework for Fuzzy Spectral Clustering**  
Jesse S. Ghashti, Warren Hare, John R.J. Thompson

---

## Requirements

- R (version ≥ 4.0 recommended)
- All required R packages are loaded inside `Functions.R`  
  (and will automatically install any missing packages with a prompt)

---

## Reproducing the Paper Results

Clone the repository and ensure all files are in the same directory.

Open R in the repository root and run

```r
# Section Analysis 1 (Simulated Datasets)
# NOTE: This script runs ONE of the six simulated datasets at a time
source("SA1-SampleAnalysis.R")

# Section Analysis 2 (Progressive Cluster Overlap)
source("SA2-Complete.R")

# Section Analysis 3 (Progressive Noise Dimensions)
source("SA3-Complete.R")

# Real-world datasets (29 methods, 10 UCI datasets)
source("RealData.R")

# Supplementary material experiments
source("SupplementaryAnalysis.R")
```

---

## Changing the Simulated Dataset in SA1

`SA1-SampleAnalysis.R` is configured to run one simulated dataset at a time.

To switch datasets, edit the dataset call (documented in the script, approximately line 89).

Example (current configuration – Spiral dataset):

```r
data <- genSpiralData(nPts = numPts, seed = 44869410 + iter)
```

Available alternatives:

```r
data <- genHypeData(nPts = numPts, seed = 44869410 + iter)      # Hypercube
data <- genGaussianData(nPts = numPts, seed = 44869410 + iter)  # Gaussian
data <- genWormData(nPts = numPts, seed = 44869410 + iter)      # Worm
data <- genRingData(nPts = numPts, seed = 44869410 + iter)      # Ring
data <- genSpiralData(nPts = numPts, seed = 44869410 + iter)    # Spiral
data <- genPieData(nPts = numPts, seed = 44869410 + iter)       # Pie
```

---

## Expected Outputs

Each experiment produces `.rds` files containing:

- Accuracy and FARI scores
- Runtime measurements
- Entropy and XB indices
- Periodic checkpoint files (saved every 20–50 iterations)

Output locations are defined directly within each script.


## Repository Structure

### Core Files (Required)

- `Functions.R`  
  Helper functions and required packages, including the FARI metric  
  (Andrews, Browne & Hvingelby, 2022; Andrews, 2013 GitHub)

- `datasets.rds`  
  Ten real datasets used in the analysis (from the UCI ML Repository)

- `SA1-Datasets.R`  
  Functions for generating the six simulated datasets used in SA1


### Method Implementations

- `VWKFC.R`  
  Viewpoint-Based Kernel Fuzzy Clustering  
  Tang et al., IEEE TETCI, 2023

- `FCMM.R`  
  Fuzzy C-Multiple-Means  
  Yang et al., IEEE GRSL, 2023

- `FSSc.R`  
  Fuzzy Subspace Clustering  
  Borgelt, 2009

- `FEc.R`  
  Fuzzy Entropy Clustering  
  Tran & Wagner, IEEE FUZZ, 2000

- `SNNs.R`  
  Shared Nearest Neighbors Spectral Clustering  
  Ye & Sakurai, IJCNN, 2015

- `LDAs.R`  
  Local Density Adaptive Similarity  
  Zhang et al., Pattern Recognition Letters, 2011


### Experiment Scripts

- `SA1-SampleAnalysis.R`  
  Reproduces Sensitivity Analysis 1 (simulated datasets)

- `SA2-Complete.R`  
  Reproduces Sensitivity Analysis 2 (progressive cluster overlap)

- `SA3-Complete.R`  
  Reproduces Sensitivity Analysis 3 (progressive noise dimensions)

- `RealData.R`  
  29 methods evaluated on 10 real datasets

- `SupplementaryAnalysis.R`  
  Supplementary experiments from the paper


## Notes on Reproducibility

- All experiments set fixed random seeds, see in-script notes on configurations.
- Scripts source all required dependencies internally
- Results are deterministic given the same R version and package set


## References

### Agreement Metrics and Evaluation

Andrews, J. L., Browne, R., & Hvingelby, C. D. (2022).  
**On assessments of agreement between fuzzy partitions.**  
*Journal of Classification*, 39, 326–342.  
https://doi.org/10.1007/s00357-021-09416-4

Andrews, J. L. (2013).  
**FARI: Fuzzy Adjusted Rand Index.**  
GitHub repository: https://github.com/its-likeli-jeff/FARI

### Fuzzy and Spectral Clustering Methods

Tang, Y., Pan, Z., Pedrycz, W., Ren, F., & Song, X. (2023).  
**Viewpoint-Based Kernel Fuzzy Clustering With Weight Information Granules.**  
*IEEE Transactions on Emerging Topics in Computational Intelligence*, 7(2), 342–356.  
https://doi.org/10.1109/TETCI.2022.3201620

Yang, X., Zhu, M., Sun, B., Wang, Z., & Nie, F. (2023).  
**Fuzzy C-Multiple-Means Clustering for Hyperspectral Image.**  
*IEEE Geoscience and Remote Sensing Letters*, 20, 1–5, Article 5503205.  
https://doi.org/10.1109/LGRS.2023.3246633

Borgelt, C. (2009).  
**Fuzzy Subspace Clustering.**  
In A. Fink, B. Lausen, W. Seidel, & A. Ultsch (Eds.),  
*Advances in Data Analysis, Data Handling and Business Intelligence*  
(Studies in Classification, Data Analysis, and Knowledge Organization).  
Springer, Berlin, Heidelberg.  
https://doi.org/10.1007/978-3-642-01044-6_8

Tran, D., & Wagner, M. (2000).  
**Fuzzy entropy clustering.**  
In *Proceedings of the Ninth IEEE International Conference on Fuzzy Systems (FUZZ-IEEE 2000)*  
(pp. 152–157, Vol. 1). San Antonio, TX, USA.  
https://doi.org/10.1109/FUZZY.2000.838650

Ye, X., & Sakurai, T. (2015).  
**Spectral clustering using robust similarity measure based on closeness of shared nearest neighbors.**  
In *Proceedings of the International Joint Conference on Neural Networks (IJCNN 2015)*  
(pp. 1–8). Killarney, Ireland.  
https://doi.org/10.1109/IJCNN.2015.7280495

Zhang, X., Li, J., & Yu, H. (2011).  
**Local density adaptive similarity measurement for spectral clustering.**  
*Pattern Recognition Letters*, 32(2), 352–358.  
https://doi.org/10.1016/j.patrec.2010.09.014


### Data Sources

Dua, D., & Graff, C. (2019).  
**UCI Machine Learning Repository.**  
University of California, Irvine.  
https://archive.ics.uci.edu/ml


## Contact

For questions or issues, please contact:

**Jesse Ghashti**  
jesse.ghashti@ubc.ca


## License

This code is provided for research purposes only.  
Please cite the paper if you use this code in your work.
