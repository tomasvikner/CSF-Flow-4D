# Semi-Automatic Analysis of 4D CSF Flow MRI Data

MATLAB-based GUI tool for **interactive and semi-automatic analysis of cerebrospinal fluid (CSF) flow** using **4D flow MRI** data.  
Allows visualization, local segmentation, and flow direction estimation in CSF spaces using either **automated (PCA-based)** or **manual** specification of flow-direction.

---

- **GUI-based workflow** for 4D CSF flow MRI analysis
- Supports **NIfTI (.nii)** image loading (preliminary version)
- Two operating modes:
  - **Full-volume mode:** Loads complete 4D flow MRI data (CSF_GUI.m; large memory footprint, not up-to-date)
  - **Brain-masked mode:** Loads only flow velocities inside a predefined brain mask (CSF_GUI_brain.m) — faster and more memory-efficient
- **Flow direction estimation**
  - Click on any middle-row panel (ax(2,1), ax(2,2), or ax(2,3)) to select a voxel  
  - The tool computes the **principal component direction (PCA)** of local flow to estimate the main CSF flow direction  
  - If the PCA-derived direction appears anatomically inconsistent, switch to **Manual Direction Mode** and define the direction interactively
- **Local segmentation**
  - Combines **T2 CUBE intensity** (CSF-sensitive) with **CSF velocity standard deviation (SD)** to facilitate local CSF segmentation 
  - Adjustable **local threshold** and **plane size** via dropdown controls for fine-tuning local segmentation
- **Interactive visualization**
  - Multi-panel display with linked zooming and panning between key views  
  - Displays the selected local plane and segmentation overlay in real-time
- **Loading and saving**
  - By default, files are **loaded from and written to remote servers** unless a local path is explicitly specified in the code or GUI
 
---
vikner@wisc.edu 
