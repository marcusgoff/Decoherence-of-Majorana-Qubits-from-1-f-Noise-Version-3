# Decoherence of Majorana Qubits from 1/f Noise

Code for all numerical calculations in *A. Alase, M. C. Goffage, M. C. Cassidy, S. N. Coppersmith (2025), Decoherence of Majorana Qubits from 1/f Noise, arXiv preprint: arXiv:2506.22394.* 

# Implementation of QPP_Library for: Decoherence in Majorana Qubits by 1/f Noise

This repository contains the code used to generate all numerical results presented in our paper *A. Alase, M. C. Goffage, M. C. Cassidy, S. N. Coppersmith, Decoherence in Majorana Qubits by 1/f Noise (2025)*. This code uses the covariance matrix method to calculate the time evolution of the tetron qubit, composed of two Kitaev chains, in the presence of a two-level-fluctuator (TLF) and returns the probability of exciting a quasiparticle pair in a single Kitaev chain of the qubit. See the Methods section and Sec. 3 of the Supplementary Information of our paper for further details on the calculations.  

For a more in-depth overview of our numerical package, called **QPP_Library**, please see Appendix B of our earlier paper, *M. C. Goffage, A. Alase, M. C. Cassidy, S. N. Coppersmith, Leakage at zero temperature from changes in chemical potential in Majorana qubits, arXiv:2504.17485 (2025)*, which employs the same library as the submitted paper.  

Note that since our numerical results are statistical averages over multiple two-level-fluctuator (TLF) noise realisations, we expect minor differences in the plots between consecutive runs.  

---

## Environment and Dependencies
Requires MATLAB 2024a or more recent.  

---

## How to Run the Code
1. Click the green **"Code"** drop-down menu at the top of this page.  
2. Select **"Download ZIP"**.  
3. Unzip the contents into a local directory.  
4. Open **`Run_all_Figures.m`** in MATLAB. Alternatively, open "Code and Figures" and navigate to the directory corresponding to the figure of interest, e.g. **`Figure 2`**, and open the .m file in that directory, e.g. "Run_Figure2b.m". 
5. In MATLAB, go to the **Editor** tab, select the **Editor** tab and then select **Run**.  
6. If prompted to *Change Folder* to the current directory, choose **Change Folder** (highlighted option).  

---

## Program Details – Code Demo
**`Run_Code_Demo.m`** in the directory **`Code Demo`** runs a simple demo of our code calculating $P_{QPP}$ versus time for a 3 µm nanowire.

**`Run_Code_Demo.m`** run time (tested on a MacBook Pro): ~4 minutes.

---

## Program Details – Run All Figures
**`Run_all_Figures.m`** generates all numerical results and corresponding figures presented in the paper. The results are saved in the subdirectories corresponding to those figures.

**`Run_all_Figures.m`** run time: >72 hours.  

You may also run individual figures by opening the **`Code and Figures`** directory and then chosing the subdirectory corresponding to the figure or interest. Then open the figure run script of interest in MATLAB (for example, `Run_Figure2b.m`). You may then run the code for just that figure by selecting **Run** in MATLAB. Note that `Run_FigureS1.m` requires you to first execute `Run_AddFig.m` located at **`Code and Figures/figures for supplemental/Additional Figures/Run_AddFig.m`**.

---

## Output
The code and results pertaining to each figure are located in a subdirectory reserved for that figure. Please navigate to the appropriate subdirectory to view the results, e.g. **`Code and Figures/Figure 2`** contains the code and output generated for Figure 2 in the paper. The Matlab code generates .mat files and .csv files for the data generated and .fig, .png. and .svg files for the figures. 

Note that figure 2b as formatted for the paper was generated using a Mathematica file, which has been included in **`Code and Figures/Figure 2/figure 2B numerical data of rate versus Gamma`** for convenience, however the data is generated using the Matlab code in **`Code and Figures/Figure 2/Run_Figure2b.m`**.
---

## Troubleshooting
If you run into any issues or have any questions, please contact m.goffage@unsw.edu.au.  
