# GCNEB
Grand-Canonical Nudged Elastic Band 

Approach to run grand-canonical DFT calculations of transition states to obtain activation barriers as a function of applied bias at a solvated surface interface. Useful for obtaining reaction kinetics of electrocatalysts. 

Dependencies:
1. JDFTx: https://jdftx.org/CompilingBasic.html 
2. Python 3
3. ASE: https://wiki.fysik.dtu.dk/ase/install.html


Setup: 
1. Clone Github Repo to computer
2. Setup & source .jdftx.bashrc in home directory based on TEMPLATE.jdftx.bashrc
3. Add GBRV_v1.5/ folder to jdftx/pseudopotentials/ (provided in pseudos)
4. Follow tutorials to run molecule calculations
5. Follow tutorial to run bulk Au calculation
6. Follow tutorial to run solvated / biased Au(111) surface calculation
7. Follow tutorials to run solvated / biased Au(111) + CO adsorbate calculations
8. Follow GCNEB tutorial starting from converged initial and final state calculations (step 7)
9. Study electrochemical systems with state-of-the-art methods and power level > 9000

If you use these scripts in your own work, please cite: 

https://pubs.acs.org/doi/full/10.1021/jacs.2c03661

Singstock, N.; Musgrave, C. How the Bio-Inspired Fe2Mo6S8 Chevrel Breaks Electrocatalytic Nitrogen Reduction Scaling Relations. JACS. 2022.
