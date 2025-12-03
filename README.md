**"Thermodynamics of consciousness: A non-invasive perturbational framework"**

*Authors: Tomas Berjaga-Buisan, Juan Manuel Monti, Martina Cortada, Michele A. Colombo, Sebastian M. Geli, Gianluca Gaglioti, Simone Sarasso, Morten L. Kringelbach, Maurizio Corbetta, Maria V. Sanchez-Vives, Marcello Massimini, Yonatan Sanz Perl, Gustavo Deco*

### 🧠 Overview

This repository accompanies the manuscript. Researchers are encouraged to contact the corresponding author for clarification or assistance in applying the methods. Please cite if you use this code:

> **Berjaga-Buisan, T., et al. (2025). Thermodynamics of consciousness: A non-invasive perturbational framework. [Preprint / DOI link to be added]**

The quest for reliable and objective measures of consciousness is critical in basic and clinical neuroscience. Across species, the Perturbational Complexity Index (PCI) has emerged as a robust empirical marker by directly perturbing the brain, yet its underlying **principles of physics** remain unclear. Here, we bridge this gap by introducing a non-invasive framework based on **generative whole-brain models of non-equilibrium brain dynamics**. Using these models, we identified **violations of the Fluctuation-Dissipation Theorem (FDT) in humans and rodents across wakefulness, anesthesia, and disorders of consciousness.**  

Mirroring PCI, we found **decreased FDT violations in unresponsive disorders of consciousness and anesthesia compared to conscious conditions**. This reveals a close link between PCI and non-equilibrium dynamics in spontaneous brain signals, **grounding PCI in fundamental principles of physics.**

Overall, this framework offers complementary, non-invasive, model-based avenues for understanding consciousness and for developing objective tools to assess its loss and recovery in health and disease. 

## 🧰 Requirements

- MATLAB R2022b or later
- Required Toolboxes: Signal Processing, Statistics, and Parallel Computing
- Computation setup:
  
                - Langevin model fitting and supplementary thermodynamic metrics (Entropy Production, Asymmetry, Irreprocity) can be computed on a standard personal computer.
  
                - FDT violation analysis requires HPC computation with SLURM due to the large-scale matrix operations and heavy computational load.
  
                - A SLURM script is provided to handle random initializations and job management.
  
                - A local (lightweight) version of the scripts is included for demonstration only — it is not suitable for real analyses.
  
                - Synthetic data are provided to allow testing of the pipeline without access to empirical datasets.
  
                - (Optional): Python ≥ 3.10 with NumPy, SciPy, and Matplotlib for auxiliary visualization and analysis scripts.
  
The FDT violations simulation pipeline runs on the institutional HPC cluster of MELIS, using Slurm for job scheduling and system-specific MPI and library modules.
Due to these dependencies, reproducing the complete run outside this environment is not straightforward. Researchers are encouraged to contact the corresponding author for clarification or assistance in applying the methods.

**NOTE**: A faster formulation of the FDT violation metric is currently under development to enable faster computation on standard personal computers. Please contact the corresponding authors before using the code, as this implementation may be available upon request at the time of contact. All results from the paper were obtained using the numerical implementation included here.

 **DISCLAIMER:** USE THE ORDER STATED IN THE REPOSITORY (GEC --> Analytical --> FDT Violations): THESE ARE THE FIRST STEP REQUIRED TO ENSURE EVERYTHING WORKS CORRECTLY!

All code and analysis pipelines were created by **Tomas Berjaga-Buisan & Juan Manuel-Monti**.

This work was conducted within the Computational Neuroscience Research Group at Universitat Pompeu Fabra (Barcelona, Spain) in collaboration with CONICET (Universidad Nacional de Rosario & Universidad de San Andrés) , Institut d’Investigacions Biomèdiques August Pi i Sunyer (IDIBAPS), Università degli Studi di Milano, University of Oxford, Aarhus University, University of Padova, IRCCS and Institució Catalana de Recerca i Estudis Avançats (ICREA)

Supervised by Prof. Dr. Gustavo Deco.


### 📫 Contact

**Tomas Berjaga-Buisan**

PhD Student – Computational Neuroscience Research Group

Universitat Pompeu Fabra (Barcelona, Spain)

📧 Email: tomas.berjaga@upf.edu

🌐 GitHub: https://github.com/tomasberjagabuisan

