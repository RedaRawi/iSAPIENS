# iSAPIENS
An in silico tool to optimize protein interaction energies

## Installation

iSAPIENS has been successfully tested on Linux systems


### Environment
- Install R and required libraries:
  - R
  - R libraries
    - Within R execute the following command: install.packages( "bio3d" )
    - - Within R execute the following command: install.packages( "doParallel" )
    - Download R's Interpol package (https://cran.r-project.org/src/contrib/Archive/Interpol/Interpol_1.3.1.tar.gz) and install manually:
      Within R terminal execute: install.packages("Interpol_1.3.1.tar.gz", source = TRUE, repos = NULL )
  - FoldX (https://foldxsuite.crg.eu/)
  - VMD (https://www.ks.uiuc.edu/Research/vmd/)
    - We used VMD version "VMD for LINUXAMD64, version 1.9.3"
  - NAMD (https://www.ks.uiuc.edu/Research/namd/)
    - We used NAMD version "NAMD_2.13_Linux-x86_64-multicore"
    - Download required NAMD CHARMM parameter and topology files; in particular
      -  top_all36_prot.rtf
      -  top_all36_carb.rtf
      -  toppar_water_ions_namd.str
      -  par_all36_prot_cu_patch.prm
      -  par_all36_carb.prm
      -  toppar_water_ions_namd.str
  -  YASARA (https://www.yasara.org/)


## Run
iSAPIENS can be run in the command line

### Multiple input arguments are required to run iSAPIENS
- Antigen - Antibody complex
  - Antigen (chain A)
  - Antibody (chain H and chain L)
- Number of CPUs
- Path to FoldX executable
- Path to VMD executable
- Path to NAMD executable
- NAMD topology files separated by "___"
- NAMD parameter files separated by "___"
- Path to YASARA executable
- Path to YASARA homology modeling script (attached here: YASARA_EM.mcr)
- Cutoff (We used 12Å)
- Option (1: Interface, 2: Stability, 3: Interface + Stability)
- Interface Chain(s) 1: We used "A"
- Interface Chain(s) 2: We used "H_L"
- Stability Chain(s) 1: We used "H_L"
- Stability Chain(s) 2: We used "H_L"

### Execute in command line

Need to generate mutation file using FoldX syntax (see attached example: "individual_list.txt")

- Option 1: Antigen-antibody interface protein interaction energies
  - R --no-save < iSAPIENS_v1.0.0.R File.pdb nCPU PathToFoldx PathToVMD PathToNAMD top_all36_prot.rtf___top_all36_carb.rtf___toppar_water_ions_namd.str par_all36_prot_cu_patch.prm___par_all36_carb.prm___toppar_water_ions_namd.str PathToYASARA YASARA_EM.mcr 12 1 A H_L
- Option 2: Antibody stability and inter-chain interaction energies
  - R --no-save < iSAPIENS_v1.0.0.R File.pdb nCPU PathToFoldx PathToVMD PathToNAMD top_all36_prot.rtf___top_all36_carb.rtf___toppar_water_ions_namd.str par_all36_prot_cu_patch.prm___par_all36_carb.prm___toppar_water_ions_namd.str PathToYASARA YASARA_EM.mcr 12 2 H_L H_L
- Option 3: Option 1 and 2
  - R --no-save < iSAPIENS_v1.0.0.R File.pdb nCPU PathToFoldx PathToVMD PathToNAMD top_all36_prot.rtf___top_all36_carb.rtf___toppar_water_ions_namd.str par_all36_prot_cu_patch.prm___par_all36_carb.prm___toppar_water_ions_namd.str PathToYASARA YASARA_EM.mcr 12 3 A H_L H_L H_L


 Examples (use files provided in repository):
 
- R --no-save < iSAPIENS_v1.0.0.R m42_HH28K_13_renum.pdb 48 PathToFoldx PathToVMD PathToNAMD top_all36_prot.rtf___top_all36_carb.rtf___toppar_water_ions_namd.str par_all36_prot_cu_patch.prm___par_all36_carb.prm___toppar_water_ions_namd.str PathToYASARA YASARA_EM.mcr 12 1 A H_L
- R --no-save < iSAPIENS_v1.0.0.R m42_HH28K_13_renum.pdb 48 PathToFoldx PathToVMD PathToNAMD top_all36_prot.rtf___top_all36_carb.rtf___toppar_water_ions_namd.str par_all36_prot_cu_patch.prm___par_all36_carb.prm___toppar_water_ions_namd.str PathToYASARA YASARA_EM.mcr 12 2 H_L H_L
- R --no-save < iSAPIENS_v1.0.0.R m42_HH28K_13_renum.pdb 48 PathToFoldx PathToVMD PathToNAMD top_all36_prot.rtf___top_all36_carb.rtf___toppar_water_ions_namd.str par_all36_prot_cu_patch.prm___par_all36_carb.prm___toppar_water_ions_namd.str PathToYASARA YASARA_EM.mcr 12 3 A H_L H_L H_L


# Reference
Mateo Reveiz et al.; In Silico Improvement of Highly Protective Anti-Malarial Antibodies, bioRxiv, 2022, https://www.biorxiv.org/content/10.1101/2022.04.08.487687v1

# Contact
Reda Rawi: reda.rawi@nih.gov

