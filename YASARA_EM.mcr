# YASARA MACRO
# TOPIC:       5. Structure prediction
# TITLE:       Building a homology model
# REQUIRES:    Structure
# AUTHOR:      Elmar Krieger
# LICENSE:     GPL
# DESCRIPTION: This macro builds a homology model using a FASTA sequence of the target, and optionally template structures and alignments. It takes some shortcuts (alignment quality, number of templates, refinement) to finish as quickly as possible

# Parameter section - adjust as needed, but NOTE that some changes only take effect
# if you start an entirely new modeling job, not if you continue an existing one. 
# =================================================================================

# You can either set the target structure by clicking on Options > Macro > Set target,
# by providing it as command line argument (see docs at Essentials > The command line),
# or by uncommenting the line below and specifying it directly.
MacroTarget='/Users/rawir/Documents/work/vrc/projects/trimer_design/ADROIT2.0/data/tmp/modelling/A.BG505.W6M.C2/A.BG505.W6M.C2_YASARA_prot_superimposed_VRC01'

# Load structure
LoadPDB '(MacroTarget).pdb'

# Prepare the structure for simulation at the chosen pH
CleanAll

# Optimize the hydrogen-bonding network
OptHydAll

# Add Caps
AddCap ACE
AddCap NH2


# Perform energy minimization
Experiment Minimization
Experiment On
Wait ExpEnd

# Save PDB
SavePDB 1,'(MacroTarget)_minimized.pdb',Format=PDB


# Exit YASARA if this macro was provided as command line argument in console mode and not included from another macro
if runWithMacro and ConsoleMode and !IndentationLevel
  Exit
if !ConsoleMode
  # Show homology modeling report
  ShowURL file://(MacroTarget).html
