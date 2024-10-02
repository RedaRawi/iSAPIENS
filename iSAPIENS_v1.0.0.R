#===============================================================================
#===============================================================================
#===============================================================================
# Libraries

library( bio3d )
library( doParallel )
library( Interpol )
data( "AAindex" )



#===============================================================================
#===============================================================================
#===============================================================================
# Functions

#===============================================================================
# Run FoldX to build the mutant model
run.FoldX <- function( exe.foldx,
                       file.pdb )
{
  command.buildmodel <- paste( exe.foldx, " ",
                               "--command=BuildModel ",
                               "--pdb=", file.pdb, " ",
                               "--mutant-file=individual_list.txt ",
                               "--numberOfRuns=1 ",
                               "> buildmodel.log",
                               sep = "" )
  system( command.buildmodel )  
  file.mutant <- paste( unlist( strsplit( file.pdb,
                                          "\\.pdb" ) ),
                        "_1.pdb",
                        sep = "" )
  
  return( file.mutant )
}
#===============================================================================

#===============================================================================
# YASARA energy minimization function
YASARA.EM <- function( wd,
                       exe.yasara,
                       file.yasara.em,
                       file.target.pdb )
{
  print( "==================================================" )
  print( "Run YASARA EM function" )
  
  #===============================================================================
  # Move to working directory
  setwd( wd )
  
  #===============================================================================
  # Copy template YASARA HM file to current directory
  file.yasara.em.local <- "YASARA_EM.mcr"
  system( paste( "cp",
                 file.yasara.em,
                 file.yasara.em.local ) )
  
  # Adjust file according to user inputs
  lines <- readLines( file.yasara.em.local )
  
  # Adjust target sequence
  lines[ 16 ] <- paste( "MacroTarget='",
                        wd, "/", unlist( strsplit( file.target.pdb,
                                                   "\\.pdb" ) ),
                        "'",
                        sep = "" )
  
  writeLines( lines,
              file.yasara.em.local )
  
  #===============================================================================
  # Start EM using YASARA
  system( paste( exe.yasara,
                 "-txt",
                 file.yasara.em.local,
                 "> EM.log" ) )
  
  #===============================================================================
  # Return target PDB file name
  return( paste( unlist( strsplit( file.target.pdb,
                                   "\\.pdb" ) ),
                 "_minimized.pdb",
                 sep = "" ) )
}
#===============================================================================

#===============================================================================
# Clean PDB
clean.PDB <- function( file.pdb,
                       type = "prot" )
{
  #===============================================================================
  # Define variables
  vec.aa1.raw <- unique( aa.table$aa1 )
  vec.aa1 <- vec.aa1.raw[ vec.aa1.raw != "X" ]
  file.pdb.clean <- paste( unlist( strsplit( file.pdb,
                                             "\\.pdb" ) ),
                           "_prot.pdb",
                           sep = "" )
  
  #===============================================================================
  # Load PDB
  pdb <- read.pdb( file.pdb )
  pdb.atom <- pdb$atom
  
  #===============================================================================
  # Clean PDB using chosen type
  if( type == "prot" )
  {
    # Extract PDB information for protein only
    df.chain.resno.insert <- data.frame( cbind( pdb.atom[ pdb.atom$elety == "CA", ]$chain,
                                                pdb.atom[ pdb.atom$elety == "CA", ]$resno,
                                                pdb.atom[ pdb.atom$elety == "CA", ]$insert ) )
    colnames( df.chain.resno.insert ) <- c( "chain",
                                            "resno",
                                            "insert" )
    df.chain.resno.insert$chain <- as.character( df.chain.resno.insert$chain )
    df.chain.resno.insert$resno <- as.numeric( as.character( df.chain.resno.insert$resno ) )
    df.chain.resno.insert$insert <- as.character( df.chain.resno.insert$insert )
    
    vec.ind.prot <- which( aa321( pdb.atom$resid ) %in% vec.aa1 )
    
    # Start saving protein only
    pdb.atom.clean <- pdb.atom[ vec.ind.prot, ]
    write.pdb( file = file.pdb.clean,
               type = pdb.atom.clean$type,
               eleno = pdb.atom.clean$eleno,
               elety = pdb.atom.clean$elety,
               alt = pdb.atom.clean$alt,
               resid = pdb.atom.clean$resid,
               chain = pdb.atom.clean$chain,
               resno = pdb.atom.clean$resno,
               insert = pdb.atom.clean$insert,
               xyz = as.vector( t( as.matrix( pdb.atom.clean[ , 9:11 ] ) ) ),
               o = pdb.atom.clean$o,
               b = pdb.atom.clean$b )
    
    return( file.pdb.clean )
  }
}
#===============================================================================

#===============================================================================
# Calculate residue-residue VdW and Elec energies
PaDRINO.energy <- function( file.pdb,
                            chain.1 = NA,
                            chain.2 = NA,
                            cutoff = 12 )
{
  #===============================================================================
  # Prepare PDB and PSF file for all pairwise comparisons
  vec.files.autopsf <- generate.PSF( file.pdb,
                                     exe.vmd,
                                     vec.files.parameters.prepare,
                                     "autopsf" )
  file.pdb.autopsf.pdb <- vec.files.autopsf[ 1 ]
  file.pdb.autopsf.psf <- vec.files.autopsf[ 2 ]
  
  # Load PDB
  pdb <- read.pdb( file.pdb.autopsf.pdb )
  pdb.atom <- pdb$atom
  
   #===============================================================================
  # Prepare residue-residue pairs
  if( is.na( chain.1 )[ 1 ] & is.na( chain.2 )[ 1 ] )
  {
    df.residue.1 <- data.frame( cbind( pdb.atom[ pdb.atom$elety == "CA", ]$chain,
                                       pdb.atom[ pdb.atom$elety == "CA", ]$resno ) )
    df.residue.2 <- data.frame( cbind( pdb.atom[ pdb.atom$elety == "CA", ]$chain,
                                       pdb.atom[ pdb.atom$elety == "CA", ]$resno ) )
  } else
  {
    df.residue.1 <- data.frame( cbind( pdb.atom[ pdb.atom$elety == "CA" &
                                                   pdb.atom$chain %in% chain.1, ]$chain,
                                       pdb.atom[ pdb.atom$elety == "CA" &
                                                   pdb.atom$chain %in% chain.1, ]$resno ) )
    df.residue.2 <- data.frame( cbind( pdb.atom[ pdb.atom$elety == "CA" &
                                                   pdb.atom$chain %in% chain.2, ]$chain,
                                       pdb.atom[ pdb.atom$elety == "CA" &
                                                   pdb.atom$chain %in% chain.2, ]$resno ) )
  }
  colnames( df.residue.1 ) <- c( "chain",
                                 "resno" )
  df.residue.1$chain <- as.character( df.residue.1$chain )
  df.residue.1$resno <- as.numeric( as.character( df.residue.1$resno ) )
  n <- nrow( df.residue.1 )
  
  colnames( df.residue.2 ) <- c( "chain",
                                 "resno" )
  df.residue.2$chain <- as.character( df.residue.2$chain )
  df.residue.2$resno <- as.numeric( as.character( df.residue.2$resno ) )
  p <- nrow( df.residue.2 )
  
  # Start running residue-residue calculations
  # system( "rm VdW* Elec* Total*" )
  
  print( paste( "Dimensions:", 
                n,
                "x",
                p ) )
  
  file.VdW <- paste( "VdW_",
                     paste( chain.1,
                            collapse = "" ),
                     "_",
                     paste( chain.2,
                            collapse = "" ),
                     ".txt",
                     sep = "" )
  file.Elec <- paste( "Elec_",
                      paste( chain.1,
                             collapse = "" ),
                      "_",
                      paste( chain.2,
                             collapse = "" ),
                      ".txt",
                      sep = "" )
  file.Total <- paste( "Total_",
                       paste( chain.1,
                              collapse = "" ),
                       "_",
                       paste( chain.2,
                              collapse = "" ),
                       ".txt",
                       sep = "" )
  
  for( i in 1:n )
  {
    print( paste( i,
                  n,
                  sep = "/" ) )
    
    #===============================================================================
    #===============================================================================
    # Calculate only for residue-residue pairs within cutoff
    pdb.atom.i <- pdb.atom[ pdb.atom$chain == df.residue.1[ i, 1 ] &
                              pdb.atom$resno == df.residue.1[ i, 2 ], 9:11 ]
    vec.j.dist <- foreach( j = 1:p, .combine = c ) %dopar%
      {
        pdb.atom.j <- pdb.atom[ pdb.atom$chain == df.residue.2[ j, 1 ] &
                                  pdb.atom$resno == df.residue.2[ j, 2 ], 9:11 ]
        min( dist.xyz( a = pdb.atom.i,
                       b = pdb.atom.j ) )
      }
    vec.ind.j <- which( vec.j.dist < cutoff )
    print( paste( "To analyse:",
                  length( vec.ind.j ) ) )
    #===============================================================================
    #===============================================================================
    
    # mat.energies <- foreach( j = 1:p, .combine = cbind ) %dopar%
    mat.energies.tmp <- foreach( j = vec.ind.j, .combine = cbind ) %dopar%
      {
        # Prepare pairwise PDB
        file.pdb.pair <- paste( unlist( strsplit( file.pdb.autopsf.pdb,
                                                  "\\.pdb" ) ),
                                "_", i, "_", j, ".pdb",
                                sep = "" )
        prepare.PDB.pairwise( pdb.atom,
                              file.pdb.pair,
                              df.residue.1[ i, ]$chain,
                              df.residue.1[ i, ]$resno,
                              df.residue.2[ j, ]$chain,
                              df.residue.2[ j, ]$resno )
        
        # Prepare PSF for one-step MD
        output.prefix <- unlist( strsplit( file.pdb.pair,
                                           "\\.pdb" ) )
        
        file.output.pair <- paste( unlist( strsplit( file.pdb.pair,
                                                     "\\.pdb" ) ),
                                   "_pairwiseMD",
                                   sep = "" )
        file.namd.pair <- paste( "md_", i, "_", j, ".conf",
                                 sep = "" )
        prepare.NAMD.config.file( file.pdb.autopsf.psf,
                                  file.pdb.autopsf.pdb,
                                  file.pdb.pair,
                                  vec.files.parameters.MD,
                                  file.output.pair,
                                  file.namd.pair,
                                  cutoff )
        file.trajectory.pair <- paste( unlist( strsplit( file.pdb.pair,
                                                         "\\.pdb" ) ),
                                       "_pairwiseMD.dcd",
                                       sep = "" )
        
        # Run one-step MD
        file.namd.pair.log <- paste( unlist( strsplit( file.namd.pair,
                                                       "\\.conf" ) ),
                                     ".log",
                                     sep = "" )
        command.md.pair <- paste( exe.namd, " ",
                                  "+p", 1, " ",
                                  file.namd.pair, " ",
                                  " > ",
                                  file.namd.pair.log,
                                  sep = "" )
        system( command.md.pair )
        
        # Load energies
        vec.lines <- readLines( file.namd.pair.log )
        ind.line <- which( grepl( "ENERGY:", 
                                  vec.lines ) )
        ind.line <- ind.line[ length( ind.line ) ]
        line.split <- unlist( strsplit( vec.lines[ ind.line ],
                                        " " ) )
        line.split <- line.split[ line.split != "" ]
        
        vec.energies <- as.numeric( c( line.split[ 7 ],
                                       line.split[ 8 ],
                                       line.split[ 12 ] ) )
        
        return( vec.energies )
      }
    
    # Fill in values to sparse matrix    
    mat.energies <- matrix( 0,
                            nrow = 3,
                            ncol = p )
    mat.energies[ , vec.ind.j ] <- mat.energies.tmp
    
    # Add calculated energies to matrices
    vec.Elec <- rep( NA,
                     p )
    write( format( mat.energies[ 1, ], scientific = FALSE ),
           file = file.Elec,
           ncolumns = p,
           append = TRUE )
    
    vec.VdW <- rep( NA,
                    p )
    write( format( mat.energies[ 2, ], scientific = FALSE ),
           file = file.VdW,
           ncolumns = p,
           append = TRUE )
    
    vec.Total <- rep( NA,
                      p )
    write( format( mat.energies[ 3, ], scientific = FALSE ),
           file = file.Total,
           ncolumns = p,
           append = TRUE )
    
    # Delete files
    prefix <- paste( unlist( strsplit( file.pdb.autopsf.pdb,
                                       "\\.pdb" ) ),
                     "_", i,
                     sep = "" )
    system( paste( "rm ",
                   prefix, "*",
                   sep = "" ) )
    system( "rm md*conf" )
    system( "rm md*log" )
  }
  
  mat.VdW <- read.table( file.VdW )
  mat.Elec <- read.table( file.Elec )
  mat.Total <- read.table( file.Total )
  vec.energy.sum <- c( sum( mat.VdW[ -c( 1, nrow( mat.VdW ) ), ] ),
                       sum( mat.Elec[ -c( 1, nrow( mat.Elec ) ), ] ),
                       sum( mat.Total[ -c( 1, nrow( mat.Total ) ), ] ) )
  
  return( c( file.pdb.autopsf.pdb,
             vec.energy.sum ) )
}
#===============================================================================

#===============================================================================
# Make PSF file
# !!! MAKE SURE NO POINTS IN FILE NAMES!!!
generate.PSF <- function( file.pdb,
                          exe.vmd,
                          vec.files.parameters.prepare,
                          output.prefix )
{
  #===============================================================================
  # Generate TCL script to execute AutoPSF
  file.tcl <- paste( output.prefix,
                     ".tcl",
                     sep = "" )
  write( "package require autopsf",
         file = file.tcl )
  write( paste( "mol new",
                file.pdb ),
         file = file.tcl,
         append = TRUE )
  write( paste( "autopsf",
                "-mol 0",
                "-top", paste( vec.files.parameters.prepare, collapse = " " ) ),
         file = file.tcl,
         append = TRUE )
  write( "exit",
         file = file.tcl,
         append = TRUE )
  
  #===============================================================================
  # Execute TCL script
  file.VMD.output <- paste( file.tcl,
                            ".out",
                            sep = "" )
  vec.lines.VMD.output <- system( paste( exe.vmd,
                                         "-dispdev text",
                                         "-e",
                                         file.tcl ),
                                  intern = TRUE )
  write( vec.lines.VMD.output,
         file.VMD.output )
  
  #===============================================================================
  # Return file names
  prefix <- unlist( strsplit( file.pdb,
                              "\\.pdb" ) )
  vec.psf.files <- c( paste( prefix, "_autopsf.pdb", sep = "" ),
                      paste( prefix, "_autopsf.psf", sep = "" ) )
  
}
#===============================================================================

#===============================================================================
# Function to prepare NAMD config file for pair of residues
prepare.NAMD.config.file <- function( file.psf,
                                      file.pdb,
                                      file.pdb.i.j,
                                      vec.files.parameters.MD,
                                      file.output.i.j,
                                      file.namd.i.j,
                                      cutoff )
{
  # Define prefix variable
  
  
  # Start writing NAMD config file
  write( paste( "structure", file.psf ),
         file = file.namd.i.j )
  write( paste( "coordinates", file.pdb ),
         file = file.namd.i.j,
         append = TRUE )
  write( "paraTypeCharmm on",
         file = file.namd.i.j,
         append = TRUE )
  for( i in 1:length( vec.files.parameters.MD ) )
  {
    write( paste( "parameters", vec.files.parameters.MD[ i ] ),
           file = file.namd.i.j,
           append = TRUE )
    
  }
  write( "numsteps 1",
         file = file.namd.i.j,
         append = TRUE )
  write( "switching on",
         # write( "switching off",
         file = file.namd.i.j,
         append = TRUE )
  write( "exclude scaled1-4",
         file = file.namd.i.j,
         append = TRUE )
  write( paste( "outputname", file.output.i.j ),
         file = file.namd.i.j,
         append = TRUE )
  write( "temperature 0",
         file = file.namd.i.j,
         append = TRUE )
  write( "COMmotion yes",
         file = file.namd.i.j,
         append = TRUE )
  write( paste( "cutoff", cutoff ),
         file = file.namd.i.j,
         append = TRUE )
  write( "dielectric 1",
         file = file.namd.i.j,
         append = TRUE )
  write( "switchdist 10.0",
         file = file.namd.i.j,
         append = TRUE )
  # write( "pairInteraction off",
  #        file = file.namd.i.j,
  #        append = TRUE )
  write( "pairInteraction on",
         file = file.namd.i.j,
         append = TRUE )
  write( "pairInteractionGroup1 1",
         file = file.namd.i.j,
         append = TRUE )
  write( paste( "pairInteractionFile", file.pdb.i.j ),
         file = file.namd.i.j,
         append = TRUE )
  write( "pairInteractionGroup2 2",
         file = file.namd.i.j,
         append = TRUE )
  write( "set ts 0",
         file = file.namd.i.j,
         append = TRUE )
  # write( paste( "coorfile open dcd", file.trajectory.i.j ),
  #        file = file.namd.i.j,
  #        append = TRUE )
  # write( "while { ![coorfile read] } {",
  #        file = file.namd.i.j,
  #        append = TRUE )
  # write( "  firstTimeStep $ts",
  #        file = file.namd.i.j,
  #        append = TRUE )
  write( "  run 0",
         file = file.namd.i.j,
         append = TRUE )
  # write( "  incr ts 1",
  #        file = file.namd.i.j,
  #        append = TRUE )
  # write( "}",
  #        file = file.namd.i.j,
  #        append = TRUE )
  # write( "coorfile close",
  #        file = file.namd.i.j,
  #        append = TRUE )
}
#===============================================================================

#===============================================================================
# Prepare pairwise PDB
prepare.PDB.pairwise <- function( pdb.atom,
                                  file.pdb.pair,
                                  chain.1 = "",
                                  resno.1 = "",
                                  chain.2 = "",
                                  resno.2 = "" )
{
  pdb.atom.new <- pdb.atom
  pdb.atom.new$b <- 0
  pdb.atom.new[ pdb.atom.new$chain == chain.1 &
                  pdb.atom.new$resno == resno.1, ]$b <- 1
  pdb.atom.new[ pdb.atom.new$chain == chain.2 &
                  pdb.atom.new$resno == resno.2, ]$b <- 2
  
  write.pdb( file = file.pdb.pair,
             type = pdb.atom.new$type,
             eleno = pdb.atom.new$eleno,
             elety = pdb.atom.new$elety,
             alt = pdb.atom.new$alt,
             resid = pdb.atom.new$resid,
             chain = pdb.atom.new$chain,
             resno = pdb.atom.new$resno,
             insert = pdb.atom.new$insert,
             xyz = as.vector( t( as.matrix( pdb.atom.new[ , 9:11 ] ) ) ),
             o = pdb.atom.new$o,
             b = pdb.atom.new$b )
  
}
#===============================================================================


#===============================================================================
#===============================================================================
#===============================================================================
# Main

#===============================================================================
# start time measurement
start.main <- proc.time()
#===============================================================================

#===============================================================================
# Command line arguments
file.pdb <- commandArgs()[ 3 ]
n.cores <- as.numeric( commandArgs()[ 4 ] )
exe.foldx <- commandArgs()[ 5 ]
exe.vmd <- commandArgs()[ 6 ]
exe.namd <- commandArgs()[ 7 ]
vec.files.parameters.prepare <- unlist( strsplit( commandArgs()[ 8 ],
                                                  "___" ) )
vec.files.parameters.MD <- unlist( strsplit( commandArgs()[ 9 ],
                                             "___" ) )
exe.yasara <- commandArgs()[ 10 ]
file.yasara.em <- commandArgs()[ 11 ]
cutoff <- as.numeric( commandArgs()[ 12 ] )

option <- as.numeric( commandArgs()[ 13 ] )

if( 1 == option )
{
  print( "PaDRINO interface only!" )

  # Chains of first interaction partner
  vec.chains.interface.1 <- unlist( strsplit( commandArgs()[ 14 ],
                                              "_" ) )
  # Chains of second interaction partner
  vec.chains.interface.2 <- unlist( strsplit( commandArgs()[ 15 ],
                                              "_" ) )
} else if( 2 == option )
{
  print( "PaDRINO stability only!" )

  # Chains of first interaction partner
  vec.chains.stability.1 <- unlist( strsplit( commandArgs()[ 14 ],
                                              "_" ) )
  # Chains of second interaction partner
  vec.chains.stability.2 <- unlist( strsplit( commandArgs()[ 15 ],
                                              "_" ) )
} else if( 3 == option )
{
  print( "PaDRINO interface and stability only!" )

  # Chains of first interaction partner
  vec.chains.interface.1 <- unlist( strsplit( commandArgs()[ 14 ],
                                              "_" ) )
  # Chains of second interaction partner
  vec.chains.interface.2 <- unlist( strsplit( commandArgs()[ 15 ],
                                              "_" ) )

  # Chains of first interaction partner
  vec.chains.stability.1 <- unlist( strsplit( commandArgs()[ 16 ],
                                              "_" ) )
  # Chains of second interaction partner
  vec.chains.stability.2 <- unlist( strsplit( commandArgs()[ 17 ],
                                              "_" ) )
} else
{
  stop( "Restart! You have to choose between options 1 2 3!" )
}


# #===============================================================================
# #===============================================================================
# #===============================================================================
# # Manual test
# file.pdb <- "m42_HH28K_13_renum.pdb"
# n.cores <- 70
# exe.foldx <- "/home/rawir/software/FoldX/foldx"
# exe.vmd <- "/usr/local/bin/vmd"
# exe.namd <- "/home/rawir/software/NAMD_2.13_Linux-x86_64-multicore/namd2"
# vec.files.parameters.prepare <- c( "/home/rawir/Documents/work/vrc/data/charmm/top_all36_prot.rtf",
#                                    "/home/rawir/Documents/work/vrc/data/charmm/top_all36_carb.rtf",
#                                    "/home/rawir/Documents/work/vrc/data/charmm/toppar_water_ions_namd.str" )
# vec.files.parameters.MD <- c( "/home/rawir/Documents/work/vrc/data/charmm/par_all36_prot_cu_patch.prm",
#                               "/home/rawir/Documents/work/vrc/data/charmm/par_all36_carb.prm",
#                               "/home/rawir/Documents/work/vrc/data/charmm/toppar_water_ions_namd.str")
# exe.yasara <- "/home/rawir/software/YASARA.app/yasara/yasara"
# file.yasara.em <- "/data2/rawir/Documents/work/vrc/papers/P3-43/code/software/dependencies/YASARA_EM.mcr"
# cutoff <- 12
# option <- 3
# vec.chains.interface.1 <- "A"
# vec.chains.interface.2 <- c( "H", "L" )
# vec.chains.stability.1 <- c( "H", "L" )
# vec.chains.stability.2 <- c( "H", "L" )
# #===============================================================================
# #===============================================================================
#===============================================================================

#===============================================================================
# Source function

#===============================================================================
# Set working directory
wd <- getwd()
setwd ( wd )

#===============================================================================
# Set variables
registerDoParallel( n.cores )
namd.cores <- n.cores


#===============================================================================
#===============================================================================
# For each mutant
# (i) Generate mutant structure
# (ii) Minimize structure
# (iii) Calculate PaDRINO energies

#===============================================================================
# (i) Generate mutant structure
file.pdb.mutant <- run.FoldX( exe.foldx,
                              file.pdb )

#===============================================================================
# (ii) Minimize structure
file.minimized <- YASARA.EM( getwd(),
                             exe.yasara,
                             file.yasara.em,
                             file.pdb.mutant )
# Clean protein (keep only protein content)
file.pdb.minimized.clean <- clean.PDB( file.minimized )
print( file.pdb.minimized.clean )


#===============================================================================
# (iii) Calculate PaDRINO energies
if( 1 == option )
{
  #==================================================
  # PaDRINO.bind (Interface)
  # "Calculate residue-residue pairwise energies"
  print( "Option 1: Calculate residue-residue pairwise energies only" )
  print( file.pdb.minimized.clean )
  print( vec.chains.interface.1 )
  print( vec.chains.interface.2 )
  print( cutoff )
  
  vec.energies.interface <- PaDRINO.energy( file.pdb.minimized.clean,
                                            chain.1 = vec.chains.interface.1,
                                            chain.2 = vec.chains.interface.2,
                                            cutoff = cutoff )
} else if( 2 == option )
{
  #==================================================
  # PaDRINO.bind (Stability)
  # "Calculate residue-residue pairwise energies"
  
  vec.energies.stability <- PaDRINO.energy( file.pdb.minimized.clean,
                                            chain.1 = vec.chains.stability.1,
                                            chain.2 = vec.chains.stability.2,
                                            cutoff = cutoff )
  
  pdb <- read.pdb( file.pdb.minimized.clean )
  pdb.atom <- pdb$atom
  n.H <- nrow( pdb.atom[ pdb.atom$elety == "CA" &
                           pdb.atom$chain == "H", ] )
  n.L <- nrow( pdb.atom[ pdb.atom$elety == "CA" &
                           pdb.atom$chain == "L", ] )
  
  mat.VdW <- read.table( "VdW_HL_HL.txt" )
  mat.VdW.H.H <- mat.VdW[ 1:n.H, 1:n.H ]
  mat.VdW.L.L <- mat.VdW[ (n.H+1):(n.H+n.L), (n.H+1):(n.H+n.L) ]
  mat.VdW.H.L <- mat.VdW[ 1:n.H, (n.H+1):(n.H+n.L) ]
  
  mat.Elec <- read.table( "Elec_HL_HL.txt" )
  mat.Elec.H.H <- mat.Elec[ 1:n.H, 1:n.H ]
  mat.Elec.L.L <- mat.Elec[ (n.H+1):(n.H+n.L), (n.H+1):(n.H+n.L) ]
  mat.Elec.H.L <- mat.Elec[ 1:n.H, (n.H+1):(n.H+n.L) ]
  
  mat.Total <- read.table( "Total_HL_HL.txt" )
  mat.Total.H.H <- mat.Total[ 1:n.H, 1:n.H ]
  mat.Total.L.L <- mat.Total[ (n.H+1):(n.H+n.L), (n.H+1):(n.H+n.L) ]
  mat.Total.H.L <- mat.Total[ 1:n.H, (n.H+1):(n.H+n.L) ]
  vec.energies.stability <- c( sum( mat.VdW.H.H ),
                               sum( mat.Elec.H.H ),
                               sum( mat.Total.H.H ),
                               sum( mat.VdW.L.L ),
                               sum( mat.Elec.L.L ),
                               sum( mat.Total.L.L ),
                               sum( mat.VdW.H.L ),
                               sum( mat.Elec.H.L ),
                               sum( mat.Total.H.L ) )
} else if( 3 == option )
{
  #==================================================
  # PaDRINO.bind (Interface+)
  # "Calculate residue-residue pairwise energies"
  
  # Interface
  vec.energies.interface <- PaDRINO.energy( file.pdb.minimized.clean,
                                            chain.1 = vec.chains.interface.1,
                                            chain.2 = vec.chains.interface.2,
                                            cutoff = cutoff )
  
  # Stability
  vec.energies.stability <- PaDRINO.energy( file.pdb.minimized.clean,
                                            chain.1 = vec.chains.stability.1,
                                            chain.2 = vec.chains.stability.2,
                                            cutoff = cutoff )
  
  pdb <- read.pdb( file.pdb.minimized.clean )
  pdb.atom <- pdb$atom
  n.H <- nrow( pdb.atom[ pdb.atom$elety == "CA" &
                           pdb.atom$chain == "H", ] )
  n.L <- nrow( pdb.atom[ pdb.atom$elety == "CA" &
                           pdb.atom$chain == "L", ] )
  
  mat.VdW <- read.table( "VdW_HL_HL.txt" )
  mat.VdW.H.H <- mat.VdW[ 1:n.H, 1:n.H ]
  mat.VdW.L.L <- mat.VdW[ (n.H+1):(n.H+n.L), (n.H+1):(n.H+n.L) ]
  mat.VdW.H.L <- mat.VdW[ 1:n.H, (n.H+1):(n.H+n.L) ]
  
  mat.Elec <- read.table( "Elec_HL_HL.txt" )
  mat.Elec.H.H <- mat.Elec[ 1:n.H, 1:n.H ]
  mat.Elec.L.L <- mat.Elec[ (n.H+1):(n.H+n.L), (n.H+1):(n.H+n.L) ]
  mat.Elec.H.L <- mat.Elec[ 1:n.H, (n.H+1):(n.H+n.L) ]
  
  mat.Total <- read.table( "Total_HL_HL.txt" )
  mat.Total.H.H <- mat.Total[ 1:n.H, 1:n.H ]
  mat.Total.L.L <- mat.Total[ (n.H+1):(n.H+n.L), (n.H+1):(n.H+n.L) ]
  mat.Total.H.L <- mat.Total[ 1:n.H, (n.H+1):(n.H+n.L) ]
  vec.energies.stability <- c( sum( mat.VdW.H.H ),
                               sum( mat.Elec.H.H ),
                               sum( mat.Total.H.H ),
                               sum( mat.VdW.L.L ),
                               sum( mat.Elec.L.L ),
                               sum( mat.Total.L.L ),
                               sum( mat.VdW.H.L ),
                               sum( mat.Elec.H.L ),
                               sum( mat.Total.H.L ) )
}

#==================================================
# Save energies (lines 1-3: VdW, Elec, Total)
if( 1 == option )
{
  write( vec.energies.interface[ 2:4 ],
         file = "Energies_option1.txt",
         ncolumns = 1 )
} else if( 2 == option )
{
  write( vec.energies.stability,
         file = "Energies_option2.txt",
         ncolumns = 1 )
} else if( 3 == option )
{
  write( c( vec.energies.interface[ 2:4 ],
            vec.energies.stability ),
         file = "Energies_option3.txt",
         ncolumns = 1 )
}





#===============================================================================
# Measure and print out time needed for the script
end.main <- proc.time()
duration.main <- end.main-start.main
print( paste( "Script duration:", round( duration.main[3] / 60, 2 ), "min") )
#===============================================================================

# Main end
#===============================================================================
#===============================================================================
#===============================================================================
