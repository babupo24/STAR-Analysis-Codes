# STAR-Analysis-Codes
This repo contains the updated macros for the dipion cross-section analysis using Run 12 pp data.

Step 1: Produce run 12 pp200 data trees
: Run 12 (picoDst) data trees needs to be produced, which can be done using macros within DataTreesCodes.
How to Run:
#load library 
starver SL21d #MuDst -> picoDst production library
cons #compiles code
sh Run12ppTreeMaker.sh #this will submit jobs for each run in the KevinRunsList.list using catlog for hpss files.
: Safer to direct outputs in the SCRATCH. 
: Outputs root trees combined for each run with codes add_trees* 
  Change input and output directories in the add_trees.sh and submit jos to add root trees with 
  star-submit add_trees.xml
  
Step 2: Produce embedding trees for 2012 pp 200 datasets
: Full sample embedding was produced by Dmitry, the same embedding sample is used for this analysis. 
: Embedding sample location:- /star/embed/embedding/pp200_production_2012/v2/
: Embedding trees can be produced using codes within EmbedTreeCodes. 
: Embedding trees are produced per run as in the data trees, even though it is produced in partonic pt bins. 
How to use codes? 
  1. Create list of MuDsts for each run and put it in ptListPerRun directory. 
  2. RunEmbeddingTreeMaker.sh prepares and submit jobs in SUMS calling RunEmbeddingTreeMaker.xml. Jobs is submitted for each run with the MuDst.root and corresponding
 pythia.root file simultaneously. pythia.root is for pythia event information.
  3. Commands
    starver SL12d_embed #embedding production library
    cons #compiles codes within StRoot
    sh RunEmbeddingTreeMaker.sh #submit jobs for all/selected runs
  4. Combine root trees for each runs using add_trees* scripts.
  
 Step 3:  
  
  
  
  
