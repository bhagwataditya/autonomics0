# Install 
# source('https://bitbucket.org/graumannlabtools/autonomics/downloads/install_autonomics.R')

# Load
library(autonomics)
library(autonomics.import)

# Create a temporary analysis dir
analysis_dir <- tempdir()
analysis_dir <- gsub('\\', '/', analysis_dir, fixed = TRUE)  # Correct for weird path notation on windows
analysis_dir <- paste0(analysis_dir, '/', 'billing_2016_stem_cell_comparison')
dir.create(analysis_dir)

# Create sample design file
# Open it with Excell
# Complete columns L, M, H using column 'Names'
# Save as tab separated txt file
sample_design_file <- paste0(analysis_dir, '/sample_design.txt')
create_sample_design_file(
   infile  = get_billing_experimental_design_template_file(), 
   outfile = sample_design_file
)
sample_design_file

# First have a look at the proteinGroups.txt file
# to get an understanding of what we are doing
(protein_groups_file <- get_billing_protein_groups_file())

# Then load the protein_groups file into an ProtSet object
pset <- load_max_quant_table(max_quant_file = protein_groups_file, sample_file = sample_design_file)

# Invert some ratios to correct for label swap
sdata(pset)
pset <- invert(pset, subgroup_levels = c('E_EM', 'E_BM', 'EM_BM'))
sdata(pset)

# Analyze protein groups
pset <- analyze_eset(pset, analysis_dir, organism = 'hsa')
analysis_dir
dir(analysis_dir)
dir(paste0(analysis_dir, '/BM_E'))
