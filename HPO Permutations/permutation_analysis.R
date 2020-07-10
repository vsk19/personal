# Conducts permutation analysis on patient Human Phenotype Ontology terms and disease affected status

library(dplyr)
library(argparser)

# Function: split up hpo terms by Patient ID
split_tibble <- function(tibble, col = 'col') tibble %>% split(., .[,col])


######## Parse user input #########
parser <- arg_parser('Conducts HPO term permutation analysis.')
parser <- add_argument(parser, '--outdir', help = 'Output directory')
parser <- add_argument(parser, '--input_dir', help = 'Input directory where hpo data is stored')
parser <- add_argument(parser, '--num_permutations', help = 'Number of permutations to conduct')

args <- parse_args(parser)
outdir <- args$outdir
indir <- args$input_dir
num_permutations <- as.numeric(args$num_permutations)

######## Build list of HPO Terms per subject #########

# all HPO terms for all Patient IDs
hpo_terms <- read.csv(paste(c(indir, 'hpo_ancestors.csv'), collapse = ''),
                      header = TRUE, stringsAsFactors = FALSE)[,1:4]
hpo_terms <- dplyr::select(hpo_terms, Patient_ID, Term, HPO_ID)

hpo_list <- split_tibble(hpo_terms, 'Patient_ID')


# Only want to use Probands with HPO terms 
probands_hpo <- read.csv(paste(c(indir, 'subj_variants_hpo.csv'), collapse = ''),
                         header = TRUE, stringsAsFactors = FALSE, na.strings = c("", " ", "NA"))
probands_hpo <- dplyr::select(probands_hpo, Patient_ID, Proband, Gene, HPO_IDs)
probands_hpo <- unique(probands_hpo)

probands_hpo <- probands_hpo[probands_hpo$Proband == "Yes",]
# probands_wo_hpo <- probands_hpo[is.na(probands_hpo$HPO_IDs),]
probands_hpo <- probands_hpo[!is.na(probands_hpo$HPO_IDs),]

# Removes genes that only occur once in the table (>= two fams don't share variants in the same gene anymore
# since we removed some probands/families if they didn't have HPO Terms).
probands_hpo <- probands_hpo[probands_hpo$Gene %in% probands_hpo$Gene[duplicated(probands_hpo$Gene)], ]


######## Set up permutation table ########

ids <- unique(probands_hpo$Patient_ID)
genes <- unique(probands_hpo$Gene)#[1:50]
probands_hpo <- probands_hpo[probands_hpo$Gene %in% genes,]

# These HPO terms will be used in the permutations
perm_hpo_terms <- unique(hpo_terms[hpo_terms$Patient_ID %in% ids, c('Term', 'HPO_ID')])

# Build permutation table. Cols are Genes and Rows are Phen IDs
mut_status <- data.frame(matrix(ncol = length(genes), nrow = length(ids)), row.names = ids)
names(mut_status) <- genes

# Put "Y" for people with a mutation in the gene and "N" for everyone else
for (id in ids) {
  curr_genes <- probands_hpo$Gene[probands_hpo$Patient_ID == id]
  if (!rlang::is_empty(curr_genes))
    mut_status[id, curr_genes] <- 'Y'
}
mut_status[is.na(mut_status)] <- 'N'


######## Permute for all HPO terms in these probands #########

pvals <- data.frame(matrix(0, nrow = length(perm_hpo_terms$HPO_ID), ncol = length(genes)), row.names = perm_hpo_terms$HPO_ID)
names(pvals) <- genes

for (hpo_id in perm_hpo_terms$HPO_ID) {
  print(hpo_id)
  
  # The first row stores the ACTUAL frequency in the dataset
  freq_tab <- data.frame(matrix(nrow = num_permutations + 1, ncol = length(genes)))
  names(freq_tab) <- genes
  
  # permute the mutation status 1000 times and see how many of these peeps have skin rash as a term
  # First iter is on the un-randomized table to get ACTUAL frequencies
  perm_tab <- mut_status
  for (i in 1: (num_permutations + 1)) {
    
    for (gene in names(perm_tab)) {
      # count number of people with the term
      
      have_mutation <- perm_tab[perm_tab[, gene] == 'Y' ,]
      num_w_mutation <- nrow(have_mutation)
      num_w_term <- 0
      
      # loop through people with the mutation
      for (id in row.names(have_mutation)) {
        
        # if person has the HPO term, increase the freq count
        if (grepl(hpo_id, hpo_list[[id]]['HPO_ID']))
          num_w_term <- num_w_term + 1
      }
      
      # Store the frequency of the gene/HPO combo
      freq_tab[i, gene] <- num_w_term/num_w_mutation
      
      # Count how many times the actual frequency > random frequency
      if (freq_tab[1, gene] > freq_tab[i, gene])
        pvals[hpo_id, gene] <- pvals[hpo_id, gene] + 1
    }
    
    # randomize the table for the next iteration
    perm_tab <- perm_tab[sample(nrow(perm_tab)),]
    rownames(perm_tab) <- ids
  }

}

# output the raw counts (how many times the actual freq > permutation freq)
write.csv(pvals, paste(c(outdir, 'raw_pval_count.csv'), collapse = ''), row.names = TRUE)

# Convert count of actual > random frequencies to a p-value
pvals[] <- lapply(pvals, function(x) 1-x/num_permutations)

write.csv(pvals, paste(c(outdir, 'permutation_test.csv'), collapse = ''), row.names = TRUE)


# Now, build list of HPO terms that have a significant association with a mutation in at least one gene
significant_term_ids <- which(apply(pvals, 1, function(r) any(r <= 0.05)))
significant_term_ids <- names(significant_terms_ids)

# Subset the permutation HPO terms to the significant ones
sig_hpo_tab <- perm_hpo_terms[perm_hpo_terms$HPO_ID %in% significant_term_ids, ] 
write.csv(sig_hpo_tab, paste(c(outdir, 'significant_hpo_terms.csv'), collapse = ''), row.names = FALSE)

