###########################
#  Eric Anderson filtered the vcf file to remove all of the read info from all of the samples, since
#  only the SNP position info was needed. This is done with vcftools (not an R script). 
#  We have vcftools on Ahi. For a Mac, download from the vcftools githup site in Terminal using:
#     git clone https://github.com/vcftools/vcftools.git

# vcftools --vcf Ppho_gtseq_q30dp10mac3.recode.vcf --out Ppho_gtseq4microhap --recode --indv Ppho5103
#   The first sample in this file is Ppho5103, so all samples other than Ppho5103 is removed
# zip the file:
#  gzip Ppho_gtseq4microhap.recode.vcf
# also zip the filtered SNP data from microhaplot:
#  gzip NS276-333_observed_filtered_haplotype.csv

library(tidyverse)
library(vcfR)

description = "Ppho_340loc_trimmed2_largefiles_q30dp10mac3_4microhap"
vcf <- read.vcfR("data/Ppho_340loc_trimmed2_largefiles_q30dp10mac3_4microhap.recode.vcf")

# get a tidy data frame of it, and also find the index of each 
# SNP within each amplicon
vtidy <- vcfR2tidy(vcf, info_only = TRUE)$fix %>%
  select(CHROM:QUAL) %>%
  arrange(CHROM, POS) %>%
  group_by(CHROM) %>%
  mutate(idx = 1:n()) %>%
  ungroup()

### Get the microhap data and explode it into SNPs
# import the compressed filtered haplotype data from microhaplot:
mh_raw <- read_csv("data/Ppho_340loc_trimmed2_largefiles_q30dp10mac3_filteredHaps.csv")

# now, we want to summarise it by locus and haplo.  For each, we want to know
# the number or distinct individuals it was seen in, and the total read depth
mh_summ <- mh_raw %>%
  group_by(locus, haplo) %>%
  summarise(nseen = n_distinct(indiv.ID),
            total_depth = sum(depth)) %>%
  ungroup()

# now we want to explode those haplotypes into SNPs so we can count up
# number of occurrences of Ns etc.  This is a job for dplyr do()

# we make a helper function first
explode <- function(hap) {
  s <- str_split(hap, "")[[1]]
  tibble(idx = 1:length(s),
         base = s
  )
}

# then explode them
mh_exploded <- mh_summ %>%
  group_by(locus, haplo, nseen, total_depth) %>%
  do(explode(.$haplo))

# now we can tally up how many Ns have been seen in each position (idx)
mh_snps_sum <- mh_exploded %>%
  group_by(locus, idx, base) %>%
  summarise(nseen = sum(nseen),
            total_depth = sum(total_depth)) %>%
  mutate(seen_fract = nseen / sum(nseen)) %>%
  ungroup()

# that last line gets the fraction of base calls at each position in these
# called haplotypes that are Ns...or something like that. it isn't clear to me
# exactly how that plays out, but I believe it is at least a reasonable measure of 
# how often N's are seen.

#Now, with that you can make informed decisions about which SNPs you want to toss.
#In some cases, it would be better to just no-call individuals that have any Ns (like, 
#if there are only 5% of base calls being N at a position.)  Meaning, you should leave some 
#sites that have a couple Ns in your data set.  But if an individual carries a haplotype at high
#frequency that has an N, you will just record that indivdidual as having missing data at that locus.

#You can look at all the Ns like this:
mh_snps_sum %>% filter(base == "N") %>% View

### Filter things out, etc.

#OK.  Now let's imagine we will toss any site in which seen_fract is greater than 5%.
#We can pick those out like this:
tossers <- mh_snps_sum %>%
filter(base == "N" & seen_fract > 0.05)
# SNPs that you would toss.  Here is the distribution of their positions
tossers %>%
  count(idx)
#plot the distribution of SNPs to remove
ggplot(tossers, aes(x = seen_fract)) +
  geom_histogram()
#We might find it interesting to plot the seen_fract against the idx (SNP position?) in the amplicon:
ggplot(tossers, aes(x = seen_fract, y = idx)) +
  geom_jitter(alpha = 0.5, colour = "blue")

#So, finally, we can toss those guys from our VCF. We need to get the actual POSes for them.
#That is easy with a join or two.
### Check that we have the right VCF, etc.
#First, let's verify that every idx from microhaplot has a corresponding 
#one from the VCF, etc.
mh_snps_sum %>%
rename(CHROM = locus) %>%
anti_join(vtidy, by = c("CHROM", "idx"))
# or vice verse:
mh_snps_sum %>%
  rename(CHROM = locus) %>%
  anti_join(vtidy, ., by = c("CHROM", "idx"))
# not sure what these do. The first didn't produce a table, just a comment:
# A tibble: 0 x 6
# ... with 6 variables: CHROM <chr>, idx <int>, base <chr>, nseen <int>, total_depth <int>, seen_fract <dbl>

# The second produced a table with 324 rows; not sure what to look for. 

### Get the POSes of loci to exclude
exclude_em <- tossers %>%
  rename(CHROM = locus) %>%
  left_join(., vtidy, by = c("CHROM", "idx"))
# Look at the table to see how many of the positions are only represented by a few samples (e.g., <20)
View(exclude_em)

# write a table of positions to be excluded:
exclude_em %>%
  select(CHROM, POS) %>%
  write.table(., file =paste(description, "_drop.txt",sep=""), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
# table includes 820 positions (out of 2479 in original vcf file ('vtidy')

#######################
# Copy the "drop.txt" file to Ahi, to folder containing the original vcf file, then run vcftools to filter:
in terminal on Ahi:
  vcftools --vcf Ppho_340loc_trimmed2_largefiles_q30dp10mac3_4microhap.recode.vcf --out Ppho_340loc_trimmed2_largefiles_q30dp10mac3_N_cleaned --recode --exclude-positions Ppho_340loc_trimmed2_largefiles_q30dp10mac3_4microhap_drop.txt 
# After filtering, kept 1659 out of a possible 2479 Sites
# output vcf file = "Ppho_N_cleaned.recode.vcf". copy back to Mac. 
#######################

#Let's just quickly count up how many amplicons were retained:
fv <- read.vcfR("Ppho_340loc_trimmed2_largefiles_q30dp10mac3_N_cleaned.recode.vcf") %>%
vcfR2tidy(., info_only = TRUE) %>%
.$fix

n_distinct(fv$CHROM)
# ouptut = 248
n_distinct(vtidy$CHROM)
# output = 333 (original number of loci in vcf file)

