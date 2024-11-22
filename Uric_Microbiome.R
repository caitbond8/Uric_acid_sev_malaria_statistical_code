library(phyloseq)
library(microViz)
library(tidyverse)


asv <- import_biom("fixed_ndi.biom") #otu and phylogeny
meta <- read_csv("NDI_samples_enroll_ext_20240619_clean_meta.csv", col_select = 2:5) #cleaned metadata from stool samples
new <- read_csv("uric_study_meta.csv") #study metadata

meta_new <- left_join(meta, new, by = "studyid") %>% as.data.frame #merge stool and study metadata
rownames(meta_new) <- meta_new$SampleID 
meta_new <- meta_new %>% select(-SampleID)
sample_meta <- sample_data(meta_new) 
asv <- merge_phyloseq(asv, sample_meta) #create phyloseq object

#clean up tax table
tax_table(asv)[, colnames(tax_table(asv))] <- gsub(tax_table(asv)[, colnames(tax_table(asv))],     pattern = "[a-z]__", replacement = "")
colnames(tax_table(asv))[1:8] <- c("Domain",  "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species", "Strain")
colnames(tax_table(asv))
asv <- tax_fix(asv, verbose = FALSE)
asv <- tax_rename(asv, rank = "Species") %>% ps_mutate(sick = group %in% c(1,2,3,4,5))

#create necessary phyloseq objects by filtering based on metadata
asv<- asv %>%
  ps_mutate(hyper_uri = case_when(uric_se1 > 7 ~ 1,
                                  uric_se1 <= 7 ~ 0),
            tff3_inj = case_when(tff3_pl1 >= 4.078 ~ 1,
                                 tff3_pl1 < 4.078 ~ 0),
            date_to_stool = (as.Date(stool_dt0, format="%d/%m/%Y") - as.Date(doe, format="%d/%m/%Y")),
            stool_day = case_when(date_to_stool == 0 ~ "day 0",
                                  date_to_stool == 1 ~ "day 1",
                                  date_to_stool == 2 ~ "day 2",
                                  date_to_stool == 3 ~ "day 3",
                                  date_to_stool == 4 ~ "day 4",
                                  date_to_stool == 5 ~ "day 5",
                                  date_to_stool == 6 ~ "day 6",
                                  date_to_stool == 7 ~ "day 7",
                                  date_to_stool > 7 ~ "past one week")) 

asv_0 <- asv %>% ps_filter(month == 0)
asv_sick <- asv_0 %>% ps_filter(sick == TRUE)


bb_models_uric <- asv_sick %>% 
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("uric_se1")
  )
bb_stats_uric <- taxatree_models2stats(bb_models_uric, param = "mu")
bb_stats_uric
bb_stats_uric %>% taxatree_stats_get()
tree_group <- bb_stats_uric %>%
  taxatree_plots(
    node_size_range = c(1, 4), palette = "Blue-Red 3"
  )
tree_group 
key <- taxatree_plotkey(
  data = bb_stats_uric,
  taxon_renamer = function(x) stringr::str_remove(x, "[PCFGS]: "),
  # 2 lines of conditions below, for filtering taxa to be labelled
  rank == "Family" | rank == "Genus" |rank == "Species", #& prevalence > 0.2,
  p.value < 0.05, 
  !grepl("Kingdom", taxon)
) +
  # add a bit more space for the longer labels by expanding the x axis
  scale_x_continuous(expand = expansion(mult = 0.2))
ggsave("16S_corncob_uric_acid_circle_tree_genus.png", tree_group, width = 13, height = 5.5, dpi = 1200, device = "png")


hyperuri_barplot <- asv_sick %>% 
  phyloseq::merge_samples(group = "hyper_uri") %>%
  comp_barplot(tax_level = "Family", bar_width = 0.8, n_taxa = 3, tax_order = c("Enterobacteriaceae", "Bacteroidaceae","Enterococcaceae"))+
  labs(x= NULL, y = NULL, title = "Hyperuri")
ggsave("family_sick_hyperuric.png", hyperuri_barplot, width = 8.3, height = 10, dpi = 1200, device = "png")

#file for intestinal damage plot
plot_data <- asv_sick %>%
  tax_fix() %>% 
  tax_transform("compositional", rank = "Genus") %>%
  tax_transform("log", zero_replace = "halfmin", chain = TRUE) %>%
  ps_get() %>%
  ps_otu2samdat() %>% 
  samdat_tbl()
plot_data <- plot_data %>% select(studyid, Escherichia, Shigella, Klebsiella, Enterobacter,Bacteroides, Enterococcus)
write_csv(plot_data, "sm_abund_uric_genus.csv")
