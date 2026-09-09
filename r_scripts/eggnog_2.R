library(dplyr)
library(tidyr)
library(readxl)
library(clusterProfiler)
library(enrichplot)

######## ORA ######## 
# query         → gene ID (used in all files)
# GOs           → comma-separated GO terms
# KEGG_Pathway  → comma-separated pathways
# read eggnog annotations output for YBM_2,1
egg_universe <- read_excel("/Users/stavrt/marmot/eggmap/YBM_2.1_emapper_v2_files/YBM_2.1.emapper.gene.annotations.xlsx")
View(egg_universe)
names(egg_universe)

prom_snow <- read_excel("/Users/stavrt/marmot/eggmap/YBM_2.1_emapper_v2_files/YBM_2.1.emapper.gene.annotations.promoters_intersect_snowmelt.xlsx")
prom_temp <- read_excel("/Users/stavrt/marmot/eggmap/YBM_2.1_emapper_v2_files/YBM_2.1.emapper.gene.annotations.promoters_intersect_temperature.xlsx")
gen_snow <- read_excel("/Users/stavrt/marmot/eggmap/YBM_2.1_emapper_v2_files/YBM_2.1.emapper.gene.annotations.genes_intersect_snowmelt.xlsx")
gen_temp <- read_excel("/Users/stavrt/marmot/eggmap/YBM_2.1_emapper_v2_files/YBM_2.1.emapper.gene.annotations.genes_intersect_temperature.xlsx")

# FOREGROUND LIST: first mention of query gene to avoid any duplicated CpGs in ORA
fore_list <- list(prom_snow = prom_snow, prom_temp = prom_temp, gen_snow  = gen_snow, gen_temp  = gen_temp)

# deduplicate >1 CpGs falling in the same gene by keeping only the FIRST mention of each genename (query)
for (nm in names(fore_list)) {
  fore_list[[nm]] <- fore_list[[nm]] %>%
    distinct(query, .keep_all = TRUE)}
# number of queries per analysis after deduplication
sapply(fore_list, function(df) length(unique(df$query)))

# UNIVERSE LIST: as a vector of queries
universe_list <- unique(egg_universe$query)
head(universe_list)

# extract GO terms from "GOs" eggnog column which has ":" after GO and is comma sep. 
term2gene_go <- egg_universe %>%
  select(query, GOs) %>%
  filter(!is.na(GOs), GOs != "", GOs != "-") %>%  # filter empty strings
  separate_rows(GOs, sep = ",") %>%
  transmute(term = GOs, gene = query) %>%
  distinct()

# extract GO terms and names (Description column) for labeling plot
term2name_go <- egg_universe %>%
  select(GOs, Description) %>%
  filter(!is.na(GOs), GOs != "", GOs != "-") %>%
  separate_rows(GOs, sep = ",") %>%
  distinct(GOs, Description) %>%
  group_by(GOs) %>%
  slice(1) %>%                    # keep first Description per GO term in case there are multiple
  ungroup() %>%
  transmute(term = GOs, name = Description)

# extract KEGG from "KEGG_Pathway" eggnog column which just comma sep.
term2gene_kg <- egg_universe %>%
  select(query, KEGG_Pathway) %>%
  filter(!is.na(KEGG_Pathway), KEGG_Pathway != "", KEGG_Pathway != "-") %>%
  separate_rows(KEGG_Pathway, sep = ",") %>%
  transmute(term = KEGG_Pathway, gene = query) %>%
  distinct()

# extract KEGG pathways and names (Description column) for labeling plot
term2name_kg <- egg_universe %>%
  select(KEGG_Pathway, Description) %>%
  filter(!is.na(KEGG_Pathway), KEGG_Pathway != "", KEGG_Pathway != "-") %>%
  separate_rows(KEGG_Pathway, sep = ",") %>%
  distinct(KEGG_Pathway, Description) %>%
  group_by(KEGG_Pathway) %>%
  ungroup() %>%
  transmute(term = KEGG_Pathway, name = Description)

# how many unique GO terms and KEGG pathways in egg universe file?
length(unique(term2gene_go$term)) # 23780
length(unique(term2gene_kg$term)) # 786

# correct term format?
head(term2gene_go$term, 10)

########  TEST ORA for GO and KEGG, using clusterProfiler sig. defaults ######## 
go_enrich   <- list()
kegg_enrich <- list()

for (nm in names(fore_list)) {
  
  fg_genes <- fore_list[[nm]]$query
  
  go_enrich[[nm]] <- enricher(
    gene          = fg_genes,
    universe      = universe_list,
    TERM2GENE     = term2gene_go,
    TERM2NAME     = term2name_go,
# enrichplots show term IDs unless you explicitly provide a TERM2NAME mapping due to custom eggnog annts
# assigned names based on the FIRST Keeping the first description per GO erm is standard and acceptable
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.20
  )
  
  kegg_enrich[[nm]] <- enricher(
    gene          = fg_genes,
    universe      = universe_list,
    TERM2GENE     = term2gene_kg,
    TERM2NAME     = term2name_kg,
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.20
  )
}

# plotting GO terms
snow_prom_GO <- dotplot(go_enrich$prom_snow, showCategory = 30) + ggtitle("Gene promoters, snowmelt")     # A)
#dotplot(go_enrich$prom_temp, showCategory = 30) + ggtitle("Gene promoters, temperature")  # NULL
#dotplot(go_enrich$gen_snow,  showCategory = 10) + ggtitle("Gene bodies, snowmelt")        # NULL
temp_bod_GO <- dotplot(go_enrich$gen_temp,  showCategory = 30) + ggtitle("Gene bodies, temperature")     # B)

# plotting kegg pathways
dotplot(kegg_enrich$prom_snow, showCategory = 10) + ggtitle("KEGG, prom_snow")
#dotplot(kegg_enrich$prom_temp, showCategory = 10) + ggtitle("KEGG ORA – prom_temp")       # NULL
#dotplot(kegg_enrich$gen_snow,  showCategory = 10) + ggtitle("KEGG ORA – gen_snow")        # NULL
dotplot(kegg_enrich$gen_temp,  showCategory = 10) + ggtitle("KEGG – gen_temp")

# save two GO plots
ggsave("/Users/stavrt/marmot/FIGURES/snow_prom_GO.png", snow_prom_GO, width = 6, height = 3, dpi = 300)
ggsave("/Users/stavrt/marmot/FIGURES/temp_bod_GO.png", temp_bod_GO, width = 6, height = 3, dpi = 300)

# print out results that were not null
go_snow <- as.data.frame(go_enrich$prom_snow)
go_snow[1, c("ID", "Count", "GeneRatio", "BgRatio", "p.adjust")]

kegg_snow <- as.data.frame(kegg_enrich$prom_snow)
kegg_snow[1, c("ID", "Count", "GeneRatio", "BgRatio", "p.adjust")]

go_temp <- as.data.frame(go_enrich$gen_temp)
go_temp[1, c("ID", "Count", "GeneRatio", "BgRatio", "p.adjust")]

kegg_temp <- as.data.frame(kegg_enrich$gen_temp)
kegg_temp[1, c("ID", "Count", "GeneRatio", "BgRatio", "p.adjust")]


# save results as a table 
print(go_table)
