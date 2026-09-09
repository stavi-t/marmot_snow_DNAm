library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)

adc_raw <- read_excel("/Users/stavrt/marmot/snowmelt/model_results/promoter_regions_characterization.xlsx", 
                     sheet = 1,range = cell_rows(1:18), na = c("nan", "NA"))
View(adc_raw)

# highlight sig. cs
cpgs <- c(2089073, 2089194)
TSS_longest <- 2088810

# first split "site" into scaffold + position
adc_clean <- adc_raw %>%
  separate(Site, into = c("scaffold", "pos"), sep = ":", convert = TRUE)

# long format and compute methylation proportion over region
adc_long <- adc_clean %>% 
    pivot_longer(
    cols = -c(scaffold, pos),
    names_to = c("sample", "measure"),
    names_pattern = "(.*)_(meth|total)_cytosine",
    values_to  = "value"
  ) %>%
  pivot_wider(
    names_from = measure,   # meth / total
    values_from = value
  ) %>%
  mutate(
    prop    = meth / total,
    focal_cpgs   = pos %in% cpgs)
  
# calculate per-site mean methylation for the smoothing curve + drop nans
adc_sites_summary <- adc_long %>%
    group_by(pos) %>%
    summarise(mean_prop = mean(prop, na.rm = TRUE),
              .groups   = "drop")
  
# plotting meth over region
adcy4 <- ggplot() +
    # all CpGs
    geom_point(
      data = adc_long %>% filter(!focal_cpgs, !is.na(prop)), # drop nans
      aes(x = pos, y = prop)) +
    # highlighted cpgs
    geom_point(
      data = adc_long %>% filter(focal_cpgs, !is.na(prop)),  # drop nans
      aes(x = pos, y = prop),
      color = "lightblue") +
    geom_smooth(
      data = adc_sites_summary,
      aes(x = pos, y = mean_prop),
      method = "loess") +
  geom_vline(xintercept = TSS_longest, linetype = "dashed", color="red") +
  labs(x = "Position on scaffold NW_023144851.1",
         y = "Proportion of methylated reads",
    ) +  theme_bw() +
  theme(panel.grid.major = element_line(color = "grey95"),
        panel.grid.minor = element_line(color = "grey100"),
        axis.text.x  = element_text(size = 16),
        axis.text.y  = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)
        )
adcy4
ggsave("/Users/stavrt/marmot/FIGURES/adcy4_v2.png", adcy4, width = 8, height = 4, dpi = 300)

adc_site_coverage <- adc_long %>%
  group_by(pos) %>%
  summarise(
    mean_site_cov = mean(total, na.rm = TRUE),
    .groups = "drop"
  )

adc_site_coverage

############################################################
slc_raw <- read_excel("/Users/stavrt/marmot/snowmelt/model_results/promoter_regions_characterization.xlsx", 
                     sheet = 2,range = cell_rows(1:9), na = c("nan", "NA"))
# highlight sig. cs
cpg <- c(9891095)
TSS <- 9889784

# first split "site" into scaffold + position and convert position "pos" to numeric for later plotting
slc_clean <- slc_raw %>%
  separate(Site, into = c("scaffold", "pos"), sep = ":", convert = TRUE)

# long format and compute methylation proportion over region
slc_long <- slc_clean %>% 
  pivot_longer(
    cols = -c(scaffold, pos),
    names_to = c("sample", "measure"),
    names_pattern = "(.*)_(meth|total)_cytosine",
    values_to  = "value"
  ) %>%
  pivot_wider(
    names_from = measure,   # meth / total
    values_from = value
  ) %>%
  mutate(
    prop    = meth / total,
    focal_cpg   = pos %in% cpg)

# calculate per-site mean methylation for the smoothing curve + drop nans
slc_sites_summary <- slc_long %>%
  group_by(pos) %>%
  summarise(mean_prop = mean(prop, na.rm = TRUE),
            .groups   = "drop")

# plotting meth over region
ggplot() +
  # all CpGs
  geom_point(
    data = slc_long %>% filter(!focal_cpg, !is.na(prop)), # drop nans
    aes(x = pos, y = prop),
  ) +
  # highlighted cpgs
  geom_point(
    data = slc_long %>% filter(focal_cpg, !is.na(prop)),  # drop nans
    aes(x = pos, y = prop),
    color = "purple"
  ) +
  geom_vline(xintercept = TSS, linetype = "dashed", color="red") +
  labs(x = "Position on scaffold NW_023144710.1 (bp)",
       y = "Proportion methylated reads",title = "CpGs sequenced upstream (~2kb) of Slc38a1",
  ) +  theme_bw() +
  theme(panel.grid.major = element_line(color = "grey95"),
        panel.grid.minor = element_line(color = "grey100")
  )

ggsave("/Users/stavrt/marmot/FIGURES/slc38a1.png", slc, width = 6, height = 4, dpi = 300)


slc_site_coverage <- slc_long %>%
  group_by(pos) %>%
  summarise(
    mean_site_cov = mean(total, na.rm = TRUE),
    .groups = "drop"
  )

slc_site_coverage

############################################################
tarm_raw <- read_excel("/Users/stavrt/marmot/snowmelt/model_results/promoter_regions_characterization.xlsx", 
                      sheet = 3,range = cell_rows(1:6), na = c("nan", "NA"))
# highlight sig. cs
cpg_2 <- c(204417)
TSS_2 <- 204169

# first split "site" into scaffold + position and convert position "pos" to numeric for later plotting
tarm_clean <- tarm_raw %>%
  separate(Site, into = c("scaffold", "pos"), sep = ":", convert = TRUE)

# long format and compute methylation proportion over region
tarm_long <- tarm_clean %>% 
  pivot_longer(
    cols = -c(scaffold, pos),
    names_to = c("sample", "measure"),
    names_pattern = "(.*)_(meth|total)_cytosine",
    values_to  = "value"
  ) %>%
  pivot_wider(
    names_from = measure,   # meth / total
    values_from = value
  ) %>%
  mutate(
    prop    = meth / total,
    focal_cpg   = pos %in% cpg_2)

# calculate per-site mean methylation for the smoothing curve + drop nans
tarm_sites_summary <- tarm_long %>%
  group_by(pos) %>%
  summarise(mean_prop = mean(prop, na.rm = TRUE),
            .groups   = "drop")

# plotting meth over region
ggplot() +
  # all CpGs
  geom_point(
    data = tarm_long %>% filter(!focal_cpg, !is.na(prop)), # drop nans
    aes(x = pos, y = prop),
  ) +
  # highlighted cpgs
  geom_point(
    data = tarm_long %>% filter(focal_cpg, !is.na(prop)),  # drop nans
    aes(x = pos, y = prop),
    color = "lightgreen"
  ) +
  geom_vline(xintercept = TSS_2, linetype = "dashed", color="red") +
  labs(x = "Position on scaffold NW_023145038.1 (bp)",
       y = "Proportion methylated reads",title = "CpGs sequenced upstream (~2kb) of Tarm1",
  ) +  theme_bw() +
  theme(panel.grid.major = element_line(color = "grey95"),
        panel.grid.minor = element_line(color = "grey100")
  )
