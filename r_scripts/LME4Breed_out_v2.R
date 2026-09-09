###############################################################################
# read output from lme4breeding models run with "575k_0.80_165_count.txt"
# pre-processed methylation freq file to remove high (>0.9) and low (<0.1) meth
# final values used in manuscript here.
###############################################################################

snow <- readRDS("/Users/stavrt/Dropbox (Princeton)/RMBL_MARMOTS/575k_output/fixed_coefs_575k_snow.rds") # model1 out
sprT <- readRDS("/Users/stavrt/Dropbox (Princeton)/RMBL_MARMOTS/575k_output/fixed_coefs_575k_sprT.rds") # model2 out

dim(snow)   # an array of arrays with 5 rows (one for each fixed and random effect), 
            # 4 columns for each output parameter (estimates/effect size, SE, z, p), 
            # and 575,639 layers/array for each CpG (model) ran

age_snow <- snow[,,"age"] # repeat for age + snow effects: extract a df with 575,639 rows and 4 columns
snowmelt_snow <- snow[,,"snow_melt"]

age_sprT <- sprT[,,"age"] # repeat for temp + snow effects: extract a df with 575,639 rows and 4 columns
temp_sprT <- sprT[,,"spr_avg_temp"]

View(age_snow)
View(snowmelt_snow)
View(age_sprT)
View(temp_sprT)

# get mean summaries
summary(age_snow)       # positive effect - earlier snowmelt = increased meth in aging cytosines
summary(age_sprT)       # positive effect - lower spring temp. = increased meth in aging cytosines
summary(snowmelt_snow)  # means now both negative for snow and temp. - trend in same direction
summary(temp_sprT) 

# rename pval column for filtering
colnames(age_snow)[colnames(age_snow) == "Pr(>|z|)"] <- "p"
colnames(snowmelt_snow)[colnames(snowmelt_snow) == "Pr(>|z|)"] <- "p"
colnames(age_sprT)[colnames(age_sprT) == "Pr(>|z|)"] <- "p"
colnames(temp_sprT)[colnames(temp_sprT) == "Pr(>|z|)"] <- "p"

# output as df
age_from_snowmelt <- as.data.frame(age_snow)
snow_from_snowmelt <- as.data.frame(snowmelt_snow)
age_from_temp <- as.data.frame(age_sprT)
sprT_from_temp <- as.data.frame(temp_sprT)

# add FDR threshold column
age_from_snowmelt$FDR <- p.adjust(age_from_snowmelt$p, method = "fdr")
snow_from_snowmelt$FDR <- p.adjust(snow_from_snowmelt$p, method = "fdr")
age_from_temp$FDR <- p.adjust(age_from_temp$p, method = "fdr")
sprT_from_temp$FDR <- p.adjust(sprT_from_temp$p, method = "fdr")

# show significance in each df for 20% FDR
age_from_snowmelt$sig_fdr20 <- age_from_snowmelt$FDR < 0.20 # 22,616 sig. sites with age (31,388 in marmot array!)
snow_from_snowmelt$sig_fdr20 <- snow_from_snowmelt$FDR < 0.20
age_from_temp$sig_fdr20 <- age_from_temp$FDR < 0.20
sprT_from_temp$sig_fdr20 <- sprT_from_temp$FDR < 0.20
View(snow_from_snowmelt)

# how many sig. sites for each predictor?
sum(age_from_snowmelt$FDR < 0.20) # 31,679 sites sig. with this threshold for single model, now 9,578
sum(age_from_temp$FDR < 0.20) # 31,679 sites sig. with this threshold for single model, now 13,038
sum(snow_from_snowmelt$FDR < 0.20) # only 20 sites sig. with this threshold for single model, now 77
sum(sprT_from_temp$FDR < 0.20) # 59 sites sig. with this threshold for single model, now 209

# final draft check: are the 2 snowmelt DMS in adcy4 also significant at 10% FDR?
snow_from_snowmelt$sig_fdr10 <- snow_from_snowmelt$FDR < 0.10
sum(snow_from_snowmelt$FDR < 0.10) # 0 sites sig. at 10% fdr

# raw p values: confirm uniform distribution for snow and temp.
hist(age_from_temp$p, breaks = 100, main = "raw p-values – age x snowmelt", xlab = "p-val")
hist(snow_from_snowmelt$p, breaks = 100, main = "raw p-values – snowmelt", xlab = "p-val")
hist(sprT_from_temp$p, breaks = 100, main = "raw p-values – spring temp.", xlab = "p-val")

# write each df for site annt - must be txt file
#out_path <- "/Users/stavrt/marmot/snowmelt/model_results/"
#write.table(age_from_snowmelt, file.path(out_path, "age_snowmelt_model_fdr20.txt"), row.names=T, sep="\t")
#write.table(snow_from_snowmelt, file.path(out_path, "snowdate_snowmelt_model_fdr20.txt"), row.names=T, sep="\t")
#write.table(age_from_temp, file.path(out_path, "age_sprT_model_fdr20.txt"), row.names=T, sep="\t")
#write.table(sprT_from_temp, file.path(out_path, "tmean_sprT_model_fdr20.txt"), row.names=T, sep="\t")

###############################################################################
# QC plotting
###############################################################################
library(ggplot2)
library(cowplot)
library(patchwork)
library(lattice)
library(ggvenn)

# make a combined plot for age-related to concat the dataframes for plotting
age_from_snowmelt$source <- "snowmelt"
age_from_temp$source <- "sprT"
age_both <- bind_rows(age_from_snowmelt, age_from_temp)

# volcano plot for beta and p values - snow x age, plotting Estimate as effect size
age_vol <- ggplot(age_both, aes(x = Estimate, y = -log10(p), color = sig_fdr20)) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Age (Estimate)", y = "-log10(p-value)") +
  scale_color_discrete(
    name = " ",
    breaks = c(TRUE, FALSE),
    labels = c("<20% FDR", "n.s.")) +
  theme_minimal()
plot(age_vol)

# volcano plot for beta and p values - snowmelt
snow_vol <- ggplot(snow_from_snowmelt, aes(x = Estimate, y = -log10(p), color = sig_fdr20)) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Snowmelt date (Estimate)", y = "-log10(p-value)") +
  scale_color_discrete(
    name = " ",
    breaks = c(TRUE, FALSE),
    labels = c("<20% FDR", "n.s.")) +
  theme_minimal()
plot(snow_vol)

# volcano plot for beta and p values - spring temp
temp_vol <- ggplot(sprT_from_temp, aes(x = Estimate, y = -log10(p), color = sig_fdr20)) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Spring temperature (Estimate)", y = "-log10(p-value)") +
  scale_color_discrete(
    name = " ",
    breaks = c(TRUE, FALSE),
    labels = c("<20% FDR", "n.s.")) +
  theme_minimal()

all_vol <- age_vol + snow_vol + temp_vol + plot_layout(ncol = 3, widths = c(1, 1, 1), guides = "collect")
plot(all_vol)
ggsave("/Users/stavrt/marmot/FIGURES/volcano_combined_3.png", all_vol, width = 12, height = 4, dpi = 300)

# QQ plot: compare -log10 raw p values compared to theoretical p values from a uniform distribution
qq_snow <- qqmath(~-log10(snow_from_snowmelt$p),
                 distribution=function(x){-log10(qunif(1-x))},
                 pch = 19,
                 xlab = "Expected -log[10](p-value)", ylab = "Observed -log[10](p-value)",
                 main = list("Quantile-quantile plot of model p-values, annual 50% snowmelt date"),
                 panel = function(x, y, ...) {
                   # Draw the QQ plot points
                   panel.qqmath(x, y, ...)
                   # Add a diagonal trendline y = x
                   panel.abline(a = 0, b = 1, col = "red", lty = 3)
                 })
plot(qq_snow)

qq_sprT <- qqmath(~-log10(sprT_from_temp$p),
                      distribution=function(x){-log10(qunif(1-x))},
                      pch = 19,
                      xlab = "Expected -log[10](p-value)", ylab = "Observed -log[10](p-value)",
                      main = list("Quantile-quantile plot of model p-values, spring daily mean temperature"),
                      panel = function(x, y, ...) {
                        # Draw the QQ plot points
                        panel.qqmath(x, y, ...)
                        # Add a diagonal trendline y = x
                        panel.abline(a = 0, b = 1, col = "red", lty = 3)
                      })
plot(qq_sprT)

qq_2 <- plot_grid(qq_snow, qq_sprT, nrow=2, rel_heights=c(1,1), labels = "AUTO")

ggsave("/Users/stavrt/marmot/FIGURES/qq_2.png", qq_2, width = 12, height = 10, dpi = 300)

# lollipop plot: confirm “significant” sites are not concentrated all on the same scaffolds 
# plot fdr20 sites, p value vs.scaffold locations to look at physical clustering

# name rownames column site id
loll_plot <- cbind(site = rownames(snow_from_snowmelt), snow_from_snowmelt)
View(loll_plot)

# split scaffold and CpG position
loll_plot <- loll_plot %>%
  separate(site, into = c("scaffold", "position"), sep = ":", convert = TRUE) %>%
  mutate(logp = -log10(p),
         sig = ifelse(sig_fdr20 == TRUE, "FDR < 0.20", "NS"))

ggplot(loll_plot, aes(x = scaffold, y = p, color = sig)) +
  geom_point() +
  geom_segment(aes(xend = scaffold, yend = 0), size = 0.3) +
  geom_point() +
  facet_wrap(~ scaffold, scales = "free_x", ncol = 1) +
  scale_color_manual(values = c("FDR < 0.20" = "red", "NS" = "grey")) +
  labs(
    x = "Scaffold",
    y = expression(-log[10](p)),
    color = "FDR 20% significance") + 
  theme_minimal()

# venn diagram of significant sites pre-annt
# combine sig. age results from snow and sprT
snow_sprt_age_sig <- rbind(
  age_from_snowmelt[age_from_snowmelt$sig_fdr20 == TRUE, ],
  age_from_temp[age_from_temp$sig_fdr20 == TRUE, ]
)
# get snow and temp. significant sites separately
snow_sig <- snow_from_snowmelt[snow_from_snowmelt$sig_fdr20 == TRUE, ]
temp_sig <- sprT_from_temp[sprT_from_temp$sig_fdr20 == TRUE, ]

# add column for rownames for plot
snow_sprt_age_sig <- cbind(site = rownames(snow_sprt_age_sig), snow_sprt_age_sig)
snow_sig <- cbind(site = rownames(snow_sig), snow_sig)
temp_sig <- cbind(site = rownames(temp_sig), temp_sig)

# extract cpg sites as vectors
snow_sprt_age_sig_sites <- as.vector(snow_sprt_age_sig$site)
snow_sig_sites <- as.vector(snow_sig$site)
temp_sig_sites <- as.vector(temp_sig$site)

# age x snow x temperature
age_snow_temp <- ggvenn(
  list(age = snow_sprt_age_sig_sites,
       snowmelt = snow_sig_sites,
       temp = temp_sig_sites),
  show_percentage = FALSE,
  fill_color = c("#00BE67", "#00A9FF", "#F8766D"),
  stroke_size = 0.4,
  set_name_size = 6)

plot(age_snow_temp)

# snow x temperature 
snow_temp <- ggvenn(
  list(
    snowmelt = snow_sig_sites,
    temp = temp_sig_sites),
  show_percentage = FALSE,
  fill_color = c("#00A9FF", "#F8766D"),
  stroke_size = 0.4,
  set_name_size = 6
)
plot(snow_temp)

ggsave("/Users/stavrt/marmot/FIGURES/age_snow_temp_venndiagram.png", age_snow_temp, width = 6, height = 4, dpi = 300)
ggsave("/Users/stavrt/marmot/FIGURES/snow_temp_venndiagram.png", snow_temp, width = 6, height = 4, dpi = 300)

