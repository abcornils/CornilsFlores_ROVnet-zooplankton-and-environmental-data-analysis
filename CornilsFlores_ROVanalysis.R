## =============================================================================
## 02_run_analysis.R
##
## PS122 ROV data analysis - MOSAiC expedition
## STEP 2 of 2: Full analysis, all Tables and Figures
##
## This script ONLY requires "CornilsFlores_ROVanalysis_essential_data.RData"
## (produced by 01_prepare_data.R) - no raw data files or private helper
## scripts are needed to reproduce every result below.
##
## ------------------------------- TABLE OF CONTENTS -------------------------
##  0. Setup: libraries, load data, colour schemes
##  1. Environmental data: derive ps122_environ / ps122_environ_nMDS
##  --- FIGURE 1 --- Environmental time series + station map
##  --- FIGURE 2 --- NMDS (panel B) and PCA (panel A) of community/environment
##                    [elbow plots + dendrograms below = SUPPL. FIGURE 1]
##  --- SUPPL. TABLE 2 (part 1) --- Abundance totals per taxon/regime/depth
##  --- FIGURE 3 --- Copepod & non-copepod abundance composition
##  --- FIGURE 4 & SUPPL. FIGURE 5 --- Polar plots (panel components only;
##                    final figures were assembled manually in Inkscape)
##  --- SUPPL. FIGURE 2 --- Biovolume composition
##  --- SUPPL. TABLE 2 (part 2) --- Biovolume totals per taxon/regime/depth
##  --- SUPPL. FIGURE 3 --- Sympagic fauna (abundance + biovolume)
##  --- SUPPL. FIGURE 4 --- Calanus glacialis / Metridia longa / Paraeuchaeta
##  --- FIGURE 5 --- NBSS vs. environmental drivers
##  --- SUPPL. TABLE 7 --- NBSS slopes per sample
##  Statistics: ANOVA/post-hoc on NBSS slopes (supports Figure 5 results text)
##  --- SUPPL. TABLE 6 --- Multivariate dispersion homogeneity (PERMDISP)
##  Support calculations for Table 1: Biomass, Carbon biomass, Ingestion/Egestion
##  --- TABLE 1 --- Final summary table
##
## Note: additional supplementary tables in the manuscript (e.g. Suppl. Table 1)
## were compiled outside R and are not part of this script.
## -----------------------------------------------------------------------------

## -----------------------------------------------------------------------
## 0. Libraries
## -----------------------------------------------------------------------
library(ggOceanMaps)
library(ggplot2)
library(ggspatial)
library(lubridate)
library(patchwork)
library(readr)
library(tidyverse)
library(ggpubr)
library(dplyr)
library(vegan)
library(missMDA)
library(corrplot)
library(ggforce)
library(ggrepel)
library(emmeans)
library(scales)
# NOTE: plyr is deliberately NOT attached via library(plyr) - it is only used
# via the fully-qualified plyr::rbind.fill() calls in the Suppl. Table 2
# sections below. Attaching plyr after dplyr/tidyverse would mask dplyr's
# summarise()/mutate()/n() etc. and break downstream code (e.g. the NBSS
# calc_nbss() function, which relies on dplyr's n() inside summarise()).

## -----------------------------------------------------------------------
## Load the essential data prepared by 01_prepare_data.R
## -----------------------------------------------------------------------
#setwd("your/local/output/directory")
setwd("/Users/acornils/Desktop/MOSAiC/MOSAiC_ROV_Netze/2025_ROV_submission/2026_2nd_revision/Data for Github/")
load("CornilsFlores_ROVanalysis_essential_data.RData")

## -----------------------------------------------------------------------
## Reproducibility: fix the random seed ONCE, at the top, for the whole script
## -----------------------------------------------------------------------
## Several steps below draw random numbers and will give different results on
## every run unless the seed is fixed here: metaMDS() (random starting
## configurations), envfit() and adonis2()/PERMANOVA (permutation tests),
## betadisper()/permutest() (permutation tests), the NBSS slope bootstrap in
## Figure 5, and the kmeans() elbow plots. Run this script top-to-bottom in a
## fresh R session (do not skip/re-order chunks) to get identical output
## every time.
set.seed(42)

## -----------------------------------------------------------------------
## Colour schemes (used throughout all figures below)
## -----------------------------------------------------------------------
# regime colors
regime_colors1 <- c(
  "WIN" = "#56B4E9",
  "11" = "#56B4E9", "12" = "#56B4E9", "1" = "#56B4E9", "2" = "#56B4E9",
  "3" = "#009E73",
  "SPR" = "#009E73", "4" = "#009E73", "5" = "#009E73",
  "SUM" = "#F0E442", "6" = "#F0E442", "7" = "#F0E442",
  "AUT" = "#E69F00", "8" = "#E69F00", "9" = "#E69F00"
)

regime_base_colors <- c(
  "WIN" = "#56B4E9",
  "SPR" = "#009E73",
  "SUM" = "#F0E442",
  "AUT" = "#E69F00"
)

win_stages <- c("11", "12", "1", "2")
spr_stages <- c("3", "4", "5")
sum_stages <- c("6", "7")
aut_stages <- c("8", "9")

win_colors <- colorRampPalette(c("#FFFFFF", regime_base_colors["WIN"]))(length(win_stages) + 1)[-1]
spr_colors <- colorRampPalette(c("#FFFFFF", regime_base_colors["SPR"]))(length(spr_stages) + 1)[-1]
sum_colors <- colorRampPalette(c("#FFFFFF", regime_base_colors["SUM"]))(length(sum_stages) + 1)[-1]
aut_colors <- colorRampPalette(c("#FFFFFF", regime_base_colors["AUT"]))(length(aut_stages) + 1)[-1]

names(win_colors) <- win_stages
names(spr_colors) <- spr_stages
names(sum_colors) <- sum_stages
names(aut_colors) <- aut_stages

regime_colors2 <- c(regime_base_colors, win_colors, spr_colors, sum_colors, aut_colors)

## -----------------------------------------------------------------------
## Helper for boxplots: suppress the box for groups with < 3 observations
## -----------------------------------------------------------------------
## A boxplot summarizing fewer than 3 points is not meaningful. For every
## boxplot below, geom_boxplot() is given a filtered copy of the data (only
## groups with >= MIN_N_BOXPLOT observations) via this helper, while
## geom_jitter() keeps using the full, unfiltered data - so individual points
## are always shown, and a box is only drawn on top of them when there are
## enough points to summarize.
MIN_N_BOXPLOT <- 3

filter_min_n <- function(df, group_vars, min_n = MIN_N_BOXPLOT) {
  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::filter(dplyr::n() >= min_n) %>%
    dplyr::ungroup()
}

## =============================================================================
## 1. Environmental data: derive ps122_environ / ps122_environ_nMDS
##    (feeds Figure 2 NMDS/PCA, Figure 5 Table S1, Table 1)
## =============================================================================
ps122_environ <- merged_final

ps122_environ$sample <- paste(ps122_environ$sample, ps122_environ$depth_real, sep = "_")

ps122_environ <- column_to_rownames(ps122_environ, 'sample')
ps122_environ$site <- rownames(ps122_environ)
ps122_environ$Region <- ifelse(ps122_environ$Julian_Day > 67 & ps122_environ$Julian_Day < 100, "GR",
                         ifelse(ps122_environ$Julian_Day > 100 & ps122_environ$Julian_Day < 131, "NB",
                         ifelse(ps122_environ$Julian_Day > 178 & ps122_environ$Julian_Day < 194, "YK",
                         ifelse(ps122_environ$Julian_Day > 200 & ps122_environ$Julian_Day < 209, "FS", "AB"))))
ps122_environ$Regime <- ifelse(ps122_environ$Julian_Day > 66 & ps122_environ$Julian_Day < 140, "SPR",
                         ifelse(ps122_environ$Julian_Day > 140 & ps122_environ$Julian_Day < 210, "SUM",
                         ifelse(ps122_environ$Julian_Day > 210 & ps122_environ$Julian_Day < 260, "AUT", "WIN")))

ps122_environ$month <- month(ps122_environ$date)
ps122_environ <- ps122_environ %>% filter(!month == 6) %>% filter(!month == 7)

ps122_environ_nMDS <- select(ps122_environ, "Julian_Day", "Integrated_POC", "Integrated_Chla",
                              "Chla_ice_mg_m2", "mean_draft", "PAR", "Light",
                              "SST", "SSS", "MLD")

ps122_environ_nMDS <- ps122_environ_nMDS %>%
  mutate(
    Integrated_POC_log  = log1p(Integrated_POC),
    Integrated_Chla_log = log1p(Integrated_Chla),
    Chla_ice_mg_m2_log  = log1p(Chla_ice_mg_m2),
    MLD_log        = log1p(MLD),
    mean_draft_log = log1p(mean_draft),
    PAR_log   = log1p(PAR),
    Light_log = log1p(Light)
  )

ps122_environ_nMDS <- select(ps122_environ_nMDS,
                              -c("Integrated_POC", "Integrated_Chla", "Chla_ice_mg_m2", "MLD", "mean_draft", "PAR", "Light"))

# correlation matrix among environmental predictors (diagnostic, not a numbered figure)
ps122_environ_nMDS_cor <- cor(ps122_environ_nMDS, use = "pairwise.complete.obs")
ps122_environ_nMDS_cor[upper.tri(ps122_environ_nMDS_cor)] <- NA

ps122_environ_nMDS <- select(ps122_environ_nMDS,
                              c("Julian_Day", "MLD_log", "Integrated_Chla_log",
                                "Chla_ice_mg_m2_log", "mean_draft_log",
                                "PAR_log", "Integrated_POC_log"))

## =============================================================================
## ===========================    FIGURE 1    =================================
## Environmental time series (SST, SSS, MLD, PAR/light, sea-ice draft, chl a,
## POC) with station map (panel A).
## =============================================================================

## --- Figure 1A: station map ---------------------------------------------------
meta_schulz_filtered <- meta_schulz1 |>
  dplyr::filter(!is.na(Regime) & Regime %in% names(regime_colors1))

dt <- expand.grid(lon = c(-30, 35, 100), lat = c(78, 84, 90))

map1 <- cowplot::plot_grid(
  basemap(dt, expand.factor = 1.1, bathy.style = "rcb") +
    geom_spatial_path(
      data = meta_schulz,
      aes(x = longitude, y = latitude),
      color = "black"
    ) +
    geom_spatial_point(
      data = meta_schulz_filtered,
      aes(x = longitude, y = latitude, color = Regime)
    ) +
    scale_color_manual(values = regime_colors1) +
    theme(axis.title = element_blank(), legend.position = "none")
)
map1

## --- Figure 1B (SST, SSS) and 1D (MLD) -----------------------------------------
ps122_rov_sst <- ggplot(meta_schulz1, aes(Date, SST)) +
  geom_path(color = "blue") +
  scale_x_date(date_breaks = "1 month",
               limits = ymd(c("2019/08/15", "2020/05/14")), date_labels = "%b") +
  annotate("rect", xmin = as.Date("2019-10-01"), xmax = as.Date("2019-11-01"),
           ymin = -Inf, ymax = Inf, alpha = 1, fill = "#F5F5F5") +
  geom_vline(aes(xintercept = as.numeric(as.Date("2020-03-18"))), linetype = 2, col = "black", alpha = 0.5) +
  geom_vline(aes(xintercept = as.numeric(as.Date("2020-04-09"))), linetype = 2, col = "black", alpha = 0.5) +
  xlab("") +
  ylab("T (°C)") +
  theme_bw() +
  theme(panel.grid = element_blank())
ps122_rov_sst

ps122_rov_sss <- ggplot(meta_schulz1, aes(Date, SSS)) +
  geom_path(color = "red") +
  annotate("rect", xmin = as.Date("2019-10-01"), xmax = as.Date("2019-11-01"),
           ymin = -Inf, ymax = Inf, alpha = 1, fill = "#F5F5F5") +
  scale_x_date(date_breaks = "1 month", limits = ymd(c("2019/08/15", "2020/05/14")),
               date_labels = "%b") +
  geom_vline(aes(xintercept = as.numeric(as.Date("2020-03-18"))), linetype = 2, col = "black", alpha = 0.5) +
  geom_vline(aes(xintercept = as.numeric(as.Date("2020-04-09"))), linetype = 2, col = "black", alpha = 0.5) +
  xlab("") +
  ylab("Salinity") +
  theme_bw() +
  theme(panel.grid = element_blank())
ps122_rov_sss

ps122_rov_mld <- ggplot(meta_schulz1, aes(Date, MLD)) +
  geom_path(color = "black") +
  annotate("rect", xmin = as.Date("2019-10-01"), xmax = as.Date("2019-11-01"),
           ymin = -Inf, ymax = Inf, alpha = 1, fill = "#F5F5F5") +
  scale_x_date(date_breaks = "1 month", limits = ymd(c("2019/08/15", "2020/05/14")),
               date_labels = "%b") +
  scale_y_reverse() +
  geom_vline(aes(xintercept = as.numeric(as.Date("2020-03-18"))), linetype = 2, col = "black", alpha = 0.5) +
  geom_vline(aes(xintercept = as.numeric(as.Date("2020-04-09"))), linetype = 2, col = "black", alpha = 0.5) +
  xlab("") +
  ylab("MLD (m)") +
  theme_bw() +
  theme(panel.grid = element_blank())
ps122_rov_mld

## --- Figure 1C: PAR / under-ice light -----------------------------------------
light_colors <- c(Light = "yellow", PAR = "skyblue")

ps122_rov_par <-
  ggplot(combined_data_0m, aes(as.Date(Date), mean_value, color = source)) +
  geom_point() +
  geom_smooth(method = "loess") +
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), width = 0.2) +
  annotate("rect", xmin = as.Date("2019-10-01"), xmax = as.Date("2019-11-01"),
           ymin = -Inf, ymax = Inf, alpha = 1, fill = "#F5F5F5") +
  scale_x_date(date_breaks = "1 month", limits = ymd(c("2019/08/15", "2020/05/14")), date_labels = "%b") +
  scale_color_manual(values = light_colors) +
  geom_vline(xintercept = as.numeric(as.Date("2020-03-18")), linetype = 2, col = "black", alpha = 0.5) +
  geom_vline(xintercept = as.numeric(as.Date("2020-04-09")), linetype = 2, col = "black", alpha = 0.5) +
  xlab("") +
  ylab(expression("PAR (W" ~ m^{-2} ~ ")")) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none")
ps122_rov_par
#ggsave("ps122_rov_par.pdf", device="pdf", scale = 1, dpi=200, limitsize = TRUE, bg = NULL)

## --- Figure 1: sea-ice draft ---------------------------------------------------
ps122_rov_ice <-
  ggplot(summary_ice_10m, aes(as.Date(Date), mean_value)) +
  geom_point(color = "grey") +
  geom_smooth(method = "loess", color = "grey", se = TRUE) +
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), width = 0.2, color = "grey") +
  annotate("rect", xmin = as.Date("2019-10-01"), xmax = as.Date("2019-11-01"),
           ymin = -Inf, ymax = Inf, alpha = 1, fill = "#F5F5F5") +
  scale_x_date(date_breaks = "1 month", limits = ymd(c("2019/08/15", "2020/05/14")), date_labels = "%b") +
  geom_vline(xintercept = as.numeric(as.Date("2020-03-18")), linetype = 2, col = "black", alpha = 0.5) +
  geom_vline(xintercept = as.numeric(as.Date("2020-04-09")), linetype = 2, col = "black", alpha = 0.5) +
  xlab("") +
  ylab(expression("Draft (m)")) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none")
ps122_rov_ice
#ggsave("ps122_rov_ice.pdf", device="pdf", scale = 1.2, dpi=200, limitsize = TRUE, bg = NULL)

## --- Figure 1: chlorophyll a, water column --------------------------------------
ps122_rov_chla <-
  ggplot(chla_integrated, aes(as.Date(date), Integrated_Chla)) +
  geom_point(color = "#E69F00") +
  geom_smooth(method = "loess", color = "#E69F00") +
  annotate("rect", xmin = as.Date("2019-10-01"), xmax = as.Date("2019-11-01"),
           ymin = -Inf, ymax = Inf, alpha = 1, fill = "#F5F5F5") +
  scale_x_date(date_breaks = "1 month", limits = ymd(c("2019/08/15", "2020/05/14")), date_labels = "%b") +
  geom_vline(xintercept = as.numeric(as.Date("2020-03-18")), linetype = 2, col = "black", alpha = 0.5) +
  geom_vline(xintercept = as.numeric(as.Date("2020-04-09")), linetype = 2, col = "black", alpha = 0.5) +
  xlab("") +
  ylab(expression("Chl a (mg" ~ m^{-2} ~ ")")) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none")
ps122_rov_chla
#ggsave("ps122_rov_chla.pdf", device="pdf", scale = 1, dpi=200, limitsize = TRUE, bg = NULL)

## --- Figure 1: chlorophyll a, sea ice ---------------------------------------------
ps122_rov_chla_ice <-
  ggplot(rov_chla_ice_sel, aes(as.Date(Date), Chla_ice_mg_m2)) +
  geom_point(color = "#0072B2") +
  geom_smooth(method = "loess", color = "#0072B2") +
  annotate("rect", xmin = as.Date("2019-10-01"), xmax = as.Date("2019-11-01"),
           ymin = -Inf, ymax = Inf, alpha = 1, fill = "#F5F5F5") +
  scale_x_date(date_breaks = "1 month", limits = ymd(c("2019/08/15", "2020/05/15")), date_labels = "%b") +
  geom_vline(xintercept = as.numeric(as.Date("2020-03-18")), linetype = 2, col = "black", alpha = 0.5) +
  geom_vline(xintercept = as.numeric(as.Date("2020-04-09")), linetype = 2, col = "black", alpha = 0.5) +
  xlab("") +
  ylab(expression("Chl a (ice) mg" ~ m^{-2} ~ ")")) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none")
ps122_rov_chla_ice
#ggsave("ps122_rov_chla_ice.pdf", device="pdf", scale = 1, dpi=200, limitsize = TRUE, bg = NULL)

## --- Figure 1: POC, water column -----------------------------------------------
ps122_rov_poc <-
  ggplot(poc_integrated, aes(as.Date(date), Integrated_POC_mg)) +
  geom_point(color = "#8B4513") +
  geom_smooth(method = "loess", color = "#8B4513") +
  annotate("rect", xmin = as.Date("2019-10-01"), xmax = as.Date("2019-11-01"),
           ymin = -Inf, ymax = Inf, alpha = 1, fill = "#F5F5F5") +
  scale_x_date(date_breaks = "1 month", limits = ymd(c("2019/08/15", "2020/05/14")), date_labels = "%b") +
  geom_vline(xintercept = as.numeric(as.Date("2020-03-18")), linetype = 2, col = "black", alpha = 0.5) +
  geom_vline(xintercept = as.numeric(as.Date("2020-04-09")), linetype = 2, col = "black", alpha = 0.5) +
  xlab("") +
  ylab(expression("POC (mg" ~ m^{-2} ~ ")")) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none")
ps122_rov_poc
#ggsave("ps122_rov_poc.pdf", device="pdf", scale = 1, dpi=200, limitsize = TRUE, bg = NULL)

## --- Figure 1: final assembly ---------------------------------------------------
patch <- ps122_rov_sst + ps122_rov_sss + ps122_rov_mld + ps122_rov_par +
  ps122_rov_chla + ps122_rov_poc + ps122_rov_chla_ice + ps122_rov_ice +
  plot_layout(ncol = 2)

# FIGURE 1 (final):
map1 + patch + plot_layout(ncol = 1, heights = c(1, 2)) + plot_annotation(tag_levels = 'A')

ggsave("Figure_1.pdf", device="pdf", scale = 2, dpi=300, limitsize = TRUE, bg = NULL, width = 4, height = 5)
#Graphs are arranged in InkScape

## =============================================================================
## ===========================    FIGURE 2    =================================
## NMDS (community, panel B) and PCA (environment, panel A) ordinations.
## The elbow plots (WCSS) and cluster dendrograms below are SUPPL. FIGURE 1.
## =============================================================================

## --- NMDS ordination on zooplankton abundance (feeds Figure 2, panel B) -------
ps122_ROV_abundance_nmds <- ps122_ROV_abundance %>%
  filter(!str_detect(depth_stratum, "20 - 100 m")) %>%
  filter(!str_detect(sample, "tx0129")) %>%
  filter(!str_detect(object_annotation_category, "Unidentified Copepoda|Maxillopoda"))

ps122_ROV_abundance_nmds$month <- month(ps122_ROV_abundance_nmds$Date)
ps122_ROV_abundance_nmds <- ps122_ROV_abundance_nmds %>% filter(!month == 6) %>% filter(!month == 7)

ps122_ROV_abundance_nmds$Taxon <- ifelse(is.na(ps122_ROV_abundance_nmds$Life_stage),
                                          ps122_ROV_abundance_nmds$Species,
                                          paste(ps122_ROV_abundance_nmds$Species, ps122_ROV_abundance_nmds$Life_stage, sep = "_"))

ps122_ROV_abundance_nmds_aggr <- aggregate(abundance_m3 ~ sample + depth_real + Species, data = ps122_ROV_abundance_nmds, FUN = sum)

wide <- spread(ps122_ROV_abundance_nmds_aggr, Species, abundance_m3)

wide[is.na(wide)] <- 0
wide$sample <- paste(wide$sample, wide$depth_real, sep = "_")
wide <- wide %>% dplyr::select(-depth_real)
wide10 <- column_to_rownames(wide, 'sample')
wide10log <- log1p(wide10)

wide10log$Ctenophora <- wide10log$Ctenophora + wide10log$`Beroe cucumis` + wide10log$`Mertensia ovum`
wide10log$Eusirus <- wide10log$Eusirus + wide10log$`Eusirus holmii`
wide10log$Paraeuchaeta <- wide10log$Paraeuchaeta + wide10log$`Paraeuchaeta glacialis`
wide10log$Heterorhabdidae <- wide10log$Heterorhabdidae + wide10log$Heterorhabdus
wide10log$`Eukrohnia hamata` <- wide10log$Chaetognatha + wide10log$`Eukrohnia hamata`
wide10log$Spinocalanus <- wide10log$Spinocalanus + wide10log$`Spinocalanus longicornis`
wide10log$Cnidaria <- wide10log$Cnidaria + wide10log$Botrynema + wide10log$Sminthea
wide10log$`Scaphocalanus magnus` <- wide10log$`Scaphocalanus magnus` + wide10log$`Scaphocalanus magnus`
wide10log$Gastropoda <- wide10log$Mollusca + wide10log$Gastropoda

# remove singular occurrences and higher "non-unique" taxa
wide10log <- dplyr::select(wide10log, -c(`Beroe cucumis`, `Mertensia ovum`, `Eusirus holmii`, `Paraeuchaeta glacialis`,
                                          Chaetognatha, `Spinocalanus longicornis`, Amphipoda, Appendicularia,
                                          `Boreogadus saida`, Calanoida, Calanus, Heterorhabdus, Themisto, Botrynema,
                                          Onisimus, Sminthea, Unidentified, Mollusca, Copepoda, egg))

set.seed(42) # metaMDS() uses random starting configurations - fix seed for reproducibility
nmds_ROV <- metaMDS(comm = wide10log, distance = "bray", k = 3)
nmds_ROV

stress_value_ROV <- round(nmds_ROV$stress, 3)

data.scores_Abund <- as.data.frame(scores(nmds_ROV, "sites"))
data.scores_Abund$site <- rownames(data.scores_Abund)
data.scores_Abund <- merge(data.scores_Abund, ps122_environ, by = "site", by.y = "site")
data.scores_Abund$month <- month(data.scores_Abund$date)
data.scores_Abund <- column_to_rownames(data.scores_Abund, "site")
numeric_scores <- data.scores_Abund %>% dplyr::select(where(is.numeric))

dist_matrix_Abund <- dist(numeric_scores[, 1:3])

cluster_result_Abund <- hclust(dist_matrix_Abund, method = "ward.D")

# dendrogram (diagnostic)
#pdf("Cluster_nMDS.pdf", width = 7, height = 4)
plot(cluster_result_Abund)
rect.hclust(cluster_result_Abund, k = 3, border = "red")
#dev.off()
data.scores_Abund$Cluster <- factor(cutree(cluster_result_Abund, k = 3))

# PERMANOVA
data.scores_AB <- data.scores_Abund %>% dplyr::select(-c(NMDS1, NMDS2, NMDS3))
adonis_Abund <- adonis2(vegdist(wide10log, method = "bray") ~ Regime, data = data.scores_AB)
print(adonis_Abund)

# BETADISPER
bd_Abund <- betadisper(vegdist(sqrt(wide10log), method = "bray"), data.scores_AB$Regime)
anova(bd_Abund)

# species scores
species_scores_Abund <- as.data.frame(scores(nmds_ROV, display = "species"))
species_scores_Abund$species <- rownames(species_scores_Abund)
species_scores_Abund$magnitude <- sqrt(rowSums(species_scores_Abund[, c("NMDS1", "NMDS2")]^2))
top_species_Abund <- species_scores_Abund[order(-species_scores_Abund$magnitude), ][1:15, ]

# environmental factors (only significant)
en_Abund <- envfit(nmds_ROV, ps122_environ_nMDS, permutations = 999, na.rm = TRUE)
en_coord_cont_Abund <- as.data.frame(scores(en_Abund, "vectors")) * ordiArrowMul(en_Abund, fill = 0.07)
pvals_Abund <- as.data.frame(en_Abund$vectors$pvals)
en_coord_cont_Abund$pval <- pvals_Abund[, 1]
en_coord_cont_Abund <- en_coord_cont_Abund[en_coord_cont_Abund$pval < 0.01, ]
en_coord_cont_Abund$Variable <- rownames(en_coord_cont_Abund)

# draft NMDS plot (superseded by nmds_rov below, which is the Figure 2 panel B version)
nmds_AB <- ggplot(data = data.scores_Abund, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = factor(month, levels = c("11", "12", "1", "2", "3", "4", "5", "6", "7", "8", "9")), shape = depth_real), size = 3, ) +
  scale_color_manual(values = regime_colors2, aesthetics = c("color")) +
  geom_segment(data = en_coord_cont_Abund, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               size = 1, alpha = 0.8, colour = "blue") +
  geom_text(data = en_coord_cont_Abund, aes(x = NMDS1, y = NMDS2),
            label = rownames(en_coord_cont_Abund), colour = "black", fontface = "bold") +
  geom_text_repel(data = top_species_Abund, aes(x = NMDS1, y = NMDS2, label = species),
                  colour = "black", fontface = "italic", size = 3,
                  max.overlaps = Inf, box.padding = 0) + labs(color = "Regime (Month)") +
  theme_minimal() +
  theme(
    panel.border = element_rect(fill = NA, colour = "grey30"),
    panel.grid = element_blank(),
    axis.text = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.ticks = element_line(colour = "grey30"),
    legend.key = element_blank()
  ) + annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("Stress = ", stress_value_ROV),
           size = 4, fontface = "italic") +
  geom_mark_ellipse(aes(group = Cluster), alpha = 0.1, color = "black")
nmds_AB
#ggsave("nMDS_ROV_env_m.pdf", device="pdf", scale=1.2, dpi=300, limitsize = TRUE, bg = NULL)

## --- Elbow plot (WCSS) for cluster number selection - abundance NMDS ---------
## (part of SUPPL. FIGURE 1)
wcss <- numeric()
for (k in 1:10) {
  set.seed(42)
  kmeans_model <- kmeans(data.scores_Abund[, 1:3], centers = k, nstart = 25)
  wcss[k] <- kmeans_model$tot.withinss
}
#pdf("Elbow_nMDS.pdf", width = 7, height = 4)
plot(1:10, wcss, type = "b", pch = 19,
     xlab = "Number of clusters (k)",
     ylab = "Within-cluster sum of squares (WCSS)",
     main = "Elbow Plot")
#dev.off()

# Look for the "kink" in the curve: the point at which the reduction in WCSS is
# no longer significant enough to justify adding more clusters.

## --- Figure 2, panel B: final NMDS plot ---------------------------------------
library(vegan)
data.scores_AB <- dplyr::select(data.scores_Abund, -NMDS1, -NMDS2, -NMDS3)

# PERMANOVA with both predictors
adonis_Abund_combined <- adonis2(vegdist(wide10log, method = "bray") ~ Region + Regime, data = data.scores_AB)
print(adonis_Abund_combined)

# main effects and interactions of Region and Regime (not depth_stratum)
adonis2(vegdist(wide10log, method = "bray") ~ Region + Regime + Region:Regime, data = data.scores_AB, by = "terms")

# 3 main variables + interactions
# requires enough replicates in all combinations (each Region x Regime x depth
# should have >= 2 samples, otherwise unstable)
adonis2(vegdist(wide10log, method = "bray") ~ Region * Regime * depth_stratum,
        data = data.scores_AB,
        by = "terms")

# Final PERMANOVA
# interactions could not be used, insufficient samples
adonis2(
  vegdist(wide10log, method = "bray") ~ Regime * depth_stratum,
  data = data.scores_AB,
  by = "terms"
)

# betadisper, to check whether dispersion is homogeneous
vegan::betadisper(vegdist(wide10log, "bray"), data.scores_AB$Region)
vegan::betadisper(vegdist(wide10log, "bray"), data.scores_AB$Regime)
vegan::betadisper(vegdist(wide10log, "bray"), data.scores_AB$depth_stratum)

data.scores_Abund$depth_fill <- dplyr::recode(data.scores_Abund$depth_stratum,
                                       "0 m" = "IWI",
                                       "10 m" = "10 m")
data.scores_Abund$depth_fill <- factor(data.scores_Abund$depth_fill, levels = c("IWI", "10 m"))

depth_fill_color <- c("IWI" = "white", "10 m" = "grey")

data.scores_Abund$Region <- factor(data.scores_Abund$Region)

# FIGURE 2, panel B (final):
nmds_rov <- ggplot(data = data.scores_Abund, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(shape = Region, fill = depth_fill, color = Regime),
             size = 4, stroke = 1.2) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
  scale_fill_manual(values = depth_fill_color, name = "Depth stratum") +
  scale_color_manual(values = regime_colors2, name = "Regime") +
  geom_segment(data = en_coord_cont_Abund,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")), size = 1, color = "black") +
  geom_text(data = en_coord_cont_Abund,
            aes(x = NMDS1, y = NMDS2, label = Variable),
            colour = "black", fontface = "bold", size = 3.5) +
  geom_text_repel(data = top_species_Abund,
                  aes(x = NMDS1, y = NMDS2, label = species),
                  colour = "grey30", fontface = "italic", size = 3) +
  geom_mark_ellipse(aes(group = Cluster), color = "black") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("Stress = ", stress_value_ROV),
           size = 4, fontface = "italic") +
  labs(fill = "Depth stratum", color = "Regime", shape = "Region") +
  coord_fixed(ratio = 0.755) +
  theme_minimal() +
  theme(
    panel.border = element_rect(fill = NA, colour = "black"),
    panel.grid = element_blank(),
    axis.text = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.ticks = element_line(colour = "black"),
    legend.key = element_blank()
  )
nmds_rov
#ggsave("nMDS_ROV_env_m_Species_new.pdf", device="pdf", scale=2, dpi=300, limitsize = TRUE, bg = NULL)

## --- PCA of environmental data (feeds Figure 2, panel A) ----------------------
pca_data <- ps122_environ_nMDS

pca_data <- pca_data %>%
  dplyr::distinct(Julian_Day, .keep_all = TRUE)

pca_data <- pca_data %>% select(-Julian_Day)

nb_comp <- estim_ncpPCA(pca_data, scale = TRUE)$ncp
imputed_data <- imputePCA(pca_data, ncp = nb_comp, scale = TRUE)$completeObs
pca_result <- prcomp(imputed_data, center = TRUE, scale. = TRUE)
scores_pca <- as.data.frame(pca_result$x[, 1:2])

ps122_environ_unique <- ps122_environ %>%
  dplyr::distinct(date, .keep_all = TRUE)

scores_pca$Station <- ps122_environ_unique$site
scores_pca$Regime <- ps122_environ_unique$Regime
scores_pca$Region <- ps122_environ_unique$Region
scores_pca$depth <- ps122_environ_unique$depth_real
scores_pca$month <- month(ps122_environ_unique$date)

dist_matrix <- dist(scores_pca[, 1:2])
cluster_result <- hclust(dist_matrix, method = "ward.D2")
scores_pca$Cluster <- factor(cutree(cluster_result, k = 3))

# dendrogram (part of SUPPL. FIGURE 1)
# FIX: the original script called dev.off() here without a matching
# pdf()/png() call first (unlike the nMDS dendrogram above), which would
# error or close an unrelated graphics device. Added a matching pdf() call
# for consistency with the nMDS dendrogram above.
#pdf("Cluster_PCA.pdf", width = 7, height = 4)
plot(cluster_result)
rect.hclust(cluster_result, k = 3, border = "red")
#dev.off()

loadings <- as.data.frame(pca_result$rotation[, 1:2])
loadings$Variable <- rownames(loadings)
loadings_scaled <- loadings
loadings_scaled[, 1:2] <- loadings_scaled[, 1:2] * 5

expl_var <- round(summary(pca_result)$importance[2, 1:2] * 100, 1)
pc1_lab <- paste0("PC1 (", expl_var[1], "%)")
pc2_lab <- paste0("PC2 (", expl_var[2], "%)")

#draft PCA plot (superseded by the Figure 2 panel A version below)
pca_env <- ggplot(scores_pca, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = factor(month, levels = c("11", "12", "1", "2", "3", "4", "5", "6", "7", "8", "9")), shape = depth), size = 3) +
  scale_color_manual(values = regime_colors2, aesthetics = "color") +
  geom_segment(data = loadings_scaled,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), color = "blue") +
  geom_text(data = loadings_scaled,
            aes(x = PC1 * 1.1, y = PC2 * 1.1, label = Variable),
            color = "black", size = 4) +
  geom_mark_ellipse(aes(group = Cluster), alpha = 0.1, color = "black") +
  labs(x = pc1_lab, y = pc2_lab) + labs(color = "Month") +
  theme_minimal() +
  theme(
    panel.border = element_rect(fill = NA, colour = "grey30"),
      panel.grid = element_blank(),
    axis.text = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.ticks = element_line(colour = "grey30"),
    legend.key = element_blank()
  )
pca_env
#ggsave("PCA_ROV_env_m.png", device="png", scale=2, dpi=1000, limitsize = TRUE, bg = NULL)
#ggsave("PCA_ROV_env_m.pdf", device="pdf", scale=1.2, dpi=300, limitsize = TRUE, bg = NULL)

## --- Elbow plot (WCSS) for cluster number selection - environmental PCA ------
## (part of SUPPL. FIGURE 1)
wcss <- numeric()
for (k in 1:10) {
  set.seed(42)
  kmeans_model <- kmeans(scores_pca[, 1:2], centers = k, nstart = 25)
  wcss[k] <- kmeans_model$tot.withinss
}
#pdf("Elbow_PCA.pdf", width = 7, height = 4)
plot(1:10, wcss, type = "b", pch = 19,
     xlab = "Number of clusters (k)",
     ylab = "Within-cluster sum of squares (WCSS)",
     main = "Elbow Plot")
#dev.off()

## --- Figure 2, panel A: final PCA plot -----------------------------------------
env_scaled <- scale(pca_data)

adonis2(vegdist(env_scaled, method = "euclidean", na.rm = TRUE) ~ Regime + Region, data = ps122_environ_unique)
adonis2(vegdist(env_scaled, method = "euclidean", na.rm = TRUE) ~ Region * Regime * depth_stratum,
        data = ps122_environ_unique,
        by = "terms")

scores_pca$Region <- factor(scores_pca$Region)

# FIGURE 2, panel A (final):
pca_env <- ggplot(scores_pca, aes(x = PC1, y = PC2)) +
  geom_point(aes(
    shape = Region,
    color = Regime
  ), size = 4, stroke = 1.2) +
  scale_color_manual(values = regime_colors2, name = "Regime") +
  scale_shape_manual(values = c(21, 22, 23, 24, 25), name = "Region") +
  geom_segment(data = loadings_scaled,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", linewidth = 1) +
  geom_text(data = loadings_scaled,
            aes(x = PC1 * 1.1, y = PC2 * 1.1, label = Variable),
            color = "black", size = 4, fontface = "bold") +
  geom_mark_ellipse(aes(group = Cluster), alpha = 0.1, color = "black") +
  labs(
    x = pc1_lab,
    y = pc2_lab,
    color = "Regime",
    shape = "Region"
  ) +
  coord_fixed(ratio = 1.35) +
  theme_minimal() +
  theme(
    panel.border = element_rect(fill = NA, colour = "grey30"),
    panel.grid = element_blank(),
    axis.text = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.ticks = element_line(colour = "grey30"),
    legend.key = element_rect(fill = "white")
  )
pca_env
#ggsave("PCA_ROV_env_m_3_new.pdf", device="pdf", scale=2, dpi=300, limitsize = TRUE, bg = NULL)

# NOTE: pca_env (panel A) and nmds_rov (panel B) are produced separately here
# and combined manually into Figure 2 in Inkscape (as with all other figures).

## =============================================================================
## ===================    SUPPL. TABLE 2 (part 1: abundance)    ===============
## =============================================================================
ROV_total <- ps122_ROV_abundance %>%
  filter(!str_detect(depth_stratum, "Ridge")) %>%
  filter(!str_detect(depth_stratum, "-")) %>%
  mutate(month = lubridate::month(Date)) %>%
  mutate(year = lubridate::year(Date))
yy <- merged_final %>% dplyr::select(Julian_Day, sample)
ROV_total <- merge(ROV_total, yy, by = "sample")
ROV_total$Regime <- ifelse(ROV_total$Julian_Day > 66 & ROV_total$Julian_Day < 140, "SPR",
                     ifelse(ROV_total$Julian_Day > 140 & ROV_total$Julian_Day < 210, "SUM",
                     ifelse(ROV_total$Julian_Day > 210 & ROV_total$Julian_Day < 260, "AUT", "WIN")))

ROV_total_woNaupl <- ROV_total %>%
  filter(str_detect(object_annotation_category, "nauplii<Copepoda|Foraminifera"))

ROV_total_woNaupl_agg <- aggregate(abundance_m3 ~ sample + depth_stratum + Regime + object_annotation_category, data = ROV_total_woNaupl, FUN = sum)

ROV_total_agg <- aggregate(abundance_m3 ~ sample + depth_stratum + Regime, data = ROV_total, FUN = sum)
ROV_total_agg$object_annotation_category <- "Total"

ROV_total_rbind <- plyr::rbind.fill(ROV_total_agg, ROV_total_woNaupl_agg)

ROV_total_rbind_w <- spread(ROV_total_rbind, object_annotation_category, abundance_m3)
#write.table(ROV_total_rbind_w, "ROV_total_abundance.txt", row.names=FALSE, sep="\t")

## =============================================================================
## ===========================    FIGURE 3    =================================
## Copepod & non-copepod abundance composition (boxplots + relative composition)
## =============================================================================

## --- Figure 3: Copepoda -----------------------------------------------------
ROV_copepoda <- ps122_ROV_abundance %>% filter(!str_detect(depth_stratum, "Ridge")) %>% filter(!str_detect(depth_stratum, "-")) %>% filter(str_detect(Subclass, "Copepoda")) %>%
  mutate(month = lubridate::month(Date)) %>% mutate(year = lubridate::year(Date))
yy <- merged_final %>% dplyr::select(Julian_Day, sample)
ROV_copepoda <- merge(ROV_copepoda, yy, by = "sample")
ROV_copepoda$Regime <- ifelse(ROV_copepoda$Julian_Day > 66 & ROV_copepoda$Julian_Day < 140, "SPR", ifelse(ROV_copepoda$Julian_Day > 140 & ROV_copepoda$Julian_Day < 210, "SUM", ifelse(ROV_copepoda$Julian_Day > 210 & ROV_copepoda$Julian_Day < 260, "AUT", "WIN")))

ROV_copepoda_agg <- aggregate(abundance_m3 ~ sample + depth_stratum + Regime, data = ROV_copepoda, FUN = sum)

ROV_copepoda_agg$Regime <- factor(ROV_copepoda_agg$Regime, levels = c("AUT", "WIN", "SPR"))

ps122_rov_cop_box <- ggplot(ROV_copepoda_agg, aes(x = Regime, y = abundance_m3, fill = Regime)) +
  geom_boxplot(data = filter_min_n(ROV_copepoda_agg, c("Regime", "depth_stratum")), outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, shape = 21, color = "black") +
  facet_wrap(~ depth_stratum, scales = "fixed") +
  scale_fill_manual(values = regime_colors1) +
  labs(
    x = "Regime",
    y = expression(Abundance~(ind.~m^{-3}))
  ) +
  theme_bw() +
  theme(
      panel.grid = element_blank(),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )
ps122_rov_cop_box
#ggsave("boxplot_ROV_abund_wN.png", device="png", scale=2, dpi=1000, limitsize = TRUE, bg = NULL)
#ggsave("boxplot_ROV_abund_wN.pdf", device="pdf", scale=2, dpi=300, limitsize = TRUE, bg = NULL)

family_colors <- c("Aetideidae" = "#D55E00",
"Calanidae" = "#0072B2",
"Clausocalanidae" = "#c4dfe6",
"Cyclopinidae" = "#20948b",
"Metridinidae" = "#56B4E9",
"Oithonidae" = "#E69F00",
"Oncaeidae" = "#c05805",
"Copepoda, other" = "#cccccc",
"Copepoda, unidentified" = "darkgrey",
"Copepoda, nauplii" = "black"
)

ROV_copepoda_fam <- ROV_copepoda
ROV_copepoda_fam$Family[grepl("nauplii<Copepoda", ROV_copepoda_fam$object_annotation_category)] <- "Copepoda, nauplii"
ROV_copepoda_fam$Family[grepl("Unidentified Copepoda", ROV_copepoda_fam$object_annotation_category)] <- "Copepoda, unidentified"
ROV_copepoda_fam$Family[grepl("Copepoda<Maxillopoda", ROV_copepoda_fam$object_annotation_category)] <- "Copepoda, other"
ROV_copepoda_fam$Family[grepl("Calanoida", ROV_copepoda_fam$Family)] <- "Copepoda, other"
ROV_copepoda_fam$Family[grepl("Euchaetidae", ROV_copepoda_fam$Family)] <- "Copepoda, other"
ROV_copepoda_fam$Family[grepl("Heterorhabdidae", ROV_copepoda_fam$Family)] <- "Copepoda, other"
ROV_copepoda_fam$Family[grepl("Lucicutiidae", ROV_copepoda_fam$Family)] <- "Copepoda, other"
ROV_copepoda_fam$Family[grepl("Scolecitrichidae", ROV_copepoda_fam$Family)] <- "Copepoda, other"
ROV_copepoda_fam$Family[grepl("Spinocalanidae", ROV_copepoda_fam$Family)] <- "Copepoda, other"
ROV_copepoda_fam$Family[grepl("Ectinosomatidae", ROV_copepoda_fam$Family)] <- "Copepoda, other"
ROV_copepoda_fam$Family[grepl("Harpacticoida", ROV_copepoda_fam$Family)] <- "Copepoda, other"

ROV_copepoda_fam_agg <- aggregate(abundance_m3 ~ depth_stratum + Regime + Family, data = ROV_copepoda_fam, FUN = sum)

ROV_copepoda_fam_agg$Family <- factor(ROV_copepoda_fam_agg$Family, levels = c("Aetideidae", "Calanidae", "Clausocalanidae", "Metridinidae", "Cyclopinidae", "Oithonidae", "Oncaeidae", "Copepoda, other", "Copepoda, unidentified", "Copepoda, nauplii"))

ROV_copepoda_fam_agg$Regime <- factor(ROV_copepoda_fam_agg$Regime, levels = c("AUT", "WIN", "SPR"))

ps122_rov_cop_comp <-
  ggplot(ROV_copepoda_fam_agg, aes(x = Regime, y = abundance_m3, fill = Family)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = family_colors) +
  facet_wrap(~ depth_stratum, scale = "fixed") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Regime",
    y = "Abundance (%)"
  ) +
  theme_bw() +
  theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold")
  )
ps122_rov_cop_comp
#ggsave("comp_ROV_abund_wN.png", device="png", scale=2, dpi=1000, limitsize = TRUE, bg = NULL)
#ggsave("comp_ROV_abund_wN.pdf", device="pdf", scale=1, dpi=300, limitsize = TRUE, bg = NULL)

## --- (kept from original, not used in the manuscript) -------------------------
## Total abundance of copepods per month
if (FALSE) {
  ROV_copepoda_fam_seasonal <- aggregate(abundance_m3 ~ station + depth_stratum + Date + time + Regime, data = ROV_copepoda_fam, FUN = sum)
  ROV_copepoda_fam_seasonal <- ROV_copepoda_fam_seasonal %>%
    mutate(Date = gsub("2020-09-05", "2019-09-05", Date)) %>%
    mutate(Date = gsub("2020-09-12", "2019-09-12", Date)) %>%
    mutate(Date = gsub("2020-08-29", "2019-08-29", Date))

  ps122_abund <-
    ggplot(ROV_copepoda_fam_seasonal, aes(as.Date(Date), log1p(abundance_m3))) +
    geom_point(aes(colour = Regime, size = 3)) +
    annotate("rect", xmin = as.Date("2019-10-01"), xmax = as.Date("2019-11-01"),
             ymin = -Inf, ymax = Inf, alpha = 1, fill = "#F5F5F5") +
    scale_x_date(date_breaks = "1 month", limits = ymd(c("2019/08/15", "2020/05/14")), date_labels = "%b") +
    scale_color_manual(values = regime_colors2, name = "Regime") +
    geom_vline(xintercept = as.numeric(as.Date("2020-03-18")), linetype = 2, col = "black", alpha = 0.5) +
    geom_vline(xintercept = as.numeric(as.Date("2020-04-09")), linetype = 2, col = "black", alpha = 0.5) +
    xlab("") +
    ylab(expression("Abundance (ind." ~ m^{-3} ~ ")")) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "right")
  ps122_abund
}

## =============================================================================
## ============    FIGURE 4 & SUPPL. FIGURE 5 (polar plots)    ================
## Panel components only - the final figures were assembled manually in
## Inkscape from these panels (per the corresponding author).
## =============================================================================
ROV_all <- ps122_ROV_abundance %>% filter(!str_detect(depth_stratum, "Ridge")) %>% filter(!str_detect(depth_stratum, "-")) %>%
  mutate(month = lubridate::month(Date)) %>% mutate(year = lubridate::year(Date))
yy <- merged_final %>% dplyr::select(Julian_Day, sample)
ROV_all <- merge(ROV_all, yy, by = "sample")
ROV_all$Regime <- ifelse(ROV_all$Julian_Day > 66 & ROV_all$Julian_Day < 140, "SPR", ifelse(ROV_all$Julian_Day > 140 & ROV_all$Julian_Day < 210, "SUM", ifelse(ROV_all$Julian_Day > 210 & ROV_all$Julian_Day < 260, "AUT", "WIN")))

ROV_all$Family[grepl("nauplii<Copepoda", ROV_all$object_annotation_category)] <- "Nauplii"
ROV_all$Species[grepl("nauplii<Copepoda", ROV_all$object_annotation_category)] <- "Nauplii"

ROV_copepoda_x_agg <- aggregate(abundance_m3 ~ sample + Date + depth_stratum + Regime + Family, data = ROV_all, FUN = sum)

ROV_copepoda_x_agg <- ROV_copepoda_x_agg %>%
  mutate(Date = gsub("2020-09-05", "2019-09-05", Date)) %>%
  mutate(Date = gsub("2020-09-12", "2019-09-12", Date)) %>%
  mutate(Date = gsub("2020-08-29", "2019-08-29", Date))

ROV_copepoda_x_agg_1 <- ROV_copepoda_x_agg %>%
  dplyr::mutate(
    abundance_rel = abundance_m3 / max(abundance_m3, na.rm = TRUE),
    abundance_rel_percent = 100 * abundance_rel,
    abundance_max_dev = (abundance_m3 / max(abundance_m3, na.rm = TRUE) - 1) * 100
  ) %>%
  ungroup()

ROV_copepoda_x_agg_p1 <- ROV_copepoda_x_agg_1 %>% filter(!str_detect(Family, "Amphipoda|Calanoida|Copepoda|Appendicularia|egg|Foraminifera|Calanidae|Oithonidae|Clausocalanidae|Oncaeidae|Cyclopinidae|Oikopleuridae|Calliopiidae|Nauplii|Metridinidae|Chaetognatha|Eukrohniidae|Ectinosomatidae|Ostracoda"))

ROV_copepoda_x_agg_p1$Family <- ifelse(ROV_copepoda_x_agg_p1$Family == "Beroidae", "Ctenophora",
                             ifelse(ROV_copepoda_x_agg_p1$Family == "Mertensiidae", "Ctenophora",
                             ifelse(ROV_copepoda_x_agg_p1$Family == "Rhopalonematidae", "Cnidaria",
                             ifelse(ROV_copepoda_x_agg_p1$Family == "Halicreatidae", "Cnidaria",
                             ifelse(ROV_copepoda_x_agg_p1$Family == "Typhloscolecidae", "Polychaeta",
                             ifelse(ROV_copepoda_x_agg_p1$Family == "Gastropoda", "Mollusca",
                             ROV_copepoda_x_agg_p1$Family))))))

ROV_copepoda_x_agg_p1 <- ROV_copepoda_x_agg_p1 %>%
  dplyr::group_by(Family) %>%
  dplyr::mutate(
    abundance_rel = abundance_m3 / max(abundance_m3, na.rm = TRUE),
    abundance_rel_percent = 100 * abundance_rel,
    abundance_max_dev = (abundance_m3 / max(abundance_m3, na.rm = TRUE) - 1) * 100
  ) %>%
  ungroup()

# Figure 4 / Suppl. Figure 5, panel: by higher taxon (Family)
p1 <- ROV_copepoda_x_agg_p1 %>%
  ggplot(aes(x = as.Date(Date), y = abundance_rel_percent * 0.9, fill = Regime)) +
  geom_col(position = "identity", width = 10, alpha = 0.7) +
  annotate("rect", xmin = as.Date("2019-10-01"), xmax = as.Date("2019-11-01"),
           ymin = 0, ymax = 100 * 0.9, alpha = 1, fill = "#F5F5F5") +
  annotate("rect", xmin = as.Date("2020-06-01"), xmax = as.Date("2020-07-31"),
           ymin = 0, ymax = 100 * 0.9, alpha = 1, fill = "#F5F5F5") +
  coord_polar(theta = "x") +
  scale_y_continuous(
    trans = "sqrt",
    breaks = c(25 * 0.9, 50 * 0.9, 75 * 0.9, 100 * 0.9),
    labels = c("25", "50", "75", "100")
  ) +
  scale_fill_manual(values = regime_colors2, name = "Regime") +
  scale_x_date(
    date_breaks = "3 month",
    date_minor_breaks = "1 month",
    date_labels = "%b",
    limits = as.Date(c("2019-08-01", "2020-07-31"))
  ) +
  facet_wrap(Family ~ depth_stratum, ncol = 12, drop = FALSE) +
  theme_minimal() + ylab("") + xlab("")
p1

ROV_selected_x_agg_2 <- ROV_all %>% filter(str_detect(Species, "Calanus|Microcalanus|Oithona|Oncaeidae|Metridia|Ostracoda|Eukrohnia|Chaetognatha|Ectinosomatidae|Cyclopina|Oikopleuridae|Apherusa|Nauplii|Pseudocalanus"))

ROV_selected_x_agg_2$Species <- ifelse(ROV_selected_x_agg_2$Species == "Chaetognatha", "Eukrohnia hamata", ROV_selected_x_agg_2$Species)

ROV_selected_x_agg_p2 <- aggregate(abundance_m3 ~ sample + Date + depth_stratum + Regime + Species, data = ROV_selected_x_agg_2, FUN = sum)

ROV_selected_x_agg_p2 <- ROV_selected_x_agg_p2 %>%
  mutate(Date = gsub("2020-09-05", "2019-09-05", Date)) %>%
  mutate(Date = gsub("2020-09-12", "2019-09-12", Date)) %>%
  mutate(Date = gsub("2020-08-29", "2019-08-29", Date))

ROV_selected_x_agg_p2 <- ROV_selected_x_agg_p2 %>%
  dplyr::group_by(Species) %>%
  dplyr::mutate(
    abundance_rel = abundance_m3 / max(abundance_m3, na.rm = TRUE),
    abundance_rel_percent = 100 * abundance_rel,
    abundance_max_dev = (abundance_m3 / max(abundance_m3, na.rm = TRUE) - 1) * 100
  ) %>%
  ungroup()

# Figure 4 / Suppl. Figure 5, panel: selected species
p2 <- ROV_selected_x_agg_p2 %>%
  ggplot(aes(x = as.Date(Date), y = 0.9 * abundance_rel_percent, fill = Regime)) +
  geom_col(position = "identity", width = 10, alpha = 0.7) +
  annotate("rect", xmin = as.Date("2019-10-01"), xmax = as.Date("2019-11-01"),
           ymin = 0, ymax = 100 * 0.9, alpha = 1, fill = "#F5F5F5") +
  annotate("rect", xmin = as.Date("2020-06-01"), xmax = as.Date("2020-07-31"),
           ymin = 0, ymax = 100 * 0.9, alpha = 1, fill = "#F5F5F5") +
  coord_polar(theta = "x") +
  scale_y_continuous(
    trans = "sqrt",
    breaks = c(25 * 0.9, 50 * 0.9, 75 * 0.9, 100 * 0.9),
    labels = c("25", "50", "75", "100")
  ) +
  scale_fill_manual(values = regime_colors2, name = "Regime") +
  scale_x_date(
    date_breaks = "3 month",
    date_minor_breaks = "1 month",
    date_labels = "%b",
    limits = as.Date(c("2019-08-01", "2020-07-31"))
  ) +
  facet_wrap(Species ~ depth_stratum, ncol = 12, drop = FALSE) +
  theme_minimal() + ylab("") + xlab("")
p2
# p1 and p2 both contribute panels to Figure 4 and Suppl. Figure 5, which are
# then assembled manually in Inkscape.

## --- Figure 3: non-Copepoda -----------------------------------------------------
ROV_others <- ps122_ROV_abundance %>%
  filter(!str_detect(depth_stratum, "Ridge")) %>%
  filter(!str_detect(depth_stratum, "-")) %>%
  filter(!str_detect(Subclass, "Copepoda")) %>%
  filter(!str_detect(Family, "Foraminifera")) %>%
  filter(!str_detect(object_annotation_category, "egg")) %>%
  mutate(month = lubridate::month(Date)) %>%
  mutate(year = lubridate::year(Date))
ROV_others <- merge(ROV_others, yy, by = "sample")
ROV_others$Regime <- ifelse(ROV_others$Julian_Day > 66 & ROV_others$Julian_Day < 140, "SPR", ifelse(ROV_others$Julian_Day > 140 & ROV_others$Julian_Day < 210, "SUM", ifelse(ROV_others$Julian_Day > 210 & ROV_others$Julian_Day < 260, "AUT", "WIN")))

ROV_others_agg <- aggregate(abundance_m3 ~ sample + depth_stratum + Regime, data = ROV_others, FUN = sum)

ROV_others_agg$Regime <- factor(ROV_others_agg$Regime, levels = c("AUT", "WIN", "SPR"))

ps122_rov_nc_box <- ggplot(ROV_others_agg, aes(x = Regime, y = abundance_m3, fill = Regime)) +
  geom_boxplot(data = filter_min_n(ROV_others_agg, c("Regime", "depth_stratum")), outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, shape = 21, color = "black") +
  facet_wrap(~ depth_stratum, scales = "fixed") +
  scale_fill_manual(values = regime_colors1) +
  labs(
    x = "Regime",
    y = expression(Abundance~(ind.~m^{-3}))
  ) +
  theme_bw() +
  theme(
      panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position = "none"
  )
ps122_rov_nc_box

other_colors <- c("Amphipoda" = "#D55E00",
"Annelida" = "#0072B2",
"Appendicularia" = "#c4dfe6",
"Chaetognatha" = "#20948b",
"Cnidaria" = "#1e1f26",
"Euphausiacea" = "#56B4E9",
"Isopoda" = "#E69F00",
"Ostracoda" = "#c05805",
"Gastropoda" = "#c05",
"Zooplankton, other" = "#cccccc"
)

ROV_others_fam <- ROV_others
ROV_others_fam$Family[grepl("Amphipoda", ROV_others_fam$Order)] <- "Amphipoda"
ROV_others_fam$Family[grepl("Annelida", ROV_others_fam$Phylum)] <- "Annelida"
ROV_others_fam$Family[grepl("Appendicularia", ROV_others_fam$Class)] <- "Appendicularia"
ROV_others_fam$Family[grepl("Chaetognatha", ROV_others_fam$Phylum)] <- "Chaetognatha"
ROV_others_fam$Family[grepl("Cnidaria", ROV_others_fam$Phylum)] <- "Cnidaria"
ROV_others_fam$Family[grepl("Euphausiacea", ROV_others_fam$Order)] <- "Euphausiacea"
ROV_others_fam$Family[grepl("Isopoda", ROV_others_fam$Order)] <- "Isopoda"
ROV_others_fam$Family[grepl("Ostracoda", ROV_others_fam$Order)] <- "Ostracoda"
ROV_others_fam$Family[grepl("Gastropoda", ROV_others_fam$Order)] <- "Gastropoda"
ROV_others_fam$Family[grepl("Mollusca", ROV_others_fam$Order)] <- "Gastropoda"
ROV_others_fam$Family[grepl("Gadidae", ROV_others_fam$Family)] <- "Zooplankton, other"
ROV_others_fam$Family[grepl("Ctenophora", ROV_others_fam$Phylum)] <- "Zooplankton, other"
ROV_others_fam$Family[grepl("Unidentified", ROV_others_fam$Phylum)] <- "Zooplankton, other"

ROV_others_fam_agg <- aggregate(abundance_m3 ~ depth_stratum + Regime + Family, data = ROV_others_fam, FUN = sum)

ROV_others_fam_agg$Regime <- factor(ROV_others_fam_agg$Regime, levels = c("AUT", "WIN", "SPR"))

ps122_rov_nc_comp <-
  ggplot(ROV_others_fam_agg, aes(x = Regime, y = abundance_m3, fill = Family)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = other_colors) +
  facet_wrap(~ depth_stratum, scale = "fixed") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Regime",
    y = "Abundance (%)"
  ) +
  theme_bw() +
  theme(
      panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank()
  )
ps122_rov_nc_comp
#ggsave("comp_ROV_abund_other.png", device="png", scale=2, dpi=1000, limitsize = TRUE, bg = NULL)
#ggsave("comp_ROV_abund_other.pdf", device="pdf", scale=2, dpi=300, limitsize = TRUE, bg = NULL)

## --- Figure 3: final assembly ---------------------------------------------------
# FIGURE 3 (final):
ps122_rov_cop_box + ps122_rov_cop_comp + ps122_rov_nc_box + ps122_rov_nc_comp +
  plot_annotation(tag_levels = ("A")) +
  plot_layout(ncol = 2, guides = "collect", axis_titles = "collect_y", axes = "collect")

## =============================================================================
## ===========================    SUPPL. FIGURE 2    ============================
## Biovolume composition (Copepoda + non-Copepoda)
## =============================================================================
ROV_copepoda_b <- ps122_ROV_biovolume %>%
  filter(!str_detect(depth_stratum, "Ridge")) %>%
  filter(!str_detect(depth_stratum, "-")) %>%
  filter(str_detect(Subclass, "Copepoda")) %>%
  mutate(month = lubridate::month(Date)) %>%
  mutate(year = lubridate::year(Date))
ROV_copepoda_b <- merge(ROV_copepoda_b, yy, by = "sample")
ROV_copepoda_b$Regime <- ifelse(ROV_copepoda_b$Julian_Day > 66 & ROV_copepoda_b$Julian_Day < 140, "SPR", ifelse(ROV_copepoda_b$Julian_Day > 140 & ROV_copepoda_b$Julian_Day < 210, "SUM", ifelse(ROV_copepoda_b$Julian_Day > 210 & ROV_copepoda_b$Julian_Day < 260, "AUT", "WIN")))

ROV_copepoda_b_agg <- aggregate(biovol_maj_mm3_m3 ~ sample + depth_stratum + Regime, data = ROV_copepoda_b, FUN = sum)

ROV_copepoda_b_agg$Regime <- factor(ROV_copepoda_b_agg$Regime, levels = c("AUT", "WIN", "SPR"))

ps122_rov_cop_box_b <- ggplot(ROV_copepoda_b_agg, aes(x = Regime, y = biovol_maj_mm3_m3, fill = Regime)) +
  geom_boxplot(data = filter_min_n(ROV_copepoda_b_agg, c("Regime", "depth_stratum")), outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, shape = 21, color = "black") +
  facet_wrap(~ depth_stratum, scales = "fixed") +
  scale_fill_manual(values = regime_colors1) +
  labs(
    x = "Regime",
    y = expression(Biovolume~(mm^3~m^{-3}))
  ) +
  theme_bw() +
  theme(
      panel.grid = element_blank(),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )
ps122_rov_cop_box_b

ROV_copepoda_fam_b <- ROV_copepoda_b
ROV_copepoda_fam_b$Family[grepl("nauplii<Copepoda", ROV_copepoda_fam_b$object_annotation_category)] <- "Copepoda, nauplii"
ROV_copepoda_fam_b$Family[grepl("Unidentified Copepoda", ROV_copepoda_fam_b$object_annotation_category)] <- "Copepoda, unidentified"
ROV_copepoda_fam_b$Family[grepl("Copepoda<Maxillopoda", ROV_copepoda_fam_b$object_annotation_category)] <- "Copepoda, other"
ROV_copepoda_fam_b$Family[grepl("Calanoida", ROV_copepoda_fam_b$Family)] <- "Copepoda, other"
ROV_copepoda_fam_b$Family[grepl("Euchaetidae", ROV_copepoda_fam_b$Family)] <- "Copepoda, other"
ROV_copepoda_fam_b$Family[grepl("Heterorhabdidae", ROV_copepoda_fam_b$Family)] <- "Copepoda, other"
ROV_copepoda_fam_b$Family[grepl("Lucicutiidae", ROV_copepoda_fam_b$Family)] <- "Copepoda, other"
ROV_copepoda_fam_b$Family[grepl("Scolecitrichidae", ROV_copepoda_fam_b$Family)] <- "Copepoda, other"
ROV_copepoda_fam_b$Family[grepl("Spinocalanidae", ROV_copepoda_fam_b$Family)] <- "Copepoda, other"
ROV_copepoda_fam_b$Family[grepl("Ectinosomatidae", ROV_copepoda_fam_b$Family)] <- "Copepoda, other"
ROV_copepoda_fam_b$Family[grepl("Harpacticoida", ROV_copepoda_fam_b$Family)] <- "Copepoda, other"

ROV_copepoda_fam_agg_b <- aggregate(biovol_maj_mm3_m3 ~ depth_stratum + Regime + Family, data = ROV_copepoda_fam_b, FUN = sum)

ROV_copepoda_fam_agg_b$Family <- factor(ROV_copepoda_fam_agg_b$Family, levels = c("Aetideidae", "Calanidae", "Clausocalanidae", "Metridinidae", "Cyclopinidae", "Oithonidae", "Oncaeidae", "Copepoda, other", "Copepoda, unidentified", "Copepoda, nauplii"))

ROV_copepoda_fam_agg_b$Regime <- factor(ROV_copepoda_fam_agg_b$Regime, levels = c("AUT", "WIN", "SPR"))

ps122_rov_cop_comp_b <-
  ggplot(ROV_copepoda_fam_agg_b, aes(x = Regime, y = biovol_maj_mm3_m3, fill = Family)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = family_colors) +
  facet_wrap(~ depth_stratum, scale = "fixed") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Regime",
    y = "Biovolume (%)"
  ) +
  theme_bw() +
  theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold")
  )
ps122_rov_cop_comp_b
#ggsave("comp_ROV_abund_wN.png", device="png", scale=2, dpi=1000, limitsize = TRUE, bg = NULL)
#ggsave("comp_ROV_abund_wN.pdf", device="pdf", scale=2, dpi=300, limitsize = TRUE, bg = NULL)

## =============================================================================
## ===================    SUPPL. TABLE 2 (part 2: biovolume)    ===============
## =============================================================================
ROV_total_b <- ps122_ROV_biovolume %>%
  filter(!str_detect(depth_stratum, "Ridge")) %>%
  filter(!str_detect(depth_stratum, "-")) %>%
  mutate(month = lubridate::month(Date)) %>%
  mutate(year = lubridate::year(Date))
yy <- merged_final %>% dplyr::select(Julian_Day, sample)
ROV_total_b <- merge(ROV_total_b, yy, by = "sample")
ROV_total_b$Regime <- ifelse(ROV_total_b$Julian_Day > 66 & ROV_total_b$Julian_Day < 140, "SPR", ifelse(ROV_total_b$Julian_Day > 140 & ROV_total_b$Julian_Day < 210, "SUM", ifelse(ROV_total_b$Julian_Day > 210 & ROV_total_b$Julian_Day < 260, "AUT", "WIN")))

ROV_total_b_woNaupl <- ROV_total_b %>%
  filter(str_detect(object_annotation_category, "nauplii<Copepoda|Foraminifera"))

ROV_total_b_woNaupl_agg <- aggregate(biovol_maj_mm3_m3 ~ sample + depth_stratum + Regime + object_annotation_category, data = ROV_total_b_woNaupl, FUN = sum)

ROV_total_b_agg <- aggregate(biovol_maj_mm3_m3 ~ sample + depth_stratum + Regime, data = ROV_total_b, FUN = sum)
ROV_total_b_agg$object_annotation_category <- "Total"

ROV_total_b_rbind <- plyr::rbind.fill(ROV_total_b_agg, ROV_total_b_woNaupl_agg)
ROV_total_b_rbind_w <- spread(ROV_total_b_rbind, object_annotation_category, biovol_maj_mm3_m3)

## --- (kept from original, not used in the manuscript) -------------------------
## Total biovolume per month
if (FALSE) {
  ROV_copepoda_fam_seasonal_b <- aggregate(biovol_maj_mm3_m3 ~ station + depth_stratum + Date + time + Regime, data = ROV_copepoda_b, FUN = sum)
  ROV_copepoda_fam_seasonal_b <- ROV_copepoda_fam_seasonal_b %>%
    mutate(Date = gsub("2020-09-05", "2019-09-05", Date)) %>%
    mutate(Date = gsub("2020-09-12", "2019-09-12", Date)) %>%
    mutate(Date = gsub("2020-08-29", "2019-08-29", Date))

  ps122_biov <-
    ggplot(ROV_copepoda_fam_seasonal_b, aes(as.Date(Date), log1p(biovol_maj_mm3_m3))) +
    geom_point(aes(colour = Regime, shape = depth_stratum, size = 3)) +
    annotate("rect", xmin = as.Date("2019-10-01"), xmax = as.Date("2019-11-01"),
             ymin = -Inf, ymax = Inf, alpha = 1, fill = "#F5F5F5") +
    scale_x_date(date_breaks = "1 month", limits = ymd(c("2019/08/15", "2020/05/14")), date_labels = "%b") +
    scale_color_manual(values = regime_colors2, name = "Regime") +
    geom_vline(xintercept = as.numeric(as.Date("2020-03-18")), linetype = 2, col = "black", alpha = 0.5) +
    geom_vline(xintercept = as.numeric(as.Date("2020-04-09")), linetype = 2, col = "black", alpha = 0.5) +
    xlab("") +
    ylab(expression("Biovolume (mm3" ~ m^{-3} ~ ")")) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "right")
  ps122_biov
}

## --- non-Copepoda biovolume (panel for the composition figure below) ----------
ROV_others_b <- ps122_ROV_biovolume %>%
  filter(!str_detect(depth_stratum, "Ridge")) %>%
  filter(!str_detect(depth_stratum, "-")) %>%
  filter(!str_detect(Subclass, "Copepoda")) %>%
  filter(!str_detect(Family, "Foraminifera")) %>%
  filter(!str_detect(object_annotation_category, "egg")) %>%
  mutate(month = lubridate::month(Date)) %>%
  mutate(year = lubridate::year(Date))

ROV_others_b <- merge(ROV_others_b, yy, by = "sample")

ROV_others_b$Regime <- ifelse(ROV_others_b$Julian_Day > 66 & ROV_others_b$Julian_Day < 140, "SPR", ifelse(ROV_others_b$Julian_Day > 140 & ROV_others_b$Julian_Day < 210, "SUM", ifelse(ROV_others_b$Julian_Day > 210 & ROV_others_b$Julian_Day < 260, "AUT", "WIN")))

ROV_others_agg_b <- aggregate(biovol_maj_mm3_m3 ~ sample + depth_stratum + Regime, data = ROV_others_b, FUN = sum)

ROV_others_agg_b$Regime <- factor(ROV_others_agg_b$Regime, levels = c("AUT", "WIN", "SPR"))

ps122_rov_nc_box_b <- ggplot(ROV_others_agg_b, aes(x = Regime, y = biovol_maj_mm3_m3, fill = Regime)) +
  geom_boxplot(data = filter_min_n(ROV_others_agg_b, c("Regime", "depth_stratum")), outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, shape = 21, color = "black") +
  facet_wrap(~ depth_stratum, scales = "fixed") +
  scale_fill_manual(values = regime_colors1) +
  labs(
    x = "Regime",
    y = expression(Biovolume~(mm^3~m^{-3}))
  ) +
  theme_bw() +
  theme(
      panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position = "none"
  )
ps122_rov_nc_box_b

ROV_others_fam_b <- ROV_others_b
ROV_others_fam_b$Family[grepl("Amphipoda", ROV_others_fam_b$Order)] <- "Amphipoda"
ROV_others_fam_b$Family[grepl("Annelida", ROV_others_fam_b$Phylum)] <- "Annelida"
ROV_others_fam_b$Family[grepl("Appendicularia", ROV_others_fam_b$Class)] <- "Appendicularia"
ROV_others_fam_b$Family[grepl("Chaetognatha", ROV_others_fam_b$Phylum)] <- "Chaetognatha"
ROV_others_fam_b$Family[grepl("Cnidaria", ROV_others_fam_b$Phylum)] <- "Cnidaria"
ROV_others_fam_b$Family[grepl("Euphausiacea", ROV_others_fam_b$Order)] <- "Euphausiacea"
ROV_others_fam_b$Family[grepl("Isopoda", ROV_others_fam_b$Order)] <- "Isopoda"
ROV_others_fam_b$Family[grepl("Ostracoda", ROV_others_fam_b$Order)] <- "Ostracoda"
ROV_others_fam_b$Family[grepl("Gastropoda", ROV_others_fam_b$Order)] <- "Gastropoda"
ROV_others_fam_b$Family[grepl("Mollusca", ROV_others_fam_b$Order)] <- "Gastropoda"
ROV_others_fam_b$Family[grepl("Gadidae", ROV_others_fam_b$Family)] <- "Zooplankton, other"
ROV_others_fam_b$Family[grepl("Ctenophora", ROV_others_fam_b$Phylum)] <- "Zooplankton, other"
ROV_others_fam_b$Family[grepl("Unidentified", ROV_others_fam_b$Phylum)] <- "Zooplankton, other"

ROV_others_fam_agg_b <- aggregate(biovol_maj_mm3_m3 ~ depth_stratum + Regime + Family, data = ROV_others_fam_b, FUN = sum)

ROV_others_fam_agg_b$Regime <- factor(ROV_others_fam_agg_b$Regime, levels = c("AUT", "WIN", "SPR"))

ps122_rov_nc_comp_b <-
  ggplot(ROV_others_fam_agg_b, aes(x = Regime, y = biovol_maj_mm3_m3, fill = Family)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = other_colors) +
  facet_wrap(~ depth_stratum, scale = "fixed") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Regime",
    y = "Biovolume (%)"
  ) +
  theme_bw() +
  theme(
      panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank()
  )
ps122_rov_nc_comp_b

## --- Biovolume composition figure: final assembly ------------------------------
# SUPPL. FIGURE 2 (final):
ps122_rov_cop_box_b + ps122_rov_cop_comp_b + ps122_rov_nc_box_b + ps122_rov_nc_comp_b +
  plot_annotation(tag_levels = ("A")) +
  plot_layout(ncol = 2, guides = "collect", axis_titles = "collect_y", axes = "collect")

## =============================================================================
## ===========================    SUPPL. FIGURE 3    ===========================
## Sympagic fauna (abundance + biovolume)
## =============================================================================
ROV_sympagic <- ps122_ROV_abundance %>% filter(!str_detect(depth_stratum, "Ridge")) %>% filter(!str_detect(depth_stratum, "-")) %>% filter(str_detect(object_annotation_category, "Apherusa|Eusirus|Onisimus|Cyclopina|Harpacticoida|Ectinosomatidae|Jaschnovia|Polychaeta")) %>%
  mutate(month = lubridate::month(Date)) %>% mutate(year = lubridate::year(Date))
ROV_sympagic <- merge(ROV_sympagic, yy, by = "sample")
ROV_sympagic$Regime <- ifelse(ROV_sympagic$Julian_Day > 66 & ROV_sympagic$Julian_Day < 140, "SPR", ifelse(ROV_sympagic$Julian_Day > 140 & ROV_sympagic$Julian_Day < 210, "SUM", ifelse(ROV_sympagic$Julian_Day > 210 & ROV_sympagic$Julian_Day < 260, "AUT", "WIN")))

ROV_sympagic_b <- ps122_ROV_biovolume %>% filter(!str_detect(depth_stratum, "Ridge")) %>% filter(!str_detect(depth_stratum, "-")) %>% filter(str_detect(object_annotation_category, "Apherusa|Eusirus|Onisimus|Cyclopina|Harpacticoida|Ectinosomatidae|Jaschnovia|Polychaeta")) %>% filter(!str_detect(object_annotation_category, "part<Polychaeta")) %>%
  mutate(month = lubridate::month(Date)) %>% mutate(year = lubridate::year(Date))
ROV_sympagic_b <- merge(ROV_sympagic_b, yy, by = "sample")
ROV_sympagic_b$Regime <- ifelse(ROV_sympagic_b$Julian_Day > 66 & ROV_sympagic_b$Julian_Day < 140, "SPR", ifelse(ROV_sympagic_b$Julian_Day > 140 & ROV_sympagic_b$Julian_Day < 210, "SUM", ifelse(ROV_sympagic_b$Julian_Day > 210 & ROV_sympagic_b$Julian_Day < 260, "AUT", "WIN")))

ROV_sympagic_agg <- aggregate(abundance_m3 ~ sample + depth_stratum + Regime, data = ROV_sympagic, FUN = sum)
ROV_sympagic_agg$Regime <- factor(ROV_sympagic_agg$Regime, levels = c("AUT", "WIN", "SPR"))

ROV_sympagic_agg_b <- aggregate(biovol_maj_mm3_m3 ~ sample + depth_stratum + Regime, data = ROV_sympagic_b, FUN = sum)
ROV_sympagic_agg_b$Regime <- factor(ROV_sympagic_agg_b$Regime, levels = c("AUT", "WIN", "SPR"))

ps122_rov_symp_box <- ggplot(ROV_sympagic_agg, aes(x = Regime, y = abundance_m3, fill = Regime)) +
  geom_boxplot(data = filter_min_n(ROV_sympagic_agg, c("Regime", "depth_stratum")), outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, shape = 21, color = "black") +
  facet_wrap(~ depth_stratum, scales = "fixed") +
  scale_fill_manual(values = regime_colors1) +
  labs(
    x = "Regime",
    y = expression(Abundance~(ind.~m^{-3}))
  ) +
  theme_bw() +
  theme(
      panel.grid = element_blank(),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )
ps122_rov_symp_box

ps122_rov_symp_box_b <- ggplot(ROV_sympagic_agg_b, aes(x = Regime, y = biovol_maj_mm3_m3, fill = Regime)) +
  geom_boxplot(data = filter_min_n(ROV_sympagic_agg_b, c("Regime", "depth_stratum")), outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, shape = 21, color = "black") +
  facet_wrap(~ depth_stratum, scales = "fixed") +
  scale_fill_manual(values = regime_colors1) +
  labs(
    x = "Regime",
    y = expression(Biovolume~(mm^3 ~m^{-3}))
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position = "none"
  )
ps122_rov_symp_box_b

symp_colors <- c("Apherusa glacialis" = "#D55E00",
"Cyclopina spp." = "#0072B2",
"Ectinosomatidae" = "#c4dfe6",
"Harpacticoida, other" = "#20948b",
"Jaschnovia spp." = "#56B4E9",
"Eusirus spp." = "#E69F00",
"Onisimus spp." = "#c05805",
"Polychaeta" = "#1e1f26"
)

ROV_sympagic_a <- ROV_sympagic
ROV_sympagic_a$object_annotation_category[grepl("Eusirus", ROV_sympagic_a$Genus)] <- "Eusirus spp."
ROV_sympagic_a$object_annotation_category[grepl("Onisimus", ROV_sympagic_a$object_annotation_category)] <- "Onisimus spp."
ROV_sympagic_a$object_annotation_category[grepl("Cyclopina", ROV_sympagic_a$object_annotation_category)] <- "Cyclopina spp."
ROV_sympagic_a$object_annotation_category[grepl("Harpacticoida", ROV_sympagic_a$object_annotation_category)] <- "Harpacticoida, other"
ROV_sympagic_a$object_annotation_category[grepl("Jaschnovia", ROV_sympagic_a$object_annotation_category)] <- "Jaschnovia spp."

ROV_sympagic_ab <- ROV_sympagic_b
ROV_sympagic_ab$object_annotation_category[grepl("Eusirus", ROV_sympagic_ab$Genus)] <- "Eusirus spp."
ROV_sympagic_ab$object_annotation_category[grepl("Onisimus", ROV_sympagic_ab$object_annotation_category)] <- "Onisimus spp."
ROV_sympagic_ab$object_annotation_category[grepl("Cyclopina", ROV_sympagic_ab$object_annotation_category)] <- "Cyclopina spp."
ROV_sympagic_ab$object_annotation_category[grepl("Harpacticoida", ROV_sympagic_ab$object_annotation_category)] <- "Harpacticoida, other"
ROV_sympagic_ab$object_annotation_category[grepl("Jaschnovia", ROV_sympagic_ab$object_annotation_category)] <- "Jaschnovia spp."

ROV_sympagic_a_agg <- aggregate(abundance_m3 ~ depth_stratum + Regime + object_annotation_category, data = ROV_sympagic_a, FUN = sum)
ROV_sympagic_a_agg$object_annotation_category <- factor(ROV_sympagic_a_agg$object_annotation_category, levels = c("Apherusa glacialis", "Eusirus spp.", "Onisimus spp.", "Cyclopina spp.", "Ectinosomatidae", "Harpacticoida, other", "Jaschnovia spp.", "Polychaeta"))
ROV_sympagic_a_agg$Regime <- factor(ROV_sympagic_a_agg$Regime, levels = c("AUT", "WIN", "SPR"))

ROV_sympagic_a_agg_b <- aggregate(biovol_maj_mm3_m3 ~ depth_stratum + Regime + object_annotation_category, data = ROV_sympagic_ab, FUN = sum)
ROV_sympagic_a_agg_b$object_annotation_category <- factor(ROV_sympagic_a_agg_b$object_annotation_category, levels = c("Apherusa glacialis", "Eusirus spp.", "Onisimus spp.", "Cyclopina spp.", "Ectinosomatidae", "Harpacticoida, other", "Jaschnovia spp.", "Polychaeta"))
ROV_sympagic_a_agg_b$Regime <- factor(ROV_sympagic_a_agg_b$Regime, levels = c("AUT", "WIN", "SPR"))

ps122_rov_symp_comp <-
  ggplot(ROV_sympagic_a_agg, aes(x = Regime, y = abundance_m3, fill = object_annotation_category)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = symp_colors) +
  facet_wrap(~ depth_stratum, scale = "fixed") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Regime",
    y = "Abundance (%)"
  ) +
  theme_bw() +
  theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold")
  )
ps122_rov_symp_comp

ps122_rov_symp_comp_b <-
  ggplot(ROV_sympagic_a_agg_b, aes(x = Regime, y = biovol_maj_mm3_m3, fill = object_annotation_category)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = symp_colors) +
  facet_wrap(~ depth_stratum, scale = "fixed") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Regime",
    y = "Biovolume (%)"
  ) +
  theme_bw() +
  theme(
      panel.grid = element_blank(),
   strip.background = element_blank(),
    strip.text = element_blank()
  )
ps122_rov_symp_comp_b

# SUPPL. FIGURE 3 (final):
ps122_rov_symp_box + ps122_rov_symp_comp + ps122_rov_symp_box_b + ps122_rov_symp_comp_b +
  plot_annotation(tag_levels = ("A")) +
  plot_layout(ncol = 2, guides = "collect", axis_titles = "collect_y", axes = "collect")

## =============================================================================
## ===========================    SUPPL. FIGURE 4    ===========================
## Calanus glacialis, Metridia longa, Paraeuchaeta (abundance + biovolume by
## life stage)
## =============================================================================
calanus_Palette <- c("CIstage" = "#D55E00",
"CIIstage" = "#0072B2",
"CIIIstage" = "#20948b",
"CIVstage" = "#E69F00",
"CVstage" = "#1e1f26",
"female" = "#56B4E9",
"male" = "#c4dfe6",
"CI-CV" = "#c05805")

ROV_species <- ps122_ROV_abundance %>% mutate(Species = gsub("Paraeuchaeta glacialis", "Paraeuchaeta", Species)) %>% filter(!str_detect(depth_stratum, "Ridge")) %>% filter(!str_detect(depth_stratum, "-")) %>% filter(str_detect(object_annotation_category, "Calanus glacialis|Metridia longa|Paraeuchaeta")) %>%
  mutate(month = lubridate::month(Date)) %>% mutate(year = lubridate::year(Date))
ROV_species <- merge(ROV_species, yy, by = "sample")
ROV_species$Regime <- ifelse(ROV_species$Julian_Day > 66 & ROV_species$Julian_Day < 140, "SPR", ifelse(ROV_species$Julian_Day > 140 & ROV_species$Julian_Day < 210, "SUM", ifelse(ROV_species$Julian_Day > 210 & ROV_species$Julian_Day < 260, "AUT", "WIN")))

ROV_species_b <- ps122_ROV_biovolume %>% mutate(Species = gsub("Paraeuchaeta glacialis", "Paraeuchaeta", Species)) %>% filter(!str_detect(depth_stratum, "Ridge")) %>% filter(!str_detect(depth_stratum, "-")) %>% filter(str_detect(object_annotation_category, "Calanus glacialis|Metridia longa|Paraeuchaeta")) %>%
  mutate(month = lubridate::month(Date)) %>% mutate(year = lubridate::year(Date))
ROV_species_b <- merge(ROV_species_b, yy, by = "sample")
ROV_species_b$Regime <- ifelse(ROV_species_b$Julian_Day > 66 & ROV_species_b$Julian_Day < 140, "SPR", ifelse(ROV_species_b$Julian_Day > 140 & ROV_species_b$Julian_Day < 210, "SUM", ifelse(ROV_species_b$Julian_Day > 210 & ROV_species_b$Julian_Day < 260, "AUT", "WIN")))

ROV_species_agg <- aggregate(abundance_m3 ~ sample + depth_stratum + Species + Regime, data = ROV_species, FUN = sum)
ROV_species_agg$Regime <- factor(ROV_species_agg$Regime, levels = c("AUT", "WIN", "SPR"))

ROV_species_agg_b <- aggregate(biovol_maj_mm3_m3 ~ sample + depth_stratum + Species + Regime, data = ROV_species_b, FUN = sum)
ROV_species_agg_b$Regime <- factor(ROV_species_agg_b$Regime, levels = c("AUT", "WIN", "SPR"))

ps122_rov_sp_box <- ggplot(ROV_species_agg, aes(x = Regime, y = abundance_m3, fill = Regime)) +
  geom_boxplot(data = filter_min_n(ROV_species_agg, c("Regime", "Species", "depth_stratum")), outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, shape = 21, color = "black") +
  facet_grid(Species ~ depth_stratum, scale = "free") +
  scale_fill_manual(values = regime_colors1) +
  labs(
    x = "Regime",
    y = expression(Abundance~(ind.~m^{-3}))
  ) +
  theme_bw() +
  theme(
      panel.grid = element_blank(),
    strip.background = element_rect(fill = "lightgray"),
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_blank(),
    legend.position = "none"
  )
ps122_rov_sp_box

ps122_rov_sp_box_b <- ggplot(ROV_species_agg_b, aes(x = Regime, y = biovol_maj_mm3_m3, fill = Regime)) +
  geom_boxplot(data = filter_min_n(ROV_species_agg_b, c("Regime", "Species", "depth_stratum")), outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, shape = 21, color = "black") +
  facet_grid(Species ~ depth_stratum, scale = "free") +
  scale_fill_manual(values = regime_colors1) +
  labs(
    x = "Regime",
    y = expression(Biovolume~(mm^3 ~m^{-3}))
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    strip.text.y = element_blank(),
    legend.position = "none"
  )
ps122_rov_sp_box_b

ROV_species_a_agg <- aggregate(abundance_m3 ~ depth_stratum + Regime + Species + Life_stage, data = ROV_species, FUN = sum)
ROV_species_a_agg$Life_stage <- factor(ROV_species_a_agg$Life_stage, levels = c("CIstage", "CIIstage", "CIIIstage", "CIVstage", "CVstage", "female", "male"))

ROV_species_a_agg$Regime <- factor(ROV_species_a_agg$Regime, levels = c("AUT", "WIN", "SPR"))

ROV_species_a_agg_b <- aggregate(biovol_maj_mm3_m3 ~ depth_stratum + Regime + Species + Life_stage, data = ROV_species_b, FUN = sum)

ROV_species_a_agg_b$Life_stage <- factor(ROV_species_a_agg_b$Life_stage, levels = c("CIstage", "CIIstage", "CIIIstage", "CIVstage", "CVstage", "female", "male"))

ROV_species_a_agg_b$Regime <- factor(ROV_species_a_agg_b$Regime, levels = c("AUT", "WIN", "SPR"))

ps122_rov_sp_comp <-
  ggplot(ROV_species_a_agg, aes(x = Regime, y = abundance_m3, fill = Life_stage)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = calanus_Palette) +
  facet_grid(Species ~ depth_stratum, scale = "fixed") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Regime",
    y = "Abundance (%)"
  ) +
  theme_bw() +
  theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold"),
   strip.background = element_rect(fill = "lightgray")
  )
ps122_rov_sp_comp

ps122_rov_sp_comp_b <-
  ggplot(ROV_species_a_agg_b, aes(x = Regime, y = biovol_maj_mm3_m3, fill = Life_stage)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = calanus_Palette) +
  facet_grid(Species ~ depth_stratum, scale = "fixed") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Regime",
    y = "Biovolume (%)"
  ) +
  theme_bw() +
  theme(
      panel.grid = element_blank(),
      strip.text.y = element_text(face = "bold"),
      strip.background = element_rect(fill = "lightgray"),
    strip.text.x = element_blank()
  )
ps122_rov_sp_comp_b

# SUPPL. FIGURE 4 (final):
ps122_rov_sp_box + ps122_rov_sp_comp + ps122_rov_sp_box_b + ps122_rov_sp_comp_b +
  plot_annotation(tag_levels = ("A")) +
  plot_layout(ncol = 2, guides = "collect", axis_titles = "collect_y", axes = "collect")

## =============================================================================
## ===========================    FIGURE 5    ==================================
## Normalized Biomass Size Spectra (NBSS) vs. environmental drivers.
## Also produces SUPPL. TABLE 7 (NBSS slopes per sample).
## =============================================================================
ps122_ROV_biovolume_size <- ROV_biovolume_final %>%
  dplyr::select(sample, object_lat, object_lon, Date, time, depth_real, depth_stratum,
                Leg_text, Month, Phylum, Subphylum, Class, Subclass, Order, Family,
                Genus, Species, Life_stage, object_annotation_category,
                sample_tot_vol, acq_sub_part, abundance_m3, spher_vol_maj, biovol_maj_mm3_m3) %>%
  filter(!is.na(biovol_maj_mm3_m3)) %>%
  mutate(biovol_ind_mm3_m3 = spher_vol_maj / sample_tot_vol)

biovol_ind_size <- ps122_ROV_biovolume_size %>%
  dplyr::mutate(acq_sub_part = as.integer(acq_sub_part)) %>%
  dplyr::slice(rep(1:dplyr::n(), acq_sub_part)) %>%
  dplyr::filter(!str_detect(sample, "tx0129"),
         !str_detect(object_annotation_category, "head|part|trunk"))

ps122_environ_NBSS <- merged_final %>%
  mutate(Regime = case_when(
    Julian_Day > 66  & Julian_Day < 140 ~ "SPR",
    Julian_Day >= 140 & Julian_Day < 210 ~ "SUM",
    Julian_Day >= 210 & Julian_Day < 260 ~ "AUT",
    TRUE ~ "WIN"),
    Region = case_when(
      Julian_Day > 67  & Julian_Day < 100 ~ "GR",
      Julian_Day >= 100 & Julian_Day < 131 ~ "NB",
      Julian_Day >= 178 & Julian_Day < 194 ~ "YP",
      Julian_Day >= 200 & Julian_Day < 209 ~ "FS",
      TRUE ~ "AB"))

## NOTE: all dplyr verbs inside this function are explicitly namespaced
## (dplyr::group_by, dplyr::summarise, dplyr::n(), dplyr::mutate, dplyr::lag).
## This is deliberate: if any other attached package on your system defines
## its own group_by()/summarise()/n() (e.g. plyr, Hmisc, or similar), an
## unqualified call would silently resolve to the wrong function and n()
## would fail with "Must only be used inside data-masking verbs...". Fully
## qualifying every call here makes calc_nbss() work regardless of what else
## is attached in your R session.
calc_nbss <- function(df, noBins = 22) {
  upperBound  <- 1.25 * max(ceiling(df$biovol_ind_mm3_m3), na.rm = TRUE)
  smallestBin <- upperBound / sum(2^(0:(noBins - 1)))
  bins <- tibble::tibble(
    rightBorder = smallestBin * cumsum(2^(0:(noBins - 1)))
  ) %>%
    dplyr::mutate(interval   = paste(dplyr::lag(rightBorder, default = 0), rightBorder, sep = "-"),
           binCenter = (dplyr::lag(rightBorder, default = 0) + rightBorder) / 2)

  splitpersample <- split(df, df$sample)

  dplyr::bind_rows(lapply(splitpersample, function(x) {
    x %>%
      dplyr::mutate(binID = findInterval(biovol_ind_mm3_m3, bins$rightBorder) + 1) %>%
      dplyr::group_by(binID) %>%
      dplyr::summarise(totalBiovol = sum(biovol_ind_mm3_m3),
                elements    = dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(interval          = bins$interval[binID],
             binCenter         = bins$binCenter[binID],
             rightBorder       = bins$rightBorder[binID],
             norm_biovol       = totalBiovol / rightBorder,
             log2_norm_biovol  = log2(norm_biovol),
             ind.biomass       = log2(binCenter),
             sample            = unique(x$sample))
  }))
}

argh_0 <- calc_nbss(biovol_ind_size) %>%
  dplyr::filter(ind.biomass > -15)

ROV_MN_NBSS <- argh_0 %>%
  left_join(ps122_environ_NBSS, by = "sample") %>%
  filter(!str_detect(Regime, "SUM"))

ROV_MN_NBSS$Regime <- factor(ROV_MN_NBSS$Regime, levels = c("AUT", "WIN", "SPR"))

## --- NBSS slope fit per sample -> SUPPL. TABLE 7 -------------------------------
formulas <- ROV_MN_NBSS %>%
  group_by(sample) %>%
  nest() %>%
  mutate(
    n = map_int(data, nrow),
    var_ind = map_dbl(data, ~ var(.x$ind.biomass, na.rm = TRUE)),
    model = map(data, ~ {
      d <- .x
      if (nrow(d) >= 3 && !is.na(var(d$ind.biomass)) && var(d$ind.biomass, na.rm = TRUE) > 0) {
        lm(log2_norm_biovol ~ ind.biomass, data = d)
      } else {
        NULL
      }
    }),
    intercept = map_dbl(model, ~ if (!is.null(.x)) coef(.x)[1] else NA_real_),
    slope     = map_dbl(model, ~ if (!is.null(.x)) coef(.x)[2] else NA_real_),
    r_squared = map_dbl(model, ~ if (!is.null(.x)) summary(.x)$r.squared else NA_real_),
    formula   = if_else(!is.na(intercept) & !is.na(slope),
                        paste0("y = ", round(intercept, 2), " + ", round(slope, 2), " x, R2 = ", round(r_squared, 2)),
                        NA_character_)
  ) %>%
  select(sample, n, var_ind, intercept, slope, r_squared, formula)

ps122_environ_unique <- ps122_environ_NBSS %>% distinct(sample, .keep_all = TRUE)
slope_compar <- formulas %>% left_join(ps122_environ_unique, by = "sample")
slope_compar$Regime <- factor(slope_compar$Regime, levels = c("AUT", "WIN", "SPR"))

# SUPPL. TABLE 7 (final):
table_s1 <- slope_compar %>%
  group_by(sample, Regime, Region, depth_stratum) %>%
  summarise(slope = mean(slope, na.rm = TRUE),
            intercept = mean(intercept, na.rm = TRUE),
            r_squared = mean(r_squared, na.rm = TRUE),
            n_bins = mean(n, na.rm = TRUE), .groups = "drop")

#write.csv(table_s1, "Table_S1_NBSS_slopes.csv", row.names = FALSE)

## --- Bootstrap of NBSS slopes per regime (Figure 5) ----------------------------
fit_slope <- function(df) lm(log2_norm_biovol ~ ind.biomass, data = df)$coefficients[2]

bootstrap_slope <- function(df, nboot = 2000) {
  samples <- unique(df$sample)
  boot_slopes <- replicate(nboot, {
    samp <- sample(samples, size = length(samples), replace = TRUE)
    pooled <- df %>% filter(sample %in% samp)
    fit_slope(pooled)
  })
  tibble(
    slope_mean = mean(boot_slopes, na.rm = TRUE),
    slope_lo   = quantile(boot_slopes, 0.025, na.rm = TRUE),
    slope_hi   = quantile(boot_slopes, 0.975, na.rm = TRUE),
    slopes     = list(boot_slopes)
  )
}

set.seed(42) # bootstrap_slope() uses sample() - fix seed for reproducibility
boot_results <- ROV_MN_NBSS %>%
  group_by(Regime) %>%
  group_modify(~ bootstrap_slope(.x)) %>%
  ungroup()

## --- Figure 5: final assembly ---------------------------------------------------
nbss_plot <- ggplot(ROV_MN_NBSS, aes(ind.biomass, log2_norm_biovol)) +
  geom_point(aes(colour = Regime, shape = Region, fill = Regime), alpha = 0.8, size = 2) +
  geom_smooth(method = "lm", se = FALSE, colour = "black") +
  facet_wrap(~ Regime, scales = "fixed") +
  scale_color_manual(values = regime_colors1) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
  scale_fill_manual(values = regime_colors1) +
  labs(x = "Individual biomass (log2)",
       y = "Normalized biovolume (log2)") +
  theme_bw(base_size = 12)

slopes_plot <- ggplot(boot_results, aes(x = Regime, y = slope_mean)) +
  geom_pointrange(aes(ymin = slope_lo, ymax = slope_hi, colour = Regime), size = 1) +
  scale_color_manual(values = regime_colors1) +
  theme_bw(base_size = 12) +
  labs(y = "NBSS slope (bootstrapped)", x = "")

# FIGURE 5 (final):
figure5 <- nbss_plot / slopes_plot + plot_layout(heights = c(3, 2), guides = "collect")
figure5
#ggsave("Figure_5.pdf", figure5, width = 6, height = 4, dpi = 300)

## =============================================================================
## Statistics: ANOVA & post-hoc tests on NBSS slopes (supports Figure 5 /
## results text - not an independent Table/Figure)
## =============================================================================
slope_spring <- slope_compar %>%
  filter(Regime == "SPR")

mod_slope_spring_depth <- lm(slope ~ Region * depth_stratum, data = slope_spring)
anova(mod_slope_spring_depth)

emmeans_spring_depth <- emmeans(mod_slope_spring_depth, pairwise ~ Region | depth_stratum, adjust = "tukey")
summary(emmeans_spring_depth$contrasts)

mod_slope_spring <- lm(slope ~ Region, data = slope_spring)
anova(mod_slope_spring)
emmeans_spring <- emmeans(mod_slope_spring, pairwise ~ Region, adjust = "tukey")
summary(emmeans_spring$contrasts)

# test for seasonal differences - Region AB only
slope_AB <- slope_compar %>%
  filter(Region == "AB")

mod_slope_AB_depth <- lm(slope ~ Regime * depth_stratum, data = slope_AB)
anova(mod_slope_AB_depth)

emmeans_AB_depth <- emmeans(mod_slope_AB_depth, pairwise ~ Regime | depth_stratum, adjust = "tukey")
summary(emmeans_AB_depth$contrasts)

## ===========================    SUPPL. TABLE 6    =============================
## Multivariate dispersion homogeneity (PERMDISP), supports the PERMANOVA
## results reported for Figure 2 (community composition).
## NOTE on supplementary numbering: several supplementary tables in the
## manuscript were compiled outside R and are not part of this script. Within
## this script: Suppl. Table 2 = abundance/biovolume totals per taxon (two
## parts, see above), Suppl. Table 6 = this dispersion homogeneity output,
## Suppl. Table 7 = NBSS slopes per sample (see Figure 5 section above).
## =============================================================================
bd_region <- vegan::betadisper(
  vegdist(wide10log, method = "bray"),
  data.scores_AB$Region
)
vegan::permutest(bd_region)

bd_regime <- vegan::betadisper(
  vegdist(wide10log, method = "bray"),
  data.scores_AB$Regime
)
vegan::permutest(bd_regime)

bd_depth <- vegan::betadisper(
  vegdist(wide10log, method = "bray"),
  data.scores_AB$depth_stratum
)
vegan::permutest(bd_depth)

## =============================================================================
## Support calculations for TABLE 1
## =============================================================================

## --- Biomass (dry weight) -------------------------------------------------------
## NOTE: ROV_BM_agg is not referenced again later in the script / does not
## feed a plot - kept for completeness, please verify if still needed.
ROV_BM <- ps122_ROV_biomass_final %>% filter(!str_detect(depth_stratum, "Ridge")) %>%
  filter(!str_detect(depth_stratum, "-")) %>%
  mutate(month = lubridate::month(Date)) %>%
  mutate(year = lubridate::year(Date))
ROV_BM <- merge(ROV_BM, yy, by = "sample")
ROV_BM$Regime <- ifelse(ROV_BM$Julian_Day > 66 & ROV_BM$Julian_Day < 140, "SPR", ifelse(ROV_BM$Julian_Day > 140 & ROV_BM$Julian_Day < 210, "SUM", ifelse(ROV_BM$Julian_Day > 210 & ROV_BM$Julian_Day < 260, "AUT", "WIN")))

ROV_BM_agg <- aggregate(DW_CF_mg_m3 ~ sample + depth_stratum + Regime, data = ROV_BM, FUN = sum)
ROV_BM_agg$Regime <- factor(ROV_BM_agg$Regime, levels = c("AUT", "WIN", "SPR"))

## --- Carbon biomass, per sample/depth (supporting/diagnostic figure) -----------
ROV_C <- ps122_ROV_carbon_final %>% filter(!str_detect(depth_stratum, "Ridge")) %>%
  filter(!str_detect(depth_stratum, "-")) %>%
  mutate(month = lubridate::month(Date)) %>%
  mutate(year = lubridate::year(Date))
ROV_C <- merge(ROV_C, yy, by = "sample")
ROV_C$Regime <- ifelse(ROV_C$Julian_Day > 66 & ROV_C$Julian_Day < 140, "SPR", ifelse(ROV_C$Julian_Day > 140 & ROV_C$Julian_Day < 210, "SUM", ifelse(ROV_C$Julian_Day > 210 & ROV_C$Julian_Day < 260, "AUT", "WIN")))

ROV_C_agg <- aggregate(Carbon_CF_mg_m3 ~ sample + depth_stratum + Regime, data = ROV_C, FUN = sum)
ROV_C_agg$Regime <- factor(ROV_C_agg$Regime, levels = c("AUT", "WIN", "SPR"))

ps122_rov_c_box <- ggplot(ROV_C_agg, aes(x = Regime, y = Carbon_CF_mg_m3, fill = Regime)) +
  geom_boxplot(data = filter_min_n(ROV_C_agg, c("Regime", "depth_stratum")), outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, shape = 21, color = "black") +
  facet_wrap(~ depth_stratum, scales = "fixed") +
  scale_fill_manual(values = regime_colors1) +
  labs(x = "Regime",
    y = expression("Biomass (mg C" ~m^{-3}~")")) +
  theme_bw() +
  theme(
      panel.grid = element_blank(),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )
ps122_rov_c_box
#ggsave("PS122_ROV_Carbon_2.pdf", device="pdf", scale = 2, dpi=300, limitsize = TRUE, bg = NULL, width=5, height = 5)

## --- Carbon biomass, depth-integrated per event -> feeds TABLE 1 --------------
ROV_C_agg_int <- aggregate(Carbon_CF_mg_m3 ~ station + time + sample + depth_real + depth_stratum + Regime + Date, data = ROV_C, FUN = sum)

# manual sample -> ROV-dive-event lookup table (used here and in the
# Ingestion/Egestion section below)
sample_event_df <- data.frame(
  sample = c("tx0025", "tx0047", "tx0048", "tx0050", "tx0051", "tx0053", "tx0055", "tx0057", "tx0059",
             "tx0084", "tx0085", "tx0100", "tx0113", "tx0114", "tx0127", "tx0128", "tx0144", "tx0145",
             "tx0156", "tx0157", "tx0169", "tx0170", "tx0178", "tx0179", "tx0194", "tx0195", "tx0201",
             "tx0208", "tx0209", "tx0215", "tx0216", "tx0223", "tx0224", "tx0237", "tx0238", "tx0240",
             "tx0241", "tx0254", "tx0255", "tx0279", "tx0281", "tx0301", "tx0302", "tx0325", "tx0326",
             "tx0328", "tx0329", "tx0331", "tx0332", "tx0334", "tx0335", "tx0337", "tx0341", "tx0368",
             "tx0369", "tx0388", "tx0393", "tx0410", "tx0422", "tx0436", "tx0437", "tx0458"),
  event = c(1, 2, 2, 3, 3, 4, 4, 5, 5,
            6, 6, 7, 8, 8, 9, 9, 10, 10,
            11, 11, 12, 12, 13, 13, 14, 14, 15,
            16, 16, 17, 17, 18, 18, 19, 19, 20,
            20, 21, 21, 22, 22, 23, 23, 24, 24,
            25, 25, 26, 26, 27, 27, 28, 28, 29,
            29, 30, 30, 31, 32, 33, 33, 32)
)

ROV_C_agg_int <- merge(ROV_C_agg_int, sample_event_df, by = "sample") %>%
  mutate(depth_real = as.numeric(depth_real))

mean_0m_per_regime <- ROV_C_agg_int %>%
  mutate(depth_real = as.numeric(depth_real)) %>%
  filter(depth_real == 0) %>%
  dplyr::group_by(Regime) %>%
  dplyr::summarise(mean_Carbon_0m = mean(Carbon_CF_mg_m3, na.rm = TRUE), .groups = "drop")

events_missing_0m <- ROV_C_agg_int %>%
  group_by(event) %>%
  filter(!any(depth_real == 0)) %>%
  distinct(event, Regime) %>%
  ungroup()

imputed_0m_rows <- events_missing_0m %>%
  left_join(mean_0m_per_regime, by = "Regime") %>%
  mutate(
    depth_real = 0,
    Carbon_CF_mg_m3 = mean_Carbon_0m
  ) %>%
  select(-mean_Carbon_0m)

ROV_C_complete <- bind_rows(ROV_C_agg_int, imputed_0m_rows)

integrated_carbon <- function(df) {
  df <- df %>%
    mutate(depth_real = as.numeric(depth_real)) %>%
    arrange(depth_real)

  integrate_value <- sum(8 * (head(df$Carbon_CF_mg_m3, -1) + tail(df$Carbon_CF_mg_m3, -1)) / 2)

  return(data.frame(
    event = unique(df$event),
    Regime = unique(df$Regime),
    integrated_carbon = integrate_value
  ))
}

# carbon_integrated is re-used further below in the Table 1 section
carbon_integrated <- ROV_C_complete %>%
  group_by(event) %>%
  group_split() %>%
  map_df(integrated_carbon)

average_carbon_regime <- carbon_integrated %>%
  dplyr::group_by(Regime) %>%
  dplyr::summarise(median_int_Carbon_mg_C_m2 = median(integrated_carbon, na.rm = TRUE),
                   min_int_Carbon_mg_C_m2 = min(integrated_carbon, na.rm = TRUE),
                   max_int_Carbon_mg_C_m2 = max(integrated_carbon, na.rm = TRUE), .groups = "drop")

ps122_rov_int_box <- ggplot(carbon_integrated, aes(x = Regime, y = integrated_carbon, fill = Regime)) +
   geom_boxplot(data = filter_min_n(carbon_integrated, c("Regime")), outlier.shape = NA) +
   geom_jitter(width = 0.2, alpha = 0.6, size = 2, shape = 21, color = "black") +
   scale_fill_manual(values = regime_colors1) +
   labs(x = "Regime",
     y = expression("Biomass (mg C" ~m^{-2}~")")) +
   theme_bw() +
   theme(
       panel.grid = element_blank(),
     strip.background = element_rect(fill = "lightgray"),
     strip.text = element_text(face = "bold"),
    legend.position = "none"
  )
ps122_rov_int_box

## --- Ingestion / Egestion (carbon demand), feeds TABLE 1 -----------------------
## IEB - Ingestion/Egestion Balance calculated after Bode et al. 2017 /
## Ikeda 2014 (respiration, gut fluorescence method - "GM")
ROV_Ingestion_IEB <- ps122_ROV_all_final %>% filter(!str_detect(depth_stratum, "Ridge")) %>%
  filter(!str_detect(depth_stratum, "-")) %>%
  mutate(month = lubridate::month(Date)) %>%
  mutate(year = lubridate::year(Date))
ROV_Ingestion_IEB <- merge(ROV_Ingestion_IEB, yy, by = "sample")
ROV_Ingestion_IEB$Regime <- ifelse(ROV_Ingestion_IEB$Julian_Day > 66 & ROV_Ingestion_IEB$Julian_Day < 140, "SPR", ifelse(ROV_Ingestion_IEB$Julian_Day > 140 & ROV_Ingestion_IEB$Julian_Day < 210, "SUM", ifelse(ROV_Ingestion_IEB$Julian_Day > 210 & ROV_Ingestion_IEB$Julian_Day < 260, "AUT", "WIN")))

ROV_Ingestion_IEB_agg <- aggregate(Ingestion_IEB_mug_C_m3_d_GM ~ station + time + sample + depth_real + depth_stratum + Regime + Date, data = ROV_Ingestion_IEB, FUN = sum)

ROV_Ingestion_IEB_agg$Regime <- factor(ROV_Ingestion_IEB_agg$Regime, levels = c("WIN", "SPR", "SUM", "AUT"))

ROV_Ingestion_IEB_agg_int <- merge(ROV_Ingestion_IEB_agg, sample_event_df, by = "sample") %>%
  mutate(depth_real = as.numeric(depth_real))

mean_0m_demand_per_regime <- ROV_Ingestion_IEB_agg_int %>%
  mutate(depth_real = as.numeric(depth_real)) %>%
  filter(depth_real == 0) %>%
  dplyr::group_by(Regime) %>%
  dplyr::summarise(mean_demand_0m = mean(Ingestion_IEB_mug_C_m3_d_GM, na.rm = TRUE), .groups = "drop")

events_demand_missing_0m <- ROV_Ingestion_IEB_agg_int %>%
  group_by(event) %>%
  filter(!any(depth_real == 0)) %>%
  distinct(event, Regime) %>%
  ungroup()

imputed_demand_0m_rows <- events_demand_missing_0m %>%
  left_join(mean_0m_demand_per_regime, by = "Regime") %>%
  mutate(
    depth_real = 0,
    Ingestion_IEB_mug_C_m3_d_GM = mean_demand_0m
  ) %>%
  select(-mean_demand_0m)

ROV_demand_complete <- bind_rows(ROV_Ingestion_IEB_agg_int, imputed_demand_0m_rows)

integrated_demand <- function(df) {
  df <- df %>%
    mutate(depth_real = as.numeric(depth_real)) %>%
    arrange(depth_real)

  integrate_demand_value <- sum(8 * (head(df$Ingestion_IEB_mug_C_m3_d_GM, -1) + tail(df$Ingestion_IEB_mug_C_m3_d_GM, -1)) / 2)

  return(data.frame(
    event = unique(df$event),
    Regime = unique(df$Regime),
    integrated_demand = integrate_demand_value
  ))
}

# demand_integrated is re-used further below in the Table 1 section
demand_integrated <- ROV_demand_complete %>%
  group_by(event) %>%
  group_split() %>%
  map_df(integrated_demand)

demand_integrated$Ingestion_IEB_mg_C_m2_d_GM <- demand_integrated$integrated_demand / 1000

average_demand_regime <- demand_integrated %>%
  dplyr::group_by(Regime) %>%
  dplyr::summarise(median_int_Carbon_demand_mg_C_m2_d = median(Ingestion_IEB_mg_C_m2_d_GM, na.rm = TRUE),
                   min_int_Carbon_demand_mg_C_m2_d = min(Ingestion_IEB_mg_C_m2_d_GM, na.rm = TRUE),
                   max_int_Carbon_demand_mg_C_m2_d = max(Ingestion_IEB_mg_C_m2_d_GM, na.rm = TRUE),
                   .groups = "drop")

ps122_rov_dem_b_box <- ggplot(demand_integrated, aes(x = Regime, y = Ingestion_IEB_mg_C_m2_d_GM, fill = Regime)) +
  geom_boxplot(data = filter_min_n(demand_integrated, c("Regime")), outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, shape = 21, color = "black") +
  scale_fill_manual(values = regime_colors1) +
  labs(x = "Regime",
       y = expression("Carbon demand (mg C" ~ m^-2 ~ d^-1 ~ ")")) +
  geom_text(label = "calculated after Ikeda 2014 / Bode et al. 2018",
            size = 4, y = 60, x = 2) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 10)) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )
ps122_rov_dem_b_box

## =============================================================================
## ===========================    TABLE 1    ====================================
## Final summary table: POC standing stock, integrated carbon biomass, carbon
## demand, turnover time and food supply ratio per regime.
## =============================================================================
merged_poc$Regime <- ifelse(merged_poc$Date < "2019-11-01", "AUT",
                                ifelse(merged_poc$Date > "2020-03-01", "SPR", "WIN"))

merged_poc_x <- merge(merged_poc, sample_event_df, by = "sample") %>%
  mutate(depth_real = as.numeric(depth_real))

merged_poc_x <- merged_poc_x %>%
  dplyr::distinct(event, .keep_all = TRUE)

average_POC_regime <- merged_poc_x %>%
  dplyr::group_by(Regime) %>%
  dplyr::summarise(median_int_POC_0_20m_mg_C_m2 = median(Integrated_POC_mg, na.rm = TRUE),
                   min_int_POC_0_20m_mg_C_m2 = min(Integrated_POC_mg, na.rm = TRUE),
                   max_int_POC_0_20m_mg_C_m2 = max(Integrated_POC_mg, na.rm = TRUE),
                   .groups = "drop")

table1_gesamt <- merged_poc_x %>%
  left_join(carbon_integrated, by = "event") %>%
  left_join(demand_integrated, by = "event")

# FIX: the original script referenced the undefined/typo'd column
# "Ingestion_IEE_mg_C_m2_d"; corrected to "Ingestion_IEB_mg_C_m2_d_GM"
# (confirmed).
table1_gesamt$turnover <- table1_gesamt$Ingestion_IEB_mg_C_m2_d_GM / (table1_gesamt$Integrated_POC_mg / 2)
table1_gesamt$food_supply <- (table1_gesamt$Integrated_POC_mg / 2) / table1_gesamt$integrated_carbon

# TABLE 1 (final):
average_table1_gesamt <- table1_gesamt %>%
  dplyr::group_by(Regime) %>%
  dplyr::summarise(median_int_POC_0_20m_mg_C_m2 = median(Integrated_POC_mg, na.rm = TRUE),
                   min_int_POC_0_20m_mg_C_m2 = min(Integrated_POC_mg, na.rm = TRUE),
                   max_int_POC_0_20m_mg_C_m2 = max(Integrated_POC_mg, na.rm = TRUE),
                   median_int_Carbon_mg_C_m2 = median(integrated_carbon, na.rm = TRUE),
                   min_int_Carbon_mg_C_m2 = min(integrated_carbon, na.rm = TRUE),
                   max_int_Carbon_mg_C_m2 = max(integrated_carbon, na.rm = TRUE),
                   median_int_Carbon_demand_mg_C_m2_d = median(Ingestion_IEB_mg_C_m2_d_GM, na.rm = TRUE),
                   min_int_Carbon_demand_mg_C_m2_d = min(Ingestion_IEB_mg_C_m2_d_GM, na.rm = TRUE),
                   max_int_Carbon_demand_mg_C_m2_d = max(Ingestion_IEB_mg_C_m2_d_GM, na.rm = TRUE),
                   median_int_turnover_d = median(turnover, na.rm = TRUE),
                   min_int_turnover_d = min(turnover, na.rm = TRUE),
                   max_int_turnover_d = max(turnover, na.rm = TRUE),
                   median_int_food_supply = median(food_supply, na.rm = TRUE),
                   min_int_food_supply = min(food_supply, na.rm = TRUE),
                   max_int_food_supply = max(food_supply, na.rm = TRUE),
                   .groups = "drop")

t_average_table1_gesamt <- t(average_table1_gesamt)
t_average_table1_gesamt
#write.csv(average_table1_gesamt, "Table_1.csv", row.names = FALSE)
