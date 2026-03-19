#### UPSETR COMPARAISON DES DEG PANREGIONS = DEG PAR REGION POUR VOIR LES INTERACTIONS ENTRE LES DEG ####

#### LIBRARIES ####
library(dplyr)
library(UpSetR)

#### PATHS #### 
## INPUT ##
# A modifié et mettre le ficheir avec tous les tp pour panregion 
panregion_deg.path <- "~/CPID_multiregion/data/2__differential_expression_analysis/panregion/deseq2/design_reg_group/annotation_deg_rg/panregion_deg_rg_tp1.rds"
deglist.path <- "~/CPID_multiregion/data/2__differential_expression_analysis/deglist.Rdata"

## OUTPUT ## 
upset.path <- "~/CPID_multiregion/graphs_results/panregion/deseq2_reg_group/upset_plots/upset_deg_panregion_vs_sr_tp1.png"

#### DEG DF #### 
### PANREGION ###
panregion_deg <- readRDS(panregion_deg.path)

##LABEL UP AND DOWN ##
panregion_up <-  panregion_deg %>% filter(diffexpressed == "UP")
panregion_down <-  panregion_deg %>% filter(diffexpressed == "DOWN")

### REGIONS (SINGLE REGIONS) ### 
load(deglist.path)

##LABEL UP AND DOWN ##
BLA_up   <- deg_BLA_tp1 %>% filter(diffexpressed == "UP")
BLA_down <- deg_BLA_tp1 %>% filter(diffexpressed == "DOWN")
DRN_up   <- deg_DRN_tp1 %>% filter(diffexpressed == "UP")
DRN_down <- deg_DRN_tp1 %>% filter(diffexpressed == "DOWN")
Hb_up    <-  deg_Hb_tp1 %>% filter(diffexpressed == "UP")
Hb_down  <-  deg_Hb_tp1 %>% filter(diffexpressed == "DOWN")
Ins_up   <- deg_Ins_tp1 %>% filter(diffexpressed == "UP")
Ins_down <- deg_Ins_tp1 %>% filter(diffexpressed == "DOWN")
NAc_up   <- deg_NAc_tp1 %>% filter(diffexpressed == "UP")
NAc_down <- deg_NAc_tp1 %>% filter(diffexpressed == "DOWN")
VTA_up   <- deg_VTA_tp1 %>% filter(diffexpressed == "UP")
VTA_down <- deg_VTA_tp1 %>% filter(diffexpressed == "DOWN")


#### LIST UPSET ####
listInput <- list(panregion_up   = panregion_up$label,
                  panregion_down = panregion_down$label,
                  BLA_up   = BLA_up$label,
                  BLA_down = BLA_down$label,
                  DRN_up   = DRN_up$label,
                  DRN_down = DRN_down$label,
                  Hb_up    = Hb_up$label,
                  Hb_down  = Hb_down$label,
                  Ins_up   = Ins_up$label,
                  Ins_down = Ins_down$label,
                  NAc_up   = NAc_up$label,
                  NAc_down = NAc_down$label,
                  VTA_up   = VTA_up$label,
                  VTA_down = VTA_down$label
)

#### UPSET PLOT ####
upset_plot <- (upset(fromList(listInput), 
                     mainbar.y.label = "DEGs Intersections", 
                     sets.x.label = "Number of DEGs", order.by="freq", nsets = 14))


#### SAVE ####
png(upset.path, width = 3000, height = 2400, res = 300)
upset_plot
dev.off() 