
################################################################################
# Analysis: Myotis velifer landscape genetics - least cost transect analysis
#
# Script title: Create pairwise connections between genetic sample locations
#
# Script author: Eric Palm
#
################################################################################

# Load required packages
sapply(c("readr",  "dplyr", "sf"),
       require, character.only = T)

gd <- readr::read_csv("data/GeneticDist.csv") 
xy <- readr::read_csv("data/MYVECoordinates.csv") %>% 
  dplyr::filter(Seq_ID != "310308Mv33_S145")

gd_xy <-
  gd %>% 
  dplyr::select(id_1 = Indv1,
                id_2 = Indv2,
                gd = GenDisPC12) %>% 
  dplyr::inner_join(., xy %>% 
                      dplyr::select(id_1 = Seq_ID, 
                                    x_1 = x, 
                                    y_1 = y)) %>% 
  dplyr::inner_join(., xy %>% 
                      dplyr::select(id_2 = Seq_ID, 
                                    x_2 = x, 
                                    y_2 = y)) %>% 
  dplyr::mutate(euc_geog = sqrt((x_2 - x_1)^2 + (y_2 - y_1)^2)/1000)

# Save file
saveRDS(gd_xy, "lcta/data/samples_pairwise.rds")
