###### Import data from PFHpA Project for extracting TK data and import to App ####
### excel sheet from PFHpA/output/data/dose_to_serum_data.xlsx

library(tidyverse)

getwd()
tk_params_init <- readxl::read_excel("Additional files/Datasets/general TK/dose_to_serum_data.xlsx",
                                     sheet = "Matched TK Data")

tk_params <- tk_params_init %>% 
  rename(species = species_name) %>% 
  select(chem, sex, species, 
         #strain,
         #route,
         CLTot_L_kg_day, VDss_L_kg, HLe_invivo, kabs) %>% 
  #select(-c(Study, contains("mismatch"),
   #         contains("pubmed"), contains("author_year"))) %>% 
  mutate(across(c(chem, sex, species#, 
                  #strain, 
                  #route
                  ),  as.factor)) %>% 
  rename_with(str_to_title) %>% 
  rename(PFAS = Chem,
        Clearance_L_per_kg_d = Cltot_l_kg_day,
        Volume_of_Distribution_L_per_kg = Vdss_l_kg,
        Half_Life_hr = Hle_invivo,
        Absorption_Coefficient_unitless = Kabs
        ) %>%
  # mutate(Route = case_when(
  #   Route == "dermal" ~ "Dermal",
  #   T ~ Route
  # )) %>% 
  mutate(Sex = case_when(
    Sex == "F" ~ "Female",
    Sex == "M" ~ "Male"
  )) %>% 
  mutate(PFAS = case_when(
    PFAS == "PFHPA" ~ "PFHpA",
    PFAS == "PFHXA" ~ "PFHxA",
    PFAS == "PFHXS" ~ "PFHxS",
    PFAS == "PFPRA" ~ "PFPrA",
    T ~ PFAS
  )) 


saveRDS(tk_params, "Additional files/Datasets/general TK/tk_params.rds")
