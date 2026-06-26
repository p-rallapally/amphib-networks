
library(tidyverse)
library(igraph)
library(ggraph)
library(janitor)
library(lubridate)
library(vegan)
library(knitr)
library(dplyr)

#lists column names so I know what to join by

files <- list.files("Desktop/briggs/networks/kacie_csvs", pattern="\\.csv$", full.names=TRUE)

setNames(
  lapply(files, \(f) names(read.csv(f, nrows=0))),
  basename(files)
)

#stage 1: read + standardize data

amphib_dissect <- read_csv("amphib_dissect.csv") %>%
  clean_names()

amphib_parasite <- read_csv("amphib_parasite.csv") %>%
  clean_names()

bd <- read_csv("bd.csv") %>%
  clean_names()

lut_amphib_parasite <- read_csv("lut_amphib_parasite.csv") %>%
  clean_names()

lut_freeliving <- read_csv("lut_freelving.csv") %>%
  clean_names()

netting_info <- read_csv("netting_infor.csv") %>%
  clean_names()

site_info <- read_csv("Site_information_.csv") %>%
  clean_names()

survey_spp <- read_csv("survey_spp.csv") %>%
  clean_names()

transect_sppsum_meta <- read_csv("transect_SppSUM__META_.csv") %>%
  clean_names()

visual_spp <- read_csv("visual_spp.csv") %>%
  clean_names()

wetland_info <- read_csv("wetland_infor.csv") %>%
  clean_names()

#structure check

list(
  amphib_dissect = dim(amphib_dissect),
  amphib_parasite = dim(amphib_parasite),
  bd = dim(bd),
  site_info = dim(site_info),
  survey_spp = dim(survey_spp),
  wetland_info = dim(wetland_info)
)

#checking clean names
names(amphib_dissect)
names(amphib_parasite)
names(bd)
names(site_info)
names(survey_spp)
names(wetland_info)

#date/year standardization

amphib_dissect <- amphib_dissect %>%
  mutate(
    collect_date = as_date(collect_date),
    year = year(collect_date)
  )

amphib_parasite <- amphib_parasite %>%
  mutate(
    date = as_date(date),
    year = year(date)
  )

bd <- bd %>%
  mutate(
    date = as_date(date),
    year = as.integer(year)
  )

survey_spp <- survey_spp %>%
  mutate(
    date = as_date(date),
    year = as.integer(year)
  )

wetland_info <- wetland_info %>%
  mutate(
    date = as_date(date),
    year = as.integer(year)
  )

netting_info <- netting_info %>%
  mutate(
    date = as_date(date),
    year = as.integer(year)
  )

#sanity check key IDs

amphib_dissect %>%
  summarize(
    n_rows = n(),
    n_dissect_codes = n_distinct(dissect_code),
    n_sites = n_distinct(site_code),
    n_species = n_distinct(spp_code),
    min_year = min(year, na.rm = TRUE),
    max_year = max(year, na.rm = TRUE)
  )

amphib_parasite %>%
  summarize(
    n_rows = n(),
    n_dissect_codes = n_distinct(dissect_code),
    n_parasites = n_distinct(parasite_name),
    min_year = min(year, na.rm = TRUE),
    max_year = max(year, na.rm = TRUE)
  )

site_info %>%
  count(longevity, sort = TRUE)

wetland_info %>%
  summarize(
    n_rows = n(),
    n_sites = n_distinct(site),
    min_year = min(year, na.rm = TRUE),
    max_year = max(year, na.rm = TRUE)
  )








