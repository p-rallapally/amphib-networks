
library(tidyverse)
library(igraph)
library(ggraph)
library(janitor)
library(lubridate)
library(vegan)
library(knitr)
library(dplyr)
library(igraph)
library(ggraph)
library(stringr)
library(lubridate)


#lists column names so I know what to join by

files <- list.files("Desktop/briggs/networks/kacie_csvs", pattern="\\.csv$", full.names=TRUE)

setNames(
  lapply(files, \(f) names(read.csv(f, nrows=0))),
  basename(files)
)

#stage 0: fixing the dates in some of the datasets bc the years were scuffed (likely bc excel)

#hopefully the collection code is more accurate 
amphib_dissect <- amphib_dissect %>%
  mutate(
    collect_date = str_extract(collect_code, "\\d{8}"),
    collect_date = ymd(collect_date),
    year = year(collect_date)
  )

amphib_parasite <- amphib_parasite %>%
  mutate(
    date = str_extract(collect_code, "\\d{8}"),
    date = ymd(date),
    year = year(date)
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


#stage 2: construct master host–parasite interaction table

#cleaning disssection data

dissect_clean <- amphib_dissect %>%
  select(
    dissect_code,
    site_code,
    spp_code,
    collect_date,
    year,
    parasites
  ) %>%
  filter(
    !is.na(dissect_code),
    !is.na(site_code),
    !is.na(spp_code),
    !is.na(year)
  )

#cleaning parasite
parasite_clean <- amphib_parasite %>%
  select(
    dissect_code,
    parasite_name,
    parasite_cnt,
    presence
  ) %>%
  mutate(
    parasite_name = str_trim(parasite_name)
  ) %>%
  filter(
    !is.na(dissect_code),
    !is.na(parasite_name),
    parasite_name != ""
  )

#joining parasite obs to host metadata
hp_raw <- parasite_clean %>%
  left_join(
    dissect_clean,
    by = "dissect_code"
  )

hp_raw <- hp_raw %>% #removing incomplete obs
  filter(
    !is.na(site_code),
    !is.na(spp_code),
    !is.na(year)
  )

#convert to binary presence
hp_binary <- hp_raw %>%
mutate(
  interaction = case_when(
    !is.na(parasite_cnt) & parasite_cnt > 0 ~ 1L,
    presence %in% c(TRUE, "Y", "Yes", 1) ~ 1L,
    TRUE ~ 0L
  )
)


#collapse to unique interactions (same host with the same site/year are collapsed into one)
hp_edges <- hp_binary %>%
  group_by(
    site_code,
    year,
    spp_code,
    parasite_name
  ) %>%
  summarize(
    interaction = max(interaction),
    n_hosts_examined = n_distinct(dissect_code),
    total_parasites = sum(parasite_cnt, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(interaction == 1)


#inc. pond longevity
hp_edges <- hp_edges %>%
  left_join(
    site_info %>%
      select(site_code, longevity),
    by = "site_code"
  )

#sanity checking
cat("\nUnique interactions:", nrow(hp_edges), "\n")

cat("\nHost species:", n_distinct(hp_edges$spp_code), "\n")

cat("Parasite taxa:", n_distinct(hp_edges$parasite_name), "\n")

cat("Sites:", n_distinct(hp_edges$site_code), "\n")

cat("Years:", n_distinct(hp_edges$year), "\n")#thank god they're all reasonable

# interactions per year
hp_edges %>%
  count(year)

# interactions by longevity class
hp_edges %>%
  count(longevity)

# host richness by site
hp_edges %>%
  distinct(site_code, spp_code) %>%
  count(site_code, name = "host_richness")

# parasite richness by site
hp_edges %>%
  distinct(site_code, parasite_name) %>%
  count(site_code, name = "parasite_richness")


#stage 3: build metaweb

#creating binary edge list
meta_edges <- hp_edges %>%
  distinct(
    spp_code,
    parasite_name
  ) %>%
  mutate(
    from = paste0("host_", spp_code),
    to   = paste0("parasite_", parasite_name),
    weight = 1
  )

#creating node table
host_nodes <- meta_edges %>%
  distinct(from) %>%
  transmute(
    name = from,
    type = "Host",
    label = str_remove(from, "^host_")
  )

parasite_nodes <- meta_edges %>%
  distinct(to) %>%
  transmute(
    name = to,
    type = "Parasite",
    label = str_remove(to, "^parasite_")
  )

nodes <- bind_rows(
  host_nodes,
  parasite_nodes
)

#building graph + diagnostics
g_meta <- graph_from_data_frame(
  d = meta_edges %>%
    select(from, to, weight),
  vertices = nodes,
  directed = FALSE
)

cat("Nodes:", vcount(g_meta), "\n")
cat("Edges:", ecount(g_meta), "\n\n")

table(V(g_meta)$type)

is_bipartite(g_meta)

components(g_meta)$no


#degree dist
node_degree <- tibble(
  node = V(g_meta)$name,
  type = V(g_meta)$type,
  degree = degree(g_meta)
)

host_degree <- node_degree %>%
  filter(type == "Host") %>%
  arrange(desc(degree))

parasite_degree <- node_degree %>%
  filter(type == "Parasite") %>%
  arrange(desc(degree))

host_degree

parasite_degree


#viz 
set.seed(6262026)

V(g_meta)$bipartite_type <- V(g_meta)$type == "Host"

table(V(g_meta)$type, V(g_meta)$bipartite_type)

ggraph(g_meta, layout = "bipartite", types = V(g_meta)$bipartite_type) +
  
  geom_edge_link(
    alpha = 0.4
  ) +
  
  geom_node_point(
    aes(color = type),
    size = 4
  ) +
  
  theme_void() +
  
  labs(
    title = "Binary Host–Parasite Metaweb"
  )

#centrality
meta_centrality <- tibble(
  node = V(g_meta)$name,
  type = V(g_meta)$type,
  degree = degree(g_meta),
  betweenness = betweenness(
    g_meta,
    normalized = TRUE
  ),
  eigen = eigen_centrality(g_meta)$vector
)

meta_centrality %>%
  arrange(desc(eigen))


#host richness
host_richness <- hp_edges %>%
  distinct(
    spp_code,
    parasite_name
  ) %>%
  count(
    spp_code,
    name = "parasite_richness"
  ) 

host_richness %>%
  arrange(desc(parasite_richness))


#parasite host breadth
parasite_breadth <- hp_edges %>%
  distinct(
    spp_code,
    parasite_name
  ) %>%
  count(
    parasite_name,
    name = "host_breadth"
  ) 

parasite_breadth %>%
  arrange(desc(host_breadth))







