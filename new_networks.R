
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


#stage 4: network-level metrics for subsets

#function to build graph from edge table
build_hp_graph <- function(edge_df) {
  
  edges <- edge_df %>%
    distinct(spp_code, parasite_name) %>%
    mutate(
      from = paste0("host_", spp_code),
      to   = paste0("parasite_", parasite_name),
      weight = 1
    ) %>%
    select(from, to, weight)
  
  nodes <- tibble(
    name = unique(c(edges$from, edges$to))
  ) %>%
    mutate(
      node_type = case_when(
        str_detect(name, "^host_") ~ "Host",
        str_detect(name, "^parasite_") ~ "Parasite",
        TRUE ~ NA_character_
      ),
      bipartite_type = node_type == "Host"
    )
  
  graph_from_data_frame(
    d = edges,
    vertices = nodes,
    directed = FALSE
  )
}

#function that calculates network-level metrics + summarizes graph
summarize_hp_graph <- function(g, subset_type, subset_value) {
  
  node_tbl <- tibble(
    node = V(g)$name,
    node_type = V(g)$node_type
  )
  
  n_hosts <- node_tbl %>%
    filter(node_type == "Host") %>%
    nrow()
  
  n_parasites <- node_tbl %>%
    filter(node_type == "Parasite") %>%
    nrow()
  
  n_edges <- ecount(g)
  
  possible_edges <- n_hosts * n_parasites
  
  connectance <- ifelse(
    possible_edges > 0,
    n_edges / possible_edges,
    NA_real_
  )
  
  comps <- components(g)
  
  giant_component_size <- max(comps$csize)
  
  giant_component_fraction <- giant_component_size / vcount(g)
  
  modularity_value <- ifelse(
    ecount(g) > 0 && vcount(g) > 2,
    modularity(cluster_louvain(g)),
    NA_real_
  )
  
  tibble(
    subset_type = subset_type,
    subset_value = subset_value,
    n_hosts = n_hosts,
    n_parasites = n_parasites,
    n_nodes = vcount(g),
    n_edges = n_edges,
    connectance = connectance,
    mean_degree = mean(degree(g)),
    modularity = modularity_value,
    n_components = comps$no,
    giant_component_size = giant_component_size,
    giant_component_fraction = giant_component_fraction
  )
}

#pooled metaweb summary
g_meta <- build_hp_graph(hp_edges)

summary_meta <- summarize_hp_graph(
  g = g_meta,
  subset_type = "pooled",
  subset_value = "all_years"
)

summary_meta


#year-specific summaries 
yearly_network_summary <- hp_edges %>%
  group_split(year) %>%
  map_dfr(function(df_year) {
    
    yr <- unique(df_year$year)
    
    g_year <- build_hp_graph(df_year)
    
    summarize_hp_graph(
      g = g_year,
      subset_type = "year",
      subset_value = as.character(yr)
    )
  })

yearly_network_summary


#focusing on longevity
longevity_network_summary <- hp_edges %>%
  filter(!is.na(longevity), longevity != "") %>%
  group_split(longevity) %>%
  map_dfr(function(df_longevity) {
    
    lon <- unique(df_longevity$longevity)
    
    g_lon <- build_hp_graph(df_longevity)
    
    summarize_hp_graph(
      g = g_lon,
      subset_type = "longevity",
      subset_value = as.character(lon)
    )
  })

longevity_network_summary

#combining summaries
network_summary <- bind_rows(
  summary_meta,
  yearly_network_summary,
  longevity_network_summary
)

network_summary

#plots

library(patchwork)

# A: network size through time
p1 <- yearly_network_summary %>%
  mutate(year = as.integer(subset_value)) %>%
  ggplot(aes(x = year, y = n_edges)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Year",
    y = "Number of interactions",
    title = "A. Network size through time"
  ) +
  theme_minimal()

# B: connectance through time
p2 <- yearly_network_summary %>%
  mutate(year = as.integer(subset_value)) %>%
  ggplot(aes(x = year, y = connectance)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Year",
    y = "Connectance",
    title = "B. Connectance through time"
  ) +
  theme_minimal()

# C: modularity through time
p3 <- yearly_network_summary %>%
  mutate(year = as.integer(subset_value)) %>%
  ggplot(aes(x = year, y = modularity)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Year",
    y = "Louvain modularity",
    title = "C. Modularity through time"
  ) +
  theme_minimal()

# D: comparing by longevity
p4 <- longevity_network_summary %>%
  pivot_longer(
    cols = c(n_hosts, n_parasites, n_edges, connectance, modularity),
    names_to = "metric",
    values_to = "value"
  ) %>%
  ggplot(aes(x = subset_value, y = value)) +
  geom_col() +
  facet_wrap(~ metric, scales = "free_y") +
  labs(
    x = "Pond longevity",
    y = NULL,
    title = "D. Network structure by pond longevity"
  ) +
  theme_minimal()

# 2x2 grid
network_grid <- (p1 + p2) / (p3 + p4)

network_grid

#C: mod through time
yearly_network_summary %>%
  mutate(year = as.integer(subset_value)) %>%
  ggplot(aes(x = year, y = modularity)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Year",
    y = "Louvain modularity",
    title = "Host-parasite modularity through time"
  ) +
  theme_minimal()

#D: comparing by longevity
longevity_network_summary %>%
  pivot_longer(
    cols = c(n_hosts, n_parasites, n_edges, connectance, modularity),
    names_to = "metric",
    values_to = "value"
  ) %>%
  ggplot(aes(x = subset_value, y = value)) +
  geom_col() +
  facet_wrap(~ metric, scales = "free_y") +
  labs(
    x = "Pond longevity",
    y = NULL,
    title = "Network structure by pond longevity"
  ) +
  theme_minimal()

network_grid <- (p1 + p2) / (p3 + p4)

network_grid

#ggsave(
#  filename = "network_summary_2x2.png",
#  plot = network_grid,
#  width = 12,
#  height = 9,
#  dpi = 300
#)


#stage 5: node level metrics

#function for node centality per graph
summarize_hp_nodes <- function(g, subset_type, subset_value) {
  
  node_name <- V(g)$name
  
  label <- str_remove(node_name, "^host_|^parasite_")
  
  node_type <- V(g)$node_type
  
  tibble(
    subset_type = subset_type,
    subset_value = subset_value,
    node = node_name,
    node_type = node_type,
    label = label,
    spp_code = if_else(
      node_type == "Host",
      label,
      NA_character_
    ),
    
    parasite_name = if_else(
      node_type == "Parasite",
      label,
      NA_character_
    ),
    degree = degree(g),
    betweenness = betweenness(g, normalized = TRUE),
    closeness = closeness(g, normalized = TRUE),
    eigen = eigen_centrality(g)$vector
  )
}

#pooled node centrality
g_meta <- build_hp_graph(hp_edges)

meta_node_summary <- summarize_hp_nodes(
  g_meta,
  "pooled",
  "all_years"
)

meta_node_summary %>%
  arrange(desc(eigen))

#year-specific node centrality
yearly_node_summary <- hp_edges %>%
  group_split(year) %>%
  map_dfr(function(df_year) {
    
    yr <- unique(df_year$year)
    
    g_year <- build_hp_graph(df_year)
    
    summarize_hp_nodes(
      g = g_year,
      subset_type = "year",
      subset_value = as.character(yr)
    )
  })

yearly_node_summary

#longevity-specific node centrality
longevity_node_summary <- hp_edges %>%
  filter(!is.na(longevity), longevity != "") %>%
  group_split(longevity) %>%
  map_dfr(function(df_longevity) {
    
    lon <- unique(df_longevity$longevity)
    
    g_lon <- build_hp_graph(df_longevity)
    
    summarize_hp_nodes(
      g = g_lon,
      subset_type = "longevity",
      subset_value = as.character(lon)
    )
  })

longevity_node_summary

#combined node summaries
node_summary <- bind_rows(
  meta_node_summary,
  yearly_node_summary,
  longevity_node_summary
)

node_summary


#most central hosts + parasites
top_hosts <- meta_node_summary %>%
  filter(node_type == "Host") %>%
  arrange(desc(degree))

top_hosts

top_parasites <- meta_node_summary %>%
  filter(node_type == "Parasite") %>%
  arrange(desc(degree))

top_parasites

#plots

#top hosts by parasite richness/degree
top_hosts %>%
  ggplot(aes(x = reorder(label, degree), y = degree)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Host species",
    y = "Degree",
    title = "Host parasite richness in pooled metaweb"
  ) +
  theme_minimal()


#top parasites by centrality/degree
top_parasites %>%
  ggplot(aes(x = reorder(label, degree), y = degree)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Parasite taxon",
    y = "Degree",
    title = "Parasite host breadth in pooled metaweb"
  ) +
  theme_minimal()


#host centrality thru time 
yearly_node_summary %>%
  filter(node_type == "Host") %>%
  mutate(year = as.integer(subset_value)) %>%
  ggplot(aes(x = year,
             y = degree,
             group = label,
             color = label)) +
  geom_line(alpha = 0.8, linewidth = 1) +
  geom_point(size = 2) +
  labs(
    x = "Year",
    y = "Degree",
    color = "Host species",
    title = "Host parasite richness through time"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right"
  )


#same thing for parasites 

top_n_parasites <- 8

top_parasite_labels <- yearly_node_summary %>%
  filter(node_type == "Parasite") %>%
  group_by(label) %>%
  summarize(
    mean_degree = mean(degree, na.rm = TRUE),
    max_degree = max(degree, na.rm = TRUE),
    n_years_observed = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_degree), desc(max_degree)) %>%
  slice_head(n = top_n_parasites) %>%
  pull(label)

parasite_plot_df <- yearly_node_summary %>%
  filter(node_type == "Parasite") %>%
  mutate(
    year = as.integer(subset_value),
    highlight = if_else(label %in% top_parasite_labels, label, "Other parasites"),
    is_highlight = label %in% top_parasite_labels
  )

parasite_plot <- ggplot() +
  
  # Background lines: all non-top parasites
  geom_line(
    data = parasite_plot_df %>% filter(!is_highlight),
    aes(
      x = year,
      y = degree,
      group = label
    ),
    color = "gray75",
    alpha = 0.5,
    linewidth = 0.5
  ) +
  
  geom_point(
    data = parasite_plot_df %>% filter(!is_highlight),
    aes(
      x = year,
      y = degree,
      group = label
    ),
    color = "gray75",
    alpha = 0.5,
    size = 1
  ) +
  
  # Highlighted lines: top parasites
  geom_line(
    data = parasite_plot_df %>% filter(is_highlight),
    aes(
      x = year,
      y = degree,
      group = label,
      color = label
    ),
    linewidth = 1.1
  ) +
  
  geom_point(
    data = parasite_plot_df %>% filter(is_highlight),
    aes(
      x = year,
      y = degree,
      color = label
    ),
    size = 2
  ) +
  
  labs(
    x = "Year",
    y = "Degree",
    color = "Top parasites",
    title = "Parasite host breadth through time",
    subtitle = paste0(
      "Top ", top_n_parasites,
      " parasites by average yearly degree highlighted; all others shown in gray"
    )
  ) +
  
  theme_minimal() +
  
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8),
    plot.title = element_text(face = "bold")
  )

parasite_plot

#stage 5.5: subnetwork builder

#filtering hp_edges
filter_hp_edges <- function(data,
                            years = NULL,
                            longevity_classes = NULL,
                            sites = NULL,
                            hosts = NULL,
                            parasites = NULL) {
  
  out <- data
  
  if (!is.null(years)) {
    out <- out %>% filter(year %in% years)
  }
  
  if (!is.null(longevity_classes)) {
    out <- out %>% filter(longevity %in% longevity_classes)
  }
  
  if (!is.null(sites)) {
    out <- out %>% filter(site_code %in% sites)
  }
  
  if (!is.null(hosts)) {
    out <- out %>% filter(spp_code %in% hosts)
  }
  
  if (!is.null(parasites)) {
    out <- out %>% filter(parasite_name %in% parasites)
  }
  
  out
}


#build edge list + graph from any filtered hp_edges table
make_hp_network <- function(data,
                            subset_type = "custom",
                            subset_value = "custom") {
  
  edge_df <- data %>%
    distinct(spp_code, parasite_name) %>%
    filter(
      !is.na(spp_code),
      !is.na(parasite_name),
      spp_code != "",
      parasite_name != ""
    ) %>%
    mutate(
      from = paste0("host_", spp_code),
      to = paste0("parasite_", parasite_name),
      weight = 1
    ) %>%
    select(from, to, weight)
  
  node_df <- tibble(
    name = unique(c(edge_df$from, edge_df$to))
  ) %>%
    mutate(
      node_type = case_when(
        str_detect(name, "^host_") ~ "Host",
        str_detect(name, "^parasite_") ~ "Parasite",
        TRUE ~ NA_character_
      ),
      label = str_remove(name, "^host_|^parasite_"),
      bipartite_type = node_type == "Host"
    )
  
  g <- graph_from_data_frame(
    d = edge_df,
    vertices = node_df,
    directed = FALSE
  )
  
  list(
    graph = g,
    edges = edge_df,
    nodes = node_df,
    subset_type = subset_type,
    subset_value = subset_value
  )
}

#wrapped: filter + build network
build_subset_network <- function(data,
                                 years = NULL,
                                 longevity_classes = NULL,
                                 sites = NULL,
                                 hosts = NULL,
                                 parasites = NULL,
                                 subset_type = "custom",
                                 subset_value = "custom") {
  
  filtered_data <- filter_hp_edges(
    data = data,
    years = years,
    longevity_classes = longevity_classes,
    sites = sites,
    hosts = hosts,
    parasites = parasites
  )
  
  network <- make_hp_network(
    data = filtered_data,
    subset_type = subset_type,
    subset_value = subset_value
  )
  
  network$raw_data <- filtered_data
  
  network
}

#test
test_net <- build_subset_network(
  data = hp_edges,
  years = 2017,
  subset_type = "year",
  subset_value = "2017"
)

test_net$graph

summarize_hp_graph(
  g = test_net$graph,
  subset_type = test_net$subset_type,
  subset_value = test_net$subset_value
)


#stage 6: bullfrog presence - absence

bullfrog <- "RACA"
bullfrog_site_years <- hp_edges %>%
  distinct(site_code, year, spp_code) %>%
  group_by(site_code, year) %>%
  summarize(
    bullfrog_present = any(spp_code == bullfrog),
    .groups = "drop"
  )

bullfrog_site_years %>%
  count(bullfrog_present)

#adding bullfrog presence status to all edges 
hp_edges_bf <- hp_edges %>%
  left_join(
    bullfrog_site_years,
    by = c("site_code", "year")
  ) %>%
  mutate(
    bullfrog_present = replace_na(bullfrog_present, FALSE),
    bullfrog_status = if_else(
      bullfrog_present,
      "Bullfrog present",
      "Bullfrog absent"
    )
  )

hp_edges_bf %>%
  count(bullfrog_status)

#building present/absent networks 
bf_present_net <- build_subset_network(
  data = hp_edges_bf %>% filter(bullfrog_present),
  subset_type = "bullfrog_status",
  subset_value = "present"
)

bf_absent_net <- build_subset_network(
  data = hp_edges_bf %>% filter(!bullfrog_present),
  subset_type = "bullfrog_status",
  subset_value = "absent"
)


#network level comparison
bullfrog_network_summary <- bind_rows(
  summarize_hp_graph(
    g = bf_present_net$graph,
    subset_type = "bullfrog_status",
    subset_value = "present"
  ),
  summarize_hp_graph(
    g = bf_absent_net$graph,
    subset_type = "bullfrog_status",
    subset_value = "absent"
  )
)

bullfrog_network_summary

#node level comparison 

bullfrog_node_summary <- bind_rows(
  summarize_hp_nodes(
    g = bf_present_net$graph,
    subset_type = "bullfrog_status",
    subset_value = "present"
  ),
  summarize_hp_nodes(
    g = bf_absent_net$graph,
    subset_type = "bullfrog_status",
    subset_value = "absent"
  )
)

bullfrog_node_summary

#host centrality by bullfrog status
bullfrog_node_summary %>%
  filter(node_type == "Host") %>%
  select(subset_value, label, degree, betweenness, eigen) %>%
  arrange(subset_value, desc(degree))

#parasite comparison 
bullfrog_node_summary %>%
  filter(node_type == "Parasite") %>%
  select(subset_value, label, degree, betweenness, eigen) %>%
  arrange(subset_value, desc(degree))

#bullfrog network comparison (plot)
bullfrog_network_summary %>%
  pivot_longer(
    cols = c(n_hosts, n_parasites, n_edges, connectance, modularity),
    names_to = "metric",
    values_to = "value"
  ) %>%
  ggplot(aes(x = subset_value, y = value)) +
  geom_col() +
  facet_wrap(~ metric, scales = "free_y") +
  labs(
    x = "Bullfrog status",
    y = NULL,
    title = "Host-parasite network structure by bullfrog status"
  ) +
  theme_minimal()

#parasite degree by bullfrog status
bullfrog_node_summary %>%
  filter(node_type == "Parasite") %>%
  ggplot(aes(x = reorder(label, degree), y = degree)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ subset_value) +
  labs(
    x = "Parasite taxon",
    y = "Degree",
    title = "Parasite host breadth by bullfrog status"
  ) +
  theme_minimal()



#stage 8:  incorporating Bd data

bd_clean <- bd %>%
  transmute(
    site_code = site,
    year,
    spp_code,
    bd_positive = bd_infect,
    bd_load = ave_ze_new
  ) %>%
  mutate(
    
    # convert to logical
    bd_positive = case_when(
      bd_positive %in% c("Y", "YES", "Positive", 1, TRUE) ~ TRUE,
      TRUE ~ FALSE
    ),
    
    bd_load = as.numeric(bd_load)
    
  ) %>%
  filter(
    !is.na(site_code),
    !is.na(year),
    !is.na(spp_code)
  )

#host x site x year summaries
bd_host_site_year <- bd_clean %>%
  group_by(spp_code, site_code, year) %>%
  summarize(
    n_swabbed = n(),
    prevalence = mean(bd_positive, na.rm = TRUE),
    
    mean_load = if (all(is.na(bd_load))) {
      NA_real_
    } else {
      mean(bd_load, na.rm = TRUE)
    },
    
    median_load = if (all(is.na(bd_load))) {
      NA_real_
    } else {
      median(bd_load, na.rm = TRUE)
    },
    
    max_load = if (all(is.na(bd_load))) {
      NA_real_
    } else {
      max(bd_load, na.rm = TRUE)
    },
    
    .groups = "drop"
  )

#aggregate to species level
bd_host_summary <- bd_host_site_year %>%
  group_by(spp_code) %>%
  summarize(
    total_swabbed = sum(n_swabbed),
    prevalence = weighted.mean(prevalence, w = n_swabbed, na.rm = TRUE),
    mean_load = weighted.mean(mean_load, w = n_swabbed, na.rm = TRUE),
    median_load = median(median_load, na.rm = TRUE),
    max_load = max(max_load, na.rm = TRUE),
    .groups = "drop"
  )

#bd summaries only on pooled host nodes
host_node_bd <- meta_node_summary %>%
  filter(node_type == "Host") %>%
  left_join(bd_host_summary, by = "spp_code")

host_node_bd %>%
  select(
    label,
    degree,
    betweenness,
    eigen,
    total_swabbed,
    prevalence,
    mean_load,
    median_load,
    max_load
  )


#rank by centrality
host_node_bd %>%
  arrange(desc(degree))

host_node_bd %>%
  arrange(desc(eigen))

host_node_bd %>%
  arrange(desc(betweenness))

#scatterplots 

#deg vs prevalence
p5 <- host_node_bd %>%
  ggplot(
    aes(
      degree,
      prevalence,
      label = spp_code
    )
  ) +
  geom_point(size = 3) +
  geom_text(
    nudge_y = 0.02,
    check_overlap = TRUE
  ) +
  labs(
    x = "Host degree",
    y = "Bd prevalence",
    title = "Host centrality and Bd prevalence"
  ) +
  theme_minimal()

#eig vs prevalence
p6 <- host_node_bd %>%
  ggplot(
    aes(
      eigen,
      prevalence,
      label = spp_code
    )
  ) +
  geom_point(size = 3) +
  geom_text(
    nudge_y = 0.02,
    check_overlap = TRUE
  ) +
  labs(
    x = "Eigenvector centrality",
    y = "Bd prevalence",
    title = "Eigenvector centrality and Bd prevalence"
  ) +
  theme_minimal()


#deg vs mean load
p7 <- host_node_bd %>%
  ggplot(
    aes(
      degree,
      mean_load,
      label = spp_code
    )
  ) +
  geom_point(size = 3) +
  geom_text(
    nudge_y = 0.02,
    check_overlap = TRUE
  ) +
  scale_y_continuous(trans = "log1p") +
  labs(
    x = "Host degree",
    y = "Mean Bd load (log1p)",
    title = "Host centrality and Bd intensity"
  ) +
  theme_minimal()


#eig vs mean load
p8 <- host_node_bd %>%
  ggplot(
    aes(
      eigen,
      mean_load,
      label = spp_code
    )
  ) +
  geom_point(size = 3) +
  geom_text(
    nudge_y = 0.02,
    check_overlap = TRUE
  ) +
  scale_y_continuous(trans = "log1p") +
  labs(
    x = "Eigenvector centrality",
    y = "Mean Bd load (log1p)",
    title = "Eigenvector centrality and Bd intensity"
  ) +
  theme_minimal()

bd_overlay_grid <-
  (p5 + p6) /
  (p7 + p8)

bd_overlay_grid


#stage 9: stat time 

#building summaries + checking n
site_year_network_summary <- hp_edges_bf %>%
  group_by(site_code, year) %>%
  group_split() %>%
  map_dfr(function(df_sy) {
    
    site_i <- unique(df_sy$site_code)
    year_i <- unique(df_sy$year)
    longevity_i <- unique(df_sy$longevity)
    bullfrog_i <- unique(df_sy$bullfrog_present)
    
    g_sy <- build_hp_graph(df_sy)
    
    summarize_hp_graph(
      g = g_sy,
      subset_type = "site_year",
      subset_value = paste(site_i, year_i, sep = "_")
    ) %>%
      mutate(
        site_code = site_i,
        year = year_i,
        longevity = longevity_i[1],
        bullfrog_present = bullfrog_i[1]
      )
  }) %>%
  relocate(site_code, year, longevity, bullfrog_present)

site_year_network_summary


site_year_network_summary %>%
  summarize(
    n_site_years = n(),
    n_sites = n_distinct(site_code),
    min_year = min(year),
    max_year = max(year)
  )

site_year_network_summary %>%
  count(bullfrog_present)

site_year_network_summary %>%
  count(longevity)

#lm

lm_connectance <- lm(
  connectance ~ bullfrog_present + longevity + year,
  data = site_year_network_summary
)

summary(lm_connectance)

lm_modularity <- lm(
  modularity ~ bullfrog_present + longevity + year,
  data = site_year_network_summary
)

summary(lm_modularity)

lm_parasites <- lm(
  n_parasites ~ bullfrog_present + longevity + year,
  data = site_year_network_summary
)

summary(lm_parasites)


site_year_network_summary %>%
  ggplot(aes(x = bullfrog_present, y = connectance)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.6) +
  labs(
    x = "Bullfrog present",
    y = "Connectance",
    title = "Site-year network connectance by bullfrog presence"
  ) +
  theme_minimal()

site_year_network_summary %>%
  ggplot(aes(x = longevity, y = n_parasites)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.6) +
  labs(
    x = "Pond longevity",
    y = "Parasite richness",
    title = "Parasite richness by pond longevity"
  ) +
  theme_minimal()

#Bullfrogs are associated with more parasite taxa, lower connectance, higher modularity

# stage 10: Host projection using Jaccard similarity

#Function to build a host-host Jaccard projection

build_host_jaccard_projection <- function(edge_df) {
  
  # Binary host x parasite incidence matrix
  incidence <- edge_df %>%
    distinct(spp_code, parasite_name) %>%
    mutate(present = 1L) %>%
    pivot_wider(
      names_from = parasite_name,
      values_from = present,
      values_fill = 0
    )
  
  host_names <- incidence$spp_code
  
  host_matrix <- incidence %>%
    select(-spp_code) %>%
    as.matrix()
  
  rownames(host_matrix) <- host_names
  
  # Jaccard similarity:
  # shared parasites / total unique parasites across the pair
  jaccard_matrix <- matrix(
    NA_real_,
    nrow = nrow(host_matrix),
    ncol = nrow(host_matrix),
    dimnames = list(host_names, host_names)
  )
  
  for (i in seq_len(nrow(host_matrix))) {
    for (j in seq_len(nrow(host_matrix))) {
      
      host_i <- host_matrix[i, ] > 0
      host_j <- host_matrix[j, ] > 0
      
      intersection <- sum(host_i & host_j)
      union <- sum(host_i | host_j)
      
      jaccard_matrix[i, j] <- ifelse(
        union == 0,
        NA_real_,
        intersection / union
      )
    }
  }
  
  diag(jaccard_matrix) <- 0
  
  # Convert upper triangle to edge table
  projection_edges <- as.data.frame(as.table(jaccard_matrix)) %>%
    as_tibble() %>%
    rename(
      from = Var1,
      to = Var2,
      jaccard = Freq
    ) %>%
    mutate(
      from = as.character(from),
      to = as.character(to)
    ) %>%
    filter(
      from < to,
      !is.na(jaccard),
      jaccard > 0
    )
  
  projection_nodes <- tibble(
    name = host_names,
    parasite_richness = rowSums(host_matrix > 0)
  )
  
  g_projection <- graph_from_data_frame(
    d = projection_edges,
    vertices = projection_nodes,
    directed = FALSE
  )
  
  list(
    graph = g_projection,
    edges = projection_edges,
    nodes = projection_nodes,
    incidence_matrix = host_matrix,
    jaccard_matrix = jaccard_matrix
  )
}

#Pooled host projection
host_projection <- build_host_jaccard_projection(hp_edges)

g_host_jaccard <- host_projection$graph

host_projection$edges
host_projection$jaccard_matrix

set.seed(6262026)

ggraph(g_host_jaccard, layout = "fr") +
  geom_edge_link(
    aes(
      width = jaccard,
      alpha = jaccard
    ),
    show.legend = TRUE
  ) +
  geom_node_point(
    aes(size = parasite_richness)
  ) +
  geom_node_text(
    aes(label = name),
    repel = TRUE
  ) +
  scale_edge_width(
    range = c(0.5, 4),
    name = "Jaccard similarity"
  ) +
  scale_edge_alpha(
    range = c(0.3, 1),
    name = "Jaccard similarity"
  ) +
  labs(
    title = "Host projection based on shared parasite assemblages",
    subtitle = "Edge weights represent Jaccard similarity"
  ) +
  theme_void()


# 4. Host centrality in projected network

host_projection_centrality <- tibble(
  spp_code = V(g_host_jaccard)$name,
  
  projected_degree = degree(g_host_jaccard),
  
  strength = strength(
    g_host_jaccard,
    weights = E(g_host_jaccard)$jaccard
  ),
  
  betweenness = betweenness(
    g_host_jaccard,
    directed = FALSE,
    weights = 1 / E(g_host_jaccard)$jaccard,
    normalized = TRUE
  ),
  
  closeness = closeness(
    g_host_jaccard,
    weights = 1 / E(g_host_jaccard)$jaccard,
    normalized = TRUE
  ),
  
  eigen = eigen_centrality(
    g_host_jaccard,
    weights = E(g_host_jaccard)$jaccard
  )$vector,
  
  parasite_richness = V(g_host_jaccard)$parasite_richness
) %>%
  arrange(desc(strength))

host_projection_centrality

# 6. Year-specific projected networks

yearly_host_projection_summary <- hp_edges %>%
  group_split(year) %>%
  map_dfr(function(df_year) {
    
    yr <- unique(df_year$year)
    
    projection <- build_host_jaccard_projection(df_year)
    g <- projection$graph
    
    tibble(
      year = yr,
      n_hosts = vcount(g),
      n_host_links = ecount(g),
      
      mean_jaccard = ifelse(
        ecount(g) > 0,
        mean(E(g)$jaccard),
        NA_real_
      ),
      
      max_jaccard = ifelse(
        ecount(g) > 0,
        max(E(g)$jaccard),
        NA_real_
      ),
      
      mean_strength = ifelse(
        vcount(g) > 0,
        mean(strength(g, weights = E(g)$jaccard)),
        NA_real_
      ),
      
      density = edge_density(g)
    )
  })

yearly_host_projection_summary

# 7. Mean Jaccard similarity through time

yearly_host_projection_summary %>%
  ggplot(aes(x = year, y = mean_jaccard)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Year",
    y = "Mean host-pair Jaccard similarity",
    title = "Host parasite-sharing similarity through time"
  ) +
  theme_minimal()


# 8. Host-pair similarity trajectories
yearly_host_pair_similarity <- hp_edges %>%
  group_split(year) %>%
  map_dfr(function(df_year) {
    
    yr <- unique(df_year$year)
    
    projection <- build_host_jaccard_projection(df_year)
    
    projection$edges %>%
      mutate(
        year = yr,
        host_pair = paste(from, to, sep = " – ")
      )
  })

yearly_host_pair_similarity %>%
  ggplot(
    aes(
      x = year,
      y = jaccard,
      group = host_pair,
      color = host_pair
    )
  ) +
  geom_line() +
  geom_point() +
  labs(
    x = "Year",
    y = "Jaccard similarity",
    color = "Host pair",
    title = "Shared parasite assemblages among host pairs"
  ) +
  theme_minimal()



# stage 9.5: Rarefaction of Bd data
set.seed(6262026)

# Number of swabs per species
bd_clean %>%
  count(spp_code)

# Smallest sample size
network_hosts <- c("BUBO", "PSRE", "RACA", "TAGR", "TATO")

bd_network <- bd_clean %>%
  filter(spp_code %in% network_hosts)

bd_network %>%
  count(spp_code)

min_n <- bd_network %>%
  count(spp_code) %>%
  summarize(min_n = min(n)) %>%
  pull()

min_n

rarefy_species <- function(df_species, sample_size){
  
  sampled <- df_species %>%
    slice_sample(n = sample_size)
  
  tibble(
    
    prevalence =
      mean(sampled$bd_positive),
    
    mean_load =
      mean(sampled$bd_load, na.rm = TRUE),
    
    median_load =
      median(sampled$bd_load, na.rm = TRUE),
    
    max_load =
      max(sampled$bd_load, na.rm = TRUE)
    
  )
}

#repeat 1000
n_iter <- 1000

bd_rarefied <- map_dfr(
  1:n_iter,
  function(i){
    bd_network %>%
      group_by(spp_code) %>%
      group_modify(~rarefy_species(.x, min_n)) %>%
      ungroup() %>%
      mutate(iteration = i)
  }
  
)

bd_rarefied_summary <- bd_rarefied %>%
  
  group_by(spp_code) %>%
  
  summarize(
    
    prevalence_mean =
      mean(prevalence),
    
    prevalence_sd =
      sd(prevalence),
    
    prevalence_lower =
      quantile(prevalence, 0.025),
    
    prevalence_upper =
      quantile(prevalence, 0.975),
    
    mean_load_mean =
      mean(mean_load, na.rm = TRUE),
    
    mean_load_sd =
      sd(mean_load, na.rm = TRUE),
    
    mean_load_lower =
      quantile(mean_load, 0.025, na.rm = TRUE),
    
    mean_load_upper =
      quantile(mean_load, 0.975, na.rm = TRUE),
    
    median_load_mean =
      mean(median_load, na.rm = TRUE),
    
    max_load_mean =
      mean(max_load, na.rm = TRUE),
    .groups = "drop"
  )

bd_rarefied_summary

bd_comparison <- bd_host_summary %>%
  
  select(
    spp_code,
    total_swabbed,
    prevalence,
    mean_load
  ) %>%
  
  left_join(
    bd_rarefied_summary,
    by = "spp_code"
  )

bd_comparison

ggplot(
  bd_rarefied_summary,
  aes(
    x = spp_code,
    y = prevalence_mean
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = prevalence_lower,
      ymax = prevalence_upper
    ),
    
    width = 0.15
    
  ) +
  
  labs(
    x = "Host species",
    y = "Rarefied Bd prevalence",
    title = paste(
      "Bd prevalence after rarefaction (n =",
      min_n,
      "swabs/species)"
    )
  ) +

  theme_minimal()

ggplot(
  bd_rarefied_summary,
  aes(
    x = spp_code,
    y = mean_load_mean
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = mean_load_lower,
      ymax = mean_load_upper
    ),
    width = 0.15
  ) +
  
  labs(
    x = "Host species",
    y = "Rarefied mean Bd load",
    title = paste(
      "Bd load after rarefaction (n =",
      min_n,
      "swabs/species)"
    )
  ) +
  theme_minimal()




host_projection_bd <- host_projection_centrality %>%
  left_join(
    bd_rarefied_summary,
    by = "spp_code"
  )

host_projection_bd

host_projection_bd %>%
  ggplot(
    aes(
      x = strength,
      y = prevalence_mean,
      label = spp_code
    )
  ) +
  geom_point(size = 3) +
  geom_text(
    nudge_y = 0.02,
    check_overlap = TRUE
  ) +
  geom_errorbar(
    aes(
      ymin = prevalence_lower,
      ymax = prevalence_upper
    ),
    width = 0.03
  ) +
  labs(
    x = "Projected host strength",
    y = "Rarefied Bd prevalence",
    title = "Parasite-sharing connectivity and standardized Bd prevalence"
  ) +
  theme_minimal()

host_projection_bd %>%
  ggplot(
    aes(
      x = strength,
      y = mean_load_mean,
      label = spp_code
    )
  ) +
  geom_point(size = 3) +
  geom_text(
    nudge_y = 0.02,
    check_overlap = TRUE
  ) +
  geom_errorbar(
    aes(
      ymin = mean_load_lower,
      ymax = mean_load_upper
    ),
    width = 0.03
  ) +
  scale_y_continuous(trans = "log1p") +
  labs(
    x = "Projected host strength",
    y = "Rarefied mean Bd load",
    title = "Parasite-sharing connectivity and standardized Bd load"
  ) +
  theme_minimal()


