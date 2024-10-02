#################### Calculating CEA metrics (dalys, costs, ICERS ####################
###############################
## Calculate ICERs for the wealth quintiles ----
###############################

cost_wq = int_costs_dif_horizon %>% 
          apply(c(1:4, 7:8), sum) %>% 
          sweep(c(1:4, 6), trt_costs_averted_horizon, "+") 

daly_wq = sweep((cost_wq+1)/(cost_wq+1), c(1:4, 6), daly_averted_horizon, "*")

dimnames(daly_wq)[[2]] = gsub(" ", "", dimnames(daly_wq)[[2]])
dimnames(cost_wq)[[2]] = gsub(" ", "", dimnames(cost_wq)[[2]])

icers_wq =  (cost_wq/daly_wq) %>% 
  cubelyr::as.tbl_cube(met_name="icers") %>% as_tibble 

tmp_daly = daly %>% apply(2:6, cumsum) %>% 
  cubelyr::as.tbl_cube(met_name="dalys") %>% as_tibble %>% 
  filter(year %in% as.numeric(years_lbl[c(10,20)])) %>% 
  mutate(horizon=ifelse(year==as.numeric(years_lbl[10]),"10y","20y")) %>% 
  pivot_wider(id_cols = c("horizon", "comp", "ns_strat", "wealth_quintile", "discounting"),
              names_glue = "{strat}_DALYs",
              names_from = "strat", values_from = "dalys") %>% 
  mutate(DALYs_averted = `Status Quo_DALYs` - `MAPS add_DALYs`)

tmp_cases = (typhoid*epipar["reporting"]) %>% apply(c(1,3:6), sum) %>% 
  apply(2:5, cumsum) %>% 
  cubelyr::as.tbl_cube(met_name="cases") %>% as_tibble %>% 
  filter(year %in% as.numeric(years_lbl[c(10,20)])) %>% 
  mutate(horizon=ifelse(year==as.numeric(years_lbl[10]),"10y","20y")) %>% 
  pivot_wider(id_cols = c("horizon", "comp", "ns_strat", "wealth_quintile"),
              names_glue = "{strat}_Cases",
              names_from = "strat", values_from = "cases") %>% 
  mutate(Cases_averted = `Status Quo_Cases` - `MAPS add_Cases`)
  
tmp_deaths = deaths_agegrp %>% apply(c(1,3:6), sum) %>% 
  apply(2:5, cumsum) %>% 
  cubelyr::as.tbl_cube(met_name="deaths") %>% as_tibble %>% 
  filter(year %in% as.numeric(years_lbl[c(10,20)])) %>% 
  mutate(horizon=ifelse(year==as.numeric(years_lbl[10]),"10y","20y")) %>% 
  pivot_wider(id_cols = c("horizon", "comp", "ns_strat", "wealth_quintile"),
              names_glue = "{strat}_Deaths",
              names_from = "strat", values_from = "deaths") %>% 
  mutate(Deaths_averted = `Status Quo_Deaths` - `MAPS add_Deaths`)

# tmp_costdif = cost_wq %>% 
#   cubelyr::as.tbl_cube(met_name="cost_dif") %>% as_tibble 

tmp_cost = int_costs %>% apply(c(1:6, 9:10), sum) %>% 
    sweep(c(1:6, 8), trt_costs, "+") %>% apply(c(1,3:8), sum) %>% 
    apply(2:7, cumsum) %>% 
    cubelyr::as.tbl_cube(met_name="cost") %>% as_tibble %>% 
    filter(year %in% as.numeric(years_lbl[c(10,20)])) %>% 
    mutate(horizon=ifelse(year==as.numeric(years_lbl[10]),"10y","20y")) %>% 
    pivot_wider(id_cols = c("horizon", "comp", "ns_strat", "wealth_quintile", 
                          "maps_profile", "discounting"),
              names_glue = "{strat}_Cost",
              names_from = "strat", values_from = "cost") %>% 
  mutate(Cost_dif = `MAPS add_Cost` - `Status Quo_Cost`)

tmp_int_cost = int_costs %>% apply(c(1, 3:6, 9:10), sum) %>% 
  apply(2:7, cumsum) %>%
  cubelyr::as.tbl_cube(met_name="int_cost") %>% as_tibble %>% 
  filter(year %in% as.numeric(years_lbl[c(10,20)])) %>% 
  mutate(horizon=ifelse(year==as.numeric(years_lbl[10]),"10y","20y")) %>% 
  pivot_wider(id_cols = c("horizon", "comp", "ns_strat", "wealth_quintile", 
                          "maps_profile", "discounting"),
              names_glue = "{strat}_IntCost",
              names_from = "strat", values_from = "int_cost")

tmp_trt_cost = trt_costs %>% apply(c(1,3:7), sum) %>% 
  apply(2:6, cumsum) %>%
  cubelyr::as.tbl_cube(met_name="trt_cost") %>% as_tibble %>% 
  filter(year %in% as.numeric(years_lbl[c(10,20)])) %>% 
  mutate(horizon=ifelse(year==as.numeric(years_lbl[10]),"10y","20y")) %>% 
  pivot_wider(id_cols = c("horizon", "comp", "ns_strat", "wealth_quintile", 
                          "discounting"),
              names_glue = "{strat}_TrtCost",
              names_from = "strat", values_from = "trt_cost")

icers_dalys_costs_summary = tmp_cases %>% full_join(tmp_deaths) %>% full_join(tmp_daly) %>%
  full_join(tmp_cost) %>% full_join(tmp_int_cost) %>% full_join(tmp_trt_cost) %>%  
  full_join(icers_wq)

###############################
## Calculate ICERS for the overall country ----
###############################

cost_all = int_costs_dif_horizon %>% 
  apply(c(1:3, 7:8), sum) %>% 
  sweep(c(1:3,5), apply(trt_costs_averted_horizon, c(1:3, 5), sum), "+") 

daly_all = sweep((cost_all+1)/(cost_all+1), c(1:3, 5), 
                 apply(daly_averted_horizon, c(1:3, 5), sum), "*")

dimnames(daly_all)[[2]] = gsub(" ", "", dimnames(daly_all)[[2]])
dimnames(cost_all)[[2]] = gsub(" ", "", dimnames(cost_all)[[2]])

icers_all = (cost_all/daly_all) %>% 
  cubelyr::as.tbl_cube(met_name="icers") %>% as_tibble 

tmp_daly_all = apply(daly, c(1:4,6), sum) %>% apply(2:5, cumsum) %>% 
  cubelyr::as.tbl_cube(met_name="dalys") %>% as_tibble %>% 
  filter(year %in% as.numeric(years_lbl[c(10,20)])) %>% 
  mutate(horizon=ifelse(year==as.numeric(years_lbl[10]),"10y","20y")) %>% 
  pivot_wider(id_cols = c("horizon", "comp", "ns_strat", "discounting"),
              names_glue = "{strat}_DALYs",
              names_from = "strat", values_from = "dalys") %>% 
  mutate(DALYs_averted = `Status Quo_DALYs` - `MAPS add_DALYs`)

tmp_cases_all = (typhoid*epipar["reporting"]) %>% apply(c(1,3:5), sum) %>% 
  apply(2:4, cumsum) %>% 
  cubelyr::as.tbl_cube(met_name="cases") %>% as_tibble %>% 
  filter(year %in% as.numeric(years_lbl[c(10,20)])) %>% 
  mutate(horizon=ifelse(year==as.numeric(years_lbl[10]),"10y","20y")) %>% 
  pivot_wider(id_cols = c("horizon", "comp", "ns_strat"),
              names_glue = "{strat}_Cases",
              names_from = "strat", values_from = "cases") %>% 
  mutate(Cases_averted = `Status Quo_Cases` - `MAPS add_Cases`)

tmp_deaths_all = deaths_agegrp %>% apply(c(1,3:5), sum) %>% 
  apply(2:4, cumsum) %>% 
  cubelyr::as.tbl_cube(met_name="deaths") %>% as_tibble %>% 
  filter(year %in% as.numeric(years_lbl[c(10,20)])) %>% 
  mutate(horizon=ifelse(year==as.numeric(years_lbl[10]),"10y","20y")) %>% 
  pivot_wider(id_cols = c("horizon", "comp", "ns_strat"),
              names_glue = "{strat}_Deaths",
              names_from = "strat", values_from = "deaths") %>% 
  mutate(Deaths_averted = `Status Quo_Deaths` - `MAPS add_Deaths`)

tmp_cost_all = int_costs %>% apply(c(1:6, 9:10), sum) %>% 
  sweep(c(1:6, 8), trt_costs, "+") %>% apply(c(1, 3:5, 7:8), sum) %>% 
  apply(2:6, cumsum) %>% 
  cubelyr::as.tbl_cube(met_name="cost") %>% as_tibble %>% 
  filter(year %in% as.numeric(years_lbl[c(10,20)])) %>% 
  mutate(horizon=ifelse(year==as.numeric(years_lbl[10]),"10y","20y")) %>% 
  pivot_wider(id_cols = c("horizon", "comp", "ns_strat",  
                          "maps_profile", "discounting"),
              names_glue = "{strat}_Cost",
              names_from = "strat", values_from = "cost") %>% 
  mutate(Cost_dif = `MAPS add_Cost` - `Status Quo_Cost`)

tmp_int_cost_all = int_costs %>% apply(c(1, 3:5, 9:10), sum) %>% 
  apply(2:6, cumsum) %>% 
  cubelyr::as.tbl_cube(met_name="int_cost") %>% as_tibble %>% 
  filter(year %in% as.numeric(years_lbl[c(10,20)])) %>% 
  mutate(horizon=ifelse(year==as.numeric(years_lbl[10]),"10y","20y")) %>% 
  pivot_wider(id_cols = c("horizon", "comp", "ns_strat", 
                          "maps_profile", "discounting"),
              names_glue = "{strat}_IntCost",
              names_from = "strat", values_from = "int_cost")

tmp_trt_cost_all = trt_costs %>% apply(c(1, 3:5, 7), sum) %>% 
  apply(2:5, cumsum) %>% 
  cubelyr::as.tbl_cube(met_name="trt_cost") %>% as_tibble %>% 
  filter(year %in% as.numeric(years_lbl[c(10,20)])) %>% 
  mutate(horizon=ifelse(year==as.numeric(years_lbl[10]),"10y","20y")) %>% 
  pivot_wider(id_cols = c("horizon", "comp", "ns_strat", "discounting"),
              names_glue = "{strat}_TrtCost",
              names_from = "strat", values_from = "trt_cost")

icers_dalys_costs_all_summary = tmp_cases_all %>% full_join(tmp_deaths_all) %>% full_join(tmp_daly_all) %>%
  full_join(tmp_cost_all) %>% full_join(tmp_int_cost_all) %>% full_join(tmp_trt_cost_all) %>% 
  full_join(icers_all)
