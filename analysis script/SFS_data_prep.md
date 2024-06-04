data_prep
================
2024-06-04

# Libraries

# Repository

``` r
# set directory
#repo_root <- "/Users/eabousam" 
#repo      <- sprintf('%s/SFS-coinfection-interactions/',     repo_root)
#setwd(repo)
```

# Read Data

``` r
#full encounter sample data
full_encounter_sample = readRDS('model_data.RDS')

full_encounter_sample = full_encounter_sample %>%
  filter(organism %in% c("Rhinovirus", "SARS_CoV_2", "RSV.A","RSV.B","IAV..H1N1", "IAV..H3N2", "Adenovirus", "IBV")) %>%
  filter(site %in% c("SCAN")) %>%
  na.omit("individual") %>%
  dplyr::select(individual, organism, present, sex, agebin, date,device, site, symptoms, sample, ct) 

# Counting Co-infections
co_infections <- full_encounter_sample %>%
  group_by(individual, date) %>%
  filter(sum(present) > 1) %>%
  ungroup()
```

# Figure 1: proportion circulating SFS

``` r
# Summarize data to get frequency of presence for each organism by date
summary_data <- full_encounter_sample %>%
  group_by(date, organism) %>%
  summarise(n = sum(present),
            total_count = n()) %>%
  mutate(percentage = n / sum(n),
         positive_percentage = n / total_count)
```

    ## `summarise()` has grouped output by 'date'. You can override using the
    ## `.groups` argument.

``` r
summary_data$organism <- gsub("_", "-", summary_data$organism)

# Replace specific labels with cleaner ones
summary_data$organism <- gsub("RSV.A", "RSV A", summary_data$organism)
summary_data$organism <- gsub("RSV.B", "RSV B", summary_data$organism)
summary_data$organism <- gsub("IAV..H1N1", "Influenza A H1N1", summary_data$organism)
summary_data$organism <- gsub("IAV..H3N2", "Influenza A H3N2", summary_data$organism)
summary_data$organism <- gsub("IBV", "Influenza B", summary_data$organism)
# Plot with modified legend labels and colors
prev_plot = ggplot(summary_data, aes(x = date, y = percentage, fill = factor(organism))) +
  geom_area(stat = "smooth", method = "loess", position = "fill", 
            alpha = 0.9, span = 1/4) +
  
  scale_fill_manual(name = "Organism",  # Adjust the legend name
                    values = c(
                      "Rhinovirus" = "#355070",
                      "SARS-CoV-2" = "#F94144", 
                      "RSV A" = "#F8961E",  
                      "RSV B" = "#F9C74F",  
                      "Influenza A H1N1" = "#90BE6D",
                      "Influenza A H3N2" = "#577590",
                      "Adenovirus" = "#6d597a",
                      "Influenza B" = "#4D908E"
                    )) +
  scale_color_manual(name = "Organism",  # Adjust the legend name for color
                    values = c(
                      "Rhinovirus" = "#355070",
                      "SARS-CoV-2" = "#F94144", 
                      "RSV A" = "#F8961E", 
                      "RSV B" = "#F9C74F",  
                      "Influenza A H1N1" = "#90BE6D",
                      "Influenza A H3N2" = "#577590",
                      "Adenovirus" = "#6d597a",
                      "Influenza B" = "#4D908E"
                    )) +
  sm_corr_theme() +
  labs(
    x = "Date",
    y = "Prevalence",
    #title = "Prevalence of Respiratory Pathogens Over Time from Seattle SCAN SFS data"
  )+
  scale_y_percent(scale = 100, limits = c(0,1))
```

    ## sm_corr_theme is equivalent to sm_hvgrid.

``` r
# Make another figure for absolute counts
# Create a plot to show the proportion of positive samples
positive_plot <- ggplot(summary_data, aes(x = date, y = positive_percentage, fill = factor(organism))) +
  geom_area(stat = "smooth", method = "loess",
            alpha = 0.9, span = 1/4) +
  
  scale_fill_manual(name = "Organism",  # Adjust the legend name
                    values = c(
                      "Rhinovirus" = "#355070",
                      "SARS-CoV-2" = "#F94144", 
                      "RSV A" = "#F8961E",  
                      "RSV B" = "#F9C74F",  
                      "Influenza A H1N1" = "#90BE6D",
                      "Influenza A H3N2" = "#577590",
                      "Adenovirus" = "#6d597a",
                      "Influenza B" = "#4D908E"
                    )) +
  sm_corr_theme() +
  labs(
    x = "Date",
    y = "Proportion of Positive Samples",
    #title = "Proportion of Positive Samples Over Time"
  ) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1))
```

    ## sm_corr_theme is equivalent to sm_hvgrid.

``` r
ggsave("prev_plot_sfs.png", prev_plot, width = 12, height = 6, dpi = 300,scale = 0.75, bg = 'white')
```

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: Removed 225 rows containing non-finite values (`stat_smooth()`).

``` r
ggsave("prop_pos_plot_sfs.png", positive_plot, width = 12, height = 6, dpi = 300,scale = 0.75, bg = 'white')
```

    ## `geom_smooth()` using formula = 'y ~ x'

# Data: prepare multiple encounter data

``` r
# Specify pathogens to be investigated
spec_pathogen_data <- full_encounter_sample %>%
  filter(organism %in% c("Rhinovirus", "SARS_CoV_2", "RSV.A", "RSV.B", "IAV..H1N1", "Adenovirus","IAV..H3N2", "IBV")) %>%
  mutate(organism = case_when(
    organism == "RSV.A" ~ "RSV",
    organism == "RSV.B" ~ "RSV",
    organism == "IAV..H1N1" ~ "IAV",
    organism == "IAV..H3N2" ~ "IAV",
    TRUE ~ organism
  )) %>%
  group_by(individual, organism, date) %>%
  dplyr::select(individual, organism, present, ct, sex, agebin, date,device, site, symptoms, sample) %>%
  arrange(organism, desc(present)) %>%
  slice(1) %>%
  ungroup() 

#Check co-infections
spec_pathogen_data %>%
  group_by(individual, date) %>%
  filter(sum(present) > 2) %>%
  ungroup()
```

    ## # A tibble: 54 × 11
    ##    individual        organism present    ct sex   agebin date       device site 
    ##    <chr>             <chr>    <lgl>   <dbl> <chr> <fct>  <date>     <chr>  <chr>
    ##  1 3b0e37482b6b65c8… Adenovi… TRUE     26.4 male  1-5    2022-04-18 OpenA… SCAN 
    ##  2 3b0e37482b6b65c8… IAV      FALSE    NA   male  1-5    2022-04-18 OpenA… SCAN 
    ##  3 3b0e37482b6b65c8… IBV      FALSE    NA   male  1-5    2022-04-18 OpenA… SCAN 
    ##  4 3b0e37482b6b65c8… RSV      FALSE    NA   male  1-5    2022-04-18 OpenA… SCAN 
    ##  5 3b0e37482b6b65c8… Rhinovi… TRUE     27.7 male  1-5    2022-04-18 OpenA… SCAN 
    ##  6 3b0e37482b6b65c8… SARS_Co… TRUE      8   male  1-5    2022-04-18 Taqma… SCAN 
    ##  7 466b4af98120b149… Adenovi… TRUE     21.8 male  1-5    2021-08-30 OpenA… SCAN 
    ##  8 466b4af98120b149… IAV      FALSE    NA   male  1-5    2021-08-30 OpenA… SCAN 
    ##  9 466b4af98120b149… IBV      FALSE    NA   male  1-5    2021-08-30 OpenA… SCAN 
    ## 10 466b4af98120b149… RSV      TRUE     12.1 male  1-5    2021-08-30 OpenA… SCAN 
    ## # ℹ 44 more rows
    ## # ℹ 2 more variables: symptoms <pq__text>, sample <chr>

``` r
# Pivot wider
encounter_rv_cv_pairs <- spec_pathogen_data %>%

  filter(as.numeric(difftime(lead(date), date, units = "days")) < 365) %>%
  #Pivot the organism column wider to separate organisms
  pivot_wider(id_cols = c(individual, date, sex, agebin, symptoms, sample),
              names_from = organism, values_from = present, names_prefix = "present_"
              #,values_fill = FALSE
              ) %>%
  group_by(individual) %>%
  filter(n() > 1) 


# Convert into binary 0 and 1
encounter_rv_cv_pairs_conf = encounter_rv_cv_pairs %>%
  mutate(present_rv_n = as.integer(as.logical(present_Rhinovirus)),
         present_sars_cov_2_n = as.integer(as.logical(present_SARS_CoV_2)),
         present_rsv_n = as.integer(as.logical(present_RSV)),
         present_adv_n = as.integer(as.logical(present_Adenovirus)),
         present_iav_n = as.integer(as.logical(present_IAV)),
         present_ibv_n = as.integer(as.logical(present_IBV))) %>%
  dplyr::select(individual, date, present_rv_n, present_sars_cov_2_n, present_rsv_n, present_adv_n, present_iav_n, present_ibv_n,sex, agebin, symptoms, sample)%>%
   na.omit(c(
    "present_rv_n", "present_sars_cov_2_n", "present_rsv_n", 
    "present_iav_n", "present_adv_n", "present_ibv_n"
  ))


# adding Ct
sample_ct_summary <- spec_pathogen_data %>%
  pivot_wider(
    id_cols = c(individual, date, symptoms, sample),
    names_from = organism,
    values_from =  ct,
    names_prefix = c("ct_"),
    values_fill = list(ct = NA)
  )


# Merge ct values data
encounter_rv_cv_pairs_conf <- left_join(encounter_rv_cv_pairs_conf, sample_ct_summary)
```

    ## Joining with `by = join_by(individual, date, symptoms, sample)`

# Evaluate Missingness

``` r
# Check missingness in data

propmiss <- function(dataframe) {
  m <- sapply(dataframe, function(x) {
      data.frame(
        nmiss=sum(is.na(x)),
        n=length(x),
        propmiss=sum(is.na(x))/length(x))
  })
  d <- data.frame(t(m))
  d <- sapply(d, unlist)
  d <- as.data.frame (d)
  d$variable <- row.names (d)
  row.names(d) <- NULL
  d <- cbind(d[ncol(d)], d[-ncol(d)])
  return(d[order(d$propmiss),])
}

# Look at the amount of missingness
propmiss(encounter_rv_cv_pairs_conf)
```

    ##                variable nmiss     n  propmiss
    ## 1            individual     0 22536 0.0000000
    ## 2                  date     0 22536 0.0000000
    ## 3          present_rv_n     0 22536 0.0000000
    ## 4  present_sars_cov_2_n     0 22536 0.0000000
    ## 5         present_rsv_n     0 22536 0.0000000
    ## 6         present_adv_n     0 22536 0.0000000
    ## 7         present_iav_n     0 22536 0.0000000
    ## 8         present_ibv_n     0 22536 0.0000000
    ## 9                   sex     0 22536 0.0000000
    ## 10               agebin     0 22536 0.0000000
    ## 11             symptoms     0 22536 0.0000000
    ## 12               sample     0 22536 0.0000000
    ## 17        ct_Rhinovirus 19871 22536 0.8817448
    ## 18        ct_SARS_CoV_2 21850 22536 0.9695598
    ## 13        ct_Adenovirus 22275 22536 0.9884185
    ## 16               ct_RSV 22444 22536 0.9959176
    ## 14               ct_IAV 22519 22536 0.9992457
    ## 15               ct_IBV 22529 22536 0.9996894

# MissRanger: Imputing missing data if any

``` r
# # Exclude date, ct value from imputation
# encounter_rv_cv_pairs_conf_imp = encounter_rv_cv_pairs_conf %>%
#   #-ct_rv,-ct_cv,
#   dplyr::select(-covidshot1date, -covidshot2date, -covidshot3date )
# # Select columns to add
# encounter_rv_cv_pairs_conf_bind = encounter_rv_cv_pairs_conf %>%
#   ungroup() %>%
#   dplyr::select(ct_rv,ct_cv,covidshot1date, covidshot2date, covidshot3date)
# #Run imputation
# encounter_rv_cv_pairs_conf <- missRanger(encounter_rv_cv_pairs_conf_imp)
# encounter_rv_cv_pairs_conf = cbind(encounter_rv_cv_pairs_conf, encounter_rv_cv_pairs_conf_bind)
```

# Exporting prepped data for Analysis

``` r
encounter_rv_cv_pairs_conf
```

    ## # A tibble: 22,536 × 18
    ## # Groups:   individual [8,348]
    ##    individual         date       present_rv_n present_sars_cov_2_n present_rsv_n
    ##    <chr>              <date>            <int>                <int>         <int>
    ##  1 00053caa07b05a9e3… 2022-07-18            0                    0             0
    ##  2 00053caa07b05a9e3… 2022-07-21            0                    0             0
    ##  3 00053caa07b05a9e3… 2022-07-24            0                    0             0
    ##  4 00053caa07b05a9e3… 2022-07-27            0                    1             0
    ##  5 0006bbaa54094a386… 2020-11-21            0                    0             0
    ##  6 0006bbaa54094a386… 2021-03-08            0                    0             0
    ##  7 000a08b251c0514f5… 2021-04-30            0                    0             0
    ##  8 000a08b251c0514f5… 2021-08-21            0                    0             0
    ##  9 000a08b251c0514f5… 2021-09-13            0                    1             0
    ## 10 000a08b251c0514f5… 2021-09-21            0                    0             0
    ## # ℹ 22,526 more rows
    ## # ℹ 13 more variables: present_adv_n <int>, present_iav_n <int>,
    ## #   present_ibv_n <int>, sex <chr>, agebin <fct>, symptoms <pq__text>,
    ## #   sample <chr>, ct_Adenovirus <dbl>, ct_IAV <dbl>, ct_IBV <dbl>,
    ## #   ct_RSV <dbl>, ct_Rhinovirus <dbl>, ct_SARS_CoV_2 <dbl>

``` r
saveRDS(encounter_rv_cv_pairs_conf, "encounter_data_pairs.rds")
```

# If excluding Co-infections

``` r
# If Excluding Co-infections
encounter_rv_cv_pairs<- encounter_rv_cv_pairs %>%
  filter(!(present_SARS_CoV_2 & present_Adenovirus) &
           !(present_SARS_CoV_2 & present_RSV) &
           !(present_SARS_CoV_2 & present_IAV) &
           !(present_SARS_CoV_2 & present_IBV) &
           !(present_SARS_CoV_2 & present_Rhinovirus) &
           !(present_Adenovirus & present_RSV) &
           !(present_Adenovirus & present_IAV) &
           !(present_Adenovirus & present_IBV) &
           !(present_Adenovirus & present_Rhinovirus) &
           !(present_RSV & present_IAV) &
           !(present_RSV & present_IBV) &
           !(present_RSV & present_Rhinovirus) &
           !(present_IAV & present_IBV) &
           !(present_IAV & present_Rhinovirus) &
           !(present_IBV & present_Rhinovirus) &
           !(present_SARS_CoV_2 & present_Adenovirus & present_RSV) &
           !(present_SARS_CoV_2 & present_Adenovirus & present_IAV) &
           !(present_SARS_CoV_2 & present_Adenovirus & present_IBV) &
           !(present_SARS_CoV_2 & present_Adenovirus & present_Rhinovirus) &
           !(present_SARS_CoV_2 & present_RSV & present_IAV) &
           !(present_SARS_CoV_2 & present_RSV & present_IBV) &
           !(present_SARS_CoV_2 & present_RSV & present_Rhinovirus) &
           !(present_SARS_CoV_2 & present_IAV & present_IBV) &
           !(present_SARS_CoV_2 & present_IAV & present_Rhinovirus) &
           !(present_SARS_CoV_2 & present_IBV & present_Rhinovirus) &
           !(present_Adenovirus & present_RSV & present_IAV) &
           !(present_Adenovirus & present_RSV & present_IBV) &
           !(present_Adenovirus & present_RSV & present_Rhinovirus) &
           !(present_Adenovirus & present_IAV & present_IBV) &
           !(present_Adenovirus & present_IAV & present_Rhinovirus) &
           !(present_Adenovirus & present_IBV & present_Rhinovirus) &
           !(present_RSV & present_IAV & present_IBV) &
           !(present_RSV & present_IAV & present_Rhinovirus) &
           !(present_RSV & present_IBV & present_Rhinovirus) &
           !(present_IAV & present_IBV & present_Rhinovirus) &
           !(present_SARS_CoV_2 & present_Adenovirus & present_RSV & present_IAV) &
           !(present_SARS_CoV_2 & present_Adenovirus & present_RSV & present_IBV) &
           !(present_SARS_CoV_2 & present_Adenovirus & present_RSV & present_Rhinovirus) &
           !(present_SARS_CoV_2 & present_Adenovirus & present_IAV & present_IBV) &
           !(present_SARS_CoV_2 & present_Adenovirus & present_IAV & present_Rhinovirus) &
           !(present_SARS_CoV_2 & present_Adenovirus & present_IBV & present_Rhinovirus) &
           !(present_SARS_CoV_2 & present_RSV & present_IAV & present_IBV) &
           !(present_SARS_CoV_2 & present_RSV & present_IAV & present_Rhinovirus) &
           !(present_SARS_CoV_2 & present_RSV & present_IBV & present_Rhinovirus) &
           !(present_SARS_CoV_2 & present_IAV & present_IBV & present_Rhinovirus) &
           !(present_Adenovirus & present_RSV & present_IAV & present_IBV) &
           !(present_Adenovirus & present_RSV & present_IAV & present_Rhinovirus) &
           !(present_Adenovirus & present_RSV & present_IBV & present_Rhinovirus) &
           !(present_Adenovirus & present_IAV & present_IBV & present_Rhinovirus) &
           !(present_RSV & present_IAV & present_IBV & present_Rhinovirus))
```
