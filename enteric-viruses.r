################################################
#                                              #
#    Title: Enteric Viruses (provisional)      #
#                                              #
#    Author: Lorenzo Mella                     #
#                                              #
#    Consultancy-Style Dissertation Project    #
#                                              #
#    School of Mathematics                     #
#                                              #
#    University of Edinburgh                   #
#                                              #
#    Academic Year: 2019/20                    #
#                                              #
################################################




##############################
#  Packages and Definitions  #
##############################


library(lubridate)
library(tidyverse)
library(ggmap)
library(patchwork)
library(survival)
library(survminer)


# For logging
printf = function(frmt, ...) { cat(sprintf(frmt, ...)) }

# Quickly standardise variables
standardise = function(x) { (x - mean(x)) / sd(x) }
standardize = standardise

count_NAs = function(dframe, vars, na_encoding = NA) {
  if (is.na(na_encoding)) {
    dframe %>%
      summarise_at(vars, ~sum(is.na(.x))) %>%
      as.matrix() %>%
      t() %>%
      knitr::kable()
  } else {
    dframe %>%
      summarise_at(vars, ~sum(.x == na_encoding)) %>%
      as.matrix() %>%
      t() %>%
      knitr::kable()
   }
}

pca_transform = function(m) {
  svd_result = svd(m)
  svd_result[[2]] %*% diag(svd_result[[1]])
}

# One-liner for a deviance residual test (returns the p-value)
deviance_test = function(reg_obj) {
  pval = pchisq(sum(residuals(reg_obj, type = "deviance")^2),
                reg_obj$df.residual)
  return (printf("Deviance residuals test: p = %f\n", pval))
}

# Latitude and longitute decimal converter
# (use the sign of the degrees as N/S and E/W)
geo_degrees_to_decimal = function(lat, long) {
  lat_decimal = lat[1] + lat[2] / 60 + lat[3] / 3600
  long_decimal = long[1] + long[2] / 60 + long[3] / 3600
  return (c(lat_decimal, long_decimal))
}

# Convert from geolocations to (exact) geodetic distance
# (Arguments in deg; output in km)
haversine = function(latitude_1, longitude_1, latitude_2, longitude_2) {
  phi_1 = latitude_1 * pi / 180
  phi_2 = latitude_2 * pi / 180
  theta = (longitude_2 - longitude_1) * pi / 180
  # Squared geodetic segment on the unit-sphere
  sq_segment = sin(0.5 * (phi_2 - phi_1))^2 +
    cos(phi_1) * cos(phi_2) * sin(theta / 2)^2
  # Average Earth radius times angle between locations
  return (6.371e3 * 2 * atan2(sqrt(sq_segment), sqrt(1 - sq_segment)))
}

# Pearson phi statistic on two *binary* vectors
pearson_phi = function(var1, var2) {
  # Pearform Pearson chi^2 test
  pearson_chi2_obj = chisq.test(x = var1, y = var2)
  # Extract the chi^2 statistic
  phi_stat = as.numeric(pearson_chi2_obj$statistic)
  # Compute the sign of phi using the contingency table
  phi_mat = pearson_chi2_obj$observed
  phi_sign = sign(phi_mat[1, 1] * phi_mat[2, 2] -
                    phi_mat[1, 2] * phi_mat[2, 1])
  return (phi_sign * sqrt(phi_stat / length(var1)))
}

# Table of Pearson phi statistics
pearson_phi_cols = function(data, varnames) {
  stat_vals = matrix(0, nrow = length(varnames),
                     ncol = length(varnames))
  rownames(stat_vals) = varnames
  colnames(stat_vals) = varnames
  for (i in 1:length(varnames)) {
    for (j in 1:length(varnames)) {
         stat_vals[i, j] = pearson_phi(data %>% pull(varnames[i]),
                                       data %>% pull(varnames[j]))
    }
  }
  return (stat_vals)
}

# Multiple Fisher Exact tests
fisher_exact_test_cols = function(data, varnames1, varnames2 = NULL,
                                  simulate = FALSE) {
  if (is.null(varnames2)) {
    varnames2 = varnames1
  }
  p_vals = matrix(0, nrow = length(varnames1),
                     ncol = length(varnames2))
  rownames(p_vals) = varnames1
  colnames(p_vals) = varnames2
  for (i in 1:length(varnames1)) {
    for (j in 1:length(varnames2)) {
      p_vals[i, j] = fisher.test(data %>% pull(varnames1[i]),
                                 data %>% pull(varnames2[j]),
                                 simulate.p.value = simulate)$p.value
    }
  }
  return (p_vals)
}

hoslem_for_logreg = function(model, g = NULL) {
  if (is.null(g)) {
    g = length(coefficients(model))
  }
  ResourceSelection::hoslem.test(model$y,
                               fitted(model),
                               g = g)
}


##################################
#  Dataset Loading and Cleaning  #
##################################


### Load dataset ###
vizions_tb = read_csv("./Data for enteric virus MSc.csv")


### Clean names ###

# Clean automatically disambiguated names for clarity
vizions_tb = vizions_tb %>%
  rename(Norovirus2.Ct.value = Ct.value_1,
         Norovirus1.Ct.value = Ct.value_2,
         Aichivirus.Ct.value = Ct.value_3,
         Adenovirus.Ct.value = Ct.value_4,
         Sapovirus.Ct.value = Ct.value_5,
         Astrovirus.Ct.value = Ct.value_6)

# Name relevant variable subsets for convenience
rt_pcr_vars = vizions_tb %>% select(42:55) %>% colnames()
deep_sequencing_vars = vizions_tb %>% select(56:87) %>% colnames()
symptom_covariates = vizions_tb %>% select(18:22) %>% colnames()
social_covariates = c("KeepAnimal", "KillingAnimal",
                      "EatCookRawMeat", "ContactDiar")
date_vars = c("Date of hospital entry", "AdminDate",
              "DateOnset", "DateDischOrDeath")
common_viruses = c("Rotavirus", "Norovirus", "Kobuvirus",
                   "Mastadenovirus", "Sapovirus", "Mamastrovirus")
uncommon_viruses = c("Alphapapillomavirus", "Alphapolyomavirus",
		     "Alphatorquevirus", "Betapapillomavirus",
		     "Betapolyomavirus", "Betatorquevirus",
		     "Bocaparvovirus", "Cardiovirus", "Circovirus",
		     "Cosavirus", "Cytomegalovirus", "Enterovirus",
		     "Gammatorquevirus", "Gemycircularvirus",
		     "Gemykibivirus", "Gemykrogvirus", "Husavirus",
		     "Lymphocryptovirus", "Morbillivirus",
		     "Parechovirus", "Picobirnavirus",
		     "Porprismacovirus", "Protoparvovirus",
		     "Rubulavirus", "Salivirus", "Unclassified virus")
blood_test_vars = vizions_tb %>% select(23:28) %>% colnames()
water_source_vars = c("Tap", "Well", "Rain", "River",
                      "Pond", "Bottled")


### Geography ###

# Generate macro-regions
vizions_tb = vizions_tb %>%
  mutate(macroregion = ifelse(CentrallyCity == "Dong Thap", "South",
                              ifelse(CentrallyCity %in%
                                       c("Dak Lak", "Dak Nong",
                                         "Khanh Hoa"),
                                     "Centre-South",
                                     "Centre"))) %>%
  mutate(macroregion = factor(macroregion,
                              levels = c("South", "Centre-South",
                                         "Centre")))


# We compute the distances from the closest densely inhabited settlement
dense_settlements = c("Cai Tau Ha", "Sa Dec", "Lap Vo", "TP. Cao Lanh",
                      "Hong Ngu", "Sa Rai", "My An", "Lai Vung",
                      "Tram Chim", "Thanh Binh", "My Tho",
                      "Cam Ranh", "Nha Trang", "Dien Khanh",
                      "Ninh Hoa", "Van Gia", "Ea Drang", "Buon Ho",
                      "Quang Phu", "Buon Ma Thuot", "Phuoc An",
                      "Krong Kmar", "Gia Nghia", "Dak Mil", "Hue",
                      "Dong Ha", "Quang Tri", "Kien Giang",
                      "Dong Hoi", "Ba Don", "Hoan Lao")

dense_settlement_coors = matrix(c(10.259731, 105.870438,
                                  10.2943716, 105.7588147,
10.3610793, 105.5204952, 10.4620322, 105.6357948,
10.8083867, 105.3417921, 10.8732222, 105.4623666,
10.523953, 105.843144, 10.2880817, 105.6596409,
10.670013, 105.564687, 10.561903, 105.484114,
10.4436825, 105.6979442,
11.9166696, 109.14861, 12.2431693, 109.1898675,
12.2597102, 109.1012101, 12.4934379, 109.1270865,
12.6996814, 109.2244852,
13.2034839, 108.2090235, 12.9187424, 108.2668299,
12.8144271, 108.0790559, 12.6796827, 108.0447368,
12.7100065, 108.3138426, 12.5057661, 108.3327763,
12.0005986, 107.6960259, 12.4488897, 107.6250921,
16.4638013, 107.5821911,
16.8172845, 107.1010378, 16.7546907, 107.1895518,
17.2282422, 106.7809474, 17.4695049, 106.6177781,
17.7558619, 106.4203402, 17.5833118, 106.5338218),
byrow = TRUE, ncol = 2)


# Compute the geodetic distance of all cohort members
# from the closest dense settlement
ward_lats = vizions_tb %>% pull(LATITUDE)
ward_longs = vizions_tb %>% pull(LONGITUDE)

ward_city_distance = numeric(nrow(vizions_tb))
for (i in 1:nrow(vizions_tb)) {
  ward_city_distance[i] = min(haversine(ward_lats[i], ward_longs[i],
                                        dense_settlement_coors[, 1],
                                        dense_settlement_coors[, 2]))
}
remove(ward_lats, ward_longs)

# Add resulting distances to dataframe
vizions_tb = vizions_tb %>%
  add_column(ward_city_distance, .after = "LATITUDE")

# Define a "rural area" indicator (a distance of 8 km from city centre
# curiously balances the viruses)
vizions_tb = vizions_tb %>%
  mutate(rural_area = ward_city_distance > 8)



### Format conversions ###

# Convert deep-sequencing results into logical
vizions_tb = vizions_tb %>%
 mutate_at(deep_sequencing_vars, .funs = ~.x == 1)

# Convert date variables from string to Date objects
vizions_tb = vizions_tb %>%
  mutate_at(date_vars, .funs = ~as.Date(.x, format = "%d/%m/%Y"))

# Convert social covariates and bacterial infection into logical
# (first, check for 9. We expect only a 1-2 encoding with no 9)
vizions_tb %>% summarise_at(c(social_covariates, "IF_Bacterial"),
                              ~any(.x == 9))

vizions_tb = vizions_tb %>%
  mutate_at(c("KeepAnimal", "KillingAnimal",
              "EatCookRawMeat", "IF_Bacterial"), .funs = ~.x == 1)

# Gender as factor
vizions_tb = vizions_tb %>%
  mutate(Gender = as.factor(ifelse(Gender == 1, "Male", "Female")))

# Change the 1-2-9 notation into 1-0-NA
vizions_tb = vizions_tb %>%
  mutate(ContactDiar = replace(ContactDiar, ContactDiar == 2, 0),
         BloodStool = replace(BloodStool, BloodStool == 2, 0),
         MucoidStool = replace(MucoidStool, MucoidStool == 2, 0),
         AbdominalPain = replace(AbdominalPain, AbdominalPain == 2, 0),
         ThreeDaysFever = replace(ThreeDaysFever, ThreeDaysFever == 2, 0)) %>%
  mutate(ContactDiar = replace(ContactDiar, ContactDiar == 9, NA),
         BloodStool = replace(BloodStool, BloodStool == 9, NA),
         MucoidStool = replace(MucoidStool, MucoidStool == 9, NA),
         AbdominalPain = replace(AbdominalPain, AbdominalPain == 9, NA),
         ThreeDaysFever = replace(ThreeDaysFever, ThreeDaysFever == 9, NA)) %>%
  mutate_at(c("ContactDiar", "BloodStool", "MucoidStool",
              "AbdominalPain", "ThreeDaysFever"),
            .funs = ~as.logical(.x))

# KnownTemp as logical (TRUE-FALSE instead of 1-NA)
vizions_tb = vizions_tb %>%
  mutate(KnownTemp = !is.na(Temp))

# Water sources as logical (TRUE-FALSE instead of 1-0)
vizions_tb = vizions_tb %>%
  mutate_at(water_source_vars, ~as.logical(.x))

### Coinfection count ###

# It is, then, safe to assume that the the NAs in the is_coinf variable
# should be interpreted as zeros
vizions_tb = vizions_tb %>%
  mutate(is_coinf = ifelse(is.na(is_coinf), 0, is_coinf))

# Add a factor representing coinfections
coinf_catchall = 2L
vizions_tb = vizions_tb %>%
  mutate(coinf_class = cut(is_coinf,
                        breaks = c(-Inf, 0:(coinf_catchall - 1), Inf),
                        labels = c(as.character(0:(coinf_catchall - 1)),
                                   sprintf("%d+", coinf_catchall))))

# Add a numerical version of the factor for regression
vizions_tb = vizions_tb %>%
  mutate(saturated_coinf = ifelse(is_coinf < 2, is_coinf, 2))

# Indicators of any common and uncommon virus infection
vizions_tb = vizions_tb %>% mutate(has_common_virus = Rotavirus | Norovirus |
                        Kobuvirus | Mastadenovirus |
                        Sapovirus | Mamastrovirus,
                        has_uncommon_virus = Alphapapillomavirus |
                        Alphapolyomavirus | Alphatorquevirus |
                        Betapapillomavirus | Betapolyomavirus |
                        Betatorquevirus | Bocaparvovirus |
                        Cardiovirus | Circovirus | Cosavirus |
                        Cytomegalovirus | Enterovirus | Gammatorquevirus |
                        Gemycircularvirus | Gemykibivirus | Gemykrogvirus |
                        Husavirus | Lymphocryptovirus | Morbillivirus |
                        Parechovirus | Picobirnavirus | Porprismacovirus |
                        Protoparvovirus | Rubulavirus | Salivirus |
                        `Unclassified virus`)



### Age ###

# Partition in three age clusters
vizions_tb = vizions_tb %>%
  mutate(age_cluster = cut(Age,
                         breaks = c(-Inf, 4, 65, Inf),
                         labels = c("Infant", "5-65", "Elderly")))

# Partition in regularly spaced age groups
age_group_step = 5
age_group_breaks = seq(age_group_step,
                       max(vizions_tb$Age), by = age_group_step)

age_group_labels = sprintf("Age_%d_to_%d",
                           age_group_breaks[1:length(age_group_breaks)] -
                             age_group_step + 1,
                         age_group_breaks[1:length(age_group_breaks)])

vizions_tb = vizions_tb %>%
  mutate(age_group =
    cut(Age, breaks = c(-Inf, age_group_breaks, Inf),
        labels = c(age_group_labels,
        sprintf("%d+", age_group_breaks[length(age_group_breaks)] + 1))))


### Water sources ###

# Add columns detailing water-source class
# (encoded as both boolean columns and a single factor column)
vizions_tb = vizions_tb %>%
  mutate(clean_water = Tap | Bottled | Rain | River,
         unclean_water = Well | Pond) %>%
  mutate(only_clean_water_source = clean_water & !unclean_water,
         only_unclean_water_source = unclean_water & !clean_water,
         mixed_water_source = unclean_water & clean_water) %>%
  select(-c(clean_water, unclean_water)) %>%
  mutate(water_source = as.factor(ifelse(only_clean_water_source,
                               "Only clean",
                               "Unclean water")))


### Time and season ###

# Add length of treatment
vizions_tb = vizions_tb %>%
  mutate(length_of_treatment = as.numeric(DateDischOrDeath - AdminDate))

# Add length of disease
vizions_tb = vizions_tb %>%
  mutate(length_of_illness = as.numeric(DateDischOrDeath - DateOnset))


# Seasonality deduced from DateOnset

# Horrible non-vectorised date-to-season converter
# (I didn't implement the boundary corrections. It shouldn't matter)
season = function(dates) {
  prev_wsolstice = as.Date("2000-12-21")
  sequinox = as.Date("2001-03-20")
  ssolstice = as.Date("2001-06-21")
  fequinox = as.Date("2001-09-23")
  current_wsolstice = as.Date("2001-12-21")
  seasons = character(length(dates))
  for (i in 1:length(dates)) {
    trudeau = lubridate::year(dates[i])
    lubridate::year(prev_wsolstice) = trudeau - 1
    lubridate::year(sequinox) = trudeau
    lubridate::year(ssolstice) = trudeau
    lubridate::year(fequinox) = trudeau
    lubridate::year(current_wsolstice) = trudeau
    if (dates[i] %within% interval(prev_wsolstice, sequinox - 1)) {
      seasons[i] = "Winter"
    } else if (dates[i] %within% interval(sequinox, ssolstice - 1)) {
      seasons[i] = "Spring"
    } else if (dates[i] %within% interval(ssolstice, fequinox - 1)) {
      seasons[i] = "Summer"
    } else if (dates[i] %within% interval(fequinox, current_wsolstice - 1)) {
      seasons[i] = "Autumn"
    } else {
      seasons[i] = "Winter"
    }
  }
  return (seasons)
}

vizions_tb = vizions_tb %>% 
  add_column(SeasonOnset = factor(season(vizions_tb$DateOnset),
                                  levels = c("Spring", "Summer",
                                             "Autumn", "Winter"),
                                  ordered = F),
             .after = "DateOnset")


# For when we focus on the Dong Thap area
vizions_dong_thap = vizions_tb %>% filter(CentrallyCity == "Dong Thap")


######################
#  Data Exploration  #
######################


# Subject prevalence per province and district
vizions_tb %>% 
  group_by(CentrallyCity, ProvincialCity) %>%
  tally() %>% knitr::kable()

# Subject prevalence per observed age
vizions_tb %>% 
  group_by(Age) %>%
  tally() %>% knitr::kable()

# Counts of patients positive to each virus (deep sequencing)
vizions_tb %>% 
  select(all_of(deep_sequencing_vars)) %>%
  summarise_at(deep_sequencing_vars, sum) %>% t()

# Patients that have at least one common or uncommon virus
vizions_tb %>% 
  group_by(has_uncommon_virus | has_common_virus) %>%
  tally()

vizions_tb %>%
  filter(!has_common_virus) %>% 
  group_by(is_coinf) %>%
  tally()

# Over/underdispersion in the number of coinfections by water source
vizions_tb %>% filter(Tap) %>%
  summarise(mean(is_coinf), var(is_coinf))
vizions_tb %>% filter(Bottled) %>%
  summarise(mean(is_coinf), var(is_coinf))
vizions_tb %>% filter(River) %>%
  summarise(mean(is_coinf), var(is_coinf))
vizions_tb %>% filter(Rain) %>%
  summarise(mean(is_coinf), var(is_coinf))
vizions_tb %>% filter(Well) %>%
  summarise(mean(is_coinf), var(is_coinf))

# Over/underdispersion in the number of coinfections by age group
vizions_tb %>% filter(age_cluster == "Infant") %>%
  summarise(mean(is_coinf), var(is_coinf))
vizions_tb %>% filter(age_cluster == "5-65") %>%
  summarise(mean(is_coinf), var(is_coinf))
vizions_tb %>% filter(age_cluster == "Elderly") %>%
  summarise(mean(is_coinf), var(is_coinf))

# Over/underdispersion in the number of coinfections by gender
vizions_tb %>% filter(Gender == "Female") %>%
  summarise(mean(is_coinf), var(is_coinf))
vizions_tb %>% filter(Gender == "Male") %>%
  summarise(mean(is_coinf), var(is_coinf))

# Over/underdispersion in the number of coinfections by social covariate
vizions_tb %>% filter(KeepAnimal) %>%
  summarise(mean(is_coinf), var(is_coinf))
vizions_tb %>% filter(KillingAnimal) %>%
  summarise(mean(is_coinf), var(is_coinf))
vizions_tb %>% filter(EatCookRawMeat) %>%
  summarise(mean(is_coinf), var(is_coinf))
vizions_tb %>% filter(ContactDiar) %>%
  summarise(mean(is_coinf), var(is_coinf))

# Frequency of deep-sequenced viral infections
vizions_tb %>% group_by(age_cluster) %>% 
  summarise_at(deep_sequencing_vars, sum) %>% t()

# No NAs in the social covariates
count_NAs(vizions_tb, social_covariates)

# Confusion matrix of blood-test variables with coinfections
vizions_tb %>% select(all_of(blood_test_vars), is_coinf) %>%
  cor(use = "pairwise.complete.obs")

# Not many NAs in the blood tests
count_NAs(vizions_tb, blood_test_vars)

# Most NAs in temperatures taken correspond to ThreeDaysFever == FALSE
vizions_tb %>% select(ThreeDaysFever, Temp) %>%
  mutate(NA_temp_xor_three_days_fever = xor(ThreeDaysFever, is.na(Temp)))

# Correlation between animal-contact variables
vizions_dong_thap %>% group_by(age_cluster) %>% 
  summarise(mean(KeepAnimal * KillingAnimal), # Scalar product
            mean(KeepAnimal * EatCookRawMeat),
            mean(EatCookRawMeat * KillingAnimal),
            # Pearson correlation coefficients (better with real vars)
            cor(KeepAnimal, KillingAnimal),
            cor(KeepAnimal, EatCookRawMeat),
            cor(EatCookRawMeat, KillingAnimal),
            # Pearson 'phi' coefficient (better with factor vars)
            pearson_phi(KeepAnimal, KillingAnimal),
            pearson_phi(KeepAnimal, EatCookRawMeat),
            pearson_phi(EatCookRawMeat, KillingAnimal)) %>% t()



########################
#  Data Visualisation  #
########################


### Coinfections in detail ###
# Variables should be binary or logical
pair_count = function(data, varnames, varnames2 = NULL,
                      suppress_diag = FALSE) {
  # Set default if only one vector of variables is used
  if (is.null(varnames2)) {
    varnames2 = varnames
  }
  count_matrix = matrix(0, nrow = length(varnames),
                        ncol = length(varnames2))
  rownames(count_matrix) = varnames
  colnames(count_matrix) = varnames2
  for (var1 in varnames) {
    for (var2 in varnames2) {
      if (!suppress_diag | var2 != var1) {
        count_matrix[var1, var2] = sum(data[[var1]] & data[[var2]])
      }
    }
  }
  return (count_matrix)
}

coinfections_in_detail =
  as_tibble(pair_count(vizions_tb,                                             c(common_viruses, "has_uncommon_virus"), suppress_diag = TRUE),
                                   rownames = "virus1") %>%
  gather(key = "virus2", value = "num_cases", -virus1)

coinfections_in_detail = coinfections_in_detail %>%
  mutate(virus1 =
           factor(virus1,
                  levels = c(common_viruses, "has_uncommon_virus")),
         virus2 =
           factor(virus2,
                  levels = c("has_uncommon_virus", common_viruses[6:1])))

coinfections_in_detail %>%
  ggplot(aes(virus1, virus2, col = num_cases, fill = num_cases,
             label = num_cases)) +
  geom_tile() +
  geom_text(col = "black") +
  ggtitle("Coinfection Pairs") +
  theme_minimal() +
  scale_x_discrete(name = "",
                   labels = c(common_viruses, "Non-enteric virus")) +
  scale_y_discrete(name = "",
                   labels = c("Non-enteric virus", common_viruses[6:1])) +
  scale_fill_gradient2(name = "Coinfection\ncount",
                       low = "white",  high = "magenta") +
  scale_color_gradient2(name = "Coinfection\ncount",
                        low = "white", high = "magenta") +
  theme(axis.text.x = element_text(angle = 90, size = 14),
        axis.text.y = element_text(size = 14),
        title = element_text(size = 15))
ggsave("./figures/coinfection_pairs_matrix.pdf", width = 7, height = 6.5)


# Coinfections by gender (all ages)
vizions_dong_thap %>%
  ggplot(aes(x = is_coinf, fill = Gender)) +
  geom_histogram(position = "dodge")

# Coinfections by season
vizions_tb  %>%
  ggplot(aes(x = SeasonOnset, fill = as.factor(is_coinf))) +
  geom_bar(position = "fill") +
  facet_grid(cols = vars(age_cluster)) +
  scale_x_discrete(name = "Season") +
  scale_y_continuous(name = "") +
  scale_fill_discrete(name = "Number of infections") +
  ggtitle("Proportion of Infections per Season") +
  theme(title = element_text(size = 18),
        axis.text.x = element_text(size = 13, angle = 45,
                                   hjust = 1),
        axis.text.y = element_text(size = 13),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13))
ggsave("./figures/coinf_by_season.pdf", width = 9, height = 5)

# Coinfections by age-group
vizions_tb %>% 
  ggplot(aes(coinf_class, fill = age_cluster)) +
  geom_bar(position = "dodge2") +
  labs(x = "Number of coinfections", y = "Absolute count",
       fill = "Age group", title = "Number of Infections vs Age Group") +
    theme(title = element_text(size = 18),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13))
ggsave("./figures/coinf_age_distribution.pdf", width = 7, height = 5)

# Coinfections by water quality
vizions_tb %>% 
  ggplot(aes(coinf_class, fill = water_source)) +
  geom_bar(position = "dodge2") +
  labs(x = "Number of coinfections", y = "Absolute count",
       fill = "Water source")
ggsave("./figures/coinf_water_distribution.pdf", width = 7, height = 5)

# Coinfections by gender and new age groupings
vizions_tb %>%
  ggplot(aes(x = coinf_class, fill = Gender)) +
  geom_bar(position = "dodge2") +
  facet_grid(cols = vars(age_cluster)) +
  scale_x_discrete(name = "Number of (co)infections") +
  scale_y_continuous(name = "") +
  ggtitle("Number of Infections vs Age and Gender") +
  theme(title = element_text(size = 18),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13))
ggsave("./figures/coinf_gender_age.pdf", width = 7, height = 5)

# Proportion of coinfections classes vs Age cluster or group
vizions_tb %>% # (cluster)
  ggplot(aes(x = age_cluster, fill = coinf_class)) +
  geom_bar(position = "fill")

vizions_tb %>% # (group)
  ggplot(aes(x = age_group, fill = coinf_class)) +
  geom_bar(position = "fill")

# Plot of virus prevalence
plot_name_list = c(common_viruses, "Alphatorquevirus", "Betatorquevirus",
                   "Enterovirus", "Gammatorquevirus", "Parechovirus",
                   "Picobirnavirus", "Salivirus", "Non-enteric virus",
                   "Unclassified virus")
vizions_tb %>%
  mutate(`Non-enteric virus` = Alphapapillomavirus + Alphapolyomavirus +
                              Betapapillomavirus + Betapolyomavirus +
                              Bocaparvovirus + Cardiovirus + Circovirus +
                              Cosavirus + Cytomegalovirus +
                              Gemycircularvirus +
                              Gemykibivirus + Gemykrogvirus + Husavirus +
                              Lymphocryptovirus + Morbillivirus +
                              Porprismacovirus + Protoparvovirus +
                              Rubulavirus) %>% 
  summarise_at(plot_name_list, sum) %>%
  pivot_longer(cols = 1:15, names_to = "Virus") %>%
  mutate(Virus = factor(Virus, levels = plot_name_list)) %>% 
  ggplot(aes(x = Virus, y = value, fill = c(rep("a", 6), rep("b", 8), "c"))) +
  geom_col() +
  ggtitle(label = "Frequency of infections",
          subtitle = "(enteric and other notable viruses)") +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Number of infections") +
  scale_fill_discrete(name = "Virus type",
                      labels = c("Enteric", "Other", "Unclassified")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
ggsave(filename = "./figures/infections_per_virus.pdf",
       width = 7, height = 5)

# This show that a happy cutoff for the Infant "cluster"
# should be > 5 year
dev.new()
vizions_tb %>% filter(Age <= 30) %>% mutate(age_cluster = Age >= 1) %>% 
  ggplot(aes(x = age_cluster, fill = has_uncommon_virus)) +
  geom_bar(position = "fill")
dev.new()
vizions_tb %>% filter(Age <= 30) %>% mutate(age_cluster = Age >= 2) %>% 
  ggplot(aes(x = age_cluster, fill = has_uncommon_virus)) +
  geom_bar(position = "fill")
dev.new()
vizions_tb %>% filter(Age <= 30) %>% mutate(age_cluster = Age >= 3) %>% 
  ggplot(aes(x = age_cluster, fill = has_uncommon_virus)) +
  geom_bar(position = "fill")
dev.new()
vizions_tb %>% filter(Age <= 30) %>% mutate(age_cluster = Age >= 4) %>% 
  ggplot(aes(x = age_cluster, fill = has_uncommon_virus)) +
  geom_bar(position = "fill")
dev.new()
vizions_tb %>% filter(Age <= 30) %>% mutate(age_cluster = Age >= 5) %>% 
  ggplot(aes(x = age_cluster, fill = has_uncommon_virus)) +
  geom_bar(position = "fill")
dev.new()
vizions_tb %>% filter(Age <= 30) %>% mutate(age_cluster = Age >= 6) %>% 
  ggplot(aes(x = age_cluster, fill = has_uncommon_virus)) +
  geom_bar(position = "fill")


# Age group by coinfections (just Dong Thap)
vizions_dong_thap %>%
  mutate(age_group = cut(Age,
                         breaks = c(-Inf, 5, Inf),
                         labels = c("Infant", "Child-to-adult"))) %>% 
  ggplot(aes(is_coinf, fill = age_group)) +
  geom_bar(position = "fill")

# Coinfections by age group
vizions_tb %>%# filter(age_group %in% Adult_group) %>% 
  ggplot(aes(x = coinf_class, fill = age_group)) +
  geom_bar(position = "fill")

# Coinfections by water source type
vizions_tb %>%
  ggplot(aes(x = water_source, fill = coinf_class)) +
  geom_bar(position = "fill") +
  facet_grid(cols = vars(age_cluster))


#######################
#  Survival Analysis  #
#######################


# Different means of length of illness by gender?
vizions_tb %>%
  ggplot(aes(length_of_illness, fill = Gender)) +
  geom_histogram(position = "dodge") +
  scale_x_continuous(name = "Length of illness (days)",
                     limits = c(0, 25)) +
  scale_fill_discrete(name = "Gender") +
  scale_y_continuous(name = "Absolute frequency") +
  ggtitle("Distribution of Length of Illness by Gender")
ggsave("./figures/length_of_illness_by_gender.pdf",
       width = 7, height = 5)

survtest_gender_tb = vizions_tb %>% filter(is_coinf > 0 &
                                             length_of_illness <= 40)
surv_obj = Surv(time = survtest_gender_tb %>% pull(length_of_illness),
                event = rep(TRUE, length.out = nrow(survtest_gender_tb)))
fit_for_graph = survfit(surv_obj ~ Gender,
                        data = survtest_gender_tb)
ggsurvplot(fit = fit_for_graph, survtest_gender_tb, pval = TRUE,
           title = "Gender affects length of illness",
           font.x = 15, font.y = 15, font.tickslab = 15,
           font.legend = 15, font.title = 20,
           legend = "bottom")
survdiff(surv_obj ~ Gender, data = survtest_gender_tb)
ggsave(filename = "./figures/surv-gender.pdf",
       width = 7, height = 6)
remove(survtest_gender_tb)

# Different means of length of illness by age?
vizions_tb %>%
  ggplot(aes(length_of_illness, fill = age_cluster)) +
  geom_histogram(position = "dodge") +
  scale_x_continuous(name = "Length of illness (days)",
                     limits = c(0, 25)) +
  scale_fill_discrete(name = "Age group") +
  scale_y_continuous(name = "Absolute frequency") +
  ggtitle("Distribution of Length of Illness by Age")
ggsave("./figures/length_of_illness_by_age.pdf",
       width = 7, height = 5)

survtest_age_tb = vizions_tb %>% filter(is_coinf > 0 &
                                          length_of_illness <= 40) %>%
  select(-Age) %>% 
  rename(Age = age_cluster)
surv_obj = Surv(time = survtest_age_tb %>% pull(length_of_illness),
                event = rep(TRUE, length.out = nrow(survtest_age_tb)))
fit_for_graph = survfit(surv_obj ~ Age,
                        data = survtest_age_tb)
ggsurvplot(fit = fit_for_graph, survtest_age_tb, pval = TRUE,
           title = "Age affects length of illness",
           font.x = 15, font.y = 15, font.tickslab = 15,
           font.legend = 15, font.title = 20,
           legend = "bottom")
survdiff(surv_obj ~ Age, data = survtest_age_tb)
ggsave(filename = "./figures/surv-age.pdf",
       width = 7, height = 6)
remove(survtest_age_tb)

# Different means of length of illness by rural or urban residence?
vizions_tb %>%
  ggplot(aes(length_of_illness, fill = rural_area)) +
  geom_bar(position = "dodge")


survtest_rural_area_tb = vizions_tb %>% filter(is_coinf > 0 &
                                          length_of_illness <= 40)
surv_obj = Surv(time = survtest_rural_area_tb %>% pull(length_of_illness),
                event = rep(TRUE, length.out = nrow(survtest_rural_area_tb)))
fit_for_graph = survfit(surv_obj ~ rural_area,
                        data = survtest_rural_area_tb)
ggsurvplot(fit = fit_for_graph, survtest_rural_area_tb, pval = TRUE)
survdiff(surv_obj ~ rural_area, data = survtest_rural_area_tb)
remove(survtest_rural_area_tb)

# Different means of length of illness by macroregion?
vizions_tb %>%
  ggplot(aes(length_of_illness, fill = macroregion)) +
  geom_bar(position = "dodge")

survtest_macroregion_tb = vizions_tb %>%
  filter(is_coinf > 0 &
           length_of_illness <= 40 &
           macroregion != "Centre")
surv_obj = Surv(time = survtest_macroregion_tb %>%
                  pull(length_of_illness),
                event = rep(TRUE,
                            length.out = nrow(survtest_macroregion_tb)))
fit_for_graph = survfit(surv_obj ~ macroregion,
                        data = survtest_macroregion_tb)
ggsurvplot(fit = fit_for_graph, survtest_macroregion_tb, pval = TRUE)
survdiff(surv_obj ~ macroregion, data = survtest_macroregion_tb)
remove(survtest_macroregion_tb)

# Different means of length of illness by coinfection class?
wilcox.test(vizions_tb %>% filter(coinf_class == "1") %>%
         pull(length_of_illness),
       vizions_tb %>% filter(coinf_class == "2+") %>%
         pull(length_of_illness))

wilcox.test(vizions_tb %>% filter(coinf_class == "0") %>%
         pull(length_of_illness),
       vizions_tb %>% filter(coinf_class != "0") %>%
         pull(length_of_illness))

vizions_tb %>%
  ggplot(aes(length_of_illness, fill = coinf_class)) +
  geom_bar(position = "dodge")


survtest_coinf_tb = vizions_tb %>%
  filter(length_of_illness <= 40) %>%
  mutate(Enteric_infection = has_common_virus & (is_coinf > 1))
surv_obj = Surv(time = survtest_coinf_tb %>%
                  pull(length_of_illness),
                event = rep(TRUE,
                            length.out = nrow(survtest_coinf_tb)))
fit_for_graph = survfit(surv_obj ~ Enteric_infection,
                        data = survtest_coinf_tb)
ggsurvplot(fit = fit_for_graph, survtest_coinf_tb, pval = TRUE)
survdiff(surv_obj ~ Enteric_infection, data = survtest_coinf_tb)
remove(survtest_coinf_tb)


survtest_coinf_tb = vizions_tb %>%
  filter(is_coinf > 0 & length_of_illness <= 40 & has_common_virus) %>%
  rename(Infections = coinf_class)
surv_obj = Surv(time = survtest_coinf_tb %>%
                  pull(length_of_illness),
                event = rep(TRUE,
                            length.out = nrow(survtest_coinf_tb)))
fit_for_graph = survfit(surv_obj ~ Infections,
                        data = survtest_coinf_tb)
ggsurvplot(fit = fit_for_graph, survtest_coinf_tb, pval = TRUE,
           title = "No distinction between one or more infections",
           font.x = 15, font.y = 15, font.tickslab = 15,
           font.legend = 15, font.title = 20,
           legend = "bottom")
survdiff(surv_obj ~ Infections, data = survtest_coinf_tb)
ggsave(filename = "./figures/surv-coinf-1-2.pdf",
       width = 7, height = 6)
remove(survtest_coinf_tb)

### Length of hospitalisation vs age group ###

# Model with age and gender stratification

surv_obj = Surv(time = vizions_tb %>% pull(length_of_illness),
                event = rep(TRUE, length.out = nrow(vizions_tb)))

fit_3 = coxph(surv_obj ~ coinf_class + macroregion +
              + strata(Age, Gender),
              data = vizions_tb)

summary(fit_3) # significant non-null

zph_3 = cox.zph(fit_3)

plot(zph_3[4])


fit_3_stepwise = step(fit_3, direction = "both", trace = 0)

summary(fit_3_stepwise)

surv_obj = Surv(time = vizions_tb %>% filter(is_coinf > 0)
                %>% pull(length_of_illness),
                event = rep(TRUE, length.out = nrow(vizions_tb %>%
                                    filter(is_coinf > 0))))

summary(survfit(surv_obj ~ coinf_class + age_cluster,
         data = vizions_tb %>% filter(is_coinf > 0)))


# Model with viral interactions
surv_obj = Surv(time = vizions_dong_thap %>% pull(length_of_illness),
                event = rep(TRUE, length.out = nrow(vizions_dong_thap)))

fit_4 = coxph(surv_obj ~ (Rotavirus + Norovirus + Kobuvirus +
                            Mastadenovirus + Sapovirus +
                            Mamastrovirus) * (Salivirus +
                                                Picobirnavirus +
                                                Parechovirus +
                                                Gammatorquevirus +
                                                Enterovirus +
                                                Betatorquevirus +
                                                Alphatorquevirus) +
                strata(age_cluster, Gender),
              data = vizions_dong_thap)

summary(fit_4)


# Test by age and virus
survtest_tb = vizions_tb %>%
  filter(is_coinf > 0 & age_cluster == "Infant")
surv_obj = Surv(time = survtest_tb %>% pull(length_of_illness),
                event = rep(TRUE, length.out = nrow(survtest_tb)))
fit_for_graph = survfit(surv_obj ~ Sapovirus + strata(Gender),
                        data = survtest_tb)
ggsurvplot(fit = fit_for_graph, survtest_tb, pval = TRUE)
survdiff(surv_obj ~ Sapovirus + strata(Gender),
                        data = survtest_tb)

survtest_tb = vizions_tb %>%
  filter(is_coinf > 0 & age_cluster == "5-65")
surv_obj = Surv(time = survtest_tb %>% pull(length_of_illness),
                event = rep(TRUE, length.out = nrow(survtest_tb)))
fit_for_graph = survfit(surv_obj ~ Rotavirus + strata(Gender),
                        data = survtest_tb)
ggsurvplot(fit = fit_for_graph, survtest_tb, pval = TRUE)
survdiff(surv_obj ~ Rotavirus + strata(Gender),
                        data = survtest_tb)

survtest_tb = vizions_tb %>%
  filter(is_coinf > 0 & age_cluster == "5-65")
surv_obj = Surv(time = survtest_tb %>% pull(length_of_illness),
                event = rep(TRUE, length.out = nrow(survtest_tb)))
fit_for_graph = survfit(surv_obj ~ Norovirus + strata(Gender),
                        data = survtest_tb)
ggsurvplot(fit = fit_for_graph, survtest_tb, pval = TRUE)
survdiff(surv_obj ~ Norovirus + strata(Gender),
                        data = survtest_tb)

survtest_tb = vizions_tb %>%
  filter(is_coinf > 0 & age_cluster == "Elderly")
surv_obj = Surv(time = survtest_tb %>% pull(length_of_illness),
                event = rep(TRUE, length.out = nrow(survtest_tb)))
fit_for_graph = survfit(surv_obj ~ Sapovirus + strata(Gender),
                        data = survtest_tb)
ggsurvplot(fit = fit_for_graph, survtest_tb, pval = TRUE)
survdiff(surv_obj ~ Sapovirus + strata(Gender),
                        data = survtest_tb)


# Test by region and virus
survtest_tb = vizions_tb %>%
  filter(is_coinf > 0 & macroregion == "South")
surv_obj = Surv(time = survtest_tb %>% pull(length_of_illness),
                event = rep(TRUE, length.out = nrow(survtest_tb)))
fit_for_graph = survfit(surv_obj ~ Rotavirus + strata(Gender, age_cluster),
                        data = survtest_tb)
ggsurvplot(fit = fit_for_graph, survtest_tb, pval = TRUE)
survdiff(surv_obj ~ Rotavirus + strata(Gender),
                        data = survtest_tb)

survtest_tb = vizions_tb %>%
  filter(is_coinf > 0 & macroregion == "South")
surv_obj = Surv(time = survtest_tb %>% pull(length_of_illness),
                event = rep(TRUE, length.out = nrow(survtest_tb)))
fit_for_graph = survfit(surv_obj ~ Norovirus + strata(Gender, age_cluster),
                        data = survtest_tb)
ggsurvplot(fit = fit_for_graph, survtest_tb, pval = TRUE,
           legend = "right")
survdiff(surv_obj ~ Norovirus + strata(Gender, age_cluster),
                        data = survtest_tb)

survtest_tb = vizions_tb %>%
  filter(is_coinf > 0 & macroregion == "South")
surv_obj = Surv(time = survtest_tb %>% pull(length_of_illness),
                event = rep(TRUE, length.out = nrow(survtest_tb)))
fit_for_graph = survfit(surv_obj ~ Mastadenovirus + strata(Gender, age_cluster),
                        data = survtest_tb)
ggsurvplot(fit = fit_for_graph, survtest_tb, pval = TRUE)
survdiff(surv_obj ~ Mastadenovirus + strata(Gender, age_cluster),
                        data = survtest_tb)


# Graphs
surv_obj = Surv(time = vizions_tb %>% filter(is_coinf >= 0)
                %>% pull(length_of_illness),
                event = rep(TRUE, length.out = nrow(vizions_tb %>%
                                    filter(is_coinf >= 0))))

fit_for_graph = survfit(surv_obj ~ coinf_class,
                        data = vizions_tb %>%
                          filter(is_coinf >= 0))

ggsurvplot_facet(fit = fit_for_graph, vizions_tb %>%
                                        filter(is_coinf >= 0),
                 facet.by = "age_cluster")

# Graphs showing estimated survival functions for the three age clusters
surv_times = vizions_tb %>% filter(age_cluster == "Infant") %>%
  pull(`Length of stay`)
surv_obj = Surv(time = surv_times,
                event = rep(TRUE, length.out = length(surv_times)))

fit_for_graph = survfit(surv_obj ~ coinf_class,
              data = vizions_tb %>% filter(age_cluster == "Infant"))

ggsurvplot(fit_for_graph, pval = T, conf.int = T)

surv_times = vizions_tb %>% filter(age_cluster == "5-65") %>%
  pull(`Length of stay`)
surv_obj = Surv(time = surv_times,
                event = rep(TRUE, length.out = length(surv_times)))

fit_for_graph = survfit(surv_obj ~ coinf_class,
              data = vizions_tb %>% filter(age_cluster == "5-65"))

ggsurvplot(fit_for_graph, pval = T, conf.int = T)

surv_times = vizions_tb %>% filter(age_cluster == "Elderly") %>%
  pull(`Length of stay`)
surv_obj = Surv(time = surv_times,
                event = rep(TRUE, length.out = length(surv_times)))

fit_for_graph = survfit(surv_obj ~ coinf_class,
              data = vizions_tb %>% filter(age_cluster == "Elderly"))

ggsurvplot(fit_for_graph, pval = T, conf.int = T)



#################################################
#  Coinfection Generalised Regression Analysis  #
#################################################


# Fisher exact test to compare coinfection vs no infection between
# genders, for each age strata
vizions_dong_thap %>% filter(age_cluster == "Elderly") %>%
  group_by(is_coinf >= 2) %>%
  summarise(male = sum(Gender == "Male"),
            female = sum(Gender == "Female")) %>%
  select(-`is_coinf >= 2`) %>%
  as.matrix() %>%
  fisher.test


### Stepwise selection on the full model (4 kinds) ###

coinf_single_obj = glm(formula = is_coinf ~ age_cluster + SeasonOnset +
                         (Tap + Well + River + Rain) *
                         ward_city_distance + Gender + KeepAnimal +
                         KillingAnimal + EatCookRawMeat + ContactDiar,
                    family = poisson(log),
                    data = vizions_tb %>%
                      ## slice(train_idx) %>% 
                      mutate_at("ward_city_distance", standardise) %>%
                      # Use the complete cases only (~10 rows less)
                      filter_at(c("ContactDiar"),
                                .vars_predicate = ~!is.na(.x)))

summary(coinf_single_obj)

car::vif(coinf_single_obj)

fisher.test(vizions_tb %>% pull(Tap), vizions_tb %>% pull(Bottled))
fisher.test(vizions_tb %>% pull(Tap), vizions_tb %>% pull(River))
fisher.test(vizions_tb %>% pull(River), vizions_tb %>% pull(Bottled))

stepwise_corrected_obj = step(coinf_single_obj,
                              direction = "both", trace = 0)

summary(stepwise_corrected_obj)

deviance_test(stepwise_corrected_obj)

par(mfrow = c(2, 2))
plot(stepwise_corrected_obj)


### Cross-validation on the four candidate models ###

# age_cluster + River + Rain + Bottled + ward_city_distance + Gender +
# KillingAnimal
vizions_cv_tb = vizions_tb %>%
  mutate_at("ward_city_distance", standardise) %>% 
  filter(!is.na(ContactDiar))
test_predictions = numeric(length(vizions_cv_tb))
for (i in 1:nrow(vizions_cv_tb)) {
  trained_model = glm(is_coinf ~ age_cluster + River + Rain + Bottled +
                        ward_city_distance + Gender + KillingAnimal,
                      family = poisson,
                      data = vizions_cv_tb %>% slice(-i))
    row = vizions_cv_tb %>% slice(i) %>%
      summarise(intercept = 1, age_cluster5_65 = age_cluster == "5-65",
                age_cluster_Elderly = age_cluster == "Elderly",
                river = River, rain = Rain, bottled = Bottled,
                wc_distance = ward_city_distance,
                GenderMALE = Gender == "Male",
                killing = KillingAnimal)
  test_predictions[i] =
    exp(coefficients(trained_model) %*% t(row %>% as.matrix()))
}
cv_error = mean((vizions_cv_tb %>% pull(is_coinf) - test_predictions)^2)
remove(vizions_cv_tb, test_predictions, trained_model, row)
print(cv_error)

# age_cluster + Tap + Well + ward_city_distance + Gender + Killing
vizions_cv_tb = vizions_tb %>%
  mutate_at("ward_city_distance", standardise) %>% 
  filter(!is.na(ContactDiar))
test_predictions = numeric(length(vizions_cv_tb))
for (i in 1:nrow(vizions_cv_tb)) {
  trained_model = glm(is_coinf ~ age_cluster + Tap + Well +
                        ward_city_distance + Gender + KillingAnimal,
                      family = poisson,
                      data = vizions_cv_tb %>% slice(-i))
    row = vizions_cv_tb %>% slice(i) %>%
      summarise(intercept = 1, age_cluster5_65 = age_cluster == "5-65",
                age_cluster_Elderly = age_cluster == "Elderly",
                tap = Tap, well = Well,
                wc_distance = ward_city_distance,
                male = Gender == "Male",
                killing = KillingAnimal)
  test_predictions[i] =
    exp(coefficients(trained_model) %*% t(row %>% as.matrix()))
}
cv_error = mean((vizions_cv_tb %>% pull(is_coinf) - test_predictions)^2)
remove(vizions_cv_tb, test_predictions, trained_model, row)
print(cv_error)

# age_cluster + River + Rain + Bottled + ward_city_distance +
# Gender + KillingAnimal + Rain:ward_city_distance
vizions_cv_tb = vizions_tb %>%
  mutate_at("ward_city_distance", standardise) %>% 
  filter(!is.na(ContactDiar))
test_predictions = numeric(length(vizions_cv_tb))
for (i in 1:nrow(vizions_cv_tb)) {
  trained_model = glm(is_coinf ~ age_cluster + River + Rain + Bottled +
                        ward_city_distance + Gender + KillingAnimal +
                        Rain:ward_city_distance,
                      family = poisson,
                      data = vizions_cv_tb %>% slice(-i))
    row = vizions_cv_tb %>% slice(i) %>%
      summarise(intercept = 1, age_cluster5_65 = age_cluster == "5-65",
                age_cluster_Elderly = age_cluster == "Elderly",
                river = River, rain = Rain, bottled = Bottled,
                wc_distance = ward_city_distance,
                GenderMALE = Gender == "Male",
                killing = KillingAnimal,
                rain_ward = Rain * ward_city_distance)
  test_predictions[i] =
    exp(coefficients(trained_model) %*% t(row %>% as.matrix()))
}
cv_error = mean((vizions_cv_tb %>% pull(is_coinf) - test_predictions)^2)
remove(vizions_cv_tb, test_predictions, trained_model, row)
print(cv_error)

# age_cluster + Tap + Rain + ward_city_distance + Gender +
# KillingAnimal + Rain:ward_city_distance
vizions_cv_tb = vizions_tb %>%
  mutate_at("ward_city_distance", standardise) %>% 
  filter(!is.na(ContactDiar))
test_predictions = numeric(length(vizions_cv_tb))
for (i in 1:nrow(vizions_cv_tb)) {
  trained_model = glm(is_coinf ~ age_cluster + Tap + Rain +
                        ward_city_distance + Gender +
                        KillingAnimal + Rain:ward_city_distance,
                      family = poisson,
                      data = vizions_cv_tb %>% slice(-i))
    row = vizions_cv_tb %>% slice(i) %>%
      summarise(intercept = 1, age_cluster5_65 = age_cluster == "5-65",
                age_cluster_Elderly = age_cluster == "Elderly",
                tap = Tap, rain = Rain,
                wc_distance = ward_city_distance,
                GenderMALE = Gender == "Male",
                killing = KillingAnimal,
                rain_ward = Rain * ward_city_distance)
  test_predictions[i] =
    exp(coefficients(trained_model) %*% t(row %>% as.matrix()))
}
cv_error = mean((vizions_cv_tb %>% pull(is_coinf) - test_predictions)^2)
remove(vizions_cv_tb, test_predictions, trained_model, row)
print(cv_error)


### Winning model ###

# Fitting
winning_coinf_obj = glm(is_coinf ~ age_cluster + Gender +
                          ward_city_distance +
                          Tap + Well + KillingAnimal,
                        family = poisson,
                        data = vizions_tb %>%
                          filter(!is.na(ContactDiar)))

# Tabular information
summary(winning_coinf_obj)

exp(coefficients(winning_coinf_obj))

confint(winning_coinf_obj)

# Diagnostic plots
par(mfrow = c(2, 2))
plot(winning_coinf_obj)


# Ordinal regression
coinf_ordinal_obj = MASS::polr(formula = is_coinf ~ age_cluster +
                                 SeasonOnset + Rain + Pond +
                                 River + Bottled +
                                 ward_city_distance + Gender +
                                 KeepAnimal + KillingAnimal +
                                 EatCookRawMeat + ContactDiar,
                               data = vizions_tb %>%
                                 mutate_at("is_coinf", as.factor) %>% 
                           filter_at(c("Age", "Gender", "SeasonOnset",
                                       "ContactDiar", "KeepAnimal",
                                       "ward_city_distance",
                                       "KillingAnimal", "EatCookRawMeat",
                                       "Tap", "Well", "Rain", "Pond",
                                       "River", "Bottled"),
                                     .vars_predicate = ~!is.na(.x)))

summary(coinf_ordinal_obj)

stepwise_ordinal_obj = step(coinf_ordinal_obj, direction = "both",
                            trace = 0)

summary(stepwise_ordinal_obj)




######################################
#  Disease Severity vs Coinfections  #
######################################


# Test of dependence between single enteric infections and symptoms 
fisher_exact_test_cols(vizions_tb %>% filter(age_cluster == "Infant",
                                             Gender == "Female"),
                       c("Rotavirus", "Norovirus",
                         "Mastadenovirus", "Sapovirus",
                         "Mamastrovirus"),
                       symptom_covariates,
                       simulate = T)

fisher_exact_test_cols(vizions_tb,
                       "is_coinf",
                       symptom_covariates,
                       simulate = T)


# Logistic regression on Abdominal Pain
logistic_reg_obj = glm(formula = AbdominalPain ~ age_cluster +
                         SeasonOnset + Gender + ward_city_distance +
                         (Rotavirus + Norovirus +
                            Mastadenovirus + Sapovirus +
                            Mamastrovirus) + has_uncommon_virus,
                       family = binomial,
                       data = vizions_dong_thap %>%
                         mutate_at("ward_city_distance",
                                   standardise))

summary(logistic_reg_obj)

car::vif(logistic_reg_obj)

backward_selection_obj = step(logistic_reg_obj, direction = "both",
                              trace = 0)

summary(backward_selection_obj)

confint(backward_selection_obj)

hoslem_for_logreg(backward_selection_obj)

arm::binnedplot(fitted(backward_selection_obj), 
           residuals(backward_selection_obj, type = "response"))


# Logistic Regression on Three Days of Fever
logistic_reg_obj = glm(formula = ThreeDaysFever ~ age_cluster +
                         SeasonOnset + Gender + ward_city_distance +
                         (Rotavirus + Norovirus +
                            Mastadenovirus + Sapovirus +
                            Mamastrovirus) + has_uncommon_virus,
                       family = binomial,
                       data = vizions_dong_thap %>%
                         mutate_at("ward_city_distance",
                                   standardise))

summary(logistic_reg_obj)

car::vif(logistic_reg_obj)

backward_selection_obj = step(logistic_reg_obj, direction = "both",
                              trace = 0)

summary(backward_selection_obj)

confint(backward_selection_obj)

hoslem_for_logreg(backward_selection_obj)

arm::binnedplot(fitted(backward_selection_obj), 
           residuals(backward_selection_obj, type = "response"))




# Mucoid stool interacts with a coinf class
logistic_reg_obj = glm(formula = ThreeDaysFever ~ age_cluster + SeasonOnset +
                         Gender + ward_city_distance +
                         coinf_class,
                       family = binomial,
                       data = vizions_tb %>% mutate_at(c("is_coinf",
                                                         "ward_city_distance"),
                                                       standardise))

summary(logistic_reg_obj)
hoslem_for_logreg(logistic_reg_obj)

car::vif(logistic_reg_obj)

backward_selection_obj = step(logistic_reg_obj, direction = "both",
                              trace = 0)

summary(backward_selection_obj)
hoslem_for_logreg(backward_selection_obj)


# ThreeDaysFever interacts with a coinf class
logistic_reg_obj = glm(formula = ThreeDaysFever ~ age_cluster + SeasonOnset +
                         Gender + ward_city_distance +
                         coinf_class,
                       family = binomial,
                       data = vizions_tb %>% mutate_at(c("is_coinf",
                                                         "ward_city_distance"),
                                                       standardise))

summary(logistic_reg_obj)
deviance_test(logistic_reg_obj)

car::vif(logistic_reg_obj)

backward_selection_obj = step(logistic_reg_obj, direction = "both",
                              trace = 0)

summary(backward_selection_obj)

deviance_test(backward_selection_obj)

# BloodStool (NON STRATIFIED)
logistic_reg_obj = glm(formula = BloodStool ~ is_coinf +
                         (age_cluster + Gender + macroregion),
                       family = binomial,
                       data = vizions_tb %>% filter(is_coinf >= 1))
summary(logistic_reg_obj)
hoslem_for_logreg(logistic_reg_obj)

BIC(logistic_reg_obj)

logistic_reg_obj = glmer(formula = BloodStool ~ is_coinf +
                           (1 | macroregion) +
                         age_cluster + Gender,
                       family = binomial,
                       data = vizions_tb)
summary(logistic_reg_obj)
hoslem_for_logreg(logistic_reg_obj)

BIC(logistic_reg_obj)


# MucoidStool (NON STRATIFIED)
logistic_reg_obj = glm(formula = MucoidStool ~ age_cluster,
                       family = binomial,
                       data = vizions_tb %>% filter(is_coinf >= 1))

logistic_reg_obj = glm(formula = MucoidStool ~ is_coinf *
                         (Gender + I(age_cluster == "Infant") +
                           + macroregion),
                       family = binomial,
                       data = vizions_tb %>% filter(is_coinf >= 1))
summary(logistic_reg_obj)

hoslem_for_logreg(logistic_reg_obj)

mucoidstool_stepwise_obj = step(logistic_reg_obj, direction = "both",
                            trace = 0)

summary(mucoidstool_stepwise_obj)
hoslem_for_logreg(mucoidstool_stepwise_obj)

confint(logistic_reg_obj)

BIC(logistic_reg_obj)

random_effect_obj = glmer(formula = MucoidStool ~ is_coinf +
                            SeasonOnset +
                           (1 | macroregion) +
                         age_cluster + Gender,
                       family = binomial,
                       data = vizions_tb %>% filter(is_coinf >= 1))
summary(random_effect_obj)

BIC(random_effect_obj)


# AbdominalPain (Dong Thap)
# Infant vs noninfant is significant
fisher.test(vizions_dong_thap %>%
              mutate(is_infant = age_cluster == "Infant") %>%
              pull(is_infant),
          vizions_dong_thap %>% pull(AbdominalPain))

# 5-65 vs Elderly is not significant
fisher.test(vizions_dong_thap %>% filter(age_cluster != "Infant") %>%
              pull(age_cluster),
            vizions_dong_thap %>% filter(age_cluster != "Infant") %>%
              pull(AbdominalPain))

# Gender highly non-significant
fisher.test(vizions_dong_thap %>%
              pull(Gender),
            vizions_dong_thap %>%
              pull(AbdominalPain))

# HENCE, NO STRATIFICATION BY GENDER, ONLY INFANT VS NON-INFANT
logistic_reg_obj = glm(formula = ThreeDaysFever ~ coinf_class +
                         age_cluster + Gender,
                       family = binomial,
                       data = vizions_dong_thap)
summary(logistic_reg_obj)

hoslem_for_logreg(logistic_reg_obj)

# ThreeDaysFever (Dong Thap)
# Infant vs noninfant is significant
fisher.test(vizions_dong_thap %>%
              mutate(is_infant = age_cluster == "Infant") %>%
              pull(is_infant),
          vizions_dong_thap %>% pull(BloodStool))

# 5-65 vs Elderly is not significant
fisher.test(vizions_dong_thap %>% filter(age_cluster != "Infant") %>%
              mutate(is_elderly = age_cluster == "Elderly") %>% 
              pull(is_elderly),
            vizions_dong_thap %>% filter(age_cluster != "Infant") %>%
              pull(BloodStool))

# Gender highly non-significant
fisher.test(vizions_dong_thap %>%
              pull(Gender),
            vizions_dong_thap %>%
              pull(BloodStool))

logistic_reg_obj = glm(formula = ThreeDaysFever ~ coinf_class
                       + I(age_cluster == "Infant"),
                       family = binomial,
                       data = vizions_dong_thap)
summary(logistic_reg_obj)

hoslem_for_logreg(logistic_reg_obj)

confint(logistic_reg_obj)


# AbdominalPain (NON STRATIFIED) (FINAL - COINFECTIONS NOT SIGNIFICANT)
logistic_reg_obj = glm(formula = AbdominalPain ~ coinf_class *
                         (age_cluster + Gender) + 
                         SeasonOnset + macroregion,
                       family = binomial,
                       data = vizions_tb %>% filter(is_coinf >= 1))
summary(logistic_reg_obj)

hoslem_for_logreg(logistic_reg_obj)

abdpain_stepwise_obj = step(logistic_reg_obj, direction = "both",
                            trace = 0)

summary(abdpain_stepwise_obj)

logistic_reg_obj = glmer(formula = AbdominalPain ~ (saturated_coinf +
                         age_cluster) + Gender + SeasonOnset +
                            (1 | macroregion) + (age_cluster | macroregion),
                       family = binomial,
                       data = vizions_tb)
summary(logistic_reg_obj)
deviance_test(logistic_reg_obj)


# ThreeDaysFever (NON STRATIFIED)
logistic_reg_obj = lme4::glmer(formula = ThreeDaysFever ~
                                 saturated_coinf +
                         age_cluster + SeasonOnset +
                            + (1 | macroregion),
                       family = binomial,
                       data = vizions_tb %>% filter(has_common_virus))
summary(logistic_reg_obj)


logistic_reg_obj = glm(ThreeDaysFever ~ coinf_class +
                         (age_cluster + Gender + macroregion) +
                         SeasonOnset,
                       family = binomial,
                       data = vizions_tb %>%
                         filter(has_common_virus))
summary(logistic_reg_obj)

stepwise_obj = step(logistic_reg_obj, direction = "both", trace = 0)

summary(stepwise_obj)


# MIRRORING THE STEPWISE-SELECTED MODEL ABOVE
logistic_reg_obj = lme4::glmer(formula = ThreeDaysFever ~ is_coinf +
                         I(age_cluster == "Infant") +
                          has_uncommon_virus +
                         + (1 | macroregion),
                       family = binomial,
                       data = vizions_tb)
summary(logistic_reg_obj)


# ThreeDaysFever (FINAL)
logistic_reg_obj = glm(formula = ThreeDaysFever ~ coinf_class +
                         (age_cluster + macroregion + Gender) +
                         SeasonOnset,
family = binomial,
                       data = vizions_tb %>% filter(is_coinf >= 1))
summary(logistic_reg_obj)

hoslem_for_logreg(logistic_reg_obj)

stepwise_obj = step(logistic_reg_obj, direction = "both", trace = 0)
summary(stepwise_obj)

hoslem_for_logreg(stepwise_obj)

confint(stepwise_obj)


plot(0:5, plogis(coefficients(logistic_reg_obj)[1] +
                   (0:5) * coefficients(logistic_reg_obj)[2]),
     col = "red")
plot_points = vizions_tb %>% filter (age_cluster != "5-65" &
                                    Gender == "Male") %>%
             group_by(is_coinf) %>%
             summarise(means = mean(ThreeDaysFever)) %>%
             pull(means)
lines(0:(length(plot_points) - 1), plot_points)
remove(plot_points)

# NumberDiarEpi - graphs
vizions_tb %>% ggplot(aes(x = coinf_class, y = NumberDiarEpi,
                          fill = macroregion)) +
  geom_boxplot(position = "dodge2") +
  facet_grid(rows = vars(age_cluster))


# Over/underdispersion in the number of coinfections by age group
vizions_tb %>% filter(age_cluster == "Infant") %>%
  summarise(mean(NumberDiarEpi, na.rm = T),
            var(NumberDiarEpi, na.rm = T))
vizions_tb %>% filter(age_cluster == "5-65") %>%
  summarise(mean(NumberDiarEpi, na.rm = T),
            var(NumberDiarEpi, na.rm = T))
vizions_tb %>% filter(age_cluster == "Elderly") %>%
  summarise(mean(NumberDiarEpi, na.rm = T),
            var(NumberDiarEpi, na.rm = T))

# Over/underdispersion in the number of coinfections by gender
vizions_tb %>% filter(Gender == "Female") %>%
  summarise(mean(NumberDiarEpi, na.rm = T),
            var(NumberDiarEpi, na.rm = T))
vizions_tb %>% filter(Gender == "Male") %>%
  summarise(mean(NumberDiarEpi, na.rm = T),
            var(NumberDiarEpi, na.rm = T))

# Over/underdispersion in the number of coinfections by macro-region
vizions_tb %>% filter(macroregion == "South") %>%
  summarise(mean(NumberDiarEpi, na.rm = T),
            var(NumberDiarEpi, na.rm = T))
vizions_tb %>% filter(macroregion == "Centre-South") %>%
  summarise(mean(NumberDiarEpi, na.rm = T),
            var(NumberDiarEpi, na.rm = T))
vizions_tb %>% filter(macroregion == "Centre") %>%
  summarise(mean(NumberDiarEpi, na.rm = T),
            var(NumberDiarEpi, na.rm = T))


# NumberDiarEpi (NON STRATIFIED) (FINAL)
numberdiarepi_reg_obj = glm(formula = NumberDiarEpi ~ is_coinf *
                         (age_cluster + Gender + macroregion) +
                         SeasonOnset + has_uncommon_virus,
                         family = poisson,
                       data = vizions_tb %>%
                         filter(has_common_virus))
summary(numberdiarepi_reg_obj)

anova(numberdiarepi_reg_obj)

stepwise_obj = step(numberdiarepi_reg_obj, direction = "both", trace = 0)

summary(stepwise_obj)

confint(stepwise_obj)

deviance_test(stepwise_obj)

par(mfrow = c(2, 2))
plot(stepwise_obj)


numberdiarepi_reg_obj = glm(formula = NumberDiarEpi ~
                              is_coinf *
                              (Gender + macroregion) +
                              has_uncommon_virus,
                            family = quasipoisson,
                            data = vizions_tb %>%
                              filter(has_common_virus))
summary(numberdiarepi_reg_obj)

deviance_test(numberdiarepi_reg_obj)



###############################################
#  Risk Factors Associated to Common Viruses  #
###############################################


# Rotavirus modelling
rotavirus_reg_obj = glm(formula = Norovirus ~ age_cluster +
                          ContactDiar + KeepAnimal +
                          KillingAnimal + EatCookRawMeat +
                          ward_city_distance * (Tap + Well + Pond +
                                                  Rain + River),
                        family = binomial,
                        data = vizions_dong_thap %>%
                          mutate_at(c("ward_city_distance", "Age"),
                                    standardise) %>% 
                          filter_at("ContactDiar",
                                    .vars_predicate = ~!is.na(.x)))

summary(rotavirus_reg_obj)

car::vif(rotavirus_reg_obj)

reduced_rotavirus_reg_obj = step(rotavirus_reg_obj,
                                 direction = "both",
                                 trace = FALSE)

summary(reduced_rotavirus_reg_obj)
hoslem_for_logreg(reduced_rotavirus_reg_obj)


# Norovirus modelling
norovirus_reg_obj = glm(formula = Norovirus ~ age_cluster +
                          ContactDiar + KeepAnimal +
                          KillingAnimal + EatCookRawMeat +
                          Tap + Pond + Rain + Well + River + Bottled,
                        family = binomial,
                        data = vizions_tb %>%
                          mutate_at(c("ward_city_distance", "Age"),
                                    standardise) %>%
                          filter_at(c("age_cluster",
                                      "ContactDiar", "KeepAnimal",
                                      "KillingAnimal", "EatCookRawMeat",
                                      "Tap", "Well", "River", "Bottled"),
                                    .vars_predicate = ~!is.na(.x)))


summary(norovirus_reg_obj)

reduced_norovirus_reg_obj = MASS::stepAIC(norovirus_reg_obj,
                                          trace = FALSE)

summary(reduced_norovirus_reg_obj)
pchisq(q = reduced_norovirus_reg_obj$deviance,
       df = reduced_norovirus_reg_obj$df.residual, lower.tail=FALSE)

# Kobuvirus modelling
kobuvirus_reg_obj = glm(formula = Kobuvirus ~ age_cluster +
                          ward_city_distance +
                          ContactDiar + KeepAnimal +
                          KillingAnimal + EatCookRawMeat +
                          Tap + Pond + Rain + Well + River + Bottled,
                        family = binomial,
                        data = vizions_tb %>%
                          mutate_at(c("ward_city_distance", "Age"),
                                    standardise) %>%
                          filter_at(c("age_cluster",
                                      "ContactDiar", "KeepAnimal",
                                      "KillingAnimal", "EatCookRawMeat",
                                      "Tap", "Well", "River", "Bottled"),
                                    .vars_predicate = ~!is.na(.x)))


summary(kobuvirus_reg_obj)

reduced_kobuvirus_reg_obj = MASS::stepAIC(kobuvirus_reg_obj,
                                          trace = FALSE)

summary(reduced_kobuvirus_reg_obj)
pchisq(q = reduced_kobuvirus_reg_obj$deviance,
       df = reduced_kobuvirus_reg_obj$df.residual, lower.tail=FALSE)


# Mastadenovirus modelling
mastadenovirus_reg_obj = glm(formula = Mastadenovirus ~ age_cluster +
                          ward_city_distance +
                          ContactDiar + KeepAnimal +
                          KillingAnimal + EatCookRawMeat +
                          Tap + Pond + Rain + Well + River + Bottled,
                        family = binomial,
                        data = vizions_tb %>%
                          mutate_at(c("ward_city_distance", "Age"),
                                    standardise) %>%
                          filter_at(c("age_cluster",
                                      "ContactDiar", "KeepAnimal",
                                      "KillingAnimal", "EatCookRawMeat",
                                      "Tap", "Well", "River", "Bottled"),
                                    .vars_predicate = ~!is.na(.x)))


summary(mastadenovirus_reg_obj)

reduced_mastadenovirus_reg_obj = step(mastadenovirus_reg_obj,
                                      direction = "both",
                                      trace = FALSE)

summary(reduced_mastadenovirus_reg_obj)



# Sapovirus modelling
sapovirus_reg_obj = glm(formula = Sapovirus ~ age_cluster +
                          ward_city_distance +
                          ContactDiar + KeepAnimal +
                          KillingAnimal + EatCookRawMeat +
                          Tap + Pond + Rain + Well + River + Bottled,
                        family = binomial,
                        data = vizions_tb %>%
                          mutate_at(c("ward_city_distance", "Age"),
                                    standardise) %>%
                          filter_at(c("age_cluster",
                                      "ContactDiar", "KeepAnimal",
                                      "KillingAnimal", "EatCookRawMeat",
                                      "Tap", "Well", "River", "Bottled"),
                                    .vars_predicate = ~!is.na(.x)))


summary(sapovirus_reg_obj)

reduced_sapovirus_reg_obj = MASS::stepAIC(sapovirus_reg_obj,
                                          trace = FALSE)

summary(reduced_sapovirus_reg_obj)
pchisq(q = reduced_sapovirus_reg_obj$deviance,
       df = reduced_sapovirus_reg_obj$df.residual, lower.tail=FALSE)

# Mamastrovirus modelling
mamastrovirus_reg_obj = glm(formula = Mamastrovirus ~ age_cluster +
                          ward_city_distance +
                          ContactDiar + KeepAnimal +
                          KillingAnimal + EatCookRawMeat +
                          Tap + Pond + Rain + Well + River + Bottled,
                        family = binomial,
                        data = vizions_tb %>%
                          mutate_at(c("ward_city_distance", "Age"),
                                    standardise) %>%
                          filter_at(c("age_cluster",
                                      "ContactDiar", "KeepAnimal",
                                      "KillingAnimal", "EatCookRawMeat",
                                      "Tap", "Well", "River", "Bottled"),
                                    .vars_predicate = ~!is.na(.x)))


summary(mamastrovirus_reg_obj)

reduced_mamastrovirus_reg_obj = MASS::stepAIC(mamastrovirus_reg_obj,
                                          trace = FALSE)

summary(reduced_mamastrovirus_reg_obj)
pchisq(q = reduced_mamastrovirus_reg_obj$deviance,
       df = reduced_mamastrovirus_reg_obj$df.residual, lower.tail=FALSE)


###########################################
#  Spatial Identification of the Viruses  #
###########################################


viet_bbox = list(left = 102.170435826, bottom = 8.59975962975,
                 right = 109.33526981, top = 23.3520633001)

data_presence_bbox = list(left = 104, bottom = 9,
                          right = 109.8, top = 18.5)

dong_thap_bbox = list(left = 105.10, bottom = 10.10,
                 right = 106, top = 11)

# Common-virus infections by age cluster
vizions_tb %>% group_by(age_cluster) %>%
  summarise_at(common_viruses, sum)

# Common vs uncommon virus infections by age cluster
vizions_tb %>% group_by(age_cluster) %>%
  filter(macroregion == "South") %>%
  summarise(common_virus = sum(has_common_virus),
            uncommon_virus = sum(has_uncommon_virus),
            both_infections = sum(has_common_virus | has_uncommon_virus))

# Common viruses by district
vizions_tb %>% group_by(CentrallyCity, ProvincialCity) %>%
  summarise_at(common_viruses, sum) %>%
  as.data.frame()

# Uncommon viruses by district
vizions_tb %>% group_by(CentrallyCity, ProvincialCity) %>%
  summarise_at(uncommon_viruses, sum) %>%
  as.data.frame()


vizions_tb %>% group_by(rural_area) %>%
  summarise_at(c("has_common_virus", "has_uncommon_virus"), sum)

vizions_tb %>% filter(rural_area) %>%
  group_by(ProvincialCity) %>% 
  summarise_at(common_viruses, sum)

# Try to visualise dives in: macroregion/rural vs nonrural area/common
# vs uncommon viruses/partition of the common/frequency of coinfections/age
vizions_tb %>%
  ggplot(aes(x = macroregion, fill = has_common_virus)) +
  geom_bar(position = "fill")

# Visualise viruses per city centre in Dong Thap
dev.new()
vizions_dong_thap %>%
  filter(age_cluster != "Infant") %>% 
  group_by(ProvincialCity) %>% 
  summarise_at(common_viruses, sum) %>%
  pivot_longer(cols = 2:7, names_to = "Enteric_virus") %>% 
  ggplot(aes(x = ProvincialCity, y = value, fill = Enteric_virus)) +
  geom_col(position = "stack")

### Visualise viruses per city centre in Dong Thap ###

# Main districts (Cao Lanh)
vizions_tb %>%
  filter(macroregion == "South") %>% 
  filter(ProvincialCity %in% c("TP. Cao Lanh", "Cao Lanh")) %>% 
  group_by(ProvincialCity) %>% 
  summarise_at(common_viruses, sum) %>%
  pivot_longer(cols = 2:7, names_to = "Enteric virus") %>% 
  ggplot(aes(x = ifelse(ProvincialCity == "Cao Lanh",
                        "Cao Lahn District",
                        "Cao Lahn City"),
             y = value, fill = `Enteric virus`)) +
  geom_col(width = 0.5) +
  xlab("District") +
  ylab("Number of infections") +
    ggtitle("Enteric Infections in Cao Lanh") +
  theme(title = element_text(size = 18),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13))
ggsave(filename = "./figures/common_viruses_populous_districts.pdf",
       width = 7, height = 6)

# Rural districts
vizions_tb %>%
  filter(macroregion == "South") %>% 
  filter(!(ProvincialCity %in% c("TP. Cao Lanh", "Cao Lanh",
                                 "Thi Xa Hong Ngu"))) %>% 
  group_by(ProvincialCity) %>% 
  summarise_at(common_viruses, sum) %>%
  pivot_longer(cols = 2:7, names_to = "Enteric virus") %>%
  ggplot(aes(x = ProvincialCity, y = value, fill = `Enteric virus`)) +
  geom_col() +
  xlab("District") +
  ylab("Number of infections") +
    ggtitle("Enteric Infections in Rural Dong Thap") +
  theme(title = element_text(size = 18),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 13),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13))
ggsave(filename = "./figures/common_viruses_lesser_districts.pdf",
       width = 7, height = 6)

# Visualise viruses per rural area EXPERIMENTAL (NOT VERY MEANINGFUL)
vizions_tb %>%
  filter(macroregion == "South") %>% 
  filter(rural_area) %>% 
  group_by(ProvincialCity) %>% 
  summarise_at(common_viruses, sum) %>%
  pivot_longer(cols = 2:7, names_to = "Enteric_virus") %>% 
  ggplot(aes(x = ProvincialCity, y = value, fill = Enteric_virus)) +
  geom_col()
dev.new()
vizions_tb %>%
  filter(macroregion == "South") %>% 
  filter(!rural_area) %>% 
  group_by(ProvincialCity) %>% 
  summarise_at(common_viruses, sum) %>%
  pivot_longer(cols = 2:7, names_to = "Enteric_virus") %>% 
  ggplot(aes(x = ProvincialCity, y = value, fill = Enteric_virus)) +
  geom_col()

### Uncommon viruses ###

# Visualise viruses per city centre in Dong Thap
vizions_tb %>%
  filter(macroregion == "South") %>% 
  group_by(ProvincialCity) %>% 
  summarise_at(uncommon_viruses, sum) %>%
  pivot_longer(cols = 2:7, names_to = "Enteric_virus") %>% 
  ggplot(aes(x = ProvincialCity, y = value, fill = Enteric_virus)) +
  geom_col()

# Visualise viruses per city centre in Dong Thap
vizions_tb %>%
  filter(macroregion != "South") %>% 
  group_by(ProvincialCity) %>% 
  summarise_at(uncommon_viruses, sum) %>%
  pivot_longer(cols = 2:7, names_to = "Enteric_virus") %>% 
  ggplot(aes(x = ProvincialCity, y = value, fill = Enteric_virus)) +
  geom_col()




# In Dong Thap only
dev.new()
ggmap(get_stamenmap(bbox = c(top = dong_thap_bbox[["top"]],
                             left = dong_thap_bbox[["left"]],
                             bottom = dong_thap_bbox[["bottom"]],
                             right = dong_thap_bbox[["right"]]),
                    zoom = 11), extent = "normal") +
  geom_point(data = vizions_dong_thap %>%
               filter(Rotavirus == TRUE | Norovirus == TRUE |
                        Kobuvirus == TRUE | Mastadenovirus == TRUE |
                        Sapovirus == TRUE | Mamastrovirus == TRUE) %>% 
               filter(T), 
             aes(x = LONGITUDE, y = LATITUDE, color = Rotavirus))
dev.new()
ggmap(get_stamenmap(bbox = c(top = dong_thap_bbox[["top"]],
                             left = dong_thap_bbox[["left"]],
                             bottom = dong_thap_bbox[["bottom"]],
                             right = dong_thap_bbox[["right"]]),
                    zoom = 11), extent = "normal") +
  geom_point(data = vizions_dong_thap %>%
               filter(Rotavirus == TRUE | Norovirus == TRUE |
                        Kobuvirus == TRUE | Mastadenovirus == TRUE |
                        Sapovirus == TRUE | Mamastrovirus == TRUE) %>% 
               filter(T), 
             aes(x = LONGITUDE, y = LATITUDE, color = Norovirus))
dev.new()
ggmap(get_stamenmap(bbox = c(top = dong_thap_bbox[["top"]],
                             left = dong_thap_bbox[["left"]],
                             bottom = dong_thap_bbox[["bottom"]],
                             right = dong_thap_bbox[["right"]]),
                    zoom = 11), extent = "normal") +
  geom_point(data = vizions_dong_thap %>%
               filter(Rotavirus == TRUE | Norovirus == TRUE |
                        Kobuvirus == TRUE | Mastadenovirus == TRUE |
                        Sapovirus == TRUE | Mamastrovirus == TRUE) %>% 
               filter(T), 
             aes(x = LONGITUDE, y = LATITUDE, color = Kobuvirus))
dev.new()
ggmap(get_stamenmap(bbox = c(top = dong_thap_bbox[["top"]],
                             left = dong_thap_bbox[["left"]],
                             bottom = dong_thap_bbox[["bottom"]],
                             right = dong_thap_bbox[["right"]]),
                    zoom = 11), extent = "normal") +
  geom_point(data = vizions_dong_thap %>%
               filter(Rotavirus == TRUE | Norovirus == TRUE |
                        Kobuvirus == TRUE | Mastadenovirus == TRUE |
                        Sapovirus == TRUE | Mamastrovirus == TRUE) %>% 
               filter(T), 
             aes(x = LONGITUDE, y = LATITUDE, color = Mastadenovirus))
dev.new()
ggmap(get_stamenmap(bbox = c(top = dong_thap_bbox[["top"]],
                             left = dong_thap_bbox[["left"]],
                             bottom = dong_thap_bbox[["bottom"]],
                             right = dong_thap_bbox[["right"]]),
                    zoom = 11), extent = "normal") +
  geom_point(data = vizions_dong_thap %>%
               filter(Rotavirus == TRUE | Norovirus == TRUE |
                        Kobuvirus == TRUE | Mastadenovirus == TRUE |
                        Sapovirus == TRUE | Mamastrovirus == TRUE) %>% 
               filter(T), 
             aes(x = LONGITUDE, y = LATITUDE, color = Sapovirus))
dev.new()
ggmap(get_stamenmap(bbox = c(top = dong_thap_bbox[["top"]],
                             left = dong_thap_bbox[["left"]],
                             bottom = dong_thap_bbox[["bottom"]],
                             right = dong_thap_bbox[["right"]]),
                    zoom = 11), extent = "normal") +
  geom_point(data = vizions_dong_thap %>%
               filter(Rotavirus == TRUE | Norovirus == TRUE |
                        Kobuvirus == TRUE | Mastadenovirus == TRUE |
                        Sapovirus == TRUE | Mamastrovirus == TRUE) %>% 
               filter(T), 
             aes(x = LONGITUDE, y = LATITUDE, color = Mamastrovirus))


# Density estimation map
ggmap(get_stamenmap(bbox = c(top = dong_thap_bbox[["top"]],
                             left = dong_thap_bbox[["left"]],
                             bottom = dong_thap_bbox[["bottom"]],
                             right = dong_thap_bbox[["right"]]),
                    zoom = 11, maptype = "toner-hybrid"), extent = "normal") +
  geom_density2d(data = vizions_dong_thap %>%
               filter(Rotavirus),
               aes(x = LONGITUDE, y = LATITUDE)) +
  geom_density2d(data = vizions_dong_thap %>%
               filter(Norovirus),
             aes(x = LONGITUDE, y = LATITUDE, color = "green")) +
  geom_density2d(data = vizions_dong_thap %>%
               filter(Kobuvirus),
             aes(x = LONGITUDE, y = LATITUDE, color = "orange")) +
  geom_density2d(data = vizions_dong_thap %>%
               filter(Mastadenovirus),
             aes(x = LONGITUDE, y = LATITUDE, color = "blue")) +
  geom_density2d(data = vizions_dong_thap %>%
               filter(Sapovirus),
             aes(x = LONGITUDE, y = LATITUDE, color = "magenta")) +
  geom_density2d(data = vizions_dong_thap %>%
               filter(Mamastrovirus),
             aes(x = LONGITUDE, y = LATITUDE, color = "yellow")) +
  scale_color_discrete(name = "Enteric Virus",
    labels = c("Rotavirus", "Mastadenovirus",  "Norovirus",
               "Kobuvirus", "Sapovirus", "Mamastrovirus"))


# Patchwork of single density estimations
# THEY ARE UNSURPRISINGLY CLUSTERED ON THE BIG CITY, CAO LANH
(ggmap(get_stamenmap(bbox = c(top = dong_thap_bbox[["top"]],
                             left = dong_thap_bbox[["left"]],
                             bottom = dong_thap_bbox[["bottom"]],
                             right = dong_thap_bbox[["right"]]),
                    zoom = 10, maptype = "toner-hybrid"),
      extent = "normal") +
   geom_density2d(data = vizions_dong_thap %>% 
                    filter(Rotavirus),
                 aes(x = LONGITUDE, y = LATITUDE), color = "red")) +
(ggmap(get_stamenmap(bbox = c(top = dong_thap_bbox[["top"]],
                             left = dong_thap_bbox[["left"]],
                             bottom = dong_thap_bbox[["bottom"]],
                             right = dong_thap_bbox[["right"]]),
                    zoom = 10, maptype = "toner-hybrid"),
       extent = "normal") +
    geom_density2d(data = vizions_dong_thap %>%
               filter(Norovirus),
               aes(x = LONGITUDE, y = LATITUDE))) +
(ggmap(get_stamenmap(bbox = c(top = dong_thap_bbox[["top"]],
                             left = dong_thap_bbox[["left"]],
                             bottom = dong_thap_bbox[["bottom"]],
                             right = dong_thap_bbox[["right"]]),
                    zoom = 10, maptype = "toner-hybrid"),
       extent = "normal") +
    geom_density2d(data = vizions_dong_thap %>%
               filter(Kobuvirus),
               aes(x = LONGITUDE, y = LATITUDE))) +
(ggmap(get_stamenmap(bbox = c(top = dong_thap_bbox[["top"]],
                             left = dong_thap_bbox[["left"]],
                             bottom = dong_thap_bbox[["bottom"]],
                             right = dong_thap_bbox[["right"]]),
                    zoom = 10, maptype = "toner-hybrid"),
       extent = "normal") +
    geom_density2d(data = vizions_dong_thap %>%
               filter(Mastadenovirus),
               aes(x = LONGITUDE, y = LATITUDE))) +
(ggmap(get_stamenmap(bbox = c(top = dong_thap_bbox[["top"]],
                             left = dong_thap_bbox[["left"]],
                             bottom = dong_thap_bbox[["bottom"]],
                             right = dong_thap_bbox[["right"]]),
                    zoom = 10, maptype = "toner-hybrid"),
       extent = "normal") +
    geom_density2d(data = vizions_dong_thap %>%
               filter(Sapovirus),
               aes(x = LONGITUDE, y = LATITUDE))) +
(ggmap(get_stamenmap(bbox = c(top = dong_thap_bbox[["top"]],
                             left = dong_thap_bbox[["left"]],
                             bottom = dong_thap_bbox[["bottom"]],
                             right = dong_thap_bbox[["right"]]),
                    zoom = 10, maptype = "toner-hybrid"),
       extent = "normal") +
    geom_density2d(data = vizions_dong_thap %>%
               filter(Mamastrovirus),
               aes(x = LONGITUDE, y = LATITUDE)))  +
  geom_density2d(data = vizions_dong_thap %>%
               filter(Kobuvirus),
             aes(x = LONGITUDE, y = LATITUDE, color = "orange")) +
  geom_density2d(data = vizions_dong_thap %>%
               filter(Mastadenovirus),
             aes(x = LONGITUDE, y = LATITUDE, color = "blue")) +
  geom_density2d(data = vizions_dong_thap %>%
               filter(Sapovirus),
             aes(x = LONGITUDE, y = LATITUDE, color = "magenta")) +
  geom_density2d(data = vizions_dong_thap %>%
               filter(Mamastrovirus),
             aes(x = LONGITUDE, y = LATITUDE, color = "yellow")) +
  scale_color_discrete(name = "Enteric Virus",
    labels = c("Rotavirus", "Mastadenovirus",  "Norovirus",
                        "Kobuvirus", "Sapovirus", "Mamastrovirus"))



# Viral density estimation (UNCOMMON VIRUSES)
dev.new()
ggmap(get_stamenmap(bbox = c(top = 19, #viet_bbox[["top"]],
                             left = viet_bbox[["left"]],
                             bottom = viet_bbox[["bottom"]],
                             right = viet_bbox[["right"]]),
                    zoom = 7), extent = "normal") +
  geom_point(data = vizions_tb %>%
               filter(Alphapapillomavirus == TRUE |
                        Alphapolyomavirus == TRUE |
                        Alphatorquevirus == TRUE |
                        Betapapillomavirus == TRUE |
                        Betapolyomavirus == TRUE |
                        Betatorquevirus == TRUE |
                        Bocaparvovirus == TRUE |
                        Cardiovirus == TRUE |
                        Circovirus == TRUE |
                        Cosavirus == TRUE |
                        Cytomegalovirus == TRUE |
                        Enterovirus == TRUE |
                        Gammatorquevirus == TRUE |
                        Gemycircularvirus == TRUE |
                        Gemykibivirus == TRUE |
                        Gemykrogvirus == TRUE |
                        Husavirus == TRUE |
                        Lymphocryptovirus == TRUE |
                        Morbillivirus == TRUE |
                        Parechovirus == TRUE |
                        Picobirnavirus == TRUE |
                        Porprismacovirus == TRUE |
                        Protoparvovirus == TRUE |
                        Rubulavirus == TRUE |
                        Salivirus == TRUE |
                        `Unclassified virus` == TRUE) %>%
               filter(age_cluster == "Infant"),
             aes(x = LONGITUDE, y = LATITUDE), h = 0.05)
