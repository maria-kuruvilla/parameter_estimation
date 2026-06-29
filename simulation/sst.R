library(here)
library(ersst)
library(tidyverse)


sst_download(years = 1925:2025, months = 4:7, save.dir = here("data", "sst_ersst"),
             version = 5)

sst <- sst_load(1925:2025, 4:7, here("data", "sst_ersst"), version = 5)


# subset data 

sst_subset <- sst_subset_space(sst, 
                               lat.min = 45, 
                               lat.max = 51,
                               lon.min = 225,
                               lon.max = 237)

sst_df <- sst_dataframe(sst_subset) %>% 
  mutate(lon = ifelse(lon > 180, lon - 360, lon)) %>% 
  filter(lat == 48, lon == -126) %>% 
  group_by(year) %>% 
  summarize(spring_ersst = mean(sst))

#save
write_csv(sst_df, here("data", "sst_ersst", "ersst_spring.csv"))
