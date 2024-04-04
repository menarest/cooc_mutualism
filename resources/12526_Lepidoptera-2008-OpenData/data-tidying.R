library(tidyverse)
library(readxl)
# install.packages("writexl")
library(writexl)

lepi_12526 <- read_excel("resources/BExIS_Lepidoptera_12526/12526.xlsx", 
                         sheet = "12526_Butterflies")

colnames(lepi_12526)

EP_Plot <- lepi_12526 %>%
  select(EPID, PlotID) %>%
  unique() %>%
  arrange(EPID)

write_xlsx(EP_Plot, "resources/BExIS_Lepidoptera_12526/EP_to_Plot.xlsx")

