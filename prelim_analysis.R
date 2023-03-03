library(LifeIneq)
library(wpp2022)
library(readxl)
library(tidyverse)
library(scales)
library(janitor)
library(glue)
library(ggpubr)
library(ggrepel)
library(ggcorrplot)
Sys.setenv(LANG = "en")

# Age-specific data on mortality rates for male and female from wpp2022 package
# Age 0-100, years 1950-2100 (forecast)
data(mx1dt)

# Life table calculation from mortality rates
all_lt <- mx1dt %>%
  as_tibble() %>%
  select(name:age, mx = mxB) %>%
  distinct() %>%
  filter(year < 2023) %>% 
  group_by(name, year) %>%
  mutate(n  = c(rep(1, 100), Inf),
         ax =  n + 1 / mx - n / (1 - exp(1) ^ (-n * mx))) %>%
  group_by(name) %>%
  nest() %>% 
  mutate(data = map(data, ~ .x       %>%
                      group_by(year) %>%
                      mutate(qx = n * mx / (1 + (n - ax) * mx),
                             qx = ifelse(is.nan(qx), 1, qx),
                             px = 1 - qx, 
                             lx = head(c(1, cumprod(px)), -1) * 100000,  
                             dx = lx - c(tail(lx, -1), 0),
                             Lx = c(head(n * c(tail(lx, -1), 0), -1), tail(lx, 1) / tail(mx, 1)),
                             Tx = rev(cumsum(rev(Lx))),
                             ex = Tx / lx) %>%
                      ungroup())) %>%
  unnest(data) %>% 
  ungroup()

# GINI from https://data.worldbank.org/indicator/SI.POV.GINI
# World development indicators
GINI <- read_excel("API_SI.POV.GINI_DS2_en_excel_v2_4770509.xls", sheet = "Data", skip = 2) %>%
  select(name = `Country Name`, 5:ncol(.)) %>%
  pivot_longer(-name,
               names_to  = "year",
               values_to = "gini",
               values_drop_na = TRUE) %>%
  mutate(year = as.integer(year))

# GDP PPP from https://ourworldindata.org/grapher/gini-vs-gdp-per-capita?country=MEX~NER~LBN
gdp <- read_csv("gini-vs-gdp-per-capita.csv") %>%
  select(name = Entity, year = Year, gdp = 5) %>%
  na.omit()

# use default inequality measure
sd1 <- all_lt %>%
  inner_join(GINI, by = c("name", "year")) %>%
  inner_join(gdp,  by = c("name", "year")) %>%
  # To use data with LifeIneq it complains when we have an open interval
  filter(age < 96) %>%
  group_by(name, year) %>%
  # Calculate a lifetable column for the conditional standard deviation in lifetable ages at death
  mutate(exineq = ineq_sd(age = age,
                          dx  = dx,
                          ex  = ex,
                          ax  = ax)) %>% 
  slice(c(1)) %>% # use inequality in the first year of age as a measure 
  ungroup()

# Calculate a lifetable column for the conditional Gini coefficient of inequality in survivorship
# with age at death and remaining lif
dt2 <- all_lt %>%
  inner_join(GINI, by = c("name", "year")) %>%
  inner_join(gdp,  by = c("name", "year")) %>%
  filter(age < 96) %>%
  group_by(name, year) %>%
  mutate(exginiaad = ineq_gini(age = age,
                               dx  = dx,
                               ex  = ex,
                               ax  = ax,
                               distribution_type = "aad"),
         exginiarl = ineq_gini(age = age,
                               dx  = dx,
                               ex  = ex,
                               ax  = ax,
                               distribution_type = "rl")
  ) %>% 
  slice(c(1)) %>% 
  ungroup()

# almost identical
dt2 %>% 
  select(name, year, gini, exginiaad, exginiarl) %>%
  rename(GINI = gini) %>%
  mutate(df = abs(year - 2015)) %>% 
  group_by(name) %>% 
  filter(df == min(df)) %>%
  filter(year == max(year)) %>% 
  ungroup() %>%
  mutate(name = str_sub(name, 1, 3)) %>%
  pivot_longer(c(exginiaad, exginiarl),
               names_to = "var",
               values_to = "val",
               values_drop_na = TRUE) %>% 
  ggplot(aes(x = GINI, y = val, group = var, color = var)) + 
  geom_jitter() + 
  geom_smooth(method = "lm", se = FALSE)


# Correlation matrix for various inequality indicators
cor_mat <-
  dt2 %>%
  select(gdp, exginiarl, exginiaad, ex, gini) %>% 
  mutate(lgdp = log(gdp)) %>% 
  cor()

cor_mat[upper.tri(cor_mat,TRUE)] <- NA
ggcorrplot(cor_mat, lab = TRUE)


# conditional standard deviation vs GINI
sd1 %>% 
  select(name, year, exineq, gini) %>%
  rename(GINI = gini, `sd(0)` = exineq) %>%
  mutate(df = abs(year - 2015)) %>% 
  group_by(name) %>% 
  filter(df == min(df)) %>%
  filter(year == max(year)) %>% 
  ungroup() %>%
  mutate(name = str_sub(name, 1, 3)) %>%
  ggplot(aes(x = GINI, y = `sd(0)`)) + 
  geom_jitter() + 
  stat_regline_equation(label.y = 35) +
  stat_cor(aes(label  = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")),
           r.accuracy = 0.01,
           p.accuracy = 0.001,
           label.y    = 33,
           size       = 4) +
  scale_x_continuous(breaks = pretty_breaks(8)) + 
  scale_y_continuous( breaks = pretty_breaks(8)) +
  theme_bw() +
  expand_limits(y = 0, x = 0)+ 
  geom_smooth(method  = "lm", 
              se      = FALSE, 
              color   = "red",
              formula = y ~ x) + 
  theme(legend.position   = "none",
        axis.text         = element_text(size = 10, colour = "black"), 
        axis.title        = element_text(size = 10, colour = "black", face = "bold"), 
        strip.background  = element_blank()) + 
  geom_label_repel(aes(label = name),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = "black",
                   max.overlaps  = 60)