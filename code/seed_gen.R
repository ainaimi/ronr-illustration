pacman::p_load(
  rio,          
  here,         
  skimr,        
  tidyverse,     
  lmtest,
  sandwich,
  broom
)

thm <- theme_classic() +
  theme(
    legend.position = "top",
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )
theme_set(thm)

set.seed(321, kind = "L'Ecuyer-CMRG")

seed_base <- 1:10e6

length(unique(seed_base))

random_seeds <- sample(seed_base, size = 1.5e6, replace = F)

length(unique(random_seeds))

.Random.seed

# first five seeds:
random_seeds[1:5]

plot(ecdf(random_seeds))

hist(random_seeds)

write_csv(data.frame(random_seeds), here("data","random_seed_values.csv"))


# a <- read_csv(here("data", "random_seed_values.csv"))
# 
# min(a[1:5000,])
# max(a[1:5000,])
# median(a[1:5000,]$random_seeds)