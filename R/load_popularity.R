# generate incomplete data based on the popularity dataset

# # load complete data
# library(lme4)
# pop_complete <-
#   foreign::read.spss("../Data/popular2.sav", to.data.frame = TRUE)[, 2:7]
# 
# mod <- lmer(popular ~ sex + extrav + texp + extrav * texp + (1 | class) + (1 | texp) +  (1 | extrav), pop_complete) %>% broom.mixed::tidy()
# 
# # parameters
# N = 5
# n = 25

# # synthesize some more observations
# icc(popular ~as.factor(class), pop_complete)
# 
# pop_synth <- pop_complete %>%
#   split(.$class) %>%
#   .[1:3] %>%
#   purrr::map_dfr(., ~{
#   mice::mice(
#     .,
#     m = 19,
#     maxit = 1,
#     method = "pmm",
#     where = cbind(FALSE, matrix(TRUE, nrow(.), ncol(.)-1))
#   ) %>%
#   mice::complete("long", include = TRUE)})

# maybe just use this https://debruine.github.io/tutorials/sim-lmer.html

