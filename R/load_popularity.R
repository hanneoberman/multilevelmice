# generate incomplete data based on the popularity dataset

# load complete data
pop_complete <-
  read.spss("Data/popular2.sav", to.data.frame = TRUE)[, 2:7]

icc(popular ~as.factor(class), pop_complete)

pop_synth <- pop_complete %>% 
  split(.$class) %>% 
  .[1:3] %>% 
  purrr::map_dfr(., ~{
  mice::mice(
    .,
    m = 19,
    maxit = 1,
    method = "pmm",
    where = cbind(FALSE, matrix(TRUE, nrow(.), ncol(.)-1))
  ) %>%
  mice::complete("long", include = TRUE)})

