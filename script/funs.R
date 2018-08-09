## Useful functions

## run mrgsolve simulation
run_sim <- function(mod, cmt, dose, rate=0){
  mod %>% 
    ev(cmt = cmt, amt = dose, ii = 12, addl = 13, ss = 1, rate = rate) %>% 
    mrgsim(delta = 0.1, end = 12) %>%
    as.data.frame() %>%
    dplyr::filter(row_number() > 1)
}