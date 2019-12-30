
library(devtools)
library(tidyverse)
library(resp)
load_all()

## respec test -----

# input
dat <- read_csv('data/EL-CENTRO.csv') %>% select(2)

# params
# fs <- seq(0.01, 10, by = 0.01)
fs <- seq(0.1, 10, by = 0.1)
h <- 0.05

# slots
dis <- rep(0, length(fs))
vel <- rep(0, length(fs))
acc <- rep(0, length(fs))

# naive -----
df_naive <- tibble(dis, vel, acc)

t0 <- proc.time()
for (i in 1:length(fs)) {

    # newmark beta
    nb <- newmark_beta(
        dat = as.matrix(dat), dt = 0.02, M = as.matrix(1),
        C = as.matrix(4*pi*h*fs[i]), K = as.matrix(4*pi^2*fs[i]^2))

    df_naive$dis[i] <- max(abs(nb$relDis))
    df_naive$vel[i] <- max(abs(nb$relVel))
    df_naive$acc[i] <- max(abs(nb$absAcc))
}
print(proc.time() - t0)

# using purrr
nb_calc <- function (frq) {

    nb <- newmark_beta(
        dat = as.matrix(dat), dt = 0.02, M = as.matrix(1),
        C = as.matrix(4*pi*h*frq), K = as.matrix(4*pi^2*frq^2)) %>%
        as_tibble %>%
        select(3, 5, 6) %>%
        abs %>%
        summarise_all(max)

    return(nb)
}
system.time(
    df_purrr <- fs %>% as.list() %>% map_df(nb_calc)
)

# using furrr
library(furrr)
plan(multiprocess)
system.time(
    df_furrr <- fs %>% as.list() %>% furrr::future_map(nb_calc) %>% bind_rows()
)
