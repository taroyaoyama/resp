
oldws <- 'C:/Users/taroy/OneDrive - The University of Tokyo/ws/Wave/'
list.files(oldws)

file.copy(
    paste0(oldws,c('NewmarkBeta.cpp', 'NewmarkBeta_NL.cpp')),
    to = 'src/')

devtools::document()
devtools::load_all()

## test NewmarkBeta
h <- 0.05
w <- 2*pi
dat <- rnorm(1024)
test <- NewmarkBeta(dat, M = as.matrix(1), C = as.matrix(2*h*w), K = as.matrix(w^2), dt = 0.01)
