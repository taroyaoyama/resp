
Nigam <- function(y, h, w, dt) {

    wd  <-  w * sqrt(1 - h^2)
    E   <-  exp(-h * w * dt)
    SS  <- -h * w * sin(wd * dt) - wd * cos(wd * dt)
    CC  <- -h * w * cos(wd * dt) + wd * sin(wd * dt)
    S1  <- (E * SS + wd) / w^2
    C1  <- (E * CC + h * w) / w^2
    S2  <- (E * dt * SS + h * w * S1 + wd * C1) / w^2
    C2  <- (E * dt * CC + h * w * C1 - wd * S1) / w^2
    S3  <- dt * S1 - S2
    C3  <- dt * C1 - C2

    A12 <-  E * sin(wd * dt) / wd
    A21 <- -E * sin(wd * dt) * w^2 / wd
    A11 <-  E * (cos(wd * dt) + h * w * sin(wd * dt) / wd)
    A22 <-  E * (cos(wd * dt) - h * w * sin(wd * dt) / wd)
    B11 <- -S2 / wd / dt
    B12 <- -S3 / wd / dt
    B21 <- (h * w * S2 - wd * C2) / wd / dt
    B22 <- (h * w * S3 - wd * C3) / wd / dt

    len <- length(y)
    dis <- rep(0, len)
    vel <- rep(0, len)
    acc <- rep(0, len)

    dis[1] <- 0
    vel[1] <- -y[1] * dt
    acc[1] <- 2.0 * h * w * y[1] * dt

    for (i in 1:(len-1)) {
        dis[i+1] <- A11 * dis[i] + A12 * vel[i] + B11 * y[i] + B12 * y[i]
        vel[i+1] <- A21 * dis[i] + A22 * vel[i] + B21 * y[i] + B22 * y[i]
        acc[i+1] <- -2.0 * h * w * vel[i+1] - w^2 * dis[i+1]
    }

    tim <- seq(0.0, by = dt, length = len)
    mat <- cbind(as.matrix(tim), as.matrix(y), as.matrix(dis), as.matrix(vel), as.matrix(acc))
    colnames(mat) <- c('time', 'input.acceleration', 'rel.displacement', 'rel.velocity', 'abs.acceleration')
    return(mat)
}
