# function for simulating data from different scenarios
gen.noC <- function(n, multi = 100, Tmodel, Qmodel, T.min, Q.min, Q.max, tau,
                        beta0.T, beta0.Q, beta.T = NA, beta.Q = NA, beta.T2 = NA,
                        beta.Q2 = NA,
                        shape.T = NA, shape.Q = NA, epsT.sd = NA, epsQ.sd = NA){

    Z1 = runif(multi*n, -1, 1)
    Z2 = rbinom(multi*n, size = 1, prob = 0.5) - 0.5
    Z = cbind(Z1 = Z1, Z2 = Z2)
    Zsq = cbind(Z1 = Z1, Z2 = Z2, Z1_sq = Z1^2-1/3, Z1Z2 = Z1*Z2)
    # Zsq = cbind(Z1 = Z1, Z2 = Z2, Z1_sq = sqrt(abs(Z1))-2/3, Z1Z2 = Z1*Z2)

    id = which(Z2<0)
    nn1 = length(id)
    nn2 = multi*n - nn1
    TT = rep(NA, multi*n)
    Q2 = rep(NA, multi*n)

    if(Tmodel == 'weibull1'){
        TT = T.min + rweibull(multi*n, shape = shape.T, scale = exp(-1/shape.T * (beta0.T + Z%*%beta.T)))
    }else if(Tmodel == 'weibull2'){
        TT = T.min + rweibull(multi*n, shape = shape.T, scale = exp(-1/shape.T * (beta0.T + Zsq%*%beta.T2)))
    }else if(Tmodel == 'AFT2_weibull2'){
        TT[id] = T.min + exp(beta0.T + (Zsq%*%beta.T2)[id] + rnorm(nn1, mean = 0, sd = epsT.sd))
        TT[-id] = T.min + rweibull(nn2, shape = shape.T, scale = exp(-1/shape.T * (beta0.T + (Zsq%*%beta.T2)[-id])))
    }


    Q2.max = tau - Q.min
    Q2.min = tau - Q.max
    U2 = runif(multi*n, min = 0, max = 1)
    if(Qmodel == "Cox1"){
        Q2 = Q2.min + as.vector((Q2.max-Q2.min)*(1-U2^(exp(-Z %*% beta.Q))))
    }else if(Qmodel == "Cox2"){
        Q2 = Q2.min + as.vector((Q2.max-Q2.min)*(1-U2^(exp(-Zsq %*% beta.Q2))))
    }
    QQ = tau - Q2


    dat.full = data.frame(time = TT, Q = QQ, Z)
    dat.all = dat.full
    obs.id = which(dat.full$Q < dat.full$time)
    dat.obs = dat.full[obs.id,]
    if(length(obs.id)<n){
        stop('Truncation rate is high. Need to increase the mutiplier.')
    }
    dat = dat.obs[1:n,]
    dat.full = dat.full[(1:obs.id[n]), ]

    return(list(dat = dat, dat.full = dat.full, dat.all = dat.all))
}


# function for simulating data under left truncation and the "c2" type of right censoring (censoring always after truncation)
gen.c2 <- function(n, multi = 100, Tmodel, Qmodel, T.min, Q.min, Q.max, tau,
                  beta0.T, beta0.Q, beta.T = NA, beta.Q = NA, beta.T2 = NA, beta.Q2 = NA,
                  shape.T = NA, shape.Q = NA, epsT.sd = NA, epsQ.sd = NA,
                  C.min, shape.C, scale.C){

    Z1 = runif(multi*n, -1, 1)
    Z2 = rbinom(multi*n, size = 1, prob = 0.5) - 0.5
    Z = cbind(Z1 = Z1, Z2 = Z2)
    Zsq = cbind(Z1 = Z1, Z2 = Z2, Z1_sq = Z1^2-1/3, Z1Z2 = Z1*Z2)

    id = which(Z2<0)
    nn1 = length(id)
    nn2 = multi*n - nn1
    TT = rep(NA, multi*n)
    Q2 = rep(NA, multi*n)

    if(Tmodel == 'weibull1'){
        TT = T.min + rweibull(multi*n, shape = shape.T, scale = exp(-1/shape.T * (beta0.T + Z%*%beta.T)))
    }else if(Tmodel == 'weibull2'){
        TT = T.min + rweibull(multi*n, shape = shape.T, scale = exp(-1/shape.T * (beta0.T + Zsq%*%beta.T2)))
    }else if(Tmodel == 'AFT2_weibull2'){
        TT[id] = T.min + exp(beta0.T + (Zsq%*%beta.T2)[id] + rnorm(nn1, mean = 0, sd = epsT.sd))
        TT[-id] = T.min + rweibull(nn2, shape = shape.T, scale = exp(-1/shape.T * (beta0.T + (Zsq%*%beta.T2)[-id])))
    }

    Q2.max = tau - Q.min
    Q2.min = tau - Q.max
    U2 = runif(multi*n, min = 0, max = 1)
    if(Qmodel == "Cox1"){
        Q2 = Q2.min + as.vector((Q2.max-Q2.min)*(1-U2^(exp(-Z %*% beta.Q))))
    }else if(Qmodel == "Cox2"){
        Q2 = Q2.min + as.vector((Q2.max-Q2.min)*(1-U2^(exp(-Zsq %*% beta.Q2))))
    }
    QQ = tau - Q2

    # simulate the censoring variable
    D = C.min + rweibull(multi*n, shape = shape.C, scale = scale.C)
    C = QQ+D

    X = pmin(TT, C)
    delta = as.numeric(TT<C)

    dat.full = data.frame(TT = TT, X = X, Q = QQ, C = C, delta = delta, Z)
    dat.all = dat.full[,-1]
    obs.id = which(dat.all$Q < dat.all$X)
    dat.obs = dat.all[obs.id,]
    if(length(obs.id)<n){
        stop('Truncation rate is high. Need to increase the mutiplier.')
    }
    dat = dat.obs[1:n,]
    dat.full = dat.full[(1:obs.id[n]), ]

    return(list(dat = dat, dat.full = dat.full, dat.all = dat.all))
}



# function for simulating data under left truncation and the "c1" type of right censoring (censoring can be before truncation)
gen.c1_Xcox <- function(n, multi = 100, Tmodel, Qmodel, T.min, Q.min, Q.max, tau,
                beta0.T, beta0.Q, beta.T = NA, beta.Q = NA, beta.T2 = NA, beta.Q2 = NA,
                shape.T = NA, shape.Q = NA, epsT.sd = NA, epsQ.sd = NA,
                C.min, shape.C, scale.C, type.C = "c1"){

    Z1 = runif(multi*n, -1, 1)
    Z2 = rbinom(multi*n, size = 1, prob = 0.5) - 0.5
    Z = cbind(Z1 = Z1, Z2 = Z2)
    Zsq = cbind(Z1 = Z1, Z2 = Z2, Z1_sq = Z1^2-1/3, Z1Z2 = Z1*Z2)

    id = which(Z2<0)
    nn1 = length(id)
    nn2 = multi*n - nn1
    TT = rep(NA, multi*n)
    Q2 = rep(NA, multi*n)

    if(T.min != C.min){
        stop("T.min not equal to C.min")
    }
    if(shape.T != shape.C){
        stop("shape.T not equal to shape.C")
    }
    kk = shape.T

    # simulate the censoring variable
    if(type.C == "c1"){
        C = C.min + rweibull(multi*n, shape = shape.C, scale = scale.C)
    }else {
        stop("This can only generate data for type.C = 'c1', i.e., truncation before censoring scenario!")
    }

    # simulate T
    if(Tmodel == 'weibull1'){
        TT = T.min + rweibull(multi*n, shape = kk, scale = ( exp(beta0.T+Z%*%beta.T) - scale.C^{-kk} )^{-1/kk} )
    }else{
        stop("This function can only simulate data for Tmodel = 'weibull1'!")
    }

    # simulate Q
    Q2.max = tau - Q.min
    Q2.min = tau - Q.max
    U2 = runif(multi*n, min = 0, max = 1)
    if(Qmodel == "Cox1"){
        Q2 = Q2.min + as.vector((Q2.max-Q2.min)*(1-U2^(exp(-Z %*% beta.Q))))
    }else if(Qmodel == "Cox2"){
        Q2 = Q2.min + as.vector((Q2.max-Q2.min)*(1-U2^(exp(-Zsq %*% beta.Q2))))
    }
    QQ = tau - Q2


    X = pmin(TT, C)
    delta = as.numeric(TT<C)

    dat.full = data.frame(TT = TT, X = X, Q = QQ, C = C, delta = delta, Z)
    dat.all = dat.full[,-1]
    obs.id = which(dat.all$Q < dat.all$X)
    dat.obs = dat.all[obs.id,]
    if(length(obs.id)<n){
        stop('Truncation rate is high. Need to increase the mutiplier.')
    }
    dat = dat.obs[1:n,]
    dat.full = dat.full[(1:obs.id[n]), ]

    return(list(dat = dat, dat.full = dat.full, dat.all = dat.all))
}








