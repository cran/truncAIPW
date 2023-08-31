## functions that apply the AIPW score for left truncation to a full data estimating function


### function that compute the truncAIPW estimator
#' @title Doubly Robust Estimation under Covariate-induced Dependent Left Truncation and No Censoring
#' @description Doubly robust estimation for the mean of an arbitrarily transformed survival time under covariate-induced dependent left truncation and no right censoring.
#' @references Wang, Y., Ying, A., Xu, R. (2022) "Doubly robust estimation under covariate-induced dependent left truncation" <arXiv:2208.06836>.
#' @param dat data frame that contains the data for constructing the estimating equation.
#' @param nu transformation that defines the parameter of interest.
#' @param Fuz.mx matrix for the estimated conditional CDF of the event time given covariates. Each row corresponds to a subject, and each column corresponds to a time point. The column names of the matrix are the time points. See \code{\link{F_est}} for an example of computing this conditional CDF matrix.
#' @param Gvz.mx matrix for the estimated conditional CDF of the truncation time given covariates. Each row corresponds to a subject, and each column corresponds to a time point. The column names of the matrix are the time points. See \code{\link{G_est}} for an example of computing this conditional CDF matrix.
#' @param T.name name of the event time variable.
#' @param Q.name name of the left truncation time variable.
#' @param trim constant that is used to bound from below for the denominators involved in the computation.
#' @return \code{truncAIPW()} returns a list of estimators (`dr', `IPW.Q', `Reg.T1', `Reg.T2'), and the model-based standard errors for the `dr' and `IPW.Q' estimators.
#' \item{dr}{doubly robust estimator `dr'.}
#' \item{IPW.Q}{inverse probability of truncation weighted estimator `IPW.Q'.}
#' \item{Reg.T1}{regression based estimator `Reg.T1'.}
#' \item{Reg.T2}{regression based estimator `Reg.T2'.}
#' \item{SE_dr}{standard error of the `dr' estimator based on the efficient influence function.}
#' \item{SE_IPW.Q}{standard error of the `IPW.Q' estimator computed from the robust sandwich variance estimator assuming the truncation weights are known.}
#' @export
#' @import survival stats
#' @seealso See \code{\link{truncAIPW_cen1}}, \code{\link{truncAIPW_cen2}} for the estimations also under noninformative right censoring. See \code{\link{F_est}}, \code{\link{G_est}} for examples of computing the input matrices for the conditional CDF's.
#' @examples
#' data("simu")
#' nu <- function(t){ return(as.numeric(t>3)) }
#' u = c(min(simu$time)-1e-10, sort(simu$time), max(simu$time)+1e-10)
#' v = c(min(simu$Q)-1e-10, sort(simu$Q), max(simu$Q)+1e-10)
#' Fuz.mx = F_est(simu, simu, u, "Cox", "time", "Q", "delta", c("Z1","Z2"))
#' Gvz.mx = G_est(simu, simu, v, "Cox", "time", "Q", "delta", c("Z1","Z2"))
#'
#' est = truncAIPW(simu, nu, Fuz.mx, Gvz.mx, "time", "Q", trim = 1e-7)
#' est
truncAIPW <- function(dat, nu, Fuz.mx, Gvz.mx, T.name, Q.name, trim = 1e-7){

    efs = truncAIPW_transMean_EF(dat, nu, Fuz.mx, Gvz.mx, T.name, Q.name, trim)

    Num_AIPW = efs$Num_AIPW
    Den_AIPW = efs$Den_AIPW

    Num_IPW.Q = efs$Num_IPW.Q
    Den_IPW.Q = efs$Den_IPW.Q

    Num_Reg.T1 = efs$Num_Reg.T1
    Num_Reg.T2 = efs$Num_Reg.T2
    Den_Reg = efs$Den_Reg


    est_dr = mean(Num_AIPW)/mean(Den_AIPW)
    est_IPWQ = mean(Num_IPW.Q)/mean(Den_IPW.Q)
    est_RegT1 = mean(Num_Reg.T1)/mean(Den_Reg)
    est_RegT2 = mean(Num_Reg.T2)/mean(Den_Reg)


    n = nrow(dat)

    beta = 1/mean(Den_IPW.Q)
    IF_dr = beta*(Num_AIPW - est_dr*Den_AIPW)
    se_dr = sqrt(mean(IF_dr^2))/sqrt(n)

    SF_IPW.Q = Num_IPW.Q - est_IPWQ*Den_IPW.Q
    se_IPW.Q = sqrt(mean(SF_IPW.Q^2)/(mean(Den_IPW.Q))^2)/sqrt(n)


    return(list(dr = est_dr,
                IPW.Q = est_IPWQ,
                Reg.T1 = est_RegT1,
                Reg.T2 = est_RegT2,
                SE_dr = se_dr,
                SE_IPW.Q = se_IPW.Q))
}








## function that compute different parts involved in the estimating function (EF)
# Fuz.mx: the matrix representing the CDF of T|(A,Z) for each subject in dat
# Gvz.mx: the matrix representing the CDF of Q|(A,Z) for each subject in dat

truncAIPW_transMean_EF <- function(dat, nu, Fuz.mx, Gvz.mx,
                                      T.name, Q.name, trim = 1e-7){

    if(nrow(dat) != nrow(Fuz.mx) | nrow(dat) != nrow(Gvz.mx)){
        stop("The number of rows of dat, Fuz.mx and Gvz.mx are not the same. ")
    }

    time = dat[, T.name]
    Q = dat[, Q.name]

    u = as.numeric(colnames(Fuz.mx))  # jumps.T
    v = as.numeric(colnames(Gvz.mx))  # jumps.Q

    tau2 = max(v)+1e-10

    Gtz = CDF_eval(time, Gvz.mx)
    Gqz = CDF_eval(Q, Gvz.mx)
    Fqz = CDF_eval(Q, Fuz.mx)

    DDen1 = 1/pmax(Gtz, trim)
    DDen2 = Fqz/(pmax(Gqz, trim)*pmax(1-Fqz, trim))

    Fvz.mx = CDF_eval.mx(v, Fuz.mx)

    nn = nrow(dat)
    atrisk.mx = matrix(nrow = nn, ncol = length(v))
    for(i in 1:nn){
        atrisk.mx[i,] = (Q[i] <= v & v < time[i])    # at risk indicator for subject i at times in jumps.Q
    }
    f.mx = atrisk.mx*Fvz.mx/(pmax(Gvz.mx, trim)^2*pmax(1-Fvz.mx, trim))
    DDen3 = as.vector(int_fmx_dF(tau2, f.mx, Gvz.mx))


    NNum1 = nu(time)/pmax(Gtz, trim)

    nuu.mx = matrix(rep(nu(u), nrow(dat)), nrow = nrow(dat), byrow = TRUE)
    mqz = diag(int_fmx_dF(Q, nuu.mx, Fuz.mx))
    NNum2 = mqz/(pmax(Gqz, trim)*pmax(1-Fqz, trim))

    mvz.mx = int_fmx_dF(v, nuu.mx, Fuz.mx)
    fnu.mx = atrisk.mx*mvz.mx/(pmax(Gvz.mx, trim)^2*pmax(1-Fvz.mx, trim))
    NNum3 = as.vector(int_fmx_dF(tau2, fnu.mx, Gvz.mx))

    Num_AIPW = NNum1 + NNum2 - NNum3
    Den_AIPW = DDen1 + DDen2 - DDen3


    # For the estimating equation of the IPW and Reg.T1, Reg.T2 estimators
    DDen4 = 1/pmax(1-Fqz, trim)
    NNum4 = nu(time) + mqz/pmax(1-Fqz, trim)

    tau.Tmax = max(u) + 1
    Enutz = as.vector(int_fmx_dF(tau.Tmax, nuu.mx, Fuz.mx))
    NNum5 = Enutz/pmax(1-Fqz, trim)



    return(list(Num_AIPW = Num_AIPW, Den_AIPW = Den_AIPW,
                Num_IPW.Q = NNum1, Den_IPW.Q = DDen1,
                Num_Reg.T1 = NNum4, Num_Reg.T2 = NNum5, Den_Reg = DDen4))

}





# function that computes the AIPW weights for PS estimation
truncAIPW_weights <- function(dat, T.name, Q.name, Fuz.mx, Gvz.mx, trim = 0){

    time = dat[, T.name]
    Q = dat[, Q.name]

    # jumps.T = as.numeric(colnames(Fuz.mx))  # jumps.T
    v = as.numeric(colnames(Gvz.mx))  # jumps.Q
    tau2 = max(v)+1e-10

    Gtz = CDF_eval(time, Gvz.mx)
    Gqz = CDF_eval(Q, Gvz.mx)
    Fqz = CDF_eval(Q, Fuz.mx)

    W1 = 1/pmax(Gtz, trim)
    W2 = Fqz/(pmax(Gqz, trim)*pmax(1-Fqz, trim))

    Fvz.mx = CDF_eval.mx(v, Fuz.mx)
    nn = nrow(dat)
    atrisk.mx = matrix(nrow = nn, ncol = length(v))
    for(i in 1:nn){
        atrisk.mx[i,] = (Q[i] <= v & v < time[i])    # at risk indicator for subject i at times in jumps.Q
    }
    f.mx = atrisk.mx*Fvz.mx/(pmax(Gvz.mx, trim)^2*pmax(1-Fvz.mx, trim))
    W3 = as.vector(int_fmx_dF(tau2, f.mx, Gvz.mx))

    w_AIPW = W1 + W2 - W3
    w_IPW.Q = W1
    w_IPW.T = 1/pmax(1-Fqz, trim)

    return(list(w_AIPW = w_AIPW, w_IPW.Q = w_IPW.Q, w_IPW.T = w_IPW.T))
}


