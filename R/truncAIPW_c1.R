### function that compute the truncAIPW estimator under "c1" type of right censoring - censoring always after left truncation

#' @title Doubly Robust Estimation under Covariate-induced Dependent Left Truncation and Noninformative Right Censoring where Censoring can be before Left Truncation
#' @description Doubly robust estimation of the mean of an arbitrarily transformed survival time under covariate-induced dependent left truncation and noninformative right censoring where censoring can be before left truncation. Inverse probability of censoring weighting is used to handle the right censoring.
#' @references Wang, Y., Ying, A., Xu, R. (2022) "Doubly robust estimation under covariate-induced dependent left truncation" <arXiv:2208.06836>.
#' @param dat data frame that contains the data for constructing the estimating equation.
#' @param nu transformation that defines the parameter of interest.
#' @param Fuz.mx matrix for the estimated conditional CDF of the event time given covariates. Each row corresponds to a subject, and each column corresponds to a time point. The column names of the matrix are the time points. See \code{\link{F_est}} for an example of computing this input matrix for the conditional CDF.
#' @param Gvz.mx matrix for the estimated conditional CDF of the truncation time given covariates. Each row corresponds to a subject, and each column corresponds to a time point. The column names of the matrix are the time points. See \code{\link{G_est}} for an example of computing this input matrix for the conditional CDF.
#' @param Sc a function for the censoring survival curve \eqn{S_c(\cdot)}.
#' @param X.name name of the censored event time variable X = min(T, C).
#' @param Q.name name of the left truncation time variable.
#' @param status.name name of the event time indicator.
#' @param trim constant that is used to bound from below for the denominators involved in the computation.
#' @return \code{truncAIPW_cen1()} returns a list of estimators (`dr', `IPW.Q', `Reg.T1', `Reg.T2').
#' \item{dr}{doubly robust estimator `dr'.}
#' \item{IPW.Q}{inverse probability of truncation weighted estimator `IPW.Q'.}
#' \item{Reg.T1}{regression based estimator `Reg.T1'.}
#' \item{Reg.T2}{regression based estimator `Reg.T2'.}
#' @export
#' @import survival stats
#' @seealso See also \code{\link{truncAIPW}} for estimation under no censoring, and \code{\link{truncAIPW_cen2}} for estimation under another type of noninformative right censoring. See also \code{\link{F_est}}, \code{\link{G_est}} as examples for computing the input matrices of the conditional CDF's.
#' @examples
#' library(survival)
#' data("simu_c1")
#' simu_c1$delta.1 = 1
#'
#' nu <- function(t){ return(as.numeric(t>3)) }
#' u = c(min(simu_c1$X)-1e-10, sort(simu_c1$X), max(simu_c1$X)+1e-10)
#' v = c(min(simu_c1$Q)-1e-10, sort(simu_c1$Q), max(simu_c1$Q)+1e-10)
#'
#' Fuz.mx = F_est(simu_c1, simu_c1, u, "Cox", "X", "Q", "delta.1", c("Z1","Z2"))
#' Gvz.mx = G_est(simu_c1, simu_c1, v, "Cox", "X", "Q", "delta.1", c("Z1","Z2"))
#'
#' # KM curve for Sc
#' kmfit.C = survfit(Surv(Q, X, 1-delta)~1, data = simu_c1, type = "kaplan-meier")
#' Sc = stepfun(kmfit.C$time,  c(1, kmfit.C$surv))
#'
#' est = truncAIPW_cen1(simu_c1, nu, Fuz.mx, Gvz.mx, Sc, "X", "Q", "delta", trim = 1e-7)
#' est
truncAIPW_cen1 <- function(dat, nu, Fuz.mx, Gvz.mx, Sc, X.name, Q.name, status.name, trim = 1e-7){

    efs = c1.truncAIPW_transMean_EF(dat, nu, Fuz.mx, Gvz.mx, Sc,
                                    X.name, Q.name, status.name, trim)

    Num_AIPW = efs$Num_AIPW
    Den_AIPW = efs$Den_AIPW

    Num_IPW.Q = efs$Num_IPW.Q
    Den_IPW.Q = efs$Den_IPW.Q

    Num_Reg.T1 = efs$Num_Reg.T1
    Num_Reg.T2 = efs$Num_Reg.T2
    Den_Reg.T1 = efs$Den_Reg.T1
    Den_Reg.T2 = efs$Den_Reg.T2

    est_dr = mean(Num_AIPW)/mean(Den_AIPW)
    est_IPWQ = mean(Num_IPW.Q)/mean(Den_IPW.Q)
    est_RegT1 = mean(Num_Reg.T1)/mean(Den_Reg.T1)
    est_RegT2 = mean(Num_Reg.T2)/mean(Den_Reg.T2)


    return(list(dr = est_dr,
                IPW.Q = est_IPWQ,
                Reg.T1 = est_RegT1,
                Reg.T2 = est_RegT2))
}







## function that compute different parts involved in the estimating function (EF) under the "c1" type right censoring
# Fuz.mx: the matrix representing the CDF of T|(A,Z) for each subject in dat
# Gvz.mx: the matrix representing the CDF of Q|(A,Z) for each subject in dat
# Sc: the function for the censoring survival probability S_c()

c1.truncAIPW_transMean_EF <- function(dat, nu, Fuz.mx, Gvz.mx, Sc,
                                      X.name, Q.name, status.name,
                                      trim = 1e-7){

    if(nrow(dat) != nrow(Fuz.mx) | nrow(dat) != nrow(Gvz.mx)){
        stop("The number of rows of dat, Fuz.mx and Gvz.mx are not the same. ")
    }

    X = dat[, X.name]
    Q = dat[, Q.name]
    delta = dat[, status.name]

    u = as.numeric(colnames(Fuz.mx))  # jumps.X
    v = as.numeric(colnames(Gvz.mx))  # jumps.Q

    tau2 = max(v)+1e-10

    Gtz = CDF_eval(X, Gvz.mx)
    Gqz = CDF_eval(Q, Gvz.mx)
    Fqz = CDF_eval(Q, Fuz.mx)
    Scx = Sc(X)

    DDen1 = delta/(Scx*pmax(Gtz, trim))

    delta.mx = matrix(rep(delta, length(u)), ncol = length(u), byrow = FALSE)
    Scu.mx = matrix(rep(Sc(u), nrow(dat)), nrow = nrow(dat), byrow = TRUE)

    delta_Sc.mx = delta.mx / Scu.mx
    mqz.nu1_Sc = diag(int_fmx_dF(Q, delta_Sc.mx, Fuz.mx))
    DDen2 = mqz.nu1_Sc/(pmax(Gqz, trim)*pmax(1-Fqz, trim))

    Fvz.mx = CDF_eval.mx(v, Fuz.mx)

    nn = nrow(dat)
    atrisk.mx = matrix(nrow = nn, ncol = length(v))
    for(i in 1:nn){
        atrisk.mx[i,] = (Q[i] <= v & v < X[i])    # at risk indicator for subject i at times in jumps.Q
    }
    mvz.mx.nu1_Sc = int_fmx_dF(v, delta_Sc.mx, Fuz.mx)
    f.mx = atrisk.mx*mvz.mx.nu1_Sc/(pmax(Gvz.mx, trim)^2*pmax(1-Fvz.mx, trim))
    DDen3 = as.vector(int_fmx_dF(tau2, f.mx, Gvz.mx))


    NNum1 = delta*nu(X)/(Scx*pmax(Gtz, trim))

    nuu.mx = matrix(rep(nu(u), nrow(dat)), nrow = nrow(dat), byrow = TRUE)
    mqz.Sc = diag(int_fmx_dF(Q, nuu.mx*delta_Sc.mx, Fuz.mx))
    NNum2 = mqz.Sc/(pmax(Gqz, trim)*pmax(1-Fqz, trim))

    mvz.mx.Sc = int_fmx_dF(v, nuu.mx*delta_Sc.mx, Fuz.mx)
    fnu.mx.Sc = atrisk.mx*mvz.mx.Sc/(pmax(Gvz.mx, trim)^2*pmax(1-Fvz.mx, trim))
    NNum3 = as.vector(int_fmx_dF(tau2, fnu.mx.Sc, Gvz.mx))

    Num_AIPW = NNum1 + NNum2 - NNum3
    Den_AIPW = DDen1 + DDen2 - DDen3


    # For the estimating equation of the IPW and Reg.T1, Reg.T2 estimators
    DDen4 = delta/Scx + mqz.nu1_Sc/pmax(1-Fqz, trim)
    NNum4 = nu(X)*delta/Scx + mqz.Sc/pmax(1-Fqz, trim)

    tau.Tmax = max(u) + 1
    Enutz_Sc = as.vector(int_fmx_dF(tau.Tmax, nuu.mx*delta_Sc.mx, Fuz.mx))
    NNum5 = Enutz_Sc/pmax(1-Fqz, trim)

    Enutz_nu1_Sc = as.vector(int_fmx_dF(tau.Tmax, delta_Sc.mx, Fuz.mx))
    DDen5 = Enutz_nu1_Sc/pmax(1-Fqz, trim)


    return(list(Num_AIPW = Num_AIPW, Den_AIPW = Den_AIPW,
                Num_IPW.Q = NNum1, Den_IPW.Q = DDen1,
                Num_Reg.T1 = NNum4, Num_Reg.T2 = NNum5,
                Den_Reg.T1 = DDen4, Den_Reg.T2 = DDen5))

}
