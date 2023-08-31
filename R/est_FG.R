
### function that output the conditional CDF matrices for F and G at given time points ----------------
#' @title Estimate the Conditional CDF of the Event Time given Covariates
#' @description Estimate the conditional cumulative distribution function (CDF) of the event time given covariates evaluated at given time points. The options implemented in this function are: Cox proportional hazards regression using function \code{coxph()} from  R package `survival', and the hazard model with penalized splines using function \code{survPen()} from R package `survPen'.
#' @param dat.fit data frame that is used to fit the model for the full data conditional distribution of the event time given the covariates.
#' @param dat.est data frame that contains the subjects for which the estimated conditional CDF is computed.
#' @param time.eval vector of time points at which the conditional CDF is evaluated.
#' @param model method used to estimate the conditional CDF. The options available are "Cox" and "spline", corresponding to Cox proportional hazards regression using function \code{coxph()} from  R package `survival', and the hazard model with penalized splines using function \code{survPen()} from R package `survPen', respectively.
#' @param time.name name of the event time variable.
#' @param Q.name name of the left truncation time variable.
#' @param event.name name of the event indicator.
#' @param cov.names vector of the names of covariates.
#' @param trim constant for bounding the estimated conditional CDF from 1.
#' @param formula.survPen the formula when applying the hazard model with penalized splines implemented in \code{survPen::survPen}.
#' @return \code{F_est()} returns a matrix of the estimated conditional CDF for subjects in `\code{data.est}' evaluated at the time points in the vector `\code{time.eval}'. Each row corresponds to a subject and each column corresponds to a time point. The column names of the matrix are the times in `\code{time.eval}'.
#' @export
#' @import survival stats survPen
#' @seealso  \code{\link{G_est}}
#' @examples
#' data("simu")
#' u = c(1, 1.5, 2, 2.5, 3, 3.5, 4)
#' Fuz.mx = F_est(simu, simu[1:10,], u, "Cox", "time", "Q", "delta", c("Z1","Z2"))
F_est <- function(dat.fit, dat.est = dat.fit, time.eval, model,
                 time.name, Q.name, event.name, cov.names, trim = 0,
                 formula.survPen = NA){

    u = time.eval

    # if(model == "RF"){
    #     F_hat = F_hat.ltrcrrf(dat.fit, time.name, Q.name, event.name, cov.names,
    #                           mtry, ntree)
    #     Fuz.mx = pmin(F_hat(newdata = dat.est, time.eval = u), 1-trim)   #***
    #
    # }else
    if(model == "Cox"){
        id.cox = ((dat.fit[,time.name] - dat.fit[,Q.name]) >= 10^{-7})
        dat.cox = dat.fit[id.cox, ]
        formula.T = formula(paste("Surv(", Q.name, ",", time.name, ",",
                                  event.name, ") ~ ",
                                  paste(cov.names, collapse = "+"), collapse = ""))
        # fit Cox-PH model for T|Z
        fit.T = coxph(formula.T, data = dat.cox)
        basehaz.T = basehaz(fit.T, centered = FALSE)
        beta.T = coef(fit.T)

        Z.T = as.matrix(dat.est[,cov.names])
        Fuz.mx = pmin(Fz(u, Z.T, basehaz.T, beta.T), 1-trim)

    }else if(model == "spline"){
        time = dat.fit[ ,time.name]
        Q = dat.fit[ ,Q.name]
        delta.T = dat.fit[ ,event.name]
        n.u = length(u)
        n.est = nrow(dat.est)

        mod.T = survPen(formula.survPen, data = dat.fit, t0 = Q, t1 = time, event = delta.T)

        cov.mx = dat.est[,cov.names]
        cov.mx.rep = matrix(rep(t(cov.mx), n.u), ncol = ncol(cov.mx), byrow = TRUE)
        TZ.mx = cbind(rep(u, each = n.est), cov.mx.rep)
        colnames(TZ.mx) <- c(time.name, cov.names)
        dat.est.u = as.data.frame(TZ.mx)

        Fuz.mx = 1 - matrix(predict(mod.T, dat.est.u)$surv, ncol = n.u, byrow = FALSE)
        Fuz.mx = pmin(Fuz.mx, 1-trim)
        colnames(Fuz.mx) = u

    }else{
        stop("This T model is not implemented in this function!")
    }

    return(Fuz.mx)
}







#' @title Estimate the Conditional CDF for the Left Truncation Time given Covariates
#' @description Estimate the conditional cumulative distribution function (CDF) of the left truncation time given covariates evaluated at given time points. The options implemented in this function are: Cox proportional hazards regression using function \code{coxph()} from  R package `survival', and the hazard model with penalized splines using function \code{survPen()} from R package `survPen'.
#' @param dat.fit data frame that is used to fit the model for the full data conditional distribution of the event time given the covariates.
#' @param dat.est data frame that contains the subjects for which the estimated conditional CDF is computed.
#' @param time.eval vector of time points at which the conditional CDF is evaluated.
#' @param model method used to estimate the conditional CDF. The options available are "Cox" and "spline", corresponding to Cox proportional hazards regression using function \code{coxph()} from  R package `survival', and the hazard model with penalized splines using function \code{survPen()} from R package `survPen', respectively.
#' @param time.name name of the event time variable.
#' @param Q.name name of the left truncation time variable.
#' @param event.name name of the event indicator.
#' @param cov.names vector of the names of covariates.
#' @param trim constant for bounding the estimated conditional CDF from 0.
#' @param weights vector of case weights.
#' @param formula.survPen the formula when applying the hazard model with penalized splines implemented in \code{survPen::survPen}.
#' @return \code{G_est()} returns a matrix of the estimated conditional CDF for subjects in `\code{data.est}' evaluated at the time points in the vector `\code{time.eval}'. Each row corresponds to a subject and each column corresponds to a time point. The column names of the matrix are the times in `\code{time.eval}'.
#' @export
#' @import survival stats survPen
#' @seealso  \code{\link{F_est}}
#' @examples
#' data("simu")
#' v = c(0.5, 1, 1.5, 2, 2.5, 3)
#' Gvz.mx = G_est(simu, simu[1:10,], v, "Cox", "time", "Q", "delta", c("Z1","Z2"))
G_est <- function(dat.fit, dat.est = dat.fit, time.eval, model,
                 time.name, Q.name, event.name, cov.names, trim = 0, weights = rep(1, nrow(dat.fit)),
                 formula.survPen = NA){

    v = time.eval
    tau = max(c(dat.fit[,time.name], dat.est[,time.name])) + 1

    # if(model == "RF"){
    #
    #     if(sum(weights != 1)>0){
    #         stop("The function LTRCforests::ltrcrrf does not handle case weights.")
    #     }
    #
    #     G_hat = G_hat.ltrcrrf(dat.fit, time.name, Q.name, event.name, cov.names,
    #                           mtry, ntree, tau)
    #     Gvz.mx = pmax(G_hat(newdata = dat.est, time.eval = v), trim)   #***
    #
    # }else
    if(model == "Cox"){
        formula.Q2 = formula(paste("Surv(tau-", time.name, ", tau-", Q.name,",",
                                   event.name, ") ~ ",
                                   paste(cov.names, collapse = "+"), collapse = ""))
        id.cox = ((dat.fit[, time.name] - dat.fit[, Q.name]) >= 10^{-7})
        dat.cox = dat.fit[id.cox, ]
        # fit Cox-PH model for (tau-Q)|Z
        fit.Q2 = coxph(formula.Q2, data = dat.cox, weights = weights)

        basehaz.Q2 = basehaz(fit.Q2, centered = FALSE)
        beta.Q2 = coef(fit.Q2)

        Z.Q = as.matrix(dat.est[,cov.names])
        Gvz.mx = pmax(Gz(v, Z.Q, basehaz.Q2, beta.Q2, tau), trim)

    }else if(model == "spline"){

        if(sum(weights != 1)>0){
            stop("The function survPen::survPen does not handle case weights.")
        }

        T2 = tau- dat.fit[ ,time.name]
        Q2 = tau - dat.fit[ ,Q.name]
        delta.Q = dat.fit[ ,event.name]
        dat.fit$Q2 = Q2
        # dat.fit$T2 = T2
        # dat.fit$delta.Q = delta.Q
        n.v = length(v)
        n.est = nrow(dat.est)

        mod.Q2 = survPen(formula.survPen, data = dat.fit, t0 = T2, t1 = Q2, event = delta.Q)

        cov.mx = dat.est[,cov.names]
        cov.mx.rep = matrix(rep(t(cov.mx), n.v), ncol = ncol(cov.mx), byrow = TRUE)
        Q2Z.mx = cbind(tau - rep(v, each = n.est), cov.mx.rep)
        colnames(Q2Z.mx) <- c("Q2", cov.names)
        dat.est.v = as.data.frame(Q2Z.mx)

        Gvz.mx = matrix(predict(mod.Q2, dat.est.v)$surv, ncol = n.v, byrow = FALSE)
        Gvz.mx = pmax(Gvz.mx, trim)
        colnames(Gvz.mx) = v

    }else{
        stop("This Q model is not implemented in this function!")
    }


    return(Gvz.mx)
}




### Functions needed when using Cox models ---------------------------------------
# function for computing the baseline CDF or survival function
baseCDF.single <- function(t, basehaz){
    if(t<min(basehaz[,'time'])){
        return(0)   # check 1 or 0
    }else{
        id = which((t - basehaz[,'time'])>=0)
        index = which.min((t - basehaz[,'time'])[id])

        return(1-exp(-basehaz[index,'hazard']))
    }
}

baseCDF <- function(t, basehaz){
    cdf = sapply(t, baseCDF.single, basehaz = basehaz)
    return(cdf)
}

baseS.single <- function(t, basehaz){
    if(t<min(basehaz[,'time'])){
        return(1)   # check 1 or 0
    }else{
        id = which((t - basehaz[,'time'])>=0)
        index = which.min((t - basehaz[,'time'])[id])

        return(exp(-basehaz[index,'hazard']))
    }
}

baseS <- function(t, basehaz){
    S = sapply(t, baseS.single, basehaz = basehaz)
    return(S)
}


## function for computing G(t|z)
Gz <- function(t, Z, basehaz.Q2, beta.Q2, tau){
    if(is.matrix(Z)){
        result = matrix(nrow = nrow(Z), ncol = length(t))
        for(j in 1:length(t)){
            if(tau-t[j]<min(basehaz.Q2$time)){
                result[,j] = 1
            }else{
                result[,j] = (baseS(tau-t[j], basehaz.Q2))^exp(Z %*% beta.Q2)
            }
        }
    }else{
        result = matrix(nrow = 1, ncol = length(t))
        for(j in 1:length(t)){
            if(tau-t[j]<min(basehaz.Q2$time)){
                result[,j] = 1
            }else{
                result[,j] = (baseS(tau-t[j], basehaz.Q2))^exp(sum(Z*beta.Q2))
            }
        }
    }
    colnames(result) = t

    return(result)
}

## function for computing F(t|z)
Fz <- function(t, Z, basehaz.T, beta.T){
    if(is.matrix(Z)){
        result = matrix(nrow = nrow(Z), ncol = length(t))
        for(j in 1:length(t)){
            result[,j] = 1-(baseS(t[j], basehaz.T))^exp(Z %*% beta.T)
        }
    }else{
        result = matrix(nrow = 1, ncol = length(t))
        for(j in 1:length(t)){
            result[,j] = 1-(baseS(t[j], basehaz.T))^exp(sum(Z*beta.T))
        }
    }
    colnames(result) = t

    return(result)
}





# ### Functions needed when using LTRCforests::ltrcrrf -------------------------------------
#
# # returns a function F_hat that takes 'newdata' as the input, and outputs a matrix
# # of the estimated survival probabilities for each subject in 'newdata' (row)
# # at the times in time.eval (col).
#
# F_hat.ltrcrrf <- function(dat, time.name, Q.name, event.name, cov.names,
#                           mtry, ntree){
#
#     # Z.names = c(A.names, cov.names)
#     Z.names = cov.names
#     formula.T = formula(paste(paste("Surv(", Q.name,",", time.name, ",", event.name,") ~ ", collapse = ""),
#                               paste(Z.names, collapse = "+"), collapse = ""))
#
#     T.obj <- ltrcrrf(formula = formula.T, data = dat, mtry = mtry, ntree = ntree)
#     jumps.T = sort(dat[,time.name])
#
#     F_hat <- function(newdata, time.eval = jumps.T){ # by default, evaluated at jumps.T
#         newdata[, Q.name] = 0
#         T.pred = predictProb(object = T.obj, newdata = newdata, time.eval = time.eval)  #**
#         survProb.mx = t(T.pred$survival.probs)
#         rownames(survProb.mx) = rownames(newdata)
#         colnames(survProb.mx) = time.eval
#
#         return(1-survProb.mx)
#     }
#
#     return(F_hat)
#
# }
#
#
# # returns a function F_hat that takes 'newdata' as the input, and outputs a matrix
# # of the estimated survival probabilities for each subject in 'newdata' (row)
# # at the times in time.eval (col).
#
# G_hat.ltrcrrf <- function(dat, time.name, Q.name, event.name, cov.names,
#                           mtry, ntree, tau = max(dat[,time.name])+1){
#
#     if(tau < max(dat[,time.name])+1e-10){
#         stop("Need a larger tau.")
#     }
#
#     names = c(time.name, Q.name, event.name, cov.names)
#     if(is.na(match("Q2", names)) & is.na(match("T2", names))){
#         dat$Q2 = tau - dat[,Q.name]
#         dat$T2 = tau - dat[,time.name]
#     }else{
#         stop("The names of the variables cannot be 'T2' or 'Q2'.")
#     }
#
#     # Z.names = c(A.names, cov.names)
#     Z.names = cov.names
#     formula.Q2 = formula(paste(paste("Surv(T2, Q2,", event.name,") ~ ", collapse = ""),
#                                paste(Z.names, collapse = "+"), collapse = ""))
#
#     Q2.obj <- ltrcrrf(formula = formula.Q2, data = dat, mtry = mtry, ntree = ntree)
#     jumps.Q = sort(dat[, Q.name])
#
#     G_hat <- function(newdata, time.eval = jumps.Q){ # by default, evaluated at jumps.Q
#         time.eval2 = tau-time.eval-1e-10
#         if(min(time.eval2)<0){
#             stop("Need a larger tau.")
#         }
#
#         newdata$Q2 = tau - newdata[,Q.name]
#         newdata$T2 = 0
#
#         odr = order(time.eval2)
#         Q.pred = predictProb(object = Q2.obj, newdata = newdata,
#                              time.eval = (time.eval2)[odr])   #**
#         odr2 = match((time.eval2)[odr], time.eval2)
#         prob.mx = t(Q.pred$survival.probs)[,odr2]
#         rownames(prob.mx) = rownames(newdata)
#         colnames(prob.mx) = time.eval
#
#         return(prob.mx)
#     }
#
#     return(G_hat)
# }




