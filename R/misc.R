
## functions for computing the CDF at given times ----------------------------------

CDF_eval <- function(time.eval, CDF.mx){
    if(length(time.eval) != nrow(CDF.mx)){
        stop("The number of time points does not equal the number of subjects!")
    }

    jumps = as.numeric(colnames(CDF.mx))
    CDF.mx = cbind(0, CDF.mx)

    probs = rep(NA, length(time.eval))
    for(i in 1:length(time.eval)){
        id = findInterval(time.eval[i], c(0, jumps, Inf))
        probs[i] = CDF.mx[i,id]
    }

    return(probs)
}


CDF_eval.mx.single <- function(time.eval, CDF.mx){
    jumps = as.numeric(colnames(CDF.mx))
    CDF.mx = cbind(0, CDF.mx)
    id = findInterval(time.eval, c(0, jumps, Inf))
    probs = CDF.mx[,id]

    return(probs)
}

CDF_eval.mx <- function(time.eval, CDF.mx){
    probs = sapply(time.eval, CDF_eval.mx.single, CDF.mx = CDF.mx)

    return(probs)
}



# functions for computing L-S integrals ---------------------------------------------

int_fmx_dF <- function(v, f.mx, F.mx){

    if(mean(dim(f.mx) == dim(F.mx))<1){
        stop("The dimensions of f.mx and F.mx are not the same!")
    }

    jumps = as.numeric(colnames(F.mx))
    dF.mx = F.mx - cbind(0, F.mx[,-ncol(F.mx)])
    # dF.mx = F.mx - cbind(F.mx[,1], F.mx[,-ncol(F.mx)])
    id = findInterval(v, c(0, jumps, Inf))
    vals = cbind(0, f.mx*dF.mx)

    inte = int_sum(id, vals)

    rownames(inte) = rownames(F.mx)
    colnames(inte) = v

    return(inte)
}

# id: a single index denoting the
int_sum.single <- function(id, vals){
    if(id==1){
        inte = vals[,1]
    }else{
        inte = rowSums(vals[,(1:id)])
    }
    return(inte)
}

int_sum <- function(id, vals){
    inte = sapply(id, int_sum.single, vals = vals)
    return(inte)
}



