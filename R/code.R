.checkSquareMatrix = function(x) {
    res = checkmate::checkMatrix(x, mode = "numeric")
    if (!isTRUE(res))
        return(res)
    if (nrow(x) != ncol(x))
        return("Must be square")
    return(TRUE)
}

# For assertions:
.assertSquareMatrix = checkmate::makeAssertionFunction(.checkSquareMatrix)


Jn <- function(n = 2){
    checkmate::assertCount(n, positive = TRUE)
    .Fortran("jmatrix", n = as.integer(n), J = diag(n))$J
}

Pm <- function(m = 2, nu = 0){

    coll = checkmate::makeAssertCollection()
    checkmate::assertCount(m, positive = TRUE, add = coll)
    checkmate::assertCount(nu, add = coll)
    checkmate::assertTRUE(m >= nu, add = coll)
    checkmate::reportAssertions(coll)
    .Fortran("pmatrix", m = as.integer(m),  v = as.integer(nu), P = diag(m+1))$P

}

Qm <- function(m = 2, nu = 0, P = NULL){

    coll = checkmate::makeAssertCollection()
    if(is.null(P)){
        checkmate::assertCount(m, positive = TRUE, add = coll)
        checkmate::assertCount(nu, add = coll)
        checkmate::assertTRUE(m >= nu, add = coll)
    }else{
        .assertSquareMatrix(P, add = coll)
        #warning("P is provided! m and nu will be ignored")
    }
    checkmate::reportAssertions(coll)

    if(is.null(P)) out = .Fortran("qmatrix", m = as.integer(m), v = as.integer(nu), Q = diag(m+1))$Q
    else out = .Fortran("qm", m = nrow(P)-1L, P = P,  Q = diag(nrow(P)))$Q
    return(out)
}

Km <- function(m = 3, nu = 0, J = NULL, Q = NULL){

    coll = checkmate::makeAssertCollection()
    if(is.null(J) | is.null(Q)){
        # Q and/or J not provided:
        checkmate::assertCount(m, positive = TRUE, add = coll)
        checkmate::assertCount(nu, add = coll)
        checkmate::assertTRUE(m >= nu, add = coll)
        if(!is.null(J)) warning("Q not provided. Ignoring J")
        if(!is.null(Q)) warning("J not provided. Ignoring Q")
    }else{
        .assertSquareMatrix(J, add = coll)
        .assertSquareMatrix(Q, add = coll)
        checkmate::assertTRUE(all(dim(J) == dim(Q)), add = coll)
        #warning("J and Q are provided! m and nu will be ignored")
    }
    checkmate::reportAssertions(coll)

    if(is.null(J) | is.null(Q)) out = .Fortran("Kmatrix", m = as.integer(m), nu = as.integer(nu), K = diag(m+1))$K
    else out = .Fortran("km", m = nrow(J)-1L, J = J, Q = Q, K = diag(nrow(J)))$K
    return(out)
}

Kkronm <- function(m = 3, nu = 0, h = 0, overlap = TRUE, K = NULL){

    coll = checkmate::makeAssertCollection()
    checkmate::assertCount(h, add = coll)
    checkmate::assertLogical(overlap, len = 1, add = coll)
    if(is.null(K)){
        # K is not provided:
        checkmate::assertCount(m, positive = TRUE, add = coll)
        checkmate::assertCount(nu, add = coll)
        checkmate::assertTRUE(m >= nu, add = coll)
    }else
        .assertSquareMatrix(K, add = coll)

    checkmate::reportAssertions(coll)

    overlap = as.integer(overlap)
    if(is.null(K)){
        mm = (m+1)*(h+1) - m*h*overlap
        out = .Fortran("kkronmatrix", m = as.integer(m), h = as.integer(h),
                       nu = as.integer(nu), overlap = overlap,
                       Kkron = diag(mm*(m+1)))$Kkron
    }else{
        m = nrow(K)-1L
        mm = (m+1)*(h+1) - m*h*overlap
        out = .Fortran("kkronm", m = m, h = as.integer(h), overlap = overlap,
                       K = K, Kkron = diag(mm*(m+1)))$Kkron
    }
    return(out)
}


F2dfa <- function(y, m = 3, nu = 0, overlap = TRUE){

    coll = checkmate::makeAssertCollection()
    checkmate::assertAtomicVector(y, add = coll)
    checkmate::assertAtomicVector(m, add = coll)
    checkmate::assertNumeric(m, lower = 1, add = coll)
    checkmate::assertCount(nu, add = coll)
    checkmate::assertTRUE(min(m) >= nu, add = coll)
    checkmate::assertLogical(overlap, len = 1, add = coll)
    checkmate::reportAssertions(coll)

    n = length(y)
    lm = length(m)
    out  = .Fortran("dfadcca", lm = lm,
                    m = as.integer(m), nu = as.integer(nu),
                    overlap = as.integer(overlap),
                    n = n, y1 = y, y2 = 0,
                    f1  = 1L, fd1m = numeric(lm),
                    f2 = 0L, fd2m = 0,
                    f12 = 0L, fd12m = 0,
                    rhom = 0)
    return(out$fd1m)
}



Fdcca <- function(y1, y2, m = 3, nu = 0, overlap = TRUE){

    coll = checkmate::makeAssertCollection()
    checkmate::assertAtomicVector(y1, add = coll)
    checkmate::assertAtomicVector(y2, add = coll)
    checkmate::assertAtomicVector(m, add = coll)
    checkmate::assertNumeric(m, lower = 1, add = coll)
    checkmate::assertCount(nu, add = coll)
    checkmate::assertTRUE(min(m) >= nu, add = coll)
    checkmate::assertLogical(overlap, len = 1, add = coll)
    checkmate::reportAssertions(coll)

    n1 = length(y1)
    n2 = length(y2)
    n =  as.integer(min(n1,n2))
    if(n1 != n2) warning("y1 and y2 have different lengths.",
                         "\nlength(y1) = ", n1, " and length(y2) = ", n2,
                         ". Coercing to sample size ", n)
    lm = length(m)
    out  = .Fortran("dfadcca", lm = lm, m = as.integer(m), nu = as.integer(nu),
                    overlap = as.integer(overlap),
                    n = n, y1 = y1[1:n], y2 = y2[1:n],
                    f1  = 0L, fd1m = 0,
                    f2 = 0L, fd2m = 0,
                    f12 = 1L, fd12m = numeric(lm),
                    rhom = 0)
    return(out$fd12m)
}


rhodcca <- function(y1, y2, m = 3, nu = 0, overlap = TRUE){

    coll = checkmate::makeAssertCollection()
    checkmate::assertAtomicVector(y1, add = coll)
    checkmate::assertAtomicVector(y2, add = coll)
    checkmate::assertAtomicVector(m, add = coll)
    checkmate::assertNumeric(m, lower = 1, add = coll)
    checkmate::assertCount(nu, add = coll)
    checkmate::assertTRUE(min(m) >= nu, add = coll)
    checkmate::assertLogical(overlap, len = 1, add = coll)
    checkmate::reportAssertions(coll)

    n1 = length(y1)
    n2 = length(y2)
    n =  as.integer(min(n1,n2))
    if(n1 != n2) warning("y1 and y2 have different lengths.",
                         "\nlength(y1) = ", n1, " and length(y2) = ", n2,
                         ". Coercing to sample size ", n)
    lm = length(m)

    out  = .Fortran("dfadcca", lm = lm, m = as.integer(m), nu = as.integer(nu),
                    overlap = as.integer(overlap),
                    n = n, y1 = y1[1:n], y2 = y2[1:n],
                    f1  = 1L, fd1m = numeric(lm),
                    f2 = 1L, fd2m = numeric(lm),
                    f12 = 1L, fd12m = numeric(lm),
                    rhom = numeric(lm))
    return(list(F2dfa1 = out$fd1m,
                F2dfa2 = out$fd2m,
                Fdcca = out$fd12m,
                rhodcca = out$rhom))
}



EF2dfa <- function(m = 3, nu = 0, G, K = NULL){

    coll = checkmate::makeAssertCollection()
    checkmate::assertMatrix(G, add = coll)
    if(is.null(K)){
        checkmate::assertAtomicVector(m, add = coll)
        checkmate::assertNumeric(m, lower = 1, add = coll)
        checkmate::assertCount(nu, add = coll)
        checkmate::assertTRUE(min(m) >= nu, add = coll)
    }else{
        #warning("K is provided. Ignoring m")
        .assertSquareMatrix(K, add = coll)
        m = nrow(K)-1L
    }
    checkmate::assertTRUE(all(dim(G) >= (max(m)+1)), add = coll)
    checkmate::reportAssertions(coll)

    lm = length(m)
    if(!is.null(K)) out = .Fortran("em", m = m, K = K, G = G[1:(m+1), 1:(m+1)], E = 0)$E
    else out = .Fortran("Expectms", lm = lm, m = as.integer(m),
                                    nu = as.integer(nu),
                                    c1 = 1L, G1 = G[1:(max(m)+1), 1:(max(m)+1)],
                                    c2 = 0L, G2 = G[1,1], c12 = 0L, G12 = G[1,1],
                                    E1 = numeric(lm), E2 = 0, E12 = 0, rhos = 0)$E1
    return(out)
}


EFdcca <- function(m = 3, nu = 0, G, K = NULL){

    coll = checkmate::makeAssertCollection()
    checkmate::assertMatrix(G, add = coll)
    if(is.null(K)){
        checkmate::assertAtomicVector(m, add = coll)
        checkmate::assertNumeric(m, lower = 1, add = coll)
        checkmate::assertCount(nu, add = coll)
        checkmate::assertTRUE(min(m) >= nu, add = coll)
    }else{
        #warning("K is provided. Ignoring m")
        .assertSquareMatrix(K, add = coll)
        m = nrow(K)-1L
    }
    checkmate::assertTRUE(all(dim(G) >= (max(m)+1)), add = coll)
    checkmate::reportAssertions(coll)

    lm = length(m)
    if(!is.null(K)) out = .Fortran("em", m = m, K = K, G = G[1:(m+1), 1:(m+1)], E = 0)$E
    else out = .Fortran("Expectms", lm = length(m), m = as.integer(m),
                        v = as.integer(nu),
                        c1 = 0L, G1 = G[1,1], c2 = 0L, G2 = G[1,1],
                        c12 = 1L, G = G[1:(max(m)+1), 1:(max(m)+1)],
                        E1 = 0, E2 = 0, E12 = numeric(lm), rhos = 0)$E12
    return(out)
}



rhoE <- function(m = 3, nu = 0, G1, G2, G12, K = NULL){

    coll = checkmate::makeAssertCollection()
    checkmate::assertMatrix(G1, add = coll)
    checkmate::assertMatrix(G2, add = coll)
    checkmate::assertMatrix(G12, add = coll)
    if(is.null(K)){
        checkmate::assertAtomicVector(m, add = coll)
        checkmate::assertNumeric(m, lower = 1, add = coll)
        checkmate::assertCount(nu, add = coll)
        checkmate::assertTRUE(min(m) >= nu, add = coll)
    }else{
        #warning("K is provided. Ignoring m")
        .assertSquareMatrix(K, add = coll)
        m = nrow(K)-1L
    }
    checkmate::assertTRUE(all(dim(G1) >= (max(m)+1)), add = coll)
    checkmate::assertTRUE(all(dim(G2) >= (max(m)+1)), add = coll)
    checkmate::assertTRUE(all(dim(G12) >= (max(m)+1)), add = coll)
    checkmate::reportAssertions(coll)

    lm = length(m)
    mmax = max(m)+1
    if(!is.null(K)){
        E1 = .Fortran("em", m = m, K = K, G = G1[1:(m+1), 1:(m+1)], E = 0)$E
        E2 = .Fortran("em", m = m, K = K, G = G2[1:(m+1), 1:(m+1)], E = 0)$E
        E12 = .Fortran("em", m = m, K = K, G = G12[1:(m+1), 1:(m+1)], E = 0)$E
        Erho = E12/sqrt(E1*E2)
        out = list(E1 = E1, E2 = E2, E12 = E12, Erhos = Erho)
    }
    else  out = .Fortran("Expectms", lm = length(m), m = as.integer(m),
                   nu = as.integer(nu),
                   c1 = 1L, G1 = G1[1:mmax, 1:mmax],
                   c2 = 1L, G2 = G2[1:mmax, 1:mmax],
                   c12 = 1L, G = G12[1:mmax, 1:mmax],
                   E1 = numeric(lm), E2 = numeric(lm),
                   E12 = numeric(lm), Erhos = numeric(lm))
    return(list(EF2dfa1 = out$E1, EF2dfa2 = out$E2, EFdcca = out$E12, rhoE = out$Erhos))
}


covF2dfa <- function(m = 3, nu = 0, h = 0, overlap = TRUE, G, Cumulants = NULL){

    coll = checkmate::makeAssertCollection()
    checkmate::assertAtomicVector(m, add = coll)
    checkmate::assertCount(nu, add = coll)
    checkmate::assertTRUE(min(m) >= nu, add = coll)
    checkmate::assertAtomicVector(h, add = coll)
    checkmate::assertLogical(overlap, len = 1, add = coll)
    checkmate::assertMatrix(G, add = coll)
    checkmate::reportAssertions(coll)

    nr = as.integer(max(m)+1)
    mm = as.integer((max(m)+1)*(max(h)+1) - max(m)*max(h)*as.integer(overlap))

    checkmate::assertTRUE(nrow(G) >= nr, add = coll)
    checkmate::assertTRUE(ncol(G) >= mm, add = coll)
    if(!is.null(Cumulants)){
        checkmate::assertMatrix(Cumulants, add = coll)
        checkmate::assertTRUE(any(dim(Cumulants) >= mm*nr), add = coll)
    }
    checkmate::reportAssertions(coll)

    if(any(h < 0)){
        #warning("h must be positive. Using |h|")
        h = abs(h)}

    nrc = 0
    if(is.null(Cumulants)) Cumulants = matrix(0, ncol = 1, nrow = 1)
    else{
        if(sum(abs(Cumulants)) <= .Machine$double.eps*nrow(Cumulants)^2)
            Cumulants  = matrix(0, ncol = 1, nrow = 1)
        else nrc = mm*nr
    }

    out = .Fortran("covf2dfa", lm = length(m), m = as.integer(m),
             nu = as.integer(nu), lh = length(h), h = as.integer(h),
             overlap = as.integer(overlap),
             nr1 = nr, nc1 = mm, G1 = G[1:nr,1:mm],
             nrc = nrc, Cum = Cumulants[1:max(1,nrc),1:max(1,nrc)],
             df = 1L, cova = matrix(1, nrow = length(m), ncol = length(h)))$cova

    colnames(out) = paste("h=",h, sep = "")
    rownames(out) = paste("m=",m, sep = "")
    return(out)
}

covFdcca <- function(m = 3, nu = 0, h = 0, overlap = TRUE,
                     G1, G2, G12, Cumulants = NULL){

    coll = checkmate::makeAssertCollection()
    checkmate::assertAtomicVector(m, add = coll)
    checkmate::assertCount(nu, add = coll)
    checkmate::assertTRUE(min(m) >= nu, add = coll)
    checkmate::assertAtomicVector(h, add = coll)
    checkmate::assertLogical(overlap, len = 1, add = coll)
    checkmate::assertMatrix(G1, add = coll)
    checkmate::assertMatrix(G2, add = coll)
    checkmate::assertMatrix(G12, add = coll)
    checkmate::reportAssertions(coll)

    nr = as.integer(max(m)+1)
    mm = as.integer((max(m)+1)*(max(h)+1) - max(m)*max(h)*as.integer(overlap))

    checkmate::assertTRUE(nrow(G1) >= nr, add = coll)
    checkmate::assertTRUE(ncol(G1) >= mm, add = coll)
    checkmate::assertTRUE(nrow(G2) >= mm, add = coll)
    checkmate::assertTRUE(ncol(G2) >= nr, add = coll)
    checkmate::assertTRUE(nrow(G12) >= mm, add = coll)
    checkmate::assertTRUE(ncol(G12) >= mm, add = coll)
    if(!is.null(Cumulants)){
        checkmate::assertMatrix(Cumulants, add = coll)
        checkmate::assertTRUE(any(dim(Cumulants) >= mm*nr), add = coll)
    }
    checkmate::reportAssertions(coll)

    if(any(h < 0)){
        #warning("h must be positive. Using |h|")
        h = abs(h)}

    nrc = 0
    if(is.null(Cumulants)) Cumulants = matrix(0, ncol = 1, nrow = 1)
    else{
        if(sum(abs(Cumulants)) <= .Machine$double.eps*nrow(Cumulants)^2)
            Cumulants  = matrix(0, ncol = 1, nrow = 1)
        else nrc = mm*nr
    }

    out = .Fortran("covfdcca", lm = length(m), m = as.integer(m),
             nu = as.integer(nu), lh = length(h), h = as.integer(h),
             overlap = as.integer(overlap),
             nr1 = nr, nc1 = mm, G1 = G1[1:nr,1:mm],
             nr2 = mm, nc2 = nr, G2 = G2[1:mm,1:nr],
             nr12 = mm, nc12 = mm, G12 = G12[1:mm,1:mm],
             nrc = nrc, Cum = Cumulants[1:max(1,nrc),1:max(1,nrc)],
             df = 0L, cova = matrix(0, nrow = length(m), ncol = length(h)))$cova

    colnames(out) = paste("h=",h, sep = "")
    rownames(out) = paste("m=",m, sep = "")
    return(out)

}
