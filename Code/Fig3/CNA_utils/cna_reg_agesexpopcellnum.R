.reg <- function(q, k) {
        Xpc = U[, 1:k,drop=F] 
        # beta <- t(Xpc) %*% q
        # I changed here: infer beta adjusted for covariates
        
        mat=left_join(q,Xpc%>% as.data.frame()%>%rownames_to_column("sample_uuid"),by="sample_uuid")

        for(mmm in 1:k){
            lm_tmp=lm(paste0("score ~ PC",mmm,"+ Age + Sex + pop_cov + Cellnum_per_sample"),data=mat)
            beta_tmp=summary(lm_tmp)$coefficient[2,1,drop=F]
            if(mmm==1){beta=beta_tmp}else{beta=rbind(beta,beta_tmp)}
         }
         colnames(beta)=NULL

        qhat <- Xpc %*% beta
        return(list(qhat = qhat, beta = beta))
    }

.stats <- function(yhat, ycond, k) {
        ssefull <- crossprod(yhat - ycond)
        ssered <- crossprod(ycond)
        deltasse <-  ssered - ssefull
        f <- (deltasse / k) / (ssefull/n)
        p <- -pf(f, k, n-(1+r+k), log.p = TRUE)    
        r2 <- 1 - ssefull/ssered
        return(list(p=p, r2=r2))
    }

.minp_stats <- function(q) {
        qhats <- purrr::map(ks, function(k) .reg(q, k)$qhat)
        .tmp <- purrr::map2(qhats, ks, function(qhat, k) .stats(qhat, q$score, k))
        ps <- purrr::map_dbl(.tmp, 'p')
        r2s <- purrr::map_dbl(.tmp, 'r2')
        k_ <- which.min(ps)
        return(list(k=ks[k_], p=ps[k_], r2=r2s[k_]))
    }

