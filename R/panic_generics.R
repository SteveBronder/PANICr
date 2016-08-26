
#'@export
summary.panic <- function(object,...){
  
  if (object$func == "panic04"){
  
    # Get significant figures markings and p values
    # Note: the second Results is the one we care about for the pooled tests
    commonsig_dm <- .get_sig(common_pvals,object$pooladf$Demeaned[2])
    idiosig <- .get_sig(ehat_pvals,object$pooladf$Idiosyncratic[2])
    commonsig <- .get_sig(common_pvals,object$Common$Common)
    
    commondm_fig <- paste0("< ",commonsig_dm$p_val," ",commonsig_dm$sign)
    idio_fig <- paste0("< ", idiosig$p_val," ", idiosig$sign)
    common_fig <- paste0("< ",commonsig$p_val," ", commonsig$sign)
    
    common_ans <- data.frame( Components = paste0("Component ", 1:length(common_fig)),
                              Results = object$Common$Common,
                              "P.Results" = common_fig
                              )
    
    pooled_ans <- data.frame( Tests = c("Idiosyncratic Test", "Demeaned Test"),
                            Results = c(object$pooladf$Idiosyncratic[2], object$pooladf$Demeaned[2]),
                            "P.Results" = c(idio_fig, commondm_fig))

    
    summary_results <- list(common_tests = common_ans,
                            pooled_tests = pooled_ans,
                            func = object$func,
                            nfac = object$nfac,
                            ic = object$ic,
                            criteria = object$criteria)
    
    summary_panic <- structure(summary_results,
                               class = "summary.panic")
    return(summary_panic)
  } else if(object$func == "panic10nm"){
    
    mp_testa_p <- .get_sig(ehat_pvals,object$MP.tests$model_a)
    mp_testb_p <- .get_sig(ehat_pvals,object$MP.tests$model_b)
    mp_testa_fig <- paste0("< ", mp_testa_p$p_val," ", mp_testa_p$sign)
    mp_testb_fig <- paste0("< ", mp_testb_p$p_val," ", mp_testb_p$sign)
    
    pmsb_p <- .get_sig(ehat_pvals,object$PMSB.tests$PMSB)
    pmsb_fig <- paste0("< ", pmsb_p$p_val," ", pmsb_p$sign)
    
    mp_test_a <- data.frame(Tests = c("ta","tb"),
                            Results = object$MP.tests$model_a,
                            Signif = mp_testa_fig)
    
    mp_test_b <- data.frame(Tests = c("ta","tb"),
                            Results = object$MP.tests$model_b,
                            Signif = mp_testb_fig)
    
    pmsb_test <- data.frame(Tests = "PMSB",
                            Results = object$PMSB.tests$PMSB,
                            Signif = pmsb_fig)
    
    summary_results <- list(mp_test_a = mp_test_a,
                            mp_test_b = mp_test_b,
                            pmsb_test = pmsb_test,
                            func = object$func,
                            nfac = object$nfac,
                            ic = object$ic,
                            criteria = object$criteria)
    
    summary_panic <- structure(summary_results,
                               class = "summary.panic")
    
    return(summary_panic)
    
    
  } else if(object$func == "panic10m"){
    
    pool_test_p <- .get_sig(ehat_pvals,object$MP.tests$P)
    mp_testc_p <- .get_sig(ehat_pvals,object$MP.tests$model_c)
    pool_test_fig <- paste0("< ", pool_test_p$p_val," ", pool_test_p$sign)
    mp_testc_fig <- paste0("< ", mp_testc_p$p_val," ", mp_testc_p$sign)
    
    pmsb_p <- .get_sig(ehat_pvals,object$PMSB.tests$PMSB)
    pmsb_fig <- paste0("< ", pmsb_p$p_val," ", pmsb_p$sign)
    
    pool_test <- data.frame(Tests = c("Pa","Pb"),
                            Results = object$MP.tests$P,
                            Signif = pool_test_fig)
    
    mp_test_c <- data.frame(Tests = c("ta","tb"),
                            Results = object$MP.tests$model_c,
                            Signif = mp_testc_fig)
    
    pmsb_test <- data.frame(Tests = "PMSB",
                            Results = object$PMSB.tests$PMSB,
                            Signif = pmsb_fig)
    
    summary_results <- list(pool_test = pool_test,
                            mp_test_c = mp_test_c,
                            pmsb_test = pmsb_test,
                            func = object$func,
                            nfac = object$nfac,
                            ic = object$ic,
                            criteria = object$criteria)
    
    summary_panic <- structure(summary_results,
                               class = "summary.panic")
    
    return(summary_panic)
    }
}


#'@export
print.summary.panic <- function(x,...){
  
  if (x$func == "panic04"){
    cat("\nCommon Components\n")
    print(format(x$common_tests,trim=TRUE,justify = "left"))
    cat("\n---\n")
    cat("\nPooled Tests\n")
    print(format(x$pooled_tests, trim=TRUE,justify= "left"))
    cat("\n Signif. codes:  '***' 0.01 '**' 0.025 '*' 0.05 '.' 0.1 ' ' 1\n")
    
    cat(paste0("\n Maximum Number of Factors: ", x$nfac))
    cat(paste0("\n Estimated Number of Factors: ", x$ic))
    cat(paste0("\n Criteria Used: ", x$criteria))
    
  } else if (x$func == "panic10nm"){
    
    cat("\n Moon and Perron Tests A and B\n")
    print(format(x$mp_test_a,trim=TRUE,justify = "left"))
    cat("\n")
    print(format(x$mp_test_b,trim=TRUE,justify = "left"))
    cat("\nPMSB Test\n")
    print(format(x$pmsb_test,trim=TRUE,justify="left"))
    cat("\n Signif. codes:  '***' 0.01 '**' 0.025 '*' 0.05 '.' 0.1 ' ' 1\n")
    
    cat(paste0("\n Maximum Number of Factors: ", x$nfac))
    cat(paste0("\n Estimated Number of Factors: ", x$ic))
    cat(paste0("\n Criteria Used: " , x$criteria))
    
    }else if(x$func == "panic10m"){
      
      cat("\n Results for PANIC 2010 Tests:")
      cat("\n Pooled Tests\n")
      print(format(x$pool_test,trim=TRUE,justify = "left"))
      cat("\nMoon and Perron Test C\n")
      print(format(x$mp_test_c,trim=TRUE, justify = "left"))
      cat("\nPMSB Test\n")
      print(format(x$pmsb_test,trim=TRUE,justify="left"))
      cat("\n Signif. codes:  '***' 0.01 '**' 0.025 '*' 0.05 '.' 0.1 ' ' 1\n")
      
      cat(paste0("\n Maximum Number of Factors: ", x$nfac))
      cat(paste0("\n Estimated Number of Factors: ", x$ic))
      cat(paste0("\n Criteria Used: ", x$criteria))
      }
}


.get_sig <- function(pvals,ans){
  
  sigs <- lapply(ans,
                 function(xx) {
                   pvals[which.min(abs(xx - pvals[,1])),c("p_val","sign")]
                   })
  do.call(rbind,sigs)
}

