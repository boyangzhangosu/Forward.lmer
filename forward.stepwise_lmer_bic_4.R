forward.lmer <- function(
  start.model, blocks, max.iter=1,
  sign.test=FALSE,
  zt=FALSE,
  bic=FALSE,
  print.log=TRUE)
{
  
  # forward.lmer: a function for stepwise regression using lmer mixed effects models
  # Author list: Rense Nieuwenhuis, Ferraz Dias de Moraes, Luis Eduardo, Boyang Zhang (zhang.7077@osu.edu)
  
  # Initialysing internal variables
  log.step <- 0
  log.LL <- log.p <- log.block <- zt.temp <- log.zt <- log.bic<- NA
  model.basis <- start.model
  blocks <- blocks
  #library
  library(lme4)
  library(Matrix)
  # Maximum number of iterations cannot exceed number of blocks
  if (max.iter > length(blocks)) max.iter <- length(blocks)
  
  # Setting up the outer loop
  for(i in 1:max.iter)
  {
    print(paste("Iteration ", i, sep=""))
    
    models <- list()
    
    # Iteratively updating the model with addition of one block of variable(s)
    # Also: extracting the loglikelihood of each estimated model
    for(j in 1:length(blocks))
    {
      models[[j]] <- update(model.basis, as.formula(paste(". ~ . + ", blocks[j])))
    }
    
    LL <- unlist(lapply(models, logLik))
    
    # Ordering the models based on their loglikelihood.
    # Additional selection criteria apply
    for (j in order(LL, decreasing=TRUE))
    {
      
      ##############
      ############## Selection based on ANOVA-test
      ##############
      
      if(sign.test != FALSE)
      {
        if(anova(model.basis, models[[j]])[2,8] < sign.test)
        {
          
          model.basis <- models[[j]]
          
          # Writing the logs
          log.step <- log.step + 1
          log.block[log.step] <- blocks[j]
          log.LL[log.step] <- as.numeric(logLik(model.basis))
          log.p[log.step] <- anova(model.basis, models[[j]])[2,8]
          
          blocks <- blocks[-j]
          
          break
        }
      }
      
      ##############
      ############## Selection based significance of added variable-block
      ##############
      
      if(zt != FALSE)
      {
        
        coefs <- data.frame(coef(summary(models[[j]])))
        diff.par <- setdiff(rownames(coefs), rownames(summary(model.basis)$coefficients))
        if (length(diff.par)==0) break
        sig.par <- FALSE
        
        for (k in 1:length(diff.par))
        {
          t1 <- abs(coefs[which(rownames(coefs)==diff.par[k]),3])
          
          p1 <- 2 * (1 - pnorm(t1))
          
          if(p1 < zt)
            
          {
            sig.par <- TRUE
            zt.temp <- coefs[which(rownames(coefs)==diff.par[k]),3]
            break
          }
        }
        
        if(sig.par==TRUE)
        {
          model.basis <- models[[j]]
          
          # Writing the logs
          log.step <- log.step + 1
          log.block[log.step] <- blocks[j]
          log.LL[log.step] <- as.numeric(logLik(model.basis))
          log.zt[log.step] <- zt.temp
          blocks <- blocks[-j]
          
          break
        }
      }
      
      if(bic != FALSE)
      {
        BIC_vec_1 <- BIC(model.basis)
        lmer1 <- try(models[[j]])  
        
        
        if(class(lmer1)!="try-error")
        {  
          BIC_vec[j]<- BIC(models[[j]])
          c.1 <- BIC_vec[[j]]-BIC_vec_1
          if(c.1<0)
          { model.basis <- models[[j]]
          # Writing the logs
          log.step <- log.step + 1
          log.block[log.step] <- blocks[j]
          log.LL[log.step] <- as.numeric(logLik(model.basis))
          log.bic[log.step] <- BIC(model.basis)
          blocks <- blocks[-j]
          
          
          break
          
          }
          
          
        }
      }
    }
  }
  
  ## Create and print log
  log.df <- data.frame(log.step=1:log.step, log.block, log.LL, log.p, log.zt)
  if(print.log == TRUE) print(log.df, digits=4)
  
  ## Return the 'best' fitting model
  return(model.basis)
}


