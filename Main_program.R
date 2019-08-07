#######################################################
##===================================================##
##                   preparation                     ##
##===================================================##
#######################################################


#===================================================#
#             remove previous objects               #
#===================================================#


rm(list = ls())


#===================================================#
#               install package                     #
#===================================================#


## List of packages should be install.
list.of.packages <- c("caret", "caroline", "cowplot", "dplyr", "dplyr", "ggplot2", "gridExtra", "microbenchmark", "partykit", "pscl", "pscl", "sandwich", "statmod", "VGAM", "xtable")
## The package is already installed.
installed_pkgs = installed.packages()[, "Package"]
## Check for missing packages.
new.packages <- list.of.packages[!(list.of.packages %in% installed_pkgs)]
if (length(new.packages)) install.packages(new.packages)
if (!"countreg" %in% installed_pkgs) install.packages("countreg", repos = "http://R-Forge.R-project.org")


#===================================================#
#                  load package                     #
#===================================================#


lapply(list.of.packages, require, character.only = TRUE)
library(countreg)



#######################################################
##===================================================##
##                     Method                        ##
##===================================================##
#######################################################


backeli <- function(object, data, dist = c("poisson", "negbin", "geometric"), alpha = 0.05, trace = FALSE) {
  FUN = if (class(object) != "zeroinfl") {
    ## The hurdle model considers all zeros to be “structural zeros”.
    hurdle
  } else {
    zeroinfl
  }
  dist <- match.arg(dist)
  fit <- object
  cats = sapply(fit$model[, -1], function(x) if (is.factor(x)) 
    length(levels(x)) - 1 else 1)
  rhs1 <- attr(fit$terms$count, "term.labels")
  rhs2 <- attr(fit$terms$zero, "term.labels")
  rhs1 = rep(names(cats[rhs1]), cats[rhs1])
  rhs2 = rep(names(cats[rhs2]), cats[rhs2])
  
  nj <- length(rhs1) * length(rhs2)
  j <- 1
  if (trace) {
    cat("Initial model\n")
    print(summary(fit))
  }
  RET <- matrix(NA, nrow = nj, ncol = 3)
  colnames(RET) <- c("loglik", "BIC", "AIC")
  while (T) {
    if (trace) 
      cat("\nstep", j, "\n")
    coef <- summary(fit)$coef
    
    ## excluding intercept
    d <- dim(coef$count)[1]
    ## Find maximum p-value from count model.
    if (dist != "negbin") 
      count.pval <- coef$count[-1, 4]  
    else count.pval <- coef$count[-c(1, d), 4]
    
    if (trace) 
      cat("\ncount", count.pval, "\n")
    
    zero.pval <- coef$zero[-1, 4]  
    
    if (trace) 
      cat("\nzero", zero.pval, "\n")
    
    nc <- length(count.pval)
    nz <- length(zero.pval)
    ## Which variable has maximum p-value from count model.
    if (dist != "negbin")
      count.order <- order(coef$count[-1, 4], decreasing = TRUE)
    else count.order <- order(coef$count[-c(1, d), 4], decreasing = TRUE)
    
    zero.order <- order(coef$zero[-1, 4], decreasing = TRUE)
    
    if (trace) 
      cat("\ncount", count.order, "\n")
    
    if (trace) 
      cat("\nzero", zero.order, "\n")
    
    cats = sapply(fit$model[, -1], function(x) if (is.factor(x))
      length(levels(x)) - 1 else 1)
    rhs1 <- attr(fit$terms$count, "term.labels")
    rhs2 <- attr(fit$terms$zero, "term.labels")
    rhs1 = rep(names(cats[rhs1]), cats[rhs1])
    rhs2 = rep(names(cats[rhs2]), cats[rhs2])
    
    if (trace) 
      cat("\nr1", rhs1, "\n")
    
    if (trace) 
      cat("\nr2", rhs2, "\n")
    
    kc <- 1
    kz <- 1
    count.max <- count.pval[count.order[kc]]
    zero.max <- zero.pval[zero.order[kz]]
    if (is.na(count.max) && is.na(zero.max)) 
      break else if (is.na(zero.max)) 
        zero.max <- 0 else if (is.na(count.max)) 
          count.max <- 0
    if (count.max > zero.max) 
      if (count.max > alpha) {
        newid <- count.order[kc]
        if (dist != "negbin") 
          dropvar <- rownames(coef$count)[-1][newid] else dropvar <- rownames(coef$count)[-c(1, d)][newid]
          if (trace) 
            cat("drop variable in count component: ", rhs1[newid], "\n")
          rhs1 <- rhs1[rhs1 %out% rhs1[newid]]
          ## kc <- kc + 1
      } else break else if (zero.max > alpha) {
        newid <- zero.order[kc]
        dropvar <- rownames(coef$zero)[-1][newid]
        if (trace) 
          cat("drop variable in zero component: ", rhs2[newid], "\n")
        rhs2 <- rhs2[rhs2 %out% rhs2[newid]]
      } else break
    rhs1 = unique(rhs1)
    rhs2 = unique(rhs2)
    if (length(rhs1) == 0) 
      rhs1tmp <- 1 else {
        rhs1tmp <- rhs1[1]
        if (length(rhs1) > 1) 
          for (i in 2:length(rhs1)) rhs1tmp <- paste(rhs1tmp, "+", rhs1[i])
      }
    if (length(rhs2) == 0) 
      rhs2tmp <- 1 else {
        rhs2tmp <- rhs2[1]
        if (length(rhs2) > 1) 
          for (i in 2:length(rhs2)) rhs2tmp <- paste(rhs2tmp, "+", rhs2[i])
      }
    ## response variable
    res <- deparse(terms(fit$terms$count)[[2]])
    out <- paste(res, "~", rhs1tmp, "|", rhs2tmp)
    ## Set the environment of the formula (i.e. where should R look for variables when data aren't specified?).
    environment(out) <- parent.frame()
    ofit = fit
    fit <- try(FUN(eval(parse(text = out)), data = data, dist = dist))
    if (inherits(fit, "try-error")) {
      fit = ofit
      break
    }
    if (trace) {
      print(summary(fit))
      cat("\nloglik of model", logLik(fit), "\n")
      cat("\nBIC of model", AIC(fit, k = log(dim(data)[1])))
      cat("\nAIC of model", AIC(fit))
    }
    RET[j, 1] <- logLik(fit)
    RET[j, 2] <- AIC(fit, k = log(dim(data)[1]))
    RET[j, 3] <- AIC(fit)
    j <- j + 1
  }
  if (trace) 
    print(RET[complete.cases(RET), ])  #}
  return(fit)
}


#===================================================#
#                fit zero-inflated                  #
#===================================================#


fitzi = function(method, using_data, var_role, depth = 2) {
  class(method) = method
  UseMethod("fitzi", method)
}

## constant
fitzi.constant = function(method, using_data, var_role, depth = 2) {
  if (!is.null(var_role$OFF)) {
    var_role$OFF = var_role$OFF[as.integer(rownames(using_data))]
  }
  if ((!sum(using_data$Y == 0)) | var_role$ZERO == "none") {
    temp = glm(Y ~ 1, family = "poisson", data = using_data, offset = var_role$OFF, maxit = 100)
    temp$loglik = logLik(temp)
    temp
  } else {
    if (var_role$ZERO == "inflated") {
      zeroinfl(Y ~ 1, data = using_data, dist = "poisson", offset = var_role$OFF)
    } else {
      hurdle(Y ~ 1, data = using_data, dist = "poisson", offset = var_role$OFF)
    }
  }
}

## best simple linear
fitzi.simple = function(method, using_data, var_role, depth = 2) {
  if (is.null(var_role$BEST) | !length(var_role$BEST)) {
    using_FITZ = var_role$FITZ
  } else {
    using_FITZ = var_role$BEST
  }
  if (!is.null(var_role$OFF)) {
    var_role$OFF = var_role$OFF[as.integer(rownames(using_data))]
  }
  if (!length(using_FITZ)) {
    using_FITZ = "1"
  }
  using_FITL = var_role$FITL
  if (!length(using_FITL)) {
    using_FITL = "1"
  }
  if ((!sum(using_data$Y == 0)) | var_role$ZERO == "none") {
    formulas = as.formula(paste("Y ~ ", paste(using_FITL, collapse = "+")))
    formulas = Filter(Negate(function(x) is.null(unlist(x))), formulas)
    temp = glm(formulas, family = "poisson", data = using_data, offset = var_role$OFF, maxit = 100)
    temp$loglik = logLik(temp)
    temp
  } else {
    
    # formulas=as.formula(paste('Y ~ ', paste(using_FITL,collapse='+'),'|',paste(using_FITZ,collapse='+')))
    formulas = sapply(using_FITZ, function(x) as.formula(paste("Y ~ ", paste(using_FITL, collapse = "+"), "|", paste(x, collapse = "+"))))
    formulas = lapply(formulas, function(y) Filter(Negate(function(x) is.null(unlist(x))), y))
    if (depth == 1) {
      models = lapply(formulas, function(x) {
        tryCatch({
          if (var_role$ZERO == "inflated") {
            zeroinfl(x, data = using_data, dist = "poisson", offset = var_role$OFF)
          } else {
            hurdle(x, data = using_data, dist = "poisson", offset = var_role$OFF)
          }
        }, error = function(e) {
          NULL
        })
      })
      maxll = do.call(c, lapply(models, function(x) if (!is.null(x)) 
        logLik(x) else -Inf))
    } else {
      models = lapply(formulas, function(x) {
        if (var_role$ZERO == "inflated") {
          zeroinfl(x, data = using_data, dist = "poisson", offset = var_role$OFF)
        } else {
          hurdle(x, data = using_data, dist = "poisson", offset = var_role$OFF)
        }
      })
      maxll = do.call(c, lapply(models, function(x) logLik(x)))
    }
    temp2 = which(maxll == max(maxll))
    maxll = ifelse(length(temp2) > 1, sample(temp2, 1), temp2)
    models[[maxll]]
  }
}

## multiple linear
fitzi.multiple = function(method, using_data, var_role, depth = 2) {
  using_FITL = var_role$FITL
  if (!length(using_FITL)) {
    using_FITL = "1"
  }
  if (!is.null(var_role$OFF)) {
    var_role$OFF = var_role$OFF[as.integer(rownames(using_data))]
  }
  if ((!sum(using_data$Y == 0)) | var_role$ZERO == "none") {
    formulas = as.formula(paste("Y ~ ", paste(using_FITL, collapse = "+")))
    formulas = Filter(Negate(function(x) is.null(unlist(x))), formulas)
    temp = glm(formulas, family = "poisson", data = using_data, offset = var_role$OFF, maxit = 100)
    temp$loglik = logLik(temp)
    temp
  } else {
    using_FITZ = var_role$FITZ
    if (!length(using_FITZ)) {
      using_FITZ = "1"
    }
    formulas = as.formula(paste("Y ~ ", paste(using_FITL, collapse = "+"), "|", paste(using_FITZ, collapse = "+")))
    formulas = Filter(Negate(function(x) is.null(unlist(x))), formulas)
    if (var_role$ZERO == "inflated") {
      zeroinfl(formulas, data = using_data, dist = "poisson", offset = var_role$OFF)
    } else {
      hurdle(formulas, data = using_data, dist = "poisson", offset = var_role$OFF)
    }
  }
}

fitnb = function(method, using_data, var_role, depth = 2) {
  class(method) = method
  UseMethod("fitnb", method)
}

## constant
fitnb.constant = function(method, using_data, var_role, depth = 2) {
  if (!is.null(var_role$OFF)) {
    var_role$OFF = var_role$OFF[as.integer(rownames(using_data))]
  }
  if ((!sum(using_data$Y == 0)) | var_role$ZERO == "none") {
    if (!is.null(var_role$OFF)) {
      OFFS = var_role$OFF
    }
    formula1 = as.formula(paste("Y ~ 1", if (!is.null(var_role$OFF)) {
      paste("+", "offset(log(", OFFS, "))")
    } else {
      ""
    }))
    temp = glm.nb(formula1, data = using_data, maxit = 100)
    temp$loglik = logLik(temp)
    temp
  } else {
    if (var_role$ZERO == "inflated") {
      zeroinfl(Y ~ 1, data = using_data, dist = "negbin", offset = var_role$OFF)
    } else {
      hurdle(Y ~ 1, data = using_data, dist = "negbin", offset = var_role$OFF)
    }
  }
}

## best simple linear
fitnb.simple = function(method, using_data, var_role, depth = 2) {
  if (is.null(var_role$BEST) | !length(var_role$BEST)) {
    using_FITZ = var_role$FITZ
  } else {
    using_FITZ = var_role$BEST
  }
  if (!is.null(var_role$OFF)) {
    var_role$OFF = var_role$OFF[as.integer(rownames(using_data))]
  }
  if (!length(using_FITZ)) {
    using_FITZ = "1"
  }
  using_FITL = var_role$FITL
  if (!length(using_FITL)) {
    using_FITL = "1"
  }
  if ((!sum(using_data$Y == 0)) | var_role$ZERO == "none") {
    if (!is.null(var_role$OFF)) {
      OFFS = var_role$OFF
    }
    formulas = as.formula(paste("Y ~ ", paste(using_FITL, collapse = "+"), if (!is.null(var_role$OFF)) {
      paste("+", "offset(log(", OFFS, "))")
    } else {
      ""
    }))
    formulas = Filter(Negate(function(x) is.null(unlist(x))), formulas)
    temp = glm.nb(formulas, data = using_data, maxit = 100)
    temp$loglik = logLik(temp)
    temp
  } else {
    # formulas=as.formula(paste('Y ~ ', paste(using_FITL,collapse='+'),'|',paste(using_FITZ,collapse='+')))
    formulas = sapply(using_FITZ, function(x) as.formula(paste("Y ~ ", paste(using_FITL, collapse = "+"), "|", paste(x, collapse = "+"))))
    formulas = lapply(formulas, function(y) Filter(Negate(function(x) is.null(unlist(x))), y))
    if (depth == 1) {
      models = lapply(formulas, function(x) {
        tryCatch({
          if (var_role$ZERO == "inflated") {
            zeroinfl(x, data = using_data, dist = "negbin", offset = var_role$OFF)
          } else {
            hurdle(x, data = using_data, dist = "negbin", offset = var_role$OFF)
          }
        }, error = function(e) {
          NULL
        })
      })
      maxll = do.call(c, lapply(models, function(x) if (!is.null(x)) 
        logLik(x) else -Inf))
    } else {
      models = lapply(formulas, function(x) {
        if (var_role$ZERO == "inflated") {
          zeroinfl(x, data = using_data, dist = "negbin", offset = var_role$OFF)
        } else {
          hurdle(x, data = using_data, dist = "negbin", offset = var_role$OFF)
        }
      })
      maxll = do.call(c, lapply(models, function(x) logLik(x)))
    }
    temp2 = which(maxll == max(maxll))
    maxll = ifelse(length(temp2) > 1, sample(temp2, 1), temp2)
    models[[maxll]]
  }
}



## multiple linear
fitnb.multiple = function(method, using_data, var_role, depth = 2) {
  using_FITL = var_role$FITL
  if (!length(using_FITL)) {
    using_FITL = "1"
  }
  if (!is.null(var_role$OFF)) {
    var_role$OFF = var_role$OFF[as.integer(rownames(using_data))]
  }
  if ((!sum(using_data$Y == 0)) | var_role$ZERO == "none") {
    if (!is.null(var_role$OFF)) {
      OFFS = var_role$OFF
    }
    formulas = as.formula(paste("Y ~ ", paste(using_FITL, collapse = "+"), if (!is.null(var_role$OFF)) {
      paste("+", "offset(log(", OFFS, "))")
    } else {
      ""
    }))
    formulas = Filter(Negate(function(x) is.null(unlist(x))), formulas)
    temp = glm.nb(formulas, data = using_data, maxit = 100)
    temp$loglik = logLik(temp)
    temp
  } else {
    using_FITZ = var_role$FITZ
    if (!length(using_FITZ)) {
      using_FITZ = "1"
    }
    formulas = as.formula(paste("Y ~ ", paste(using_FITL, collapse = "+"), "|", paste(using_FITZ, collapse = "+")))
    formulas = Filter(Negate(function(x) is.null(unlist(x))), formulas)
    if (var_role$ZERO == "inflated") {
      zeroinfl(formulas, data = using_data, dist = "negbin", offset = var_role$OFF)
    } else {
      hurdle(formulas, data = using_data, dist = "negbin", offset = var_role$OFF)
    }
  }
}


fitzi.back = function(method, using_data, var_role, depth = 2) {
  using_FITL = var_role$FITL
  if (!length(using_FITL)) {
    using_FITL = "1"
  }
  if (!is.null(var_role$OFF)) {
    var_role$OFF = var_role$OFF[as.integer(rownames(using_data))]
  }
  if ((!sum(using_data$Y == 0)) | var_role$ZERO == "none") {
    formulas = as.formula(paste("Y ~ ", paste(using_FITL, collapse = "+")))
    formulas = Filter(Negate(function(x) is.null(unlist(x))), formulas)
    temp = glm(formulas, family = "poisson", data = using_data, offset = var_role$OFF, maxit = 100)
    temp$loglik = logLik(temp)
    temp
  } else {
    using_FITZ = var_role$FITZ
    if (!length(using_FITZ)) {
      using_FITZ = "1"
    }
    formulas = as.formula(paste("Y ~ ", paste(using_FITL, collapse = "+"), "|", paste(using_FITZ, collapse = "+")))
    formulas = Filter(Negate(function(x) is.null(unlist(x))), formulas)
    if (var_role$ZERO == "inflated") {
      temp = zeroinfl(formulas, data = using_data, dist = "poisson", offset = var_role$OFF)
    } else {
      temp = hurdle(formulas, data = using_data, dist = "poisson", offset = var_role$OFF)
    }
    backeli(temp, dist = "poisson", data = using_data, alpha = var_role$back_alpha)
  }
}


fitnb.back = function(method, using_data, var_role, depth = 2) {
  using_FITL = var_role$FITL
  if (!length(using_FITL)) {
    using_FITL = "1"
  }
  if (!is.null(var_role$OFF)) {
    var_role$OFF = var_role$OFF[as.integer(rownames(using_data))]
  }
  if ((!sum(using_data$Y == 0)) | var_role$ZERO == "none") {
    if (!is.null(var_role$OFF)) {
      OFFS = var_role$OFF
    }
    formulas = as.formula(paste("Y ~ ", paste(using_FITL, collapse = "+"), if (!is.null(var_role$OFF)) {
      paste("+", "offset(log(", OFFS, "))")
    } else {
      ""
    }))
    formulas = Filter(Negate(function(x) is.null(unlist(x))), formulas)
    temp = glm.nb(formulas, data = using_data, maxit = 100)
    temp$loglik = logLik(temp)
    temp
  } else {
    using_FITZ = var_role$FITZ
    if (!length(using_FITZ)) {
      using_FITZ = "1"
    }
    formulas = as.formula(paste("Y ~ ", paste(using_FITL, collapse = "+"), "|", paste(using_FITZ, collapse = "+")))
    formulas = Filter(Negate(function(x) is.null(unlist(x))), formulas)
    if (var_role$ZERO == "inflated") {
      temp = zeroinfl(formulas, data = using_data, dist = "negbin", offset = var_role$OFF)
    } else {
      temp = hurdle(formulas, data = using_data, dist = "negbin", offset = var_role$OFF)
    }
    backeli(temp, dist = "negbin", data = using_data, alpha = var_role$back_alpha)
  }
}


#===================================================#
#             calculate score residual              #
#===================================================#


altered_lambda = function(x) {
  x/(1 - exp(-x))
}


## revised qscore
qscore = function(method, using_data, var_role, depth, model) {
  ## fitm = switch(var_role$cdist, negbin = fitnb, poisson = fitzi)
  part_score = var_role$part_score
  if (method != "constant") {
    allvars = keep_var(using_data, unique(c(var_role$FITL, var_role$FITZ)))
    var_role$FITL = intersect(var_role$FITL, allvars)
    var_role$FITZ = intersect(var_role$FITZ, allvars)
  }
  if (is.null(model)) {
    modelset = model
  }
  modelhl = score_m = sv = mean_l = temp = list()
  for (i in 1:length(var_role$DEP)) {
    if (is.list(model)) {
      modelset = model[[i]]
    }
    modelt = if (is.null(modelset)) {
      data_now = subset(using_data, select = c(i, match(var_role$FITL, names(using_data))))
      names(data_now)[1] = "Y"
      modelz(method, data_now, var_role, depth)
    } else {
      modelset
    }
    score_mt = sandwich::estfun(modelt)
    ## message('check')
    modelhl[[i]] = modelt
    score_m[[i]] = score_mt
    score_m[[i]] = score_m[[i]][, grep("(Intercept)", dimnames(score_m[[i]])[[2]])]
    sv[[i]] = if (part_score %in% c("all", "count")) {
      if (!is.null(dim(score_m[[i]]))) 
        score_m[[i]][, 1] else score_m[[i]]
    } else {
      NULL
    }
    temp[[i]] = if (class(modelhl[[i]])[1] %in% c("glm", "negbin")) 
      score_m[[i]] else switch(part_score, all = score_m[[i]], count = score_m[[i]][, 1], zero = score_m[[i]][, 2])
    temp[[i]] = if (is.null(dim(temp[[i]]))) 
      as.integer(temp[[i]] > 0) else apply(ifelse(temp[[i]] > 0, "+", "-"), 1, paste0, collapse = "")
    ## if(var_role$trace) print(temp)
    mean_l[[i]] = if (class(modelhl[[i]])[1] %in% c("glm", "negbin")) {
      exp(mean(VGAM::predict(modelhl[[i]])))
    } else if (class(modelhl[[i]])[1] == "zeroinfl") {
      mean(VGAM::predict(modelhl[[i]], type = "count"))
    } else {
      altered_lambda(mean(VGAM::predict(modelhl[[i]], type = "count")))
    }
  }
  return(list(temp, mean_l, var_role, model = modelhl, sv))
}


#===================================================#
#            construct chi-square table             #
#===================================================#


qctable = function(var_type, ...) {
  UseMethod("qctable", var_type)
}

## for numerical variable
qctable.numeric = function(var_type, method, using_data, score_vector, var_role, var_name, score_diff, t_index) {
  temp = ftable(table(score_vector[t_index], cut(var_type[t_index], breaks = c(-Inf, unique(quantile(var_type[t_index], na.rm = T, probs = c(1/4, 2/4, 3/4))), Inf))))
  return(temp[, colSums(temp) != 0])
}

## for factor variable
qctable.factor = function(var_type, method, using_data, score_vector, var_role, var_name, score_diff, t_index) {
  temp = ftable(table(score_vector[t_index], var_type[t_index]))
  return(temp[, colSums(temp) != 0])
}

## for interaction
qctable_interaction = function(vars_type, method, using_data, score_vector, var_role, score_diff, t_index) {
  all_interactions <- t(combn(var_role$CUT, 2))
  all_temps = lapply(1:nrow(all_interactions), function(i) {
    temps = as.list(using_data[all_interactions[i, ]])
    table_type = do.call(sum, (lapply(temps, is.factor)))
    switch(table_type + 1, ftable(table(cut(temps[[1]][t_index], breaks = c(-Inf, unique(quantile(temps[[1]][t_index], na.rm = T)[3]), Inf)), cut(temps[[2]][t_index], breaks = c(-Inf, unique(quantile(temps[[2]][t_index], na.rm = T)[3]), Inf)), score_vector[t_index])), ftable(table(cut(temps[[which(lapply(temps, is.factor) == F)]][t_index], breaks = c(-Inf, unique(quantile(temps[[which(lapply(temps, is.factor) == F)]][t_index], na.rm = T)[3]), Inf)), temps[[which(lapply(temps, is.factor) == T)]][t_index], score_vector[t_index])), ftable(table(temps[[1]][t_index], temps[[2]][t_index], score_vector[t_index])))
  })
  all_temps = lapply(all_temps, function(i) {
    i[rowSums(i) != 0, ]
  })
  return(list(table_names = all_interactions, interaction_tables = all_temps))
}


#===================================================#
#                variable selection                 #
#===================================================#
## revised qvariable_select


qvariable_select = function(method, using_data, var_role, alpha_main = 0.05, alpha_in = 0.05, depth = 1, model = model, min_ssize = NULL) {
  gama = 1
  var_role$CUT = keep_var(using_data, var_role$CUT)
  vars_type = as.list(using_data[var_role$CUT])
  sv = qscore(method, using_data, var_role, depth, model)
  score_vector = score_diff = model_now = list()
  for (i in 1:length(var_role$DEP)) {
    score_vector[[i]] = sv[[1]][[i]]
    score_diff[[i]] = sv[[2]][[i]]
    model_now[[i]] = sv[[4]][[i]]
  }
  var_role = sv[[3]]
  ## print(table(score_vector)
  
  if (length(var_role$CUT) == 0) {
    return(list(bs = NULL, cut_variable = NULL, terminal = T, why_no_split = "No variable to use", gama = gama, cts = NULL, model = model_now))
  }
  if (!is.null(min_ssize)) {
    if (dim(using_data)[1] <= (min_ssize * 2)) {
      return(list(bs = NULL, cut_variable = NULL, terminal = T, why_no_split = paste("Observations is NOT enough (minsize = ", min_ssize, ")"), gama = gama, cts = NULL, model = model_now))
    }
  }
  t_index = if (length(var_role$DEP) == 1) {
    !logical(length(using_data$Y))
  } else {
    !logical(length(as.numeric(unlist(using_data[var_role$DEP[1]]))))
  }
  ## print(all(t_index == T))
  table_single = list()
  for (j in 1:length(var_role$DEP)) {
    table_single[[j]] = lapply(names(vars_type), function(i) {
      qctable(vars_type[[i]], method, using_data, score_vector[[j]], var_role, i, score_diff, t_index)
    })
    names(table_single[[j]]) = names(vars_type)
  }
  chiD = dfD = NULL
  for (i in 1:length(var_role$CUT)) {
    if (!is.factor(vars_type[i])) {
      CT3_Data = NULL
      for (j in 1:length(var_role$DEP)) {
        P_CT3_Data = as.vector(table_single[[j]][[i]])
        CT3_Data = c(CT3_Data, P_CT3_Data)
      }
      Quasi = rep(1, length(CT3_Data))
      Zero_Cell = which(CT3_Data == 0)
      if (length(Zero_Cell) != 0) 
        Quasi[Zero_Cell] = 0
      dcc = dim(table_single[[j]][[i]])[1]
      dcc = ifelse(is.null(dcc), 1, dcc)
      ccc = length(CT3_Data)/(dcc * length(var_role$DEP))
      CT3 = array(CT3_Data, dim = c(dcc, ccc, length(var_role$DEP)), dimnames = list(c(as.character(1:dcc)), c(as.character(1:ccc)), c(as.character(1:length(var_role$DEP)))))
      QCT3 = array(Quasi, dim = c(dcc, ccc, length(var_role$DEP)), dimnames = list(c(as.character(1:dcc)), c(as.character(1:ccc)), c(as.character(1:length(var_role$DEP)))))
      ## Loglinear Model(Full Model)
      xy.xz.yz.model = loglin(CT3, list(c(1, 2), c(1, 3), c(2, 3)), start = QCT3, print = F)  
      ## Loglinear Model(Reduced Model)
      xz.yz.model = loglin(CT3, list(c(1, 3), c(2, 3)), start = QCT3, print = F)
      chiD[i] = xz.yz.model$pearson - xy.xz.yz.model$pearson
      dfD[i] = xz.yz.model$df - xy.xz.yz.model$df
    } else {
      CT3_Data = NULL
      for (j in 1:length(var_role$DEP)) {
        P_CT3_Data = as.vector(table_single[[j]][[i]])
        CT3_Data = c(CT3_Data, P_CT3_Data)
      }
      Quasi = rep(1, length(CT3_Data))
      Zero_Cell = which(CT3_Data == 0)
      if (length(Zero_Cell) != 0) 
        Quasi[Zero_Cell] = 0
      dcn = dim(table_single[[j]][[i]])[1]
      dcn = ifelse(is.null(dcn), 1, dcn)
      cn = length(CT3_Data)/(dcn * length(var_role$DEP))
      CT3 = array(CT3_Data, dim = c(dcn, cn, length(var_role$DEP)), dimnames = list(c(as.character(1:dcn)), as.character(1:cn), c(as.character(1:length(var_role$DEP)))))
      QCT3 = array(Quasi, dim = c(dcn, cn, length(var_role$DEP)), dimnames = list(c(as.character(1:dcn)), as.character(1:cn), c(as.character(1:length(var_role$DEP)))))
      ## Loglinear Model(Full Model)
      xy.xz.yz.model = loglin(CT3, list(c(1, 2), c(1, 3), c(2, 3)), start = QCT3, print = F)
      ## Loglinear Model(Reduced Model)
      xz.yz.model = loglin(CT3, list(c(1, 3), c(2, 3)), start = QCT3, print = F)  
      chiD[i] = xz.yz.model$pearson - xy.xz.yz.model$pearson
      dfD[i] = xz.yz.model$df - xy.xz.yz.model$df
    }
  }
  stat_single = list()
  ## stat_single[[1]] = list(-0.4317894,2)
  for (i in 1:length(var_role$CUT)) {
    stat_single[[i]] = list(chiD[i], dfD[i])
  }
  stat_single = lapply(stat_single, function(i) {
    suppressWarnings(max(0, ifelse(is.nan((7/9 + sqrt(as.numeric(i[2])) * ((as.numeric(i[1])/as.numeric(i[2]))^(1/3) - 1 + (2/(9 * as.numeric(i[2])))))^3), 0, (7/9 + sqrt(as.numeric(i[2])) * ((as.numeric(i[1])/as.numeric(i[2]))^(1/3) - 1 + (2/(9 * as.numeric(i[2])))))^3)))
  })
  names(stat_single) = names(vars_type)
  stat_single = stat_single[!unlist(lapply(stat_single, is.nan))]
  
  if (!is.null(unlist(stat_single))) {
    best_single = names(stat_single[stat_single == max(unlist(stat_single))])
    if (length(best_single) > 1) 
      best_single = sample(best_single, 1, replace = F)
    
    if (stat_single[[best_single]] > qchisq(1 - (alpha_main/length(stat_single)), 1)) {
      cut_variable = best_single
      terminal_node = F
      nosplit = NULL
    } else {
      
      ## if main effect not significant
      
      if (method != "simple") {
        if (length(var_role$CUT) > 1) {
          ## interaction
          table_couple = list()
          for (i in 1:length(var_role$DEP)) {
            table_couple[[i]] = qctable_interaction(vars_type, method, using_data, score_vector[[i]], var_role, score_diff, t_index)
          }
          chiD = dfD = NULL
          tcn = nrow(t(combn(var_role$CUT, 2)))
          for (i in 1:tcn) {
            CT3_Data = NULL
            for (j in 1:length(var_role$DEP)) {
              P_CT3_Data = as.vector(table_couple[[j]]$interaction_tables[[i]])
              CT3_Data = c(CT3_Data, P_CT3_Data)
            }
            Quasi = rep(1, length(CT3_Data))
            Zero_Cell = which(CT3_Data == 0)
            if (length(Zero_Cell) != 0) 
              Quasi[Zero_Cell] = 0
            dcc = dim(table_couple[[j]]$interaction_tables[[i]])[1]
            dcc = ifelse(is.null(dcc), 1, dcc)
            ccc = length(CT3_Data)/(dcc * length(var_role$DEP))
            CT3 = array(CT3_Data, dim = c(dcc, ccc, length(var_role$DEP)), dimnames = list(c(as.character(1:dcc)), c(as.character(1:ccc)), c(as.character(1:length(var_role$DEP)))))
            QCT3 = array(Quasi, dim = c(dcc, ccc, length(var_role$DEP)), dimnames = list(c(as.character(1:dcc)), c(as.character(1:ccc)), c(as.character(1:length(var_role$DEP)))))
            ## Loglinear Model(Full Model)
            xy.xz.yz.model = loglin(CT3, list(c(1, 2), c(1, 3), c(2, 3)), start = QCT3, print = F)
            ## Loglinear Model(Reduced Model)
            xz.yz.model = loglin(CT3, list(c(1, 3), c(2, 3)), start = QCT3, print = F)
            chiD[i] = xz.yz.model$pearson - xy.xz.yz.model$pearson
            dfD[i] = xz.yz.model$df - xy.xz.yz.model$df
          }
          stat_couple = list()
          for (i in 1:tcn) {
            stat_couple[[i]] = list(chiD[i], dfD[i])
          }
          stat_couple = lapply(stat_couple, function(i) {
            suppressWarnings(max(0, ifelse(is.nan((7/9 + sqrt(as.numeric(i[2])) * ((as.numeric(i[1])/as.numeric(i[2]))^(1/3) - 1 + (2/(9 * as.numeric(i[2])))))^3), 0, (7/9 + sqrt(as.numeric(i[2])) * ((as.numeric(i[1])/as.numeric(i[2]))^(1/3) - 1 + (2/(9 * as.numeric(i[2])))))^3)))
          })
          
          table_couple[[1]]$table_names = table_couple[[1]]$table_names[!unlist(lapply(stat_couple, is.nan)), ]
          stat_couple = stat_couple[!unlist(lapply(stat_couple, is.nan))]
          single_exist = apply(matrix(apply(matrix(table_couple[[1]]$table_names, ncol = 2), 2, function(x) match(x, names(stat_single), nomatch = 0) > 0), ncol = 2), 1, sum) > 1
          table_couple[[1]]$table_names = matrix(table_couple[[1]]$table_names, ncol = 2)[single_exist, ]
          stat_couple = stat_couple[single_exist]
          if (!is.null(unlist(stat_couple))) {
            best_couple = matrix(table_couple[[1]]$table_names, ncol = 2)[stat_couple == max(unlist(stat_couple)), ]
            best_couple = matrix(best_couple, ncol = 2)
            if (nrow(best_couple) > 1) {
              bc_stat = stat_single[best_couple]
              bc_name = as.matrix(names(bc_stat[bc_stat == do.call(max, bc_stat)]))
              couples = apply(as.matrix(apply(best_couple, 2, function(x) match(x, bc_name, nomatch = 0) > 0)), 1, sum)
              if (sum(couples) > 0) 
                best_couple = best_couple[(couples > 0), ]
              if (length(best_couple) > 2) 
                best_couple = best_couple[sample(1:(length(best_couple)/2), 1, replace = F), ]
            }
            best_csingle = names(which(max(unlist(stat_single[best_couple])) == unlist(stat_single[best_couple])))
            if (length(best_csingle) > 1) 
              best_csingle = sample(best_csingle, 1, replace = F)
            MAX = !is.na(colSums(apply(matrix(table_couple[[1]]$table_names, ncol = 2), 1, function(x) match(best_couple, x))))
            if (stat_couple[MAX][[1]] > qchisq(1 - (alpha_in/length(stat_couple)), 1)) {
              cut_variable = best_csingle
              terminal_node = F
              nosplit = NULL
            } else {
              cut_variable = best_single
              terminal_node = T
              nosplit = "both_non_significant"
            }
          } else {
            cut_variable = best_single
            terminal_node = F
            nosplit = "non_couple_terms"
          }
          # length>2
        } else {
          cut_variable = best_single
          terminal_node = T
          nosplit = "single_non_significant"
        }
        # simple
      } else {
        cut_variable = best_single
        terminal_node = T
        nosplit = "single_non_significant"
      }
    }  #significant
    # nosingle
  } else {
    cut_variable = NULL
    terminal_node = T
    nosplit = "chi2_error"
  }
  if (depth == 1 && method == "multiple") {
    return(list(bs = best_single, cut_variable = cut_variable, terminal = terminal_node, why_no_split = nosplit, gama = gama, cts = NULL, model = model_now))
  }
  
  return(list(bs = best_single, cut_variable = cut_variable, terminal = terminal_node, why_no_split = nosplit, gama = gama, cts = NULL, model = model_now))
}


#===================================================#
#                  point selection                  #
#===================================================#


mean_pzero = function(p) {
  temp = table(p)
  sum(temp * as.numeric(names(temp)))/length(p)
}


v_mode <- function(x) {
  x = na.omit(x)
  t = table(x)
  m = max(t)
  out = names(t)[t == m]
  if (length(out)) 
    out = sample(out, 1)
  return(out)
}


## How is data fitted for each split
modelz <- function(method, model_data, var_role, depth) {
  if (dim(model_data)[1] == 0) {
    NULL
  } else {
    tryCatch({
      fitm = switch(var_role$cdist, negbin = fitnb, poisson = fitzi)
      allvars = keep_var(model_data, unique(c(var_role$FITL, var_role$FITZ)))
      var_role$FITL = intersect(var_role$FITL, allvars)
      var_role$FITZ = intersect(var_role$FITZ, allvars)
      fitm(method, model_data, var_role, depth)
    }, warning = function(w) {
      tryCatch({
        var_role$ZERO = "none"
        fitm(method, model_data, var_role, depth)
      }, warning = function(w) {
        NULL
      }, error = function(e) {
        NULL
      })
    }, error = function(e) {
      tryCatch({
        var_role$ZERO = "none"
        fitm(method, model_data, var_role, depth)
      }, warning = function(w) {
        NULL
      }, error = function(e) {
        NULL
      })
    })
  }
}

## Revised point.selection


point.selection <- function(method, using_data, var_role, alpha_main = 0.05, alpha_in = 0.05, depth = 1, min_size = min_size, loglik_increment = 0, model = NULL) {
  naornot = sum(is.na(using_data))
  if (naornot) {
    o_data = using_data
    comcs = complete.cases(using_data)
    using_data = using_data[comcs, ]
  }
  
  temp_var = qvariable_select(method, using_data, var_role, alpha_main, alpha_in, depth, model, min_size)
  cut_variable = temp_var$cut_variable
  terminal_node = temp_var$terminal
  
  LNA = NULL
  
  if (!is.null(cut_variable)) 
    cut_one = using_data[[cut_variable]]
  
  ## Model_now = modelz(method, using_data, var_role)
  Model_now = list()
  Model_O = 0
  for (i in 1:length(var_role$DEP)) {
    Model_now[[i]] = temp_var$model[[i]]
    Model_O = Model_O + Model_now[[i]]$loglik * (-1)
  }
  if (!terminal_node) {
    best_var = var_role$BEST
    var_role$BEST = NULL
    subset_points = if (is.factor(cut_one)) {
      length2 <- if (length(unique(cut_one))%%2 == 1) {
        floor(length(unique(cut_one))%/%2)
      } else {
        length(unique(cut_one))%/%2
      }
      Achieve <- categorical.possible.combination <- c()
      for (j in 1:length2) {
        # temp = t(combn(sort(unique(cut_one)), j))
        categorical.possible.combination <- lapply(1:(length(t(combn(unique(cut_one), j)))/j), function(i) t(combn(sort(unique(cut_one)), j))[i, ])
        Achieve <- c(Achieve, categorical.possible.combination)
      }
      Achieve
    } else {
      sort(unique(cut_one))
    }
    
    model_impurity = lapply(subset_points, function(x) {
      if (is.factor(cut_one)) {
        cut.subset.L = subset(using_data, cut_one %in% x)
        cut.subset.R = subset(using_data, cut_one %out% x)
      } else {
        cut.subset.L = subset(using_data, cut_one <= x)
        cut.subset.R = subset(using_data, cut_one > x)
      }
      work_index_L = if (dim(cut.subset.L)[1] >= min_size) {
        T
      } else {
        F
      }
      work_index_R = if (dim(cut.subset.R)[1] >= min_size) {
        T
      } else {
        F
      }
      work_index = work_index_L & work_index_R
      ## print(work_index)
      if (work_index) {
        Model_LL = Model_RR = oModel_L = oModel_R = NULL
        for (i in 1:length(var_role$DEP)) {
          Y = cut.subset.L[, i]
          Model_L = modelz(method, cbind(Y, cut.subset.L[, -(1:length(var_role$DEP))]), var_role, depth + 1)
          oModel_L[[i]] = Model_L
          L_var = if (!is.null(Model_L)) {
            attr(Model_L$terms$count, "term.labels")
          } else {
            NULL
          }
          Model_LL[i] = if (is.null(Model_L$loglik)) {
            Inf
          } else {
            Model_L$loglik * (-1)
          }
          Y = cut.subset.R[, i]
          Model_R = modelz(method, cbind(Y, cut.subset.R[, -(1:length(var_role$DEP))]), var_role, depth + 1)
          oModel_R[[i]] = Model_R
          R_var = if (!is.null(Model_R)) {
            attr(Model_R$terms$count, "term.labels")
          } else {
            NULL
          }
          Model_RR[i] = if (is.null(Model_R$loglik)) {
            Inf
          } else {
            Model_R$loglik * (-1)
          }
        }
        Model_L = sum(Model_LL)
        Model_R = sum(Model_RR)
        list(Model_O - Model_L - Model_R, L_var[1], R_var[1], oModel_L, oModel_R)
      } else {
        NULL
      }
    })
    
    max_split = unlist(lapply(model_impurity, function(x) if (is.null(x)) {
      NULL
    } else {
      x[[1]]
    }))
    max_split = suppressWarnings(max(max_split))
    temp_max = max_split
    decrease_ll = max_split > loglik_increment
    ## print(max_split) ;print(loglik_increment)
    max_split = which(unlist(lapply(model_impurity, function(x) if (is.null(x)) {
      F
    } else {
      x[[1]] == max_split
    })))
    
    if (length(max_split) != 0 & decrease_ll) {
      if (length(max_split) > 1) 
        max_split = max_split[1]
      cut_value = if (is.factor(cut_one)) {
        subset_points[[max_split]]
      } else {
        mean(as.numeric(levels(factor(cut_one))[max_split:(max_split + 1)]))
      }
      L_var = model_impurity[[max_split]][[2]]
      R_var = model_impurity[[max_split]][[3]]
      L_model = model_impurity[[max_split]][[4]]
      R_model = model_impurity[[max_split]][[5]]
      ## print(L_model) ;print(R_model)
      if (is.factor(cut_one)) {
        cut.subset.L = subset(using_data, cut_one %in% cut_value)
        cut.subset.R = subset(using_data, cut_one %out% cut_value)
      } else {
        cut.subset.L = subset(using_data, cut_one <= cut_value)
        cut.subset.R = subset(using_data, cut_one > cut_value)
      }
      Why_no_split = NULL
      
      if (naornot) {
        cut_NA = is.na(o_data[[cut_variable]])
        if (is.factor(o_data[[cut_variable]])) {
          NAset = o_data[(!comcs) & (!cut_NA), ]
          if (NROW(NAset) & (!is.null(cut_value))) {
            NAindex = (NAset[[cut_variable]] %in% cut_value)
            cut.subset.L = rbind(cut.subset.L, NAset[NAindex, ])
            cut.subset.R = rbind(cut.subset.R, NAset[!NAindex, ])
          }
          if (NROW(o_data[which(cut_NA), ])) {
            na_cut = v_mode(o_data[[cut_variable]])
            if (na_cut %in% cut_value) {
              cut.subset.L = rbind(cut.subset.L, o_data[which(cut_NA), ])
              LNA = T
            } else {
              cut.subset.R = rbind(cut.subset.R, o_data[which(cut_NA), ])
              LNA = F
            }
          }
        } else {
          NAset = o_data[(!comcs) & (!cut_NA), ]
          if (NROW(NAset) & (!is.null(cut_value))) {
            NAindex = (NAset[[cut_variable]] <= cut_value)
            cut.subset.L = rbind(cut.subset.L, NAset[NAindex, ])
            cut.subset.R = rbind(cut.subset.R, NAset[!NAindex, ])
          }
          if (NROW(o_data[which(cut_NA), ])) {
            na_cut = mean(using_data[[cut_variable]], na.rm = T)
            if (na_cut <= cut_value) {
              cut.subset.L = rbind(cut.subset.L, o_data[which(cut_NA), ])
              LNA = T
            } else {
              cut.subset.R = rbind(cut.subset.R, o_data[which(cut_NA), ])
              LNA = F
            }
          }
        }
      }
      
    } else {
      cut_value = NULL
      L_var = NULL
      R_var = NULL
      terminal_node = T
      cut.subset.L = NULL
      cut.subset.R = NULL
      L_model = NULL
      R_model = NULL
      Why_no_split = if (length(max_split) == 0) 
        "node size too small" else if (temp_max == -Inf) 
          "model can NOT be fitted on child node" else "loglik decrease too little"
    }
    
    
    
    if (naornot) {
      var_role$BEST = best_var
      Model_now = modelz(method, replace_NA(o_data), var_role, depth)
      using_data = o_data
    }
    Y_bar = H_lambda = H_pzero = H_n0 = list()
    for (i in 1:length(var_role$DEP)) {
      Y_bar[[i]] = mean(using_data[, i])
      H_lambda[[i]] = if (class(Model_now[[i]])[1] == "zeroinfl" | class(Model_now[[i]])[1] == "hurdle") {
        mean(VGAM::predict(Model_now[[i]], type = "count"))
      } else {
        mean(exp(VGAM::predict(Model_now[[i]])))
      }
      
      H_pzero[[i]] = if (class(Model_now[[i]])[1] == "zeroinfl") {
        mean_pzero(VGAM::predict(Model_now[[i]], type = "zero"))
      } else if (class(Model_now[[i]])[1] == "hurdle") {
        if (Model_now[[i]]$dist[[1]] != "negbin") {
          mean_pzero(1 - (VGAM::predict(Model_now[[i]], type = "zero") * ppois(0, lambda = VGAM::predict(Model_now[[i]], type = "count"), lower.tail = FALSE)))
        } else {
          mean_pzero(1 - (VGAM::predict(Model_now[[i]], type = "zero") * pnbinom(0, mu = VGAM::predict(Model_now[[i]], type = "count"), size = Model_now[[i]]$theta, lower.tail = FALSE)))
        }
      } else {
        NA
      }
      H_n0[[i]] = if (class(Model_now[[i]])[1] == "zeroinfl") {
        if (Model_now[[i]]$dist[[1]] != "negbin") {
          sum(VGAM::predict(Model_now[[i]], type = "zero") + (1 - VGAM::predict(Model_now[[i]], type = "zero")) * ppois(0, VGAM::predict(Model_now[[i]], type = "count")))
        } else {
          sum(VGAM::predict(Model_now[[i]], type = "zero") + (1 - VGAM::predict(Model_now[[i]], type = "zero")) * pnbinom(0, mu = VGAM::predict(Model_now[[i]], type = "count"), size = Model_now[[i]]$theta))
        }
      } else if (class(Model_now[[i]])[1] == "hurdle") {
        if (Model_now[[i]]$dist[[1]] != "negbin") {
          sum(1 - (VGAM::predict(Model_now[[i]], type = "zero") * ppois(0, lambda = VGAM::predict(Model_now[[i]], type = "count"), lower.tail = FALSE)))
        } else {
          sum(1 - (VGAM::predict(Model_now[[i]], type = "zero") * pnbinom(0, mu = VGAM::predict(Model_now[[i]], type = "count"), size = Model_now[[i]]$theta, lower.tail = FALSE)))
        }
      } else if (class(Model_now[[i]])[1] != "negbin") {
        sum(ppois(0, exp(VGAM::predict(Model_now[[i]]))))
      } else {
        sum(pnbinom(0, mu = exp(VGAM::predict(Model_now[[i]])), size = Model_now[[i]]$theta))
      }
    }
    return_value = list(Cut_variable = cut_variable, Cut_point = cut_value, Terminal = terminal_node, Why_no_split = Why_no_split, Model = Model_now, R_t = Model_O, Node_size = nrow(using_data), Y_bar = Y_bar, H_lambda = H_lambda, H_pzero = H_pzero, H_n0 = H_n0, L_data = cut.subset.L, R_data = cut.subset.R, L_model = L_model, R_model = R_model, L_NA = if (is.null(LNA)) {
      NULL
    } else {
      LNA
    }, L_var = L_var, R_var = R_var)
    
  } else {
    if (naornot) {
      Model_now = modelz(method, replace_NA(o_data), var_role, depth)
      using_data = o_data
    }
    
    Y_bar = H_lambda = H_pzero = H_n0 = list()
    for (i in 1:length(var_role$DEP)) {
      Y_bar[[i]] = mean(using_data[, i])
      H_lambda[[i]] = if (class(Model_now[[i]])[1] == "zeroinfl" | class(Model_now[[i]])[1] == "hurdle") {
        mean(VGAM::predict(Model_now[[i]], type = "count"))
      } else {
        mean(exp(VGAM::predict(Model_now[[i]])))
      }
      
      H_pzero[[i]] = if (class(Model_now[[i]])[1] == "zeroinfl") {
        mean_pzero(VGAM::predict(Model_now[[i]], type = "zero"))
      } else if (class(Model_now[[i]])[1] == "hurdle") {
        if (Model_now[[i]]$dist[[1]] != "negbin") {
          mean_pzero(1 - (VGAM::predict(Model_now[[i]], type = "zero") * ppois(0, lambda = VGAM::predict(Model_now[[i]], type = "count"), lower.tail = FALSE)))
        } else {
          mean_pzero(1 - (VGAM::predict(Model_now[[i]], type = "zero") * pnbinom(0, mu = VGAM::predict(Model_now[[i]], type = "count"), size = Model_now[[i]]$theta, lower.tail = FALSE)))
        }
      } else {
        NA
      }
      H_n0[[i]] = if (class(Model_now[[i]])[1] == "zeroinfl") {
        if (Model_now[[i]]$dist[[1]] != "negbin") {
          sum(VGAM::predict(Model_now[[i]], type = "zero") + (1 - VGAM::predict(Model_now[[i]], type = "zero")) * ppois(0, VGAM::predict(Model_now[[i]], type = "count")))
        } else {
          sum(VGAM::predict(Model_now[[i]], type = "zero") + (1 - VGAM::predict(Model_now[[i]], type = "zero")) * pnbinom(0, mu = VGAM::predict(Model_now[[i]], type = "count"), size = Model_now[[i]]$theta))
        }
      } else if (class(Model_now[[i]])[1] == "hurdle") {
        if (Model_now[[i]]$dist[[1]] != "negbin") {
          sum(1 - (VGAM::predict(Model_now[[i]], type = "zero") * ppois(0, lambda = VGAM::predict(Model_now[[i]], type = "count"), lower.tail = FALSE)))
        } else {
          sum(1 - (VGAM::predict(Model_now[[i]], type = "zero") * pnbinom(0, mu = VGAM::predict(Model_now[[i]], type = "count"), size = Model_now[[i]]$theta, lower.tail = FALSE)))
        }
      } else if (class(Model_now[[i]])[1] != "negbin") {
        sum(ppois(0, exp(VGAM::predict(Model_now[[i]]))))
      } else {
        sum(pnbinom(0, mu = exp(VGAM::predict(Model_now[[i]])), size = Model_now[[i]]$theta))
      }
    }
    
    return_value = list(Cut_variable = cut_variable, Cut_point = NULL, Terminal = terminal_node, Why_no_split = temp_var$why_no_split, Model = Model_now, R_t = Model_O, Node_size = nrow(using_data), Y_bar = Y_bar, H_lambda = H_lambda, H_pzero = H_pzero, H_n0 = H_n0, L_data = NULL, R_data = NULL, L_model = NULL, R_model = NULL, L_NA = if (is.null(LNA)) {
      NULL
    } else {
      LNA
    }, L_var = NULL, R_var = NULL)
  }
  return(return_value)
}


#===================================================#
#                  build the tree                   #
#===================================================#
## Revised create_tree


create_tree = function(method, using_data, var_role, alpha_main = 0.05, alpha_in = 0.05, depth = 1, min_size, loglik_increment = 0, model = NULL, parent_node = new.env()) {
  ##  message(' size:', min_size)
  temp = new.env(parent = parent_node)
  temp$Node_info = point.selection(method, using_data, var_role, alpha_main, alpha_in, depth, min_size, loglik_increment, model = model)
  temp$LR = NULL
  if (!temp$Node_info$Terminal) {
    LVR = RVR = var_role
    allvars = keep_var(temp$Node_info$L_data, unique(c(LVR$FITL, LVR$FITZ, LVR$CUT)))
    LVR$FITL = intersect(LVR$FITL, allvars)
    LVR$FITZ = intersect(LVR$FITZ, allvars)
    LVR$CUT = intersect(LVR$CUT, allvars)
    LVR$BEST = temp$Node_info$L_var
    temp$L = create_tree(method, temp$Node_info$L_data, LVR, alpha_main, alpha_in, depth + 1, min_size, loglik_increment, model = temp$Node_info$L_model, parent_node = temp)
    temp$L$LR = T
    allvars = keep_var(temp$Node_info$R_data, unique(c(RVR$FITL, RVR$FITZ, RVR$CUT)))
    RVR$FITL = intersect(RVR$FITL, allvars)
    RVR$FITZ = intersect(RVR$FITZ, allvars)
    RVR$CUT = intersect(RVR$CUT, allvars)
    RVR$BEST = temp$Node_info$R_var
    temp$R = create_tree(method, temp$Node_info$R_data, RVR, alpha_main, alpha_in, depth + 1, min_size, loglik_increment, model = temp$Node_info$R_model, parent_node = temp)
    temp$R$LR = F
    temp$OP = T
  } else {
    temp$L = temp$R = NULL
    temp$OP = F
  }
  temp$NP = temp$OP
  return(temp)
}


#===================================================#
#              level order of the tree              #
#===================================================#


level_order = function(root) {
  result = list()
  temp = list()
  temp[!sapply(temp, is.null)]
  temp = c(temp, root)
  while (length(temp)) {
    root = temp[[1]]
    temp = temp[-1]
    if (!is.null(root)) {
      result = c(result, root)
      if (!is.null(root$L)) 
        temp = c(temp, root$L)
      if (!is.null(root$R)) 
        temp = c(temp, root$R)
    } else break
  }
  if (is.null(result[[1]]$hash)) {
    sapply(1:length(result), function(i) result[[i]]$hash = i)
  }
  result
}


#===================================================#
#                preorder of the tree               #
#===================================================#


pre_order = function(root, k = 0) {
  result = c()
  if (!is.null(root)) {
    result = c(result, root$hash)
    cat(rep(" ", k), "Node", root$hash, ": ")
    if ((root$NP)) {
      if (!is.null(root$Node_info$L_NA)) {
        if (root$Node_info$L_NA) {
          if (is.factor(root$Node_info$Cut_point)) {
            cat("Discriminant : ", root$Node_info$Cut_variable, " belongs to ", levels(root$Node_info$Cut_point)[root$Node_info$Cut_point], " with NA values", "\n")
          } else {
            cat("Discriminant : ", root$Node_info$Cut_variable, " <= ", root$Node_info$Cut_point, " with NA values", "\n")
          }
        } else {
          if (is.factor(root$Node_info$Cut_point)) {
            cat("Discriminant : ", root$Node_info$Cut_variable, " belongs to ", levels(root$Node_info$Cut_point)[root$Node_info$Cut_point], " with no NA values", "\n")
          } else {
            cat("Discriminant : ", root$Node_info$Cut_variable, " <= ", root$Node_info$Cut_point, " with no NA values", "\n")
          }
        }
      } else {
        if (is.factor(root$Node_info$Cut_point)) {
          cat("Discriminant : ", root$Node_info$Cut_variable, " belongs to ", levels(root$Node_info$Cut_point)[root$Node_info$Cut_point], "\n")
        } else {
          cat("Discriminant : ", root$Node_info$Cut_variable, " <= ", root$Node_info$Cut_point, "\n")
        }
      }
      pre_order(root$L, k + 1)
      pre_order(root$R, k + 1)
    } else {
      cat(" ", "Size=", root$Node_info$Node_size, ", mean=", paste0(unlist(root$Node_info$Y_bar), ","), ", pzero=", paste0(unlist(root$Node_info$H_pzero), ","), ", lambda=", paste0(unlist(root$Node_info$H_lambda), ","), "\n")
    }
  }
}


#===================================================#
#               construct prune table               #
#===================================================#


prune_one = function(atree) {
  temp = level_order(atree)
  P_table = as.data.frame(matrix(0, nrow = length(temp), ncol = 9, dimnames = list(paste("Node", 1:length(temp)), c("t", "N", "S", "R", "G", "g", "L_N", "R_N", "P_N"))))
  P_table$t = sapply(1:length(temp), function(i) temp[[i]]$hash = i)
  P_table$L_N = sapply(P_table$t, function(i) if (!is.null(temp[[i]]$L)) {
    temp[[i]]$L$hash
  } else {
    0
  })
  P_table$R_N = sapply(P_table$t, function(i) if (!is.null(temp[[i]]$R)) {
    temp[[i]]$R$hash
  } else {
    0
  })
  P_table$P_N = sapply(P_table$t, function(i) if (!is.null(parent.env(temp[[i]])$hash)) {
    parent.env(temp[[i]])$hash
  } else {
    0
  })
  P_table$R = sapply(P_table$t, function(i) temp[[i]]$Node_info$R_t)
  
  for (i in rev(P_table$t)) {
    if (P_table$L_N[i] < 1) {
      P_table$N[i] = 1
      P_table$S[i] = P_table$R[i]
      P_table$G[i] = Inf
    } else {
      P_table$N[i] = P_table$N[P_table$L_N[i]] + P_table$N[P_table$R_N[i]]
      P_table$S[i] = P_table$S[P_table$L_N[i]] + P_table$S[P_table$R_N[i]]
      P_table$g[i] = (P_table$R[i] - P_table$S[i])/(P_table$N[i] - 1)
      P_table$G[i] = min(P_table$g[i], P_table$G[P_table$L_N[i]], P_table$G[P_table$R_N[i]])
    }
  }
  result1 = P_table$N
  result2 = c()
  k = 1
  alpha = 0
  epi = 0.1
  repeat {
    if (P_table$G[1] > (alpha + epi)) {
      result2 = rbind(result2, c(k, P_table$N[1], alpha, P_table$S[1]))
      alpha = P_table$G[1]
      k = k + 1
    }
    if (P_table$N[1] == 1) 
      break
    t = 1
    while (P_table$G[t] < (P_table$g[t] - epi)) {
      if (P_table$G[t] == P_table$G[P_table$L_N[t]]) {
        t = P_table$L_N[t]
      } else {
        t = P_table$R_N[t]
      }
    }
    P_table$N[t] = 1
    P_table$S[t] = P_table$R[t]
    P_table$G[t] = Inf
    while (t > 1) {
      t = P_table$P_N[t]
      P_table$N[t] = P_table$N[P_table$L_N[t]] + P_table$N[P_table$R_N[t]]
      P_table$S[t] = P_table$S[P_table$L_N[t]] + P_table$S[P_table$R_N[t]]
      P_table$g[t] = (P_table$R[t] - P_table$S[t])/(P_table$N[t] - 1)
      P_table$G[t] = min(P_table$g[t], P_table$G[P_table$L_N[t]], P_table$G[P_table$R_N[t]])
    }
  }
  result2 = data.frame(cbind(result2, sqrt(result2[, 3] * c(result2[, 3][-1], Inf))))
  names(result2) = c("k", "N_TNode", "alpha", "R_Tk", "g_alpha")
  P_table$N = result1
  return(list(prune_tree = P_table[, -5], prune_alpha = result2))
}


#===================================================#
#             K-fold cross validation               #
#===================================================#


K_fold_CV = function(K = 10, main_tree = NULL, method, using_data, var_role, alpha_main = 0.05, alpha_in = 0.05, SE = 0.5, min_size, loglik_increment = 0) {
  # message(' K size:',min_size)
  if (is.null(main_tree)) 
    main_tree = create_tree(method, using_data, var_role, alpha_main, alpha_in, depth = 1, min_size, loglik_increment, parent_node = new.env())
  result1 = prune_one(main_tree)$prune_alpha
  maxsize = max(result1$N_TNode)
  ## message('max tree size: ', maxsize)
  
  G_alpha = result1$g_alpha
  
  folds = caret::createFolds(using_data$Y, K)
  result = c()  ##;result2 = list()
  for (i in 1:K) {
    Indexes = folds[[i]]
    train_S = using_data[Indexes, ]
    learn_S = using_data[-Indexes, ]
    
    nobs = dim(learn_S)[1]
    parm_num = 2 + length(var_role$FITL) + length(var_role$FITZ)
    if (min_size <= 0) {
      learn_size = max(floor(nobs * 0.05), parm_num * 10)
    } else if (min_size < 1) {
      learn_size = floor(nobs * min_size)
    } else if (min_size >= 1) {
      learn_size = floor(min_size * ((K - 1)/K))
    }
    
    learn_tree = create_tree(method, learn_S, var_role, alpha_main, alpha_in, depth = 1, learn_size, loglik_increment, parent_node = new.env())
    ## print.ezct(learn_tree) result2 = c(result2, learn_tree)
    learn_prune = prune_one(learn_tree)
    learn_level = level_order(learn_tree)
    train_fit(learn_tree, train_S, method, var_role)
    
    temp = learn_prune$prune_alpha
    Rts = tree_Rt(learn_level, learn_prune$prune_tree)
    result = cbind(result, Rts[sapply(G_alpha, function(x) rev(which(temp$alpha <= x))[1])])
    message("cross-validation iteration number: ", i)
  }
  result1$mean_Rcv = apply(result, 1, mean)
  result1$sd_Rcv = apply(result, 1, sd)
  result1$LB_Rcv = result1$mean_Rcv - SE * result1$sd_Rcv
  result1$UB_Rcv = result1$mean_Rcv + SE * result1$sd_Rcv
  N_Tnode = result1$N_TNode[max(which(result1$mean_Rcv <= result1$UB_Rcv[which.min(result1$mean_Rcv)]))]
  message("pruned tree size: ", N_Tnode)
  list(result1, result, N_Tnode, maxsize)
}


#===================================================#
#             pruning related functions             #
#===================================================#


train_fit = function(using_node, training_data, method, var_role) {
  if (!is.null(using_node)) {
    N_D = node_discriminant(using_node)
    if (!is.null(N_D)) {
      indexs = eval(parse(text = N_D))
      indexs[is.na(indexs)] = F
      N_D_NA = node_discriminant_NA(using_node)
      if (!(is.null(N_D_NA))) {
        indexsNA = eval(parse(text = N_D_NA))
        indexs = indexsNA | indexs
      }
      using_node$train_data = training_data[indexs, ]
    } else {
      using_node$train_data = training_data
    }
    fitts = switch(var_role$cdist, negbin = modelts2, poisson = modelts)
    
    
    using_node$train_model = fitts(method, using_node$train_data, using_node, var_role)
    temp = using_node$train_model$loglik
    if (is.null(temp)) {
      using_node$train_Rt = 0
    } else {
      using_node$train_Rt = temp * (-1)
    }
    
    if (!using_node$Node_info$Terminal) {
      train_fit(using_node$L, using_node$train_data, method, var_role)
      train_fit(using_node$R, using_node$train_data, method, var_role)
    }
  }
}

node_discriminant = function(using_node) {
  temp = parent.env(using_node)
  if (length(temp)) {
    LR = ""
    if (!using_node$LR) 
      LR = "!"
    temp2 = if (is.factor(temp$Node_info$Cut_point)) {
      paste0("training_data$", temp$Node_info$Cut_variable, "%in%", "as.character( \"", temp$Node_info$Cut_point, "\")")
    } else {
      paste("training_data$", temp$Node_info$Cut_variable, "<=", temp$Node_info$Cut_point)
    }
    if (length(temp2) > 1) 
      temp2 = paste(temp2, collapse = "|")
    temp2 = paste(LR, "(", temp2, ")")
    return(temp2)
  }
  return(NULL)
}

node_discriminant_NA = function(using_node) {
  temp = parent.env(using_node)
  if (length(temp)) {
    if (!is.null(temp$Node_info$L_NA)) {
      temp2 = NULL
      if (temp$Node_info$L_NA & using_node$LR) {
        temp2 = paste("is.na(training_data$", temp$Node_info$Cut_variable, ")")
      }
      if (!(temp$Node_info$L_NA) & !(using_node$LR)) {
        temp2 = paste("is.na(training_data$", temp$Node_info$Cut_variable, ")")
      }
      return(temp2)
    }
    return(NULL)
  }
  return(NULL)
}

tree_Rt = function(tree_level, prune_tree) {
  temp = prune_tree$g
  T_N = temp == 0
  result = do.call(sum, lapply(which(T_N), function(x) tree_level[[x]]$train_Rt))
  while (sum(T_N) > 1) {
    prune_cut = which(temp == min(temp[temp > 0]) & temp != 0)
    prune_cut = prune_cut[which.max(prune_tree[prune_cut, ]$N)]
    temp[prune_cut] = 0
    childs = unlist(path_to_leaf(prune_cut, prune_tree))
    T_N[childs] = F
    temp[childs] = 0
    T_N[prune_cut] = T
    result = c(result, do.call(sum, lapply(which(T_N), function(x) tree_level[[x]]$train_Rt)))
  }
  result
}

path_to_leaf = function(hash, prune_tree) {
  if (prune_tree$L_N[hash]) 
    return(list(prune_tree$t[hash], path_to_leaf(prune_tree$L_N[hash], prune_tree), path_to_leaf(prune_tree$R_N[hash], prune_tree))) else return(prune_tree$t[hash])
}

final_prune = function(tree_level, prune_tree, k) {
  temp = prune_tree$g
  T_N = temp == 0
  tree_level
  lapply(tree_level, function(x) x$NP = x$OP)
  while (sum(T_N) != k) {
    prune_cut = which(temp == min(temp[temp > 0]) & temp != 0)
    prune_cut = prune_cut[which.max(prune_tree[prune_cut, ]$N)]
    temp[prune_cut] = 0
    childs = unlist(path_to_leaf(prune_cut, prune_tree))
    T_N[childs] = F
    temp[childs] = 0
    T_N[prune_cut] = T
    lapply(childs, function(x) tree_level[[x]]$NP = T)
    tree_level[[prune_cut]]$NP = F
  }
}

#===================================================#
#           how data fitted when pruning            #
#===================================================#


modelts <- function(method, model_data, using_node, var_role) {
  if (sum(is.na(model_data))) 
    model_data = replace_NA(model_data)
  dep_var = attr(using_node$Node_info$Model$terms$count, "term.labels")
  var_role$FITL = dep_var
  if (dim(model_data)[1] == 0) {
    NULL
  } else {
    tryCatch({
      allvars = keep_var(model_data, unique(c(var_role$FITL, var_role$FITZ)))
      var_role$FITL = intersect(var_role$FITL, allvars)
      var_role$FITZ = intersect(var_role$FITZ, allvars)
      fitzi(method, model_data, var_role)
    }, warning = function(w) {
      tryCatch({
        var_role$FITZ = NULL
        fitzi(method, model_data, var_role)
      }, warning = function(w) {
        tryCatch({
          var_role$ZERO = "none"
          fitzi(method, model_data, var_role)
        }, error = function(e) {
          NULL
        })
      }, error = function(e) {
        tryCatch({
          var_role$ZERO = "none"
          fitzi(method, model_data, var_role)
        }, error = function(e) {
          NULL
        })
      })
    }, error = function(e) {
      tryCatch({
        var_role$FITZ = NULL
        fitzi(method, model_data, var_role)
      }, error = function(e) {
        tryCatch({
          var_role$ZERO = "none"
          fitzi(method, model_data, var_role)
        }, error = function(e) {
          NULL
        })
      })
    })
    
  }
}


modelts2 <- function(method, model_data, using_node, var_role) {
  if (sum(is.na(model_data))) 
    model_data = replace_NA(model_data)
  dep_var = attr(using_node$Node_info$Model$terms$count, "term.labels")
  var_role$FITL = dep_var
  if (dim(model_data)[1] == 0) {
    NULL
  } else {
    tryCatch({
      allvars = keep_var(model_data, unique(c(var_role$FITL, var_role$FITZ)))
      var_role$FITL = intersect(var_role$FITL, allvars)
      var_role$FITZ = intersect(var_role$FITZ, allvars)
      fitnb(method, model_data, var_role)
    }, warning = function(w) {
      tryCatch({
        var_role$FITZ = NULL
        fitnb(method, model_data, var_role)
      }, warning = function(w) {
        tryCatch({
          var_role$ZERO = "none"
          fitnb(method, model_data, var_role)
        }, error = function(e) {
          NULL
        })
      }, error = function(e) {
        tryCatch({
          var_role$ZERO = "none"
          fitnb(method, model_data, var_role)
        }, error = function(e) {
          NULL
        })
      })
    }, error = function(e) {
      tryCatch({
        var_role$FITZ = NULL
        fitnb(method, model_data, var_role)
      }, error = function(e) {
        tryCatch({
          var_role$ZERO = "none"
          fitnb(method, model_data, var_role)
        }, error = function(e) {
          NULL
        })
      })
    })
    
  }
}


#===================================================#
#                   use the method                  #
#===================================================#


tree_control = function(method = c("multiple", "constant", "simple", "back"), prune_or_not = T, SE = 0.5, fold = 10, trimzero = T, alpha_main = 0.05, alpha_in = 0.05, depth = 1, part_score = c("all", "count", "zero"), cdist = c("poisson", "negbin"), min_size = -1, loglik_increment = 0, zero_type = c("inflated", "altered", "none"), CV_seed = 39, back_alpha = 0.1, trace = F) {
  list(prune_or_not = prune_or_not, SE = SE, fold = fold, alpha_main = alpha_main, alpha_in = alpha_in, trimzero = trimzero, depth = depth, part_score = match.arg(part_score), cdist = match.arg(cdist), min_size = min_size, loglik_increment = loglik_increment, method = match.arg(method), zero_type = match.arg(zero_type), CV_seed = CV_seed, back_alpha = back_alpha, trace = trace)
}


## revised coretree


coretree = function(formula, data, na.action = na.pass, offset, control = tree_control(...), ...) {
  
  if (control$CV_seed == 39) {
    set.seed(39)
  }
  
  cl = match.call()
  
  if (missing(data)) 
    data = environment(formula)
  using_data = match.call(expand.dots = FALSE)
  m = match(c("formula", "data", "na.action", "offset"), names(using_data), 0)
  using_data = using_data[c(1, m)]
  using_data$drop.unused.levels = TRUE
  using_data$na.action <- na.action
  
  formula <- Formula::as.Formula(formula)
  
  if (length(formula)[2L] < 2L) {
    parts = F
  } else {
    
    if (length(formula)[2L] > 3L) {
      formula <- Formula::Formula(formula(formula, rhs = 1L:3L))
      warning("Formula must not have more than three RHS parts")
    }
    parts = T
  }
  using_data$formula <- formula
  using_data[[1]] <- as.name("model.frame")
  using_data <- eval(using_data, parent.frame())
  logics = sapply(using_data, is.logical)
  
  if (sum(logics)) {
    logicnames = names(using_data)[logics]
    using_data[logicnames] = lapply(using_data[logicnames], as.factor)
  }
  
  if (control$trace) 
    return(using_data[1:4, ])
  
  if (!missing(offset)) {
    N = NCOL(using_data)
    OFF = model.offset(using_data)
    using_data = using_data[, -N]
  }
  
  var_role = if (identical(suppressWarnings(attr(terms(formula, data = using_data, rhs = 3L), "term.labels")), character(0))) {
    list(FITL = attr(terms(formula, data = using_data, rhs = 1L), "term.labels"),
         CUT = if (!parts) {
      attr(terms(formula, data = using_data, rhs = 1L), "term.labels")
    } else {
      attr(terms(formula, data = using_data, rhs = 2L), "term.labels")
    }, DEP = as.character(formula)[2], FITZ = if (length(formula)[2L] == 3L) {
      attr(terms(formula, data = using_data, rhs = 3L), "term.labels")
    } else {
      NULL
    }, OFF = if (missing(offset)) {
      NULL
    } else {
      OFF
    }, ZERO = control$zero_type, part_score = control$part_score, trimzero = control$trimzero, cdist = control$cdist, back_alpha = control$back_alpha)
  } else {
    list(CUT = suppressWarnings(attr(terms(formula, data = using_data, lhs = 3L), "term.labels")), 
         DEP = suppressWarnings(attr(terms(formula, data = using_data, rhs = 3L), "term.labels")), 
         FITL = suppressWarnings(attr(terms(formula, data = using_data, lhs = 3L), "term.labels")), 
         FITZ = suppressWarnings(attr(terms(formula, data = using_data, lhs = 3L), "term.labels")), 
         OFF = if (missing(offset)) {
      NULL
    } else {
      OFF
    }, ZERO = control$zero_type, part_score = control$part_score, trimzero = control$trimzero, cdist = control$cdist, back_alpha = control$back_alpha)
  }
  
  for (i in 1:length(var_role$DEP)) {
    names(using_data)[i] = "Y"
    if (sum(is.na(using_data[i]))) {
      using_data = using_data[complete.cases(using_data[i]), ]
      message("Omit the cases that have missing values on response")
    }
  }
  
  if (control$trace) 
    print(var_role[1:4])
  # if(control$method=='multiple' | control$method=='simple') {if(do.call(max,lapply(using_data[var_role$FITL],nlevels))>2) {stop('Categorical regressors are allowed with only 2-level')}}
  nobs = dim(using_data)[1]
  
  if (nobs < 1) 
    stop("empty data")
  
  if (control$part_score == "count" & !is.null(var_role$FITZ)) {
    var_role$FITZ = NULL
    message("Using score of count part does NOT allow covariates on zero")
  }
  
  if (control$part_score == "zero" & !is.null(var_role$FITL)) {
    var_role$FITL = NULL
    message("Using score of zero part does NOT allow covariates on count")
  }
  
  for (i in 1:length(var_role$DEP)) {
    if (!isTRUE(all.equal(as.vector(unlist((using_data[i]))), as.integer(unlist(round(using_data[i] + 0.001)))))) 
      stop("invalid dependent variable, non-integer values")
  }
  
  for (i in 1:length(var_role$DEP)) {
    using_data[i] <- as.integer(unlist(round(using_data[i] + 0.001)))
  }
  
  for (i in 1:length(var_role$DEP)) {
    if (any(using_data[i] < 0)) 
      stop("invalid dependent variable, negative counts")
  }
  
  node_size = control$min_size
  parm_num = 2 + length(var_role$FITL) + length(var_role$FITZ)
  if (control$min_size <= 0) {
    control$min_size = max(floor(nobs * 0.05), parm_num * 10)
  } else if (control$min_size < 1) {
    control$min_size = floor(nobs * control$min_size)
  } else if (control$min_size >= 1) {
    control$min_size = control$min_size
  }
  if (control$min_size >= nobs) {
    stop("The sample size is too small to fit the model.")
  }
  
  message("The task is in execution...")
  
  temp = create_tree(control$method, using_data, var_role, control$alpha_main, control$alpha_in, control$depth, control$min_size, control$loglik_increment, parent_node = new.env())
  message("temp done")
  tree_level = level_order(temp)
  temp$max_size = do.call(sum, lapply(tree_level, function(x) !x$OP))
  temp$final_size = temp$max_size
  message("main tree constructed")
  if (control$prune_or_not) {
    if (temp$max_size > 1) {
      message("max tree size: ", temp$max_size)
      ptable = K_fold_CV(control$fold, temp, control$method, using_data, var_role, control$alpha_main, control$alpha_in, control$SE, node_size, control$loglik_increment)
      temp$prune_table = ptable
      final_prune(tree_level, prune_one(temp)$prune_tree, ptable[[3]])
      temp$final_size = ptable[[3]]
      message("pruned subtree sequence in size: {", paste0(ptable[[1]][, 2], collapse = ", "), "}")
    } else {
      message("No split is acceptable")
    }
  }
  temp$using_data = using_data
  temp$method = control$method
  temp$offsett = cl$offset
  temp$min_size = control$min_size
  temp$var_role = var_role
  temp$CV_setting = list(control$fold, control$alpha_main, control$alpha_in, control$SE, node_size, control$loglik_increment)
  class(temp) = append(class(temp), "core")
  return(temp)
}

#===================================================#
#                  helper functions                 #
#===================================================#
## delete variable with only 1 level
keep_var = function(using_data, role) {
  keep = using_data[role] %>% lapply(unique) %>% lapply(length) > 1
  keep = names(keep[keep == T])
  if (!length(keep)) {
    return(NULL)
  } else {
    return(keep)
  }
}

## opposition of %in% #
"%out%" <- function(x, y) !x %in% y

## assign variables to each group by role vector
var_roles = function(gvar, gtype) {
  CAT = gvar[gtype %in% c("c", "k", "a", "f", "K", "A", "F")]
  CUT = gvar[gtype %in% c("s", "c", "L", "Z", "B", "K", "A", "F")]
  DEP = gvar[gtype %in% c("d")]
  NOTD = gvar[gtype %in% c("s", "c", "l", "L", "k", "K", "b", "B", "f", "F", "z", "Z", "a", "A")]
  FITL = gvar[gtype %in% c("l", "L", "k", "K", "b", "B", "f", "F")]
  FITZ = gvar[gtype %in% c("z", "Z", "a", "A", "b", "B", "f", "F")]
  OFF = gvar[gtype %in% c("o")]
  USE = gvar[gtype %out% c("x")]
  return(list(CUT = CUT, DEP = DEP, FITL = FITL, FITZ = FITZ, OFF = OFF))
}


## replace NA

replace_NA = function(NAdata) {
  temp = rownames(NAdata)
  as.data.frame(lapply(NAdata, function(x) {
    if (is.factor(x)) {
      x[which(is.na(x))] = v_mode(x)
    } else {
      x[which(is.na(x))] = mean(x, na.rm = T)
    }
    x
  }), row.names = temp)
}

train_predict = function(using_node, training_data, offsett = NULL) {
  if (!is.null(using_node)) {
    N_D = node_discriminant(using_node)
    if (!is.null(N_D)) {
      indexs = eval(parse(text = N_D))
      indexs[is.na(indexs)] = F
      N_D_NA = node_discriminant_NA(using_node)
      if (!(is.null(N_D_NA))) {
        indexsNA = eval(parse(text = N_D_NA))
        indexs = indexsNA | indexs
      }
      temp = training_data[indexs, ]
      if (!is.null(nrow(temp))) {
        if (nrow(temp) > 0) {
          using_node$Node_info$Model$call$offset = offsett
          using_node$prf = VGAM::predict(using_node$Node_info$Model, newdata = temp, type = "response")
          names(using_node$prf) = rownames(temp)
        } else {
          using_node$prf = NULL
          temp = NULL
        }
      }
    } else {
      temp = training_data
      using_node$Node_info$Model$call$offset = offsett
      using_node$prf = VGAM::predict(using_node$Node_info$Model, newdata = temp, type = "response")
      names(using_node$prf) = rownames(temp)
    }
    if (!using_node$Node_info$Terminal) {
      train_predict(using_node$L, temp, offsett)
      train_predict(using_node$R, temp, offsett)
    }
  }
}


#===================================================#
#       plot the tree by package partykit           #
#===================================================#


tree_dsc2 = function(tree_level, using_data) {
  tree_kids = function(i) {
    if (!tree_level[[i]]$NP) 
      return(NULL) else return(c(tree_level[[i]]$L$hash, tree_level[[i]]$R$hash))
  }
  
  tree_split = function(i) {
    if (!tree_level[[i]]$NP) 
      return(NULL)
    if (is.factor(tree_level[[i]]$Node_info$Cut_point)) {
      temp = tree_level[[i]]$Node_info$Cut_point
      temp3 = unique(using_data[as.numeric(names(tree_level[[i]]$Node_info$Model$weights)), ][tree_level[[i]]$Node_info$Cut_variable])
      temp2 = rep(NA, length(levels(temp)))
      temp2[match(unlist(temp3), levels(temp))] = 2
      temp2[match(temp, levels(temp))] = 1
      split = partysplit(which(names(using_data) == tree_level[[i]]$Node_info$Cut_variable), index = as.integer(temp2))
    } else {
      split = partysplit(which(names(using_data) == tree_level[[i]]$Node_info$Cut_variable), breaks = tree_level[[i]]$Node_info$Cut_point)
    }
    return(split)
  }
  
  tree_node = function(i) {
    if (is.null(tree_kids(i))) 
      return(partynode(as.integer(i), info = tree_info[[i]]))
    partynode(as.integer(i), split = tree_split(i), kids = lapply(tree_kids(i), tree_node))
  }
  
  if (tree_level[[1]]$method == "simple") {
    tree_info = lapply(tree_level, function(root) {
      best = if (!is.null(root$Node_info$Model[[1]])) {
        if (class(root$Node_info$Model[[1]])[1] == "zeroinfl" | class(root$Node_info$Model[[1]])[1] == "hurdle") {
          attr(root$Node_info$Model[[1]]$terms$zero, "term.labels")
        } else {
          attr(root$Node_info$Model[[1]]$terms, "term.labels")
        }
      } else {
        NULL
      }
      best_sign = if (!is.null(root$Node_info$Model[[1]])) {
        if (class(root$Node_info$Model[[1]])[1] == "zeroinfl" | class(root$Node_info$Model[[1]])[1] == "hurdle") {
          root$Node_info$Model[[1]]$coef$zero[best] > 0
        } else {
          root$Node_info$Model[[1]]$coef[best] > 0
        }
      } else {
        NULL
      }
      list(Size = root$Node_info$Node_size, n0 = root$Node_info$H_n0, pzero = root$Node_info$H_pzero, lambda = root$Node_info$H_lambda, missing = root$Node_info$L_NA, negbin = if (class(root$Node_info$Model[[1]])[1] %in% c("zeroinfl", "hurdle")) {
        if (root$Node_info$Model[[1]]$dist[[1]] == "negbin") root$Node_info$Model[[1]]$theta else F
      } else {
        if (class(root$Node_info$Model[[1]]$Model[[1]])[1] == "negbin") root$Node_info$Model[[1]]$theta else F
      }, best = if (!is.null(best)) {
        if (is.na(best_sign)) {
          paste0(best)
        } else {
          paste0(ifelse(best_sign, "+", "-"), best)
        }
      } else {
        NULL
      })
    })
  } else {
    tree_info = lapply(tree_level, function(root) list(Size = root$Node_info$Node_size, n0 = root$Node_info$H_n0, pzero = root$Node_info$H_pzero, lambda = root$Node_info$H_lambda, missing = root$Node_info$L_NA, negbin = if (class(root$Node_info$Model[[1]])[1] %in% c("zeroinfl", "hurdle")) {
      if (root$Node_info$Model[[1]]$dist[[1]] == "negbin") root$Node_info$Model[[1]]$theta else F
    } else {
      if (class(root$Node_info$Model[[1]])[1] == "negbin") root$Node_info$Model[[1]]$theta else F
    }))
  }
  
  party_tree = tree_node(1)
  preid = as.character(nodeids(party_tree))
  party_tree = party(party_tree, using_data[0, ])
  names(party_tree) = preid
  party_tree
}

partykit_draw = function(p_t, big = F) {
  if (big) {
    temp = list(justmin = 4)
  } else {
    temp = list()
  }
  plot(p_t, ep_args = temp, tp_args = list(align = c("left"), fill = c("yellow2", "white"), just = c("top"), width = 8, FUN = function(node) {
    temp = format(round(node$lambda[[1]], 3), digits = 3)
    temp2 = format(round(node$n0[[1]], 0), digits = 0)
    temp3 = format(round(node$pzero[[1]], 3), digits = 3)
    temp4 = format(round(node$negbin[[1]], 3), digits = 3)
    c(expression(), paste(""), if (!node$negbin) bquote(paste(atop(italic(n) == .(node$Size), paste(tilde(lambda) == .(temp))))) else bquote(paste(atop(italic(n) == .(node$Size), paste(tilde(mu) == .(temp))))), paste(""), if (node$negbin) bquote(paste(italic(theta) == .(temp4))), paste(""), if (!is.na(node$pzero)) {
      c(if (node$pzero < 0.001) {
        bquote(paste(atop(tilde(italic(p)) %~~% .(temp3), paste(hat(italic(f)[0]) == .(temp2)))))
      } else {
        bquote(paste(atop(tilde(italic(p)) == .(temp3), paste(hat(italic(f)[0]) == .(temp2)))))
      }, paste(""))
    } else {
      bquote(paste(hat(italic(f)[0]) == .(temp2)))
    })
  }))
}




partykit_draw_best = function(p_t, big = F) {
  if (big) {
    temp = list(justmin = 4)
  } else {
    temp = list()
  }
  plot(p_t, ep_args = temp, tp_args = list(align = c("left"), fill = c("yellow2", "white"), just = c("top"), width = 8, FUN = function(node) {
    temp = format(round(node$lambda[[1]], 3), digits = 3)
    temp2 = format(round(node$n0[[1]], 0), digits = 0)
    temp3 = format(round(node$pzero[[1]], 3), digits = 3)
    temp4 = format(round(node$negbin[[1]], 3), digits = 3)
    c(expression(), paste(""), if (!node$negbin) bquote(paste(atop(italic(n) == .(node$Size), paste(tilde(lambda) == .(temp))))) else bquote(paste(atop(italic(n) == .(node$Size), paste(tilde(mu) == .(temp))))), paste(""), if (node$negbin) bquote(paste(italic(theta) == .(temp4))), paste(""), if (!is.na(node$pzero)) {
      c(if (node$pzero < 0.001) {
        bquote(paste(atop(tilde(italic(p)) %~~% .(temp3), paste(hat(italic(f)[0]) == .(temp2)))))
      } else {
        bquote(paste(atop(tilde(italic(p)) == .(temp3), paste(hat(italic(f)[0]) == .(temp2)))))
      }, paste(""))
    } else {
      bquote(paste(hat(italic(f)[0]) == .(temp2)))
    }, paste(""), if (!is.null(node$best)) {
      if (nchar(node$best) > 7) {
        bquote(atop(.(abbreviate(node$best, 6)) ~ ., phantom()))
      } else {
        bquote(atop(.(node$best), phantom()))
      }
    } else {
      NULL
    })
  }))
}



partykit_draw_big = function(p_t) {
  plot(p_t, ep_args = list(justmin = 4), tp_args = list(align = c("left"), fill = c("yellow2", "white"), just = c("top"), width = 8, FUN = function(node) {
    temp = format(round(node$lambda[[1]], 3), digits = 3)
    temp2 = format(round(node$n0, 0), digits = 0)
    temp3 = format(round(node$pzero, 3), digits = 3)
    temp4 = format(round(node$negbin, 3), digits = 3)
    c(expression(), paste(""), if (!node$negbin) bquote(paste(atop(italic(n) == .(node$Size), paste(tilde(lambda) == .(temp))))) else bquote(paste(atop(italic(n) == .(node$Size), paste(tilde(mu) == .(temp))))), paste(""), if (node$negbin) bquote(paste(italic(theta) == .(temp4))), paste(""), if (!is.na(node$pzero)) {
      c(if (node$pzero < 0.001) {
        bquote(paste(atop(tilde(italic(p)) %~~% .(temp3), paste(hat(italic(f)[0]) == .(temp2)))))
      } else {
        bquote(paste(atop(tilde(italic(p)) == .(temp3), paste(hat(italic(f)[0]) == .(temp2)))))
      }, paste(""))
    } else {
      bquote(paste(hat(italic(f)[0]) == .(temp2)))
    })
  }))
}



#===================================================#
#                core    functions                  #
#===================================================#


print.core = function(object) {
  cat("\n", "Text output:", "\n\n")
  pre_order(object)
}

plot.core = function(object, big = F) {
  
  
  
  temp = level_order(object)
  if (object$method == "simple") {
    partykit_draw_best(tree_dsc2(temp, object$using_data), big)
  } else if (do.call(sum, lapply(temp, function(x) !x$NP)) <= 16) {
    partykit_draw(tree_dsc2(temp, object$using_data), big)
  } else {
    partykit_draw_big(tree_dsc2(temp, object$using_data))
  }
}

summary.core = function(object, node = NULL) {
  temp = level_order(object)
  if (!is.null(node)) {
    k = node
  } else {
    k = which(unlist(lapply(temp, function(x) !x$NP)))
  }
  temp = lapply(k, function(x) summary(temp[[x]]$Node_info$Model))
  names(temp) = paste("node", k)
  temp
}

AIC.core = function(object, node = NULL) {
  temp = level_order(object)
  if (!is.null(node)) {
    k = node
  } else {
    k = which(unlist(lapply(temp, function(x) !x$NP)))
  }
  temp = lapply(k, function(x) AIC(temp[[x]]$Node_info$Model))
  names(temp) = paste("node", k)
  do.call(sum, temp)
}

BIC.core = function(object, node = NULL) {
  temp = level_order(object)
  if (!is.null(node)) {
    k = node
  } else {
    k = which(unlist(lapply(temp, function(x) !x$NP)))
  }
  temp = lapply(k, function(x) BIC(temp[[x]]$Node_info$Model))
  names(temp) = paste("node", k)
  do.call(sum, temp)
}

predict.core = function(object, newdata) {
  if (missing(newdata)) {
    temp = do.call(rbind, lapply(level_order(object), function(x) if (!x$NP) {
      temp2 = data.frame(nodeID = x$hash, observerd = x$Node_info$Model$y, fitted = x$Node_info$Model$fitted.values)
      rownames(temp2) = rownames(x$Node_info$Model$model)
      temp2
    }))
    
  } else {
    train_predict(object, newdata, object$offsett)
    temp = do.call(rbind, lapply(level_order(object), function(x) {
      temp2 = NULL
      if (!x$NP) {
        if (!is.null(x$prf)) 
          temp2 = data.frame(nodeID = x$hash, predicted = x$prf)
      }
      x$prf = NULL
      temp2
    }))
  }
  temp[order(as.numeric(attr(temp, "row.names"))), ]
}

change_size = function(object, tree_size) {
  final_prune(level_order(object), prune_one(object)$prune_tree, tree_size)
}

qqrplot_c = function(object, node = NULL, ...) {
  if (is.null(node)) {
    temp = lapply(level_order(object), function(x) if (!x$NP) 
      x)
    nodeID = paste("node", which(sapply(temp, function(x) !is.null(x))))
    temp = temp[sapply(temp, function(x) !is.null(x))]
    temp = lapply(temp, function(x) x$Node_info$Model)
    aa = length(temp)
    nums = find_par(aa)
    par(mfrow = nums)
    aa = lapply(1:aa, function(x) {
      qqrplot_from_countreg(temp[[x]], main = nodeID[x], ...)
    })
  } else {
    temp = lapply(level_order(object), function(x) x$Node_info$Model)
    lapply(node, function(x) qqrplot_from_countreg(temp[[x]], main = paste("node", x), ...))
    invisible()
  }
}


rootogram.core = function(object, node = NULL, style = "hanging", num = NULL) {
  if (is.null(node)) {
    aa = lapply(level_order(object), function(x) if (!x$NP) 
      x)
    aa = aa[sapply(aa, function(x) !is.null(x))]
    aaa = length(aa)
    nums = find_par(aaa)
    par(mfrow = nums)
    aa = lapply(aa, function(x) {
      temp = x$Node_info$Model
      rootogram(temp, max = max(temp$y), main = paste0("node ", x$hash), xlab = object$var_role$DEP, style = style, col = "black", fill = "gray70")
      ## ,cex.lab = 1.2)
    })
  } else {
    aa = lapply(node, function(x) {
      temp = level_order(object)[[x]]
      temp2 = temp$Node_info$Model
      rootogram(temp2, max = max(temp2$y), main = paste0("node ", temp$hash), xlab = object$var_role$DEP, style = style, col = "black", fill = "gray70", cex.lab = 1.2)
    })
    
  }
}

qqrplot_from_countreg <- function(object, type = c("random", "quantile"), nsim = 1L, prob = 0.5, range = FALSE, diag = TRUE, col = "black", fill = "lightgray", xlim = NULL, ylim = NULL, main = "Q-Q residuals plot", xlab = "Theoretical quantiles", ylab = "Quantile residuals", outlier = T, ...) {
  ## compute quantile residuals
  qres <- if ("glm" %out% class(object)) {
    qresiduals(object, type = type, nsim = nsim, prob = prob)
  } else {
    sapply(1:nsim, function(x) statmod::qresiduals(object))
  }
  
  if (is.null(dim(qres))) 
    qres <- matrix(qres, ncol = 1L)
  inf_res = matrix(F, ncol = ncol(qres), nrow = nrow(qres))
  for (i in 1:ncol(qres)) {
    if (sum(is.nan(qres[, i])) > 0) {
      qres[is.nan(qres[, i]), i] = Inf
    }
    if (sum(is.infinite(qres[, i])) > 0) {
      qres_range = range(qres[is.finite(qres[, i]), i])
      inf_res[, i] = is.infinite(qres[, i])
      qres[inf_res, i] = ifelse(sign(qres[inf_res, i]) + 1, qres_range[2] + 0.5, qres_range[1] - 0.5)
    }
  }
  if (is.null(dim(qres))) 
    qres <- matrix(qres, ncol = 1L)
  ## corresponding normal quantiles
  q2q <- function(y) qnorm(ppoints(length(y)))[order(order(y))]
  qnor <- apply(qres, 2L, q2q)
  
  ## default plotting ranges
  if (is.null(xlim)) 
    xlim <- range(qnor)
  if (is.null(ylim)) 
    ylim <- range(qres)
  xylim = range(c(xlim, ylim))
  xlim = ylim = xylim
  ## set up coordinates
  plot(0, 0, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main)
  
  ## polygon for range
  if (!identical(range, FALSE)) {
    if (isTRUE(range)) 
      range <- c(0.01, 0.99)
    rg <- qresiduals(object, type = "quantile", prob = range)
    y1 <- sort(rg[, 1])
    y2 <- sort(rg[, 2])
    x <- c(q2q(y1), rev(q2q(y2)))
    y <- c(y1, rev(y2))
    y[!is.finite(y)] <- 100 * sign(y[!is.finite(y)])
    x[!is.finite(x)] <- 100 * sign(x[!is.finite(x)])
    polygon(x, y, col = fill, border = fill)
    box()
  }
  
  ## add Q-Q plot(s)
  
  for (j in 1:ncol(qres)) {
    for (i in which(!inf_res[, j])) points(qnor[i, j], qres[i, j], col = col, ...)
    for (i in which(inf_res[, j])) {
      points(qnor[i, j], qres[i, j], col = col, pch = 4)
      if (outlier) {
        text(qnor[i, j], qres[i, j], labels = i, cex = 0.7, pos = 1)
      }
    }
  }
  ## reference diagol
  if (!identical(diag, FALSE)) {
    if (isTRUE(diag)) 
      diag <- "black"
    abline(0, 1, col = diag, lty = 2)
  }
  
  ## return coordinates invisibly
  invisible(list(normal = qnor, residuals = qres))
}

table_c = function(t_tree, node = NULL) {
  aa = lapply(level_order(t_tree), function(x) if (!x$NP) 
    x)
  nodeID = paste("node", which(sapply(aa, function(x) !is.null(x))))
  aa = aa[sapply(aa, function(x) !is.null(x))]
  EE = lapply(aa, function(x) {
    temp = x$Node_info$Model
    if (class(temp)[1] == "zeroinfl" | class(temp)[1] == "hurdle") {
      colSums(VGAM::predict(temp, type = "prob", at = 0L:max(temp$y)))
    } else {
      lam = exp(VGAM::predict(temp))
      len = max(temp$y) + 1
      Rva <- matrix(NA, nrow = length(temp$y), ncol = len)
      for (i in 1:len) Rva[, i] <- dpois(i - 1, lambda = lam)
      colSums(Rva)
    }
    
  })
  OO = lapply(aa, function(x) {
    temp = x$Node_info$Model
    tabulate(temp$y + 1)
  })
  temp = lapply(1:length(EE), function(x) data.frame(OO[[x]], EE[[x]]))
  temp = lapply(temp, function(x) {
    names(x) = c("Observed", "Expected")
    x
  })
  names(temp) = nodeID
  
  if (is.null(node)) {
    temp
  } else {
    temp = lapply(node, function(x) temp[[paste("node", x)]])
    names(temp) = paste("node", node)
    temp
  }
}

plot_models = function(models, FUN, windows = F) {
  readkey <- function() {
    cat("Press [enter] to continue")
    line <- readline()
  }
  if (windows) {
    for (i in 1:length(models)) {
      cat("model ", i + 1, " :  ", names((models)[i]), "\n")
      windows()
      FUN(models[[i]])
    }
  } else {
    for (i in 1:length(models)) {
      cat("model: ", names((models)[i]), "\n")
      FUN(models[[i]])
      readkey()
    }
  }
}

find_par = function(x) {
  temp = which.min(abs(((1:10)^2 - x)))
  if (x <= temp^2) 
    c(temp, temp) else c(temp, temp + 1)
}