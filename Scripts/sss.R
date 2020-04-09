library(sets)

sss <- function(objective, constraint, initialmodel, Omega, niter){
  scoreboard = list()
  model = initialmodel
  scoreboard[[paste("Model:",toString(model))]] = objective(model)
  for(i in seq(1,niter)){
    output = nextmodel(objective,constraint,Omega,model, scoreboard)
    scoreboard = output$scoreboard
    model = output$nextmodel
  }
  
  sortedlb = scoreboard[sort(unlist(scoreboard),decreasing=TRUE,index.return=TRUE)$ix]
  
}

drop1 <- function(A){
  Map(function(x){cset_difference(A,x)},sapply(A,identity))
}

add1 <- function(Omega, A){
  Ac=sapply(cset_difference(Omega,A),identity)
  Map(function(x){cset_union(x,A)},Ac)
}

swap1 <- function(Omega,A){
  if(cset_is_empty(A)){
    list()
  }else{
    Ac=sapply(cset_difference(Omega,A),identity)
    Aless1=drop1(A)
    unlist(Map(function(y){Map(function(x){cset_union(x,y)},Ac)},Aless1),recursive=FALSE)
  }
}

neighbours <- function(Omega,A) list(dropped = drop1(A), added = add1(Omega,A), swapped = swap1(Omega,A))

nextmodel <- function(objective, constraint, Omega, A, scoreboard){
  
  #Obtain neighbourhood sets "added", "dropped" and "swapped"
  neigh=neighbours(Omega,A)
  
  #filter empty neighbourhood sets
  neigh=neigh[sapply(neigh,function(x){!length(x)==0})]
  
  #filter infeasible sets
  neigh=Map(function(i){
    neigh[[i]][sapply(neigh[[i]],constraint)]
  },names(neigh)) #feasible
  
  #filter empty neighbourhood sets
  neigh=neigh[sapply(neigh,function(x){!length(x)==0})]
  
  # Sample one model each from added,dropped,swapped so long as these are not empty
  # If all 0's, then prob should just be null, so then 
  neigh_vec = unlist(neigh, recursive=FALSE)
  # If there are no neighbors, the next model is the current model, and the scoreboard is unchanged
  if(is.null(neigh_vec)){
    output = list("nextmodel" = cset(A),
                  "scoreboard" = scoreboard)
    return(output)
  }
  
  if(parallel){
    clusterExport(cl, varlist=c("scoreboard"), envir=environment())
    probs_vec = parLapply(cl, neigh_vec, function(mod) {
      p = scoreboard[[paste("Model:", toString(mod))]]
      if(is.null(p)){p = objective(mod)}
      p
    })
  }else{
    probs_vec = lapply(neigh_vec, function(mod) {
      p = scoreboard[[paste("Model:", toString(mod))]]
      if(is.null(p)){p = objective(mod)}
      p
    })
  }

  # Add all the log-probabilities to the scoreboard
  for(i in 1:length(probs_vec)){
    scoreboard[[paste("Model:", toString(neigh_vec[[i]]))]] = probs_vec[[i]]
  }
  
  # Generate list of probabilities in the same form as the 'neigh' list
  len = length(neigh)
  len2 = lapply(neigh, length)
  len3 = c(0, cumsum(len2))
  probs = list()
  for(i in 1:len){
    probs[[i]] = probs_vec[(len3[i]+1):len3[i+1]]
    probs[[i]] = as.list(rescale_probs(unlist(probs[[i]]))) #Putting them all on a re-scaled basis, so no computational 0's
    probs[[i]] = lapply(probs[[i]], exp) # Exponentiating to turn the log-scores into probabilities
  }
  
  
  # If all models have computationally 0 probability, set to NULL - which implies equal probability for all models in that group
  probs = lapply(probs, function(x) if(sum(unlist(x)) == 0){NULL}else{x}) 
  len = length(neigh)
  candidate_numbers = unlist(Map(function(v){sample(1:length(neigh[[v]]), size = 1, prob = probs[[v]])}, 1:len), recursive=FALSE)
  candidates = unlist(Map(function(v){neigh[[v]][candidate_numbers[v]]}, 1:len), recursive=FALSE)
  candidate_probs = unlist(Map(function(v){probs[[v]][candidate_numbers[v]]}, 1:len))
  if(sum(candidate_probs) == 0){candidate_probs = NULL} #May need this, or may not if all possibilities also have 0 probs
  
  
  output = list("nextmodel" = cset(unlist(sample(candidates, size = 1, prob = candidate_probs))),
                "scoreboard" = scoreboard)
  return(output)
}

fix_inf_probs <- function(probs){
  probs[!is.finite(probs)] = 1E7
  return(probs)
}

to_dataframe <- function(list){
  data.frame("model" = sapply(list, function(x) toString(x)),
             "score" = Map(objective, list),
             "time" = time,
             "SSS_action" = )
}

