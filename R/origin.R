origin <- function(genotype, data=sys.frame(sys.parent()), pedigree,
         id, id.father, id.mother, affected, codes = 2, first = FALSE, verbose=TRUE) {
  if (missing(genotype))
    stop("Missing genotype argument with no default")
  cl <- match.call()
  g <- eval(cl$genotype, data)
  gname <- deparse(cl$genotype)
  
  ## get main vars
  ped <- if(missing(pedigree)) {
    evalq(pedigree,data)
  } else {
    eval(cl$pedigree,data)
  }
  id <- if(missing(id)) {
    evalq(id,data)
  } else {
    eval(cl$id.father,data)
  }
  f <- if(missing(id.father)) {
    evalq(id.father,data)
  } else {
    eval(cl$id.father,data)
  }
  m <- if(missing(id.mother)) {
    evalq(id.mother,data)
  } else {
    eval(cl$id.mother,data)
  }
  aff <- if(missing(affected)) {
    evalq(affected,data)
  } else {
    eval(cl$affected,data)
  }
  
  
  if(is.factor(g)) {
    g.count <- as.integer(g)-1
  } else {
    g.count <- as.integer(as.factor(g))-1
  }
  
  if(verbose) {
    cat("\nalleles coded as:\n")
    print(table(genotype=g,
                count=g.count))
  }
  
  if(length(setdiff(g.count,NA))>3) {
    stop("require exactly 2 alleles\n")
  }
  
  
  
  ## make dataframe
  if(is.factor(ped))
    ped <- unfactor(ped)
  if(is.factor(m))
    m <- unfactor(m)
  if(is.factor(f))
    f <- unfactor(f)
  if(is.factor(id))
    id <- unfactor(id)
  
  ## only affecteds and their parents
  lid <- paste(ped,id,sep=":")
  lf <- paste(ped,f,sep=":")
  lm <- paste(ped,m,sep=":")
  
  x <- data.frame(ped=ped,id=id,f=f,m=m,
                  longid=lid,longf=lf,longm=lm,
                  g=g.count,aff=(aff %in% codes),
                  stringsAsFactors=FALSE)
  x <- subset(x,f==0 | aff) # founders and affected offspring only
  x <- subset(x,longf %in% lid & longm %in% lid) # drop affected offspring with missing parents
  x <- subset(x,!is.na(aff) & !is.na(g)) # drop offpsring with unknown disease state or genotype
  
  ## add parental genotype
  x$g.f <- g.count[match(x$longf,lid)]
  x$g.m <- g.count[match(x$longm,lid)]
  x <- subset(x,!is.na(g.f) & !is.na(g.m)) # drop offpsring with missing parental genotypes
  
  ## add mating type
  tmp <- x[,c("g.f","g.m")]
  wh <- which(tmp[,1]>tmp[,2])
  tmp[wh,] <- tmp[wh,c(2,1)]
  x$mating <- apply(tmp,1,paste,collapse="")
  x$mating.ordered <- apply(x[,c("g.f","g.m")],1,paste,collapse="")
  
  
  ## want only one set of parents per family
  ## one item per ped
  x <- split(x,x$ped)
  fathers <- lapply(x,function(p) setdiff(unique(p$f),0))
  mothers <- lapply(x,function(p) setdiff(unique(p$m),0))
  nfath <- sapply(fathers,length)
  nmoth <- sapply(mothers,length)
  wh <- which(nfath>1 | nmoth>1)
  if(length(wh)) {
    cat("\n",length(wh), " family/ies found with >1 father and/or >1 mother; keeping only the first set of parents in each:\n")
    cat(names(x)[wh],sep="\n")
    for(i in wh) {
      x[[i]] <- subset(x[[i]],f==fathers[[i]][1] & m==mothers[[i]][1])
    }
  }
  
  mating <- sapply(x,function(p) { unique(p$mating)})
  mating.ordered <- sapply(x,function(p) { unique(p$mating.ordered)})    
  offspring <- lapply(x,"[[","g")
  
  ## remove misinheritances
  ## 00 = no 1 or 2 offspring
  wh <- which(mating=="00")
  if(length(wh)) {
    offspring[wh] <- lapply(offspring[wh],mydrop,c(1,2))
  }
  
  ## 22 = no 1 or 0 offspring
  wh <- which(mating=="22")
  if(length(wh)) {
    offspring[wh] <- lapply(offspring[wh],mydrop,c(1,0))
  }
  
  ## 02 = no 2 or 0 offspring
  wh <- which(mating=="02")
  if(length(wh)) {
    offspring[wh] <- lapply(offspring[wh],mydrop,c(2,0))
  }
  
  ## 01 = no 2 offspring
  wh <- which(mating=="01")
  if(length(wh)) {
    offspring[wh] <- lapply(offspring[wh],mydrop,2)
  }
  
  ## 12 = no 0 offspring
  wh <- which(mating=="12")
  if(length(wh)) {
    offspring[wh] <- lapply(offspring[wh],mydrop,0)
  }
  
  if(first) # only take first aff per fam
    offspring <- lapply(offspring,"[",1)
  
  y <- data.frame(mating=mating,mating.ordered=mating.ordered,
                  offspring=sapply(offspring,paste,collapse=""),
                  ped=sapply(x,function(p) {unique(p$ped)}),
                  stringsAsFactors=FALSE)
  y <- subset(y,!(mating %in% c("00","22","11"))) # drop uninformative mating types
  y <- subset(y,offspring!="") ## drop families with no offspring without misinheritance
  
  if(verbose) {
    cat("\n",nrow(y)," families identified for testing. # aff offspring per family:\n",sep="")
    print(table(nchar(y$offspring)))
  }
  
  
  ## drop 212/122 matings
  
  z <- aggregate(rep(1,nrow(y)),
                 by=y[,c("mating","mating.ordered","offspring")],"sum")
  z <- z[order(z$mating,z$offspring,z$mating.ordered),]
  z$Imp <- as.integer(z$mating.ordered %in% c("01","02","12"))
  
  ## beta=M+P>1
  z$Ibeta=as.integer(z$mating %in% c("02","12"))*nchar(z$offspring)
  
  ## gamma = 1 if M+P==1
  ## gamma = -1 if M+P>2
  z$Igamma <- 0
  z$Igamma[z$mating=="01"] <- nchar(z$offspring)[z$mating=="01"]
  z$Igamma[z$mating=="12"] <- -1*nchar(z$offspring)[z$mating=="12"]
  
  n1 <- nchar(z$offspring) - nchar(gsub("1","",z$offspring))
  z$Ialpha=as.integer(n1)
  z$Ialpha[!(z$mating %in% c("01","02","12"))] <- 0
  
  ## model impr + mat geno
  m3 <- glm(Imp~Ialpha+Ibeta+Igamma - 1,data=z,
            family="binomial",
            weights=x)
  
  ## model mat geno only
  m2 <- glm(Imp~Igamma+Ibeta-1,data=z,
            family="binomial",
            weights=x)
  
  ## model impr only
  m1 <- glm(Imp~Ialpha-1,data=z,
            family="binomial",
            weights=x)
  
  ## model nothing (for completeness)
  m0 <- glm(Imp~-1,data=z,
            family="binomial",
            weights=x)
  
  M <- list(M3=m3,M2=m2,M1=m1,M0=m0)
  class(M) <- c(class(M),"origin")
  
  ## print matrix comparing models
  print.origin(M)
  
  ## return models
  invisible(M)
}

################################################################################

## print

print.origin <-
function(x,...) {
  mat <- cbind(deviance=sapply(x,"[[","deviance"),
               n.param=c(3,2,1,0),
               anova.Mi.M3=sapply(x,function(mod) anova(mod,x$M3,test="Chisq")$P[2]),
               anova.Mi.M0=sapply(x,function(mod) anova(mod,x$M0,test="Chisq")$P[2]))
  rownames(mat) <- c("M3: Imprinting + maternal genotype",
                     "M2: Maternal genotype only",
                     "M1: Imprinting only",
                     "M0: Neither effect")
  cat("\nSummary of model comparison\n")
  print(mat)
}

################################################################################

## handy functions

mydrop <-
function(x,drop) {
  wh <- which(x %in% drop)
  if(length(wh)) {
    return(sort(x[-wh]))
  } else {
    return(sort(x))
  }
}

unfactor <-
function(x) {
  if(!is.factor(x))
    return(x)
  return(levels(x)[x])
}


