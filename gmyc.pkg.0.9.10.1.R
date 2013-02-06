#GMYC species delimitaion method
#written in R package style
library(ape)

#diplay a list of GMYC species membership
spec.list <- function(res,second.peak=F) {
	tr <- res$tr
	spec <- tr$tip.label
	numtip <- length(tr$tip.label)
	
	if (second.peak==F) {
	max.mrca <- res$MRCA[[which.max(res$likelihood)]] + numtip} else
	{tmp<-table(cummax(res$likelihood))
		lik.peaks<-names(tmp[tmp>20])
		peak<-which(res$likelihood==lik.peaks[(length(lik.peaks)-1)])
		max.mrca <- res$MRCA[[peak]] + numtip}
	
	numspec <- length(max.mrca)
	
	nest.tip <- function(nod, tr) {
		
		tip <- c()
		child <- tr$edge[tr$edge[,1] == nod, 2]
		
		for (ch in child) {
			if (ch <= numtip) {
				tip <- c(tip, ch)
			} else {
				tip <- c(tip, nest.tip(ch, tr))
			}
		}
		return (tip)
	}
	
	res <- c()
	
	for (i in 1:length(max.mrca)) {
		tip.name <- tr$tip.label[nest.tip(max.mrca[i], tr)]
		res <- rbind(res, cbind(i, tip.name))	
	}
	
	
	if (length(spec[-match(res[,2], spec)]) != 0) {
		numspec <- numspec + 1
		for (s in spec[-match(res[,2], spec)]) {
			res <- rbind(res, cbind(numspec, s))
			numspec <- numspec + 1
		}
	}
	
	res <- data.frame(res)
	colnames(res) <- c("GMYC_spec", "sample_name")
	
	return (res)
}
	

#model comparison between 2 results of gmyc()
compare <- function(test1,test2) {
	chi.sq<-abs(2*(max(test1$likelihood)-max(test2$likelihood)))
	dfs<- 3*(length(test2$threshold.time[[which.max(test2$likelihood)]])-1)
	cat("Comparison of single and multiple threshold GMYC\n")
	cat("\tChi.sq",chi.sq,"Degrees of freedom",dfs,"\n", sep=" ")
	cat("\tSignificance", 1-pchisq(chi.sq,dfs), sep=" ")
}

#plot summary of species delimitation
plot.gmyc <- function(res, ask=T, second.peak=F,file.name=NA,height=96) {
	#res = result of GMYC
	#plot results of GMYC analysis; likelihood, LTT plot and tree with clusters colored 
	
	if (ask) {
		par(ask=ask)
	} else {
		par(mfrow=c(3,1))
	}
	
	if (!is.na(file.name)) {pdf(file=paste(file.name,"ltt&lik.pdf",sep=""))}
	
	if (second.peak==T) {
		tmp<-table(cummax(res$likelihood))
		lik.peaks<-names(tmp[tmp>20])
		peak<-which(res$likelihood==lik.peaks[(length(lik.peaks)-1)])}
	
	if (res[["method"]] == "single") {
		#lineage through time plot with threshold time
		ltt.plot(res$tree, log="y")
		if (second.peak==F) {
		abline(v=res$threshold.time[which.max(res$likelihood)], col = "red")} else
		{abline(v=res$threshold.time[peak], col = "red")} 
		#likelihood surface
		plot(res$threshold.time, res$likelihood, type="l", xlab="Time", ylab="likelihood")
				
					if (!is.na(file.name)) {dev.off(); pdf(height=height,file=paste(file.name,"clust.pdf",sep=""))}
				
		if (second.peak==F) {
		plot.cluster(res$tree, res$MRCA[[which.max(res$likelihood)]])} else
		{plot.cluster(res$tree, res$MRCA[[peak]])}
		
					if (!is.na(file.name)) {dev.off()}
	
	}  else if (res[["method"]] == "multiple" || res[["method"]] == "exhaustive") {
		#lineage through time plot with threshold time
		ltt.plot(res$tree, log="y")
		abline(v=res$threshold.time[[which.max(res$likelihood)]], col = "red")
		plot.cluster(res$tree, res$MRCA[[which.max(res$likelihood)]])
	# } else if (res[["method"]] == "multiple") {
		# plot.cluster(res$tree, res$MRCA[[which.max(res$likelihood)]])
	}
		
	#tree with clusters
	#plot.cluster1(res$tree, res$threshold.time[which.max(res$likelihood)])
	par(ask=F)
}

#display summary
summary.gmyc <- function(res, second.peak=F) {
	#res = result of GMYC
	#display summary of GMYC; likelihood values, chi-square test, estimated parameters, etc...
	
		if (second.peak==T) {
		tmp<-table(cummax(res$likelihood))
		lik.peaks<-names(tmp[tmp>20])
		peak<-which(res$likelihood==lik.peaks[(length(lik.peaks)-1)])}

	
	cat("Result of GMYC species delimitation\n")
	cat("\n\tmethod:\t", res[["method"]], sep="")
	cat("\n\tlikelihood of null model:\t", res$likelihood[1], sep="")
		if (second.peak==F) {
			cat("\n\tmaximum likelihood of GMYC model:\t", max(res$likelihood), sep="")} else
			{cat("\n\tmaximum likelihood of GMYC model:\t", res$likelihood[peak], sep="")}
	
	#chisq test
	if (second.peak==F) {
				LR <- 2*(max(res$likelihood)-res$likelihood[1])} else
				{LR <- 2*(res$likelihood[peak]-res$likelihood[1])}
	cat("\n\tlikelihood ratio:\t", LR, sep="")

	if (res[["method"]] == "single") {
		pvalue <- 1-pchisq(LR, 3)
	}  else if (res[["method"]] == "multiple" || res[["method"]] == "exhaustive") {
		pvalue <- 1 - pchisq(LR, 3 + length(res$threshold.time[[which.max(res$likelihood)]]) - 1)
	}
	
	cat("\n\tresult of LR test:\t", pvalue, ifelse(pvalue<0.001, "***", ifelse(pvalue<0.01, "**", ifelse(pvalue<0.05, "*", "n.s."))), sep="")
	
		if (second.peak==F) {
	cat("\n\n\tnumber of ML clusters:\t", res$cluster[which.max(res$likelihood)], sep="")
	cat("\n\tnumber of ML entities:\t", res$entity[which.max(res$likelihood)], sep="")
	if (res[["method"]] == "single") {	
		cat("\n\tthreshold time:\t", res$threshold.time[which.max(res$likelihood)], "\n", sep="")
	} else if (res[["method"]] == "multiple" || res[["method"]] == "exhaustive") {
		cat("\n\tthreshold time:\t", res$threshold.time[[which.max(res$likelihood)]], "\n", sep=" ")
	}	
	cat("\n")} else
	
	{cat("\n\n\tnumber of ML clusters:\t", res$cluster[peak], sep="")
	cat("\n\tnumber of ML entities:\t", res$entity[peak], sep="")
	if (res[["method"]] == "single") {	
		cat("\n\tthreshold time:\t", res$threshold.time[peak], "\n", sep="")
	} else if (res[["method"]] == "multiple" || res[["method"]] == "exhaustive") {
		cat("\n\tthreshold time:\t", res$threshold.time[[peak]], "\n", sep=" ")
	}	
	cat("\n")}
}

#how to remove global variables...........
#01/07/08
#finish removeing glabal varibles...
#03/07/08

#plot function with multiple threshold times
plot.cluster <- function(tr, lthresh, show.tip.label=T, show.node.label=F, cex=0.5) {
	#tr = tree of ape tree class
	##thresh = a vector of threshold values
	
	numnod <- tr$Nnode
	numtip <- length(tr$tip.label)

	cdat<-array(1,2*numnod)
	ndat<-array("",numnod)
	
	bt <- -branching.times(tr)
	
	## Nested nodes
	nest.nodes<-function(tr, x ,p=0) {
		#print(paste("nest.nodes() called with x=", as.character(x)))
		numtip <- length(tr$tip.label)
		   
		nods<-array(NA,0)
		desc<-as.integer(tr$edge[,2][tr$edge[,1]==x])
		 
		if (desc[1] > numtip) {
			nods <- c(nods, desc[1], nest.nodes(tr, desc[1]))
		}
		if (desc[2] > numtip) {
			nods <- c(nods, desc[2], nest.nodes(tr, desc[2]))
		}
		   
		if (length(nods)>0) {
			return (nods) 
		} else {
			return(NULL)
		}      
	}
	
	threshold.group <- function(mrcas) {	#function to obtain multiple threshold times from a set of mrcas
		parent <- tr$edge[,1]
		child <- tr$edge[,2]
	
		thresh.group <- list()
		thresh.time <- c()
	
		mrcas <- mrcas + numtip
		k <- 1
		
		while (TRUE) {
			times <- bt[mrcas-numtip]
			thresh1.time <- min(times)
			thresh1.node <- mrcas[which.min(times)]
			
			mrcas <- mrcas[-which.min(times)]
			
			if (length(mrcas) == 0) { 
				thresh.time <- c(thresh.time, thresh1.time)	###??????????????????????? last MRCA ???
				thresh.group[[k]] <- thresh1.node
				break 
			}	
			
			member <- thresh1.node
			del <- c()
			for (i in 1:length(mrcas)) {
				#print(mrcas[i])
				#print(length(mrcas))
				par.nod <- parent[child==mrcas[i]]
				t.par <- bt[par.nod-numtip]
				#print(c(mrcas[i], par.nod, t.par, length(mrcas)))
				if (t.par < thresh1.time) {
						member <- c(member, mrcas[i])
						del <- c(del, i)
				}
				
			}
			thresh.time <- c(thresh.time, thresh1.time)
			thresh.group[[k]] <- member
			#print(thresh.time)
			k <- k+1
			
			if (length(del) != 0) {	mrcas <- mrcas[-del]}
			
			if (length(mrcas) == 0) { break }	
		}
		
		return (thresh.group)
	}
		
	group <- threshold.group(lthresh)
	

	#print(group)
	colors <- rainbow(length(group))
	
	k <- 1
	for (g in group) {
		n.col.type <- rep(0, numnod)
		for (j in 1:length(g)) {
			n.col.type[g[j]-numtip] <- 2
			n.col.type[nest.nodes(tr, g[j])-numtip] <- 1
			
		}
		cdat[match(tr$edge[,1],which((n.col.type==1)|(n.col.type==2))+numtip)>0]<-colors[k]
		#ndat[nod.type==2]<-(1:numnod)[n.col.type==2]
		k <- k + 1
	}
	
	
	plot(tr,edge.color=cdat, use.edge.length=1,show.node.label=show.node.label, show.tip.label=show.tip.label, no.margin=F, cex=cex)
}

gmyc <- function(tr, method = "single", interval = c(0, 5), quiet = FALSE) {
	#run GMYC
	#contains;
	#reading data
	#applying single threshold
	#applying multiple threshold
	
	#local.env <- globalenv()	#
	local.env <- environment()
	
	### function to read tree and data
	read.data<-function (z=1) {
		#tr<<-tt[[5]]
		#if (length(tt)>1) {tr<<-tt[[z]]} else {tr<<-tt}
	 	bt<- -branching.times(tr)
		#Set any ages after present day or below a certain limit to finite age
	 	bt[bt>-0.000001]<- -0.000001
	 	names(bt)<-NULL
	 	#sb<-sort(bt)
		
		#assing variables to parent environment
		#02/07/08...
		assign("bt", bt, envir=local.env)
		assign("sb", sort(bt), envir=local.env)
		
	 	#19/11/07
	 	##add number of tip and number of all nodes###
		#03/07/08
		##move to assign function
		assign("numnod", length(bt), envir=local.env)
		assign("numtip", length(tr$tip.label), envir=local.env)
		assign("numall", length(bt)+length(tr$tip.label), envir=local.env)
		
		internod<-sb[2:numnod]-sb[1:numnod-1]
	 	internod[numnod]<-0-sb[numnod]
		
		assign("internod", internod, envir=local.env)
		##lists of which nodes are nesting or nested within each node, and arrays with sisters and desc
			  
		#change parameters of sapply(nesting/nested.nodes) -1:numnod -> numtip+1:numall
		assign("nesting", sapply((numtip+1):numall,nesting.nodes), envir=local.env)
		assign("nested", sapply((numtip+1):numall,nest.nodes), envir=local.env)
			
		numnested<-array(NA,numnod)
			  
		numnesting<-array(NA,numnod)
		des<-matrix(NA,nrow=numnod,ncol=2)
		sis<-array(NA,numnod)
			  
		for (j in (1:numnod)) {
			numnested[j] <- length(nested[[j]])
			numnesting[j] <- length(nesting[[j]])
			des[j,]<-as.integer(tr$edge[,2])[as.integer(tr$edge[,1])== j+numtip]
			if (des[j,1] > numtip) { sis[des[j,1]-numtip] <- des[j,2]}	#??????????
			if (des[j,2] > numtip) { sis[des[j,2]-numtip] <- des[j,1]}	#??????????
		}
		
		assign("numnested", numnested, envir=local.env)
		assign("numnesting", numnesting, envir=local.env)
		
		des[des <= numtip] <- NA
		
		assign("des", des, envir=local.env)
		assign("sis", sis, envir=local.env)
		
		##ancestral nodes, column 1, of each node, column 2
			ancs<-cbind(tr$edge[pmatch((1:numnod+numtip),tr$edge[,2]),1],(1:numnod+numtip))
		##node ages for the nodes in ancs
			bt.ancs<-cbind(bt[ancs[,1]-numtip],bt[ancs[,2]-numtip])
		assign("bt.ancs", bt.ancs, envir=local.env)

				
	}

	## Nested nodes
	nest.nodes<-function(x,p=0) {
		#print(paste("nest.nodes() called with x=", as.character(x)))
		numtip <- length(tr$tip.label)
		   
		nods<-array(NA,0)
		desc<-as.integer(tr$edge[,2][tr$edge[,1]==x])
		 
		if (desc[1] > numtip) {
			nods <- c(nods, desc[1], nest.nodes(desc[1]))
		}
		if (desc[2] > numtip) {
			nods <- c(nods, desc[2], nest.nodes(desc[2]))
		}
		   
		if (length(nods)>0) {
			return (nods) 
		} else {
			return(NULL)
		}      
	}

	## Nesting nodes
	nesting.nodes<-function(x,p=0) {  
		#print(paste("nesting.nodes() called with x=", as.character(x)))
		numtip <- length(tr$tip.label)
		
		nod<-array(NA,0)
		  
		##change -1 -> numtip+1
		if (x >= numtip+2)  { 
			anc <- as.integer(tr$edge[,1][tr$edge[,2]==x])
		}  else  { 
			anc <- 1	#?????????????????????? 1 ?? 
		}	
		  
		if (anc	>= numtip+2) {
			nod <- c(nod, anc, nesting.nodes(anc)) 
		}
		if (anc == numtip+1) {
			nod <- c(nod, anc)
		}
		  
		if (length(nod)>0)  {
			return (nod)
		} else {
			return(NULL)
		}
	}
	  

	###Derived from SpeciesTestLikFunctions2005.3

	##replaced global variables with list.
	##02/07/08...
	l.prep<-function () {
		  #print ("l.prep() called")
		  assign("n", length(mrca.nodes), envir=local.env)
		  
		  i.mat<-matrix(0,ncol=numnod,nrow=(n+1))
		  s.nod<-matrix(0,ncol=numnod,nrow=(n+1))
		  
		  assign("nod", nod.type[order(bt)], envir=local.env)
		  
		  ##set of coalescent processes for each mrca node
		  for (i in (1:n)) {
			  s.nod[i,mrca.nodes[i]]<-2
			  
			  #change s.nod[i, -nested[[mrca.nodes[i]]]> s.nod[i, nested[[mrca.nodes[i]][-numtip]
			  if (!is.null(nested[[mrca.nodes[i]]])) { s.nod[i, nested[[mrca.nodes[i]]] - numtip] <- 1}
			  
			  s.nod[i,]<-s.nod[i,order(bt)]	##????
			  i.mat[i,][s.nod[i,]==2]<-2
			  i.mat[i,][s.nod[i,]==1]<-1
			  i.mat[i,]<-cumsum(i.mat[i,])
		  }
		  s.nod[s.nod==2]<-1	###???????????????
		  ##transform number of lineages to coalescent type
		  i.mat<-i.mat*(i.mat-1)
		  
		  
		  #one diversification process
		  s.nod[n+1,]<-nod==0
		  i.mat[n+1,nod==0]<-1
		  i.mat[n+1,nod==2]<--1
		  i.mat[n+1,]<-cumsum(i.mat[n+1,])+1
		  
		  assign("s.nod", s.nod, envir=local.env)
		  assign("i.mat", i.mat, envir=local.env)
		  #write(i.mat, ncolumns=numnod, "")
		  #cat(s.nod, "\n")
		  #cat(i.mat, "\n")
	}

	###LIKELIHOOD FUNCTIONS - NULL AND MINIMUM MODEL

	##1) Calc likelihood for null model, generalized Yule model with scaling power, p
	l.null<-function(p=1) {
		#print("l.null() called")
		i.div<-2:(numnod+1)
		i.div<-i.div*(i.div-1)
		
		lambda.div <- numnod / sum(internod * i.div^p)
		assign("lambda.div", lambda.div, envir=local.env)
		
		lik <- i.div^p * lambda.div * exp(-i.div^p * lambda.div * internod)
		
		return (sum(log(lik)))
	}

	#Calcul
	l.min2<-function (q=c(1, 1)) {
	  #print("l.min2() called")

	  ##prob of any event happening
		p <- c(rep(q[1],n),q[2])
		lambda <- sum(s.nod[1:n,]) / sum(i.mat[1:n,]^p[-(n+1)]%*%internod)
		lambda <- c( rep(lambda,n),sum(s.nod[n+1,]) / i.mat[n+1,]^p[n+1]%*%internod )
				
		b<-t(i.mat^p)%*%lambda
		lik<-b*exp(-b*internod)	
		lambda<-lambda
		
		assign("b", b, envir=local.env)
		assign("lik", lik, envir=local.env)
		assign("lambda", lambda, envir=local.env)
		
		return(sum(log(lik)))
	}

	gmyc.likelihood <- function() {	
		 ###APPLY SIMPLE THRESHOLD MODEL - ASSUME ALL CLUSTERS HAVE SAME BRANCHING PARAMETERS, LAMBDA AND P
		l.prep()
		
		l.results <- rep(NA, 7)
		#x<-optim(c(1, 1),l.min2,method = "BFGS",control=list(fnscale=-1))	#likelihood for GMYC model
		
		#x <- optim(temp.params, l.min2, method = "Nelder-Mead",control=list(fnscale=-1))
		
		x <- optim(c(1, 1), l.min2, method = "Nelder-Mead",control=list(fnscale=-1))
		
		#x <- optim(temp.params, l.min2, method = "L-BFGS", control=list(fnscale=-1))

		
		l.results[1:6]<-c(x$value, lambda[n+1], lambda[1], x$par[2],x$par[1],as.integer(sum(s.nod[n+1,])+1))
		l.results[7]<-n
		#results[j,]<-nod.type
		assign("temp.params",x$par,envir=local.env)
		#print(l.results)
		
		return (l.results)
	}
	
	#exhaustive search of MRCAs...03/12/08
	gmyc.exhaustive <- function() {
		
		#obtaining all combination of MRCAs
		re.combn <- function(node=numtip+1, tr=tr, se=numtip+1) {
			
			parents <- tr$edge[,1]
			children <- tr$edge[,2]

			temp.list <- list()
				
			ch <- children[which(parents==node)]
			se <- se[-(which(se==node))]
			se <- c(se, ch)
			
			temp.list <- c(temp.list, list(se))
			#print(se)
			
			if (any(ch > numtip)) {
				for (cc in ch) {
					temp.temp.list <- temp.list
					if (cc > numtip) {
						for (tl in temp.temp.list) {
							temp.list <- c(temp.list, re.combn(node=cc, tr=tr, se = tl))
						}
					}
				}
			} else {
				
			}
			return(temp.list)
		}
		
		remove.tips <- function(vec) {
			return(vec[vec > numtip])
		}

		MRCA <- lapply(re.combn(tr=tr), remove.tips)
		MRCA <- MRCA[-length(MRCA)]
		#print(MRCA)
		if (!quiet) {
			cat("exhaustive search", "\npossible combinations of MRCAs are", length(MRCA), "\n")
		}
		
		l.results<-c()#matrix(ncol=7, nrow=length(MRCA))
		
		l.mrca <- lapply(MRCA, "-", numtip)
		
		#results<-matrix(nrow=length(MRCA) ,ncol=numnod)

		
		
		assign("temp.params",c(1,1),envir=local.env)
		
		count <- 1
		for (nodes in MRCA) {
			mrca.nodes <- nodes - numtip
			
			assign("mrca.nodes", mrca.nodes, envir=local.env)
					
			#print(mrca.nodes)
			nod.type <- rep(0, numnod)
					
			#print(length(nod.type))
					
			for (j in 1:length(mrca.nodes)) {
				nod.type[mrca.nodes[j]] <- 2
				nod.type[nested[[mrca.nodes[j]]]-numtip] <- 1	#nested[[xxx]] - numtip
			}
			assign("nod.type", nod.type, envir=local.env)
			
			l.results <- rbind(l.results, gmyc.likelihood())
			if (!quiet) {
				#cat("testing nodes", mrca.nodes, "\n")
				cat(l.results[count,1], "\n")
				count <- count+1
			}
		}
		
		return (list(l.results, l.mrca))
		
	}
	
	
	gmyc.single <- function() {
		##SET UP RESULTS MATRICES
		l.results<-matrix(ncol=7,nrow=(nthresh))
		l.mrca <- list()
		results<-matrix(nrow=nthresh,ncol=numnod)

		###WORK OUT RESULTS FOR NULL MODEL - JUST ONE CLUSTER
		#lambda.div <- NULL
		#x<-optimise(l.null,c(0,5),maximum=1)		#likelihood for null model in gmyc.single()
		
		x <- optimise(l.null, interval=interval, maximum=1)
		
		l.results[1,c(6:7,1:2,4)]<-c(as.integer(1),as.integer(1),x$objective, lambda.div, x$maximum)
		
		ml <- 0
		
		if (!quiet) { cat("node\t", "T\t", "loglik\n") }

		stthresh <- 2
		while (sb[stthresh] == sb[1]) {
			stthresh <- stthresh + 1
		}

		assign("temp.params",c(1,1),envir=local.env)

		#stthresh <- 2
		for (j in (stthresh:nthresh)) {
			 
			 ##THIS FOR-LOOP SLIDES THROUGH THE TREE AND DEFINES EACH NODE AGE IN TURN AS THE THRESHOLD
			 ##FOR THE SWITCH BETWEEN CLADOGENESIS AND COALESCENCE WITHIN CLUSTERS
			threshy<-sb[j]
						
			tmp<-(bt.ancs[,1]<threshy)&(bt.ancs[,2]>=threshy)
			nod.type<-tmp+(bt>=threshy)
			mrca.nodes<-which(nod.type==2)			
			#print(nod.type)
			
			if (nod.type[1]==1) nod.type[1] <- 2
			mrca.nodes <- which(nod.type==2)
			
			assign("mrca", mrca, envir=local.env)
			assign("mrca.nodes", mrca.nodes, envir=local.env)
			assign("nod.type", nod.type, envir=local.env)
			
			l.mrca[[j]] <- mrca.nodes
			
			l.results[j,] <- gmyc.likelihood() 
			
			if (!quiet) { cat (j, threshy, l.results[j,1],"b:", b, "lambda", lambda, "l3", l.results[j,3], "l4", l.results[j,4], "l5", l.results[j,5], "\n") }

		}
		 #results[j,]<-nod.type
		return(list(l.results, l.mrca))
	}
	
	#optimizing multiple MRCAs
	gmyc.multi <- function() {
		assign("temp.params",c(1,1),envir=local.env)
		renew.mrca <- function(n) {	#optimization of MRCA nodes with dividing-fusing algoritm
			#obtain a new set of mrca nodes from the current set.
			#fusing a pair
			f.renew <- function(n) {
				parents <- tr$edge[,1]
				children <- tr$edge[,2]
				
				renewed <- list()
				
				if (length(n) > 1) {
					pair <- combn(n, 2)
					
					for (i in 1:length(n)) {
						parent.node <- parents[children==n[i]]
						sibling.nodes <- children[parents==parent.node]
						sibling <- sibling.nodes[sibling.nodes != n[i]]
						
						if (sibling <= numtip) {
							renewed <- c(renewed, list(c(n[-i], parent.node)))
						} else if (any(n == sibling)) {
							renewed <- c(renewed, list(c(n[-c(i, which(n==sibling))] , parent.node)))
						}
					}
				}
				return (unique(renewed))
			}

			#dividing a node 
			d.renew <- function(n) {
				parents <- tr$edge[,1]
				children <- tr$edge[,2]
				
				renewed <- list()
				if (length(n) > 1) {
					for (i in 1:length(n)) {
						child.nodes <- children[parents==n[i]]
						if (any(child.nodes > numtip)) {
							renewed <- c(renewed, list(c(n[-i], child.nodes[child.nodes > numtip])))
						} 
					}
				}
				return (renewed)
			}

		
			return (c(f.renew(n), d.renew(n)))
		}
		
		##rewritten to avoid redundant part of starting point setting
		##search for 3 part 11/09/08
		##revised 15/09/08
		select.start.re <- function(time, ml) {	#recurive division of starting point of d-f algoritm, avoiding trying parts with less improvement
			result <- c()
			result.mrca <- list()
			
			m <- mean(time)
			
			left <- time[time < m]
			right <- time[time >= m]	
			mid <- time[time >= mean(left) & time < mean(right)]
			
			part <- list(left, mid, right)
			
			start <- c(mean(left), m, mean(right))
			
			if (!any(is.na(start)) || length(unique(start)) == 1) {	#stop recursive division when length of starting points are less than 1 ...03/12/08
				improve <- c(FALSE , FALSE , FALSE)
				num.improve <- c(0, 0, 0)
				
				#print(start)
				
				for (i in 1:length(start)) {
					if (!quiet) {cat("start at", start[i], "\n")}
					
					temp <- renew.likelihood.from.thresh(start[i])
					lik <- temp[[1]]
					mrca <- temp[[2]]
					
					if (any(lik[,1] > ml)) {
						if (!quiet) {cat("improvement found\n")}
						result <- rbind(result, lik[which(lik[,1] > ml),])
						result.mrca <- c(result.mrca, mrca[which(lik[,1] > ml)])
						improve[i] <- TRUE
						num.improve[i] <- length(lik[lik[,1] > ml, 1])
					}
				}
				if (!quiet) {cat(num.improve, "\n")}
				
				if (!all(num.improve == 0)) {
					temp <- select.start.re(part[[which.max(num.improve)]], max(result[,1]))
					result <- rbind(result, temp[[1]])
					result.mrca <- c(result.mrca, temp[[2]])
				}
			}
			
			return(list(result, result.mrca))
			
		}
		
		
	
		
		##rewritten to avoid redundant part of starting point setting
		##search for 3 part 
		select.start <- function(time, ml) {	#recurive division of starting point of d-f algoritm
			result <- c()
			result.mrca <- list()
			
			m <- mean(time)
			
			left <- time[time < m]
			right <- time[time >= m]	
			mid <- time[time >= mean(left) & time < mean(right)]
			
			part <- list(left, mid, right)
			
			start <- c(mean(left), m, mean(right))
			#improve <- c(FALSE , FALSE , FALSE)
			
			for (i in 1:length(start)) {
			
				if (!quiet) {cat("start at", start[i], "\n") }
				
				temp <- renew.likelihood.from.thresh(start[i])
				lik <- temp[[1]]
				mrca <- temp[[2]]
				
				if (any(lik[,1] > ml)) {
					if (!quiet) { print("improvement found\n") }
					result <- rbind(result, lik[which(lik[,1] > ml),])
					result.mrca <- c(result.mrca, mrca[which(lik[,1] > ml)])
					
					if (length(part[[i]]) > 1) {
						if (!quiet) {cat("recursively trying part of", start[i], "\n")} 
						temp <- select.start(part[[i]], max(lik[,1]))
						result <- rbind(result, temp[[1]])
						result.mrca <- c(result.mrca, temp[[2]])
					}
				}
			}
			return(list(result, result.mrca))
			
		}
		
		renew.likelihood.from.thresh <- function(start) {	#find ML solution with renew.mrca and select.start
			l.results <- c()
			
			l.mrca <- list()	######add, 08/09/08
			
			#start <- mean(sb)
			#cat("start: ", start, "\n")
			max.lik <- NA
			mrca<-array(FALSE,numnod)
			mrca[2:numnod]<- (bt[as.integer(tr$edge[,1][tr$edge[,2]>numtip]) - numtip]<start)&(bt[as.integer(tr$edge[,2][tr$edge[,2]>numtip]) - numtip]>=start)	##???
						
			nod.type <- (bt>=start)+mrca

			if (nod.type[1]==1) nod.type[1] <- 2
			mrca.nodes <- which(nod.type==2)
			
			assign("mrca", mrca, envir=local.env)
			assign("mrca.nodes", mrca.nodes, envir=local.env)
			assign("nod.type", nod.type, envir=local.env)
			

			initial.mrca <- mrca.nodes
			while (TRUE) {
				
				found <- FALSE
				max.mrca <- initial.mrca
				for (nodes in renew.mrca(initial.mrca+numtip)) {
					mrca.nodes <- nodes - numtip
					if (mrca.nodes[[1]] != 1) {
						
						assign("mrca.nodes", mrca.nodes, envir=local.env)
						
						#cat(mrca.nodes, "\n")
						nod.type <- rep(0, numnod)
						
						#print(length(nod.type))
						
						for (j in 1:length(mrca.nodes)) {
							nod.type[mrca.nodes[j]] <- 2
							nod.type[nested[[mrca.nodes[j]]]-numtip] <- 1	#nested[[xxx]] - numtip
						}
						assign("nod.type", nod.type, envir=local.env)
					
						#print(c(length(nod.type), numnod))
						#print(nod.type)
						
						res <- gmyc.likelihood()
						#print(length(res))
						
						if (!is.na(max.lik)) {
							if (max.lik <= res[1]) {
								max.lik <- res[1]
								l.results <- rbind(l.results, res)
								max.mrca <- nodes - numtip
								
								l.mrca <- c(l.mrca, list(max.mrca))	######add, 08/09/08
								
								if (!quiet) {
									cat(max.lik, "\n")
								}
								#print("improved")
								found <- TRUE
							}
						} else {
							max.lik <- res[1]
							l.results <- rbind(l.results, res)
							max.mrca <- nodes - numtip
							
							l.mrca <- c(l.mrca, list(max.mrca))	######add, 08/09/08
							
							found <- TRUE
						}
					}
				}
				
				if (found) {
					initial.mrca <-  max.mrca
				} else { 
					if (!quiet) {cat("break\n")}
					break
				}
			}
			
			rownames(l.results) <- NULL
			
			return (list(l.results, l.mrca))
		
		}
		
		
		##SET UP RESULTS MATRICES
		l.results<-c()#matrix(ncol=7,nrow=(nthresh))
		#results<-matrix(nrow=nthresh,ncol=numnod)
		l.mrca <- list()
		
		
		###WORK OUT RESULTS FOR NULL MODEL - JUST ONE CLUSTER
		#lambda.div <- NULL
		#x<-optimise(l.null,c(0,5),maximum=1)	#likelihood of GMYC model in gmyc.multi()
		
		x <- optimise(l.null, interval=interval, maximum=1)
		
		l.results<- rbind(l.results, c(x$objective, lambda.div, NA, x$maximum, NA, as.integer(1),as.integer(1)))	#[1,c(6:7,1:2,4)]
		l.mrca <- c(l.mrca, list(numtip+1))
		
		if (!quiet) {
			cat("null likelihood\n")
			cat(l.results[1, 1], "\n")
		}
		
		###########################
		##get a set of renewed MRCAs for a threshold
		##10/07/08...
		##########################
		#start <- mean(sb)
		#l.results <- rbind(l.results, renew.likelihood.from.thresh(start))
		#l.results <- select.start1(sb, 0)
		
		
		
		#temp <- select.start(sb, l.results[1, 1])	#select.start(sb, 0) ????????
		
		if (!quiet) {
			cat("\nGMYC likelihood\n")
		}
		
		temp <- select.start.re(sb, l.results[1, 1])	
		
		l.results <- rbind(l.results, temp[[1]])
		l.mrca <- c(l.mrca, temp[[2]])
		
		return (list(l.results, l.mrca))
		
	}
	
	#######################
	##PROGRAM ENTRY POINT##
	#######################
	##check if the input tree is appropriate
	if (!is.ultrametric(tr)) {
		#stop("Your ultrametric tree is not ultrametric, please check")
	}
	if (!is.binary.tree(tr)) {
		stop("Your input tree is not fully bifurcating, please resolve with zero branch lengths")
	}
		
	##SET UP DATA IN FORMAT USED LATER
	read.data()
		
	ntrees<-1

	##numnod <- -min() -> numnode<-max()-numtip
	numnod <- max(as.integer(tr$edge[,1])) - numtip

	##which is last node to try for threshold model
	nthresh<-numnod
	
	#run analysis
	if (method == "single" || method == "s") {
		method <- "single"
		temp <- gmyc.single()
		l.results <- temp[[1]]
		l.mrca <- temp[[2]]
		
	} else if (method == "multiple" || method == "m") {	
		method <- "multiple"
		temp <- gmyc.multi()
		
		l.results <- temp[[1]]
		l.mrca <- temp[[2]]
	} else if (method == "exhaustive" || method == "e") {
		method <- "exhaustive"
		temp <- gmyc.exhaustive()
		l.results <- temp[[1]]
		l.mrca <- temp[[2]]
	} else {
		stop("Invalid name of optimiztion method. Only single(s) or multiple(m) and exhaustive(e) are available.")
	}
	
	#colnames(l.results)<-c("likelihood","lambda1","lambda2","p1","p2","num entities","num clusters")
	cat("\n", date(), "\n", sep="")

	#plot.cluster1(tr, sb[which.max(l.results[,1])])
	cat("finish.\n")

	#making result data structure...
	result <- list()
	result[["method"]] <- method
	result[["likelihood"]] <- l.results[,1]
	result[["parameters"]] <- l.results[,2:5]
	colnames(result[["parameters"]]) <- c("lambda.div", "lambda.coal", "p.div", "p.coal")
	result[["entity"]] <- l.results[,6]
	result[["cluster"]] <- l.results[,7]
	
	if (method == "single") {
		result[["MRCA"]] <- l.mrca
		result[["threshold.time"]] <- sb
	# } else if (method == "exhaustive") {
		# result[["MRCA"]] <- l.mrca
		# result[["threshold.time"]] <- NA
	} else if (method == "multiple" || method == "exhaustive") {
		
		reduce.threshold <- function(mrcas) {	#function to obtain multiple threshold times from a set of mrcas
 			parent <- tr$edge[,1]
			child <- tr$edge[,2]
		
			thresh.group <- list()
			thresh.time <- c()
		
			mrcas <- mrcas + numtip
			k <- 1
			
			while (TRUE) {
				times <- bt[mrcas-numtip]
				thresh1.time <- min(times)
				thresh1.node <- mrcas[which.min(times)]
				
				mrcas <- mrcas[-which.min(times)]
				
				if (length(mrcas) == 0) { 
					thresh.time <- c(thresh.time, thresh1.time)	###??????????????????????? last MRCA ???
					thresh.group[[k]] <- thresh1.node
					break 
				}	
				
				member <- thresh1.node
				del <- c()
				for (i in 1:length(mrcas)) {
					#print(mrcas[i])
					#print(length(mrcas))
					par.nod <- parent[child==mrcas[i]]
					t.par <- bt[par.nod-numtip]
					#print(c(mrcas[i], par.nod, t.par, length(mrcas)))
					if (t.par < thresh1.time) {
							member <- c(member, mrcas[i])
							del <- c(del, i)
					}
					
				}
				thresh.time <- c(thresh.time, thresh1.time)
				thresh.group[[k]] <- member
				
				k <- k+1
				
				if (length(del) != 0) {	mrcas <- mrcas[-del]}
				
				if (length(mrcas) == 0) { break }	
			}
			
			return (thresh.time)
		}
		
		result[["MRCA"]] <- l.mrca
		result[["MRCA"]] <- l.mrca
		result[["threshold.time"]] <- lapply(l.mrca, reduce.threshold)
	
		#result[["threshold.time"]] <- l.mrca
	}
	#result[["threshold.time"]] <- sb
	result[["tree"]] <- tr
	
	class(result) <- "gmyc"
	
	return (result)
}


test.tr <- read.tree(text="(((((((spec1.5:0.00051250,spec1.4:0.00051250):0.00057917,spec1.3:0.00109166):0.00243083,spec1.2:0.00352249):0.00541207,spec1.1:0.00893457):0.59264486,(spec18.3:0.01905108,(spec18.2:0.00871202,((spec18.1:0.00616976,spec18.5:0.00616976):0.00003963,spec18.4:0.00620939):0.00250263):0.01033906):0.58252835):0.35302356,(((((((spec29.4:0.11841156,spec6.2:0.11841156):0.03524040,(spec6.1:0.01532261,spec6.4:0.01532261):0.13832935):0.00230500,(spec6.3:0.01214206,spec6.5:0.01214206):0.14381489):0.03163809,((spec26.3:0.01281688,(spec26.4:0.00261800,((spec26.5:0.00004366,spec26.1:0.00004366):0.00191472,spec26.2:0.00195838):0.00065963):0.01019887):0.05207959,(((spec29.5:0.00042606,spec29.3:0.00042606):0.00842262,spec29.2:0.00884869):0.00292546,spec29.1:0.01177415):0.05312232):0.12269859):0.23610918,((((spec14.1:0.02358839,spec14.4:0.02358839):0.19823059,(spec14.3:0.22008050,((spec14.2:0.02273357,spec14.5:0.02273357):0.14713683,((((spec23.3:0.00948647,spec23.2:0.00948647):0.00129273,spec5.2:0.01077919):0.00321286,spec23.1:0.01399206):0.08719957,((spec5.1:0.00924401,spec23.5:0.00924401):0.06218262,(spec5.5:0.01499331,((spec5.4:0.00191824,spec5.3:0.00191824):0.01179373,spec23.4:0.01371197):0.00128134):0.05643333):0.02976500):0.06867872):0.05021010):0.00173847):0.03945475,((spec16.1:0.00161333,((spec16.5:0.00011218,spec16.2:0.00011218):0.00123912,spec16.4:0.00135129):0.00026204):0.00222676,spec16.3:0.00384009):0.25743363):0.03992391,(((((((spec3.2:0.00915008,spec3.5:0.00915008):0.00648934,spec4.4:0.01563942):0.00351967,spec3.1:0.01915909):0.00011569,(spec4.2:0.01026244,spec4.1:0.01026244):0.00901235):0.00083787,((spec3.3:0.00387522,spec3.4:0.00387522):0.01560206,spec4.3:0.01947728):0.00063537):0.02538681,spec4.5:0.04549946):0.04284882,((spec25.4:0.00038731,spec25.3:0.00038731):0.01176646,((spec25.5:0.00281087,spec25.1:0.00281087):0.00325763,spec25.2:0.00606851):0.00608527):0.07619455):0.21284936):0.12250655):0.18909264,(spec27.5:0.12452891,((spec27.1:0.00733008,spec27.4:0.00733008):0.01527155,(spec27.3:0.00988906,spec27.2:0.00988906):0.01271257):0.10192728):0.48826805):0.01921807,((((spec22.4:0.02967842,(spec22.2:0.01161063,spec22.3:0.01161063):0.01806779):0.04159766,spec22.1:0.07127608):0.08331611,spec22.5:0.15459219):0.16630259,(((spec19.4:0.00464504,spec19.3:0.00464504):0.02484778,spec19.5:0.02949282):0.00949410,(spec19.1:0.03892966,spec19.2:0.03892966):0.00005727):0.28190785):0.31112025):0.32258797):0.04539696,(((spec12.4:0.02417556,((spec12.1:0.00310453,(spec12.2:0.00257384,spec12.5:0.00257384):0.00053069):0.01365586,spec12.3:0.01676039):0.00741517):0.76361499,((((spec9.3:0.00547249,(spec9.5:0.00280129,(spec9.2:0.00178041,spec9.4:0.00178041):0.00102088):0.00267120):0.00068273,spec9.1:0.00615522):0.22011967,(spec21.5:0.03149763,((spec21.1:0.00350416,spec21.4:0.00350416):0.01571562,(spec21.2:0.00139018,spec21.3:0.00139018):0.01782960):0.01227785):0.19477721):0.03361161,((spec24.1:0.04534247,((spec24.5:0.01492007,(spec24.2:0.01085146,spec24.4:0.01085146):0.00406860):0.00968471,spec24.3:0.02460477):0.02073770):0.20032535,(((spec7.5:0.01225632,((spec7.1:0.00495471,spec7.2:0.00495471):0.00533212,spec7.3:0.01028683):0.00196949):0.05203172,spec7.4:0.06428804):0.13248343,((spec17.5:0.00875739,(spec17.1:0.00745892,spec17.4:0.00745892):0.00129847):0.00964934,(spec17.2:0.01495206,spec17.3:0.01495206):0.00345467):0.17836479):0.04889626):0.01421873):0.52790396):0.11321897,(((((((spec10.4:0.00191648,spec10.2:0.00191648):0.01281827,spec15.3:0.01473474):0.08756142,((spec15.1:0.01465863,(spec15.5:0.00153143,spec15.2:0.00153143):0.01312720):0.01771219,((spec10.5:0.01849071,(spec15.4:0.01536108,spec10.3:0.01536108):0.00312962):0.00323962,spec10.1:0.02173033):0.01064049):0.06992534):0.02130853,(spec2.4:0.03522054,(spec2.3:0.02253908,((spec2.5:0.00047005,spec2.2:0.00047005):0.01595391,spec2.1:0.01642395):0.00611513):0.01268146):0.08838420):0.10934199,((spec11.1:0.03444781,spec11.3:0.03444781):0.10196334,((spec11.2:0.01964460,spec11.4:0.01964460):0.00606147,spec11.5:0.02570607):0.11070509):0.09653558):0.06613387,((((spec28.2:0.02149857,(spec28.5:0.01336447,spec28.3:0.01336447):0.00813411):0.01894941,(spec28.1:0.01066503,spec28.4:0.01066503):0.02978296):0.15742330,(((spec20.1:0.04164863,spec20.3:0.04164863):0.06643827,(spec20.2:0.05174408,spec20.5:0.05174408):0.05634282):0.04739927,spec20.4:0.15548617):0.04238511):0.01445947,((spec30.3:0.00507569,(spec30.1:0.00092926,spec30.5:0.00092926):0.00414643):0.00194319,(spec30.2:0.00097875,spec30.4:0.00097875):0.00604013):0.20531192):0.08674981):0.46439420,((spec8.3:0.12432118,spec8.4:0.12432118):0.08763767,(spec8.5:0.15450705,(((spec13.1:0.02806162,spec13.5:0.02806162):0.01998279,(spec13.4:0.01187646,(spec13.2:0.00454374,spec13.3:0.00454374):0.00733272):0.03616795):0.09208632,(spec8.1:0.12742608,spec8.2:0.12742608):0.01270465):0.01437632):0.05745180):0.55151600):0.13753463):0.09899048);")

####

