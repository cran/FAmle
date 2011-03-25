metropolis <-
function(model,iter=1000,tun=2,trans.list=NULL,start=NULL,variance=NULL,prior=NULL,burn=0,uniroot.interval=c(-100,100))
{
	x <- model$x.info[,'x']
	k <- model$k
	if(!is.null(trans.list) & length(trans.list) != k)
		stop('the length of \'trans.list\' must match the number of unknown parameters!')
	else if(is.null(trans.list)) trans.list <- lapply(as.list(1:k),function(g) function(x) x)
	log.like <- function(param)
	{
		param <- as.vector(param)
		for(i in 1:k) param[i] <- trans.list[[i]](param[i])
		distr(x,model$dist,param,'d',log=TRUE)
	}
	prior.yes <- 'yes'
	if(is.null(prior))
	{
		prior <- function(x) rep(1,length(x))
		prior.yes <- 'no'
	}
	if(!is.null(trans.list))
	{
		start.trans <- sapply(as.list(1:k),function(h) uniroot(function(g)
			trans.list[[h]](g)-model$par.hat[h],uniroot.interval)$root)
		fit <- optim(start.trans,function(g) -sum(log.like(g)),hessian=TRUE)
	}
	else fit <- model$fit
	if(is.null(start)) M <- fit$par
	else M <- start
	if(is.null(variance)) V <- solve(fit$hessian)*tun
	else V <- variance*tun
	sims <- array(NaN,c(iter,k))
	sims[1,] <- M
	rate <- rep(0,iter)
	t1 <- Sys.time()
	for(i in 2:iter)
	{
		phi <- rmvnorm(1,sims[i-1,],V)
		suppressWarnings(
		ln.ratio <- sum(log.like(phi)-log.like(sims[i-1,])) + sum(log(prior(phi))-log(prior(sims[i-1,])))
		)
		if(!is.na(ln.ratio) & runif(1) < exp(ln.ratio))
		{
			sims[i,] <- phi
			rate[i] <- 1
		}
		else sims[i,] <- sims[i-1,]
	}
	t2 <- difftime(Sys.time(),t1,units='mins')
	sims.out <- sapply(as.list(1:k),function(g) trans.list[[g]](sims[,g]))
	colnames(sims.out) <- names(formals(paste('r',model$dist,sep='')))[2:(1+model$k)]
	if(burn!=0) sims.burnt <- sims.out[-c(1:burn),]
	else sims.burnt <- sims.out
	out <- list(rate=mean(rate),total.time=t2,sims.all=sims.out,sims=sims.burnt,
		input=model,iter=iter,prior=prior.yes,burn=burn)
	class(out) <- 'metropolis'
	return(out)
}

