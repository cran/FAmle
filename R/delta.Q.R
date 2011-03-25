delta.Q <-
function(p,model,ln=FALSE)
{
	k <- model$k
	cov.hat <- model$cov.hat
	prime <- sapply(as.list(1:k),function(h) Diff.3(h,model,p,ln))
	VAR <- diag(cov.hat)
	COV <- cov.hat[lower.tri(cov.hat)]
	Q.var.hat <- sum(2*COV*sapply(as.list(as.data.frame(combn(1:k,2))),function(g) prod(prime[g]))) +
		sum(VAR*prime^2)
	if(!ln) return(c(mu.hat=distr(p,model=model,type='q'),sd.hat=sqrt(Q.var.hat)))
	else return(c(mu.hat=log(distr(p,model=model,type='q')),sd.hat=sqrt(Q.var.hat)))
}

