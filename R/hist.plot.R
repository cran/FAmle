hist.plot <-
function(x,...)
{
	HIST <- hist(x$x.info[,'x'],plot=FALSE)
	find.max <- NULL
	try(find.max <- optimize(function(k) distr(k,x$dist,x$par.hat,'d'),interval=c(0,max(HIST$mids)),maximum=TRUE)$objective,silent=TRUE)
	if(find.max > max(HIST$density) & !is.null(find.max))
		hist(x$x.info[,'x'],freq=FALSE,col='gray',ylim=c(0,find.max*1.01),main='',xlab=x$data.name,border='white',cex.axis=.8)
	else hist(x$x.info[,'x'],freq=FALSE,col='gray',main='',xlab=x$data.name,border='white',cex.axis=.8)
	rug(x$x.info[,'x'],col='red',lwd=2)
	fun <- function(a) distr(a,x$dist,x$par.hat)
	curve(fun,add=TRUE,col='red')
}

