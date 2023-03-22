.minDim <- function(...){
    min(unlist(lapply(list(...), dim)))
}

softThr <- function(y, lambda=1){
	if(max(abs(y)) < lambda){
		msg <- paste0(
			"All the elements are converted to 0.\n",
			"Please specify smaller lambda value.")
		stop(msg)
	}
	sign(y) * sapply(abs(y)-lambda, function(x){max(x, 0)})
}

dummyMatrix <- function(y, center=TRUE){
    stopifnot(is.numeric(y))
    group <- unique(y)
    l <- length(group)
	out <- matrix(0, nrow=length(y), ncol=l)
	colnames(out) <- group
	for(i in seq_len(l)){        
		out[which(y == group[i]), i] <- 1
	}
    if(center){
    	out <- scale(out, center=TRUE, scale=FALSE)
    }
    out
}
