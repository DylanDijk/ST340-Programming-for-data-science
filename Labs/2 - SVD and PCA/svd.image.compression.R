# helper function to view an image
viewImage <- function(x) {
  plot(1:2,axes=FALSE,xlab="",ylab="",type='n')
  rasterImage(x,xleft=1,ybottom=1,xright=2,ytop=2)
}

image.compression <- function() {
  # choose a picture
  pic <- as.integer(readline(prompt="Choose a picture (1-5) from Gauss, Cox, von Neumann, Nightingale or checkerboard: "))
  res <- image.compress.param(pic)
  
  # plot side-by-side the original image and the singular values
  par(mfrow=c(1,2))
  viewImage(res$mtx)
  plot(res$svd$d)
  abline(h=0,lty=2)
  
  # until you choose to quit...
  k <- Inf
  while (TRUE) {
    # choose the rank of the approximation
    k <- as.integer(readline(prompt=paste("choose a number of k between 1 and ",res$p," (anything else to exit): ",sep="")))
    compressedImage <- compute.compression(k, res$p, res$mtx, res$svd)

    # view the original and compressed image side-by-side
    if (!is.null(compressedImage)) {
      viewImage(res$mtx)
      viewImage(compressedImage)
    } else {
      break;
    }
  }
}

image.compress.param <- function(pic) {
  img <- images[[pic]] 
  
  # find the size of the image
  dims <- dim(img); m <- dims[1]; n <- dims[2]
  
  if (length(dims) > 2) {
    # convert the image into greyscale
    mtx <- matrix(0,m,n)
    for (i in 1:m) {
      for (j in 1:n) {
        mtx[i,j] <- sum(img[i,j,])/3
      }
    }
  } else {
    mtx <- img
  }
  
  p <- min(m,n)
  
  # perform the decomposition
  decomposition <- svd(mtx)
  
  return(list(img=img, mtx=mtx, p=p, svd=decomposition))
}

compute.compression <- function(k, p, mtx, decomposition) {
  if (1 <= k && k <= p) {
    
    # compute the k-rank approximation
    if (k == 1) {
      approximation <- decomposition$d[1]*decomposition$u[,1]%*%t(decomposition$v[,1])
    } else {
      approximation <- decomposition$u[,1:k]%*%diag(decomposition$d[1:k])%*%t(decomposition$v[,1:k])  
    }
    
    approximation.error <- norm(mtx-approximation,type="F")
    approximation.error.theory <- sqrt(sum(decomposition$d[(k+1):p]^2))
    
    print(paste("approximation.error = ",approximation.error,sep=""))
    print(paste("approximation.error.theory = ",approximation.error.theory,sep=""))
    
    # rescale the approximation so the values of the image matrix are in [0,1]
    maxval <- max(approximation); minval <- min(approximation)
    compressedImage <- (approximation - minval)/(maxval - minval)
    return(compressedImage)
  } else {
    return(NULL)
  }
}
