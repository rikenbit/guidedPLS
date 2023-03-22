.Easy <- function(){
    # X1 (100 x 300)
    X1 <- matrix(rpois(100 * 300, lambda=1),
            nrow = 100, ncol = 300)
    X1[1:50, 1:50] <- rpois(50*50, lambda=100)
    X1[51:70, 51:100] <- rpois(20*50, lambda=200)
    X1[71:100, 151:170] <- rpois(30*20, lambda=150)
    # Y1 (100 x 50)
    Y1 <- matrix(rpois(100 * 50, lambda=1), nrow = 100, ncol = 50)
    Y1[1:50, 1:10] <- rpois(50*10, lambda=100)
    Y1[51:70, 11:25] <- rpois(20*15, lambda=200)
    Y1[71:100, 26:50] <- rpois(30*25, lambda=150)
    # Y1_dummy (100 x 5)
    Y1_dummy <- matrix(0, nrow = 100, ncol = 3)
    Y1_dummy[1:50, 1] <- 1
    Y1_dummy[51:70, 2] <- 1
    Y1_dummy[71:100, 3] <- 1
    # X2 (200 x 150)
    X2 <- matrix(rpois(200 * 150, lambda=1),
            nrow = 200, ncol = 150)
    X2[1:120, 1:20] <- rpois(120*20, lambda=100)
    X2[121:150, 21:50] <- rpois(30*30, lambda=200)
    X2[151:200, 51:70] <- rpois(50*20, lambda=150)
    # Y2 (200 x 50)
    Y2 <- matrix(rpois(200 * 50, lambda=1), nrow = 200, ncol = 50)
    Y2[1:120, 1:10] <- rpois(120*10, lambda=100)
    Y2[121:150, 11:25] <- rpois(30*15, lambda=200)
    Y2[151:200, 26:50] <- rpois(50*25, lambda=150)
    # Y2_dummy (200 x 5)
    Y2_dummy <- matrix(0, nrow = 200, ncol = 3)
    Y2_dummy[1:120, 1] <- 1
    Y2_dummy[121:150, 2] <- 1
    Y2_dummy[151:200, 3] <- 1
    # Color Vector
    col1 <- c(rep("#66C2A5", length=50),
        rep("#FC8D62", length=20),
        rep("#8DA0CB", length=30))
    col2 <- c(rep("#66C2A5", length=120),
        rep("#FC8D62", length=30),
        rep("#8DA0CB", length=50))
    # Output
    list(X1=X1, X2=X2,
        Y1=Y1, Y1_dummy=Y1_dummy,
        Y2=Y2, Y2_dummy=Y2_dummy,
        col1=col1, col2=col2)
}

.Hard <- function(){
    # X1 (100 x 300)
    X1 <- matrix(rpois(100 * 300, lambda=98),
            nrow = 100, ncol = 300)
    X1[1:50, 1:50] <- rpois(50*50, lambda=100)
    X1[51:70, 51:100] <- rpois(20*50, lambda=102)
    X1[71:100, 151:170] <- rpois(30*20, lambda=105)
    # Y1 (100 x 50)
    Y1 <- matrix(rpois(100 * 50, lambda=1), nrow = 100, ncol = 50)
    Y1[1:50, 1:10] <- rpois(50*10, lambda=100)
    Y1[51:70, 11:25] <- rpois(20*15, lambda=200)
    Y1[71:100, 26:50] <- rpois(30*25, lambda=150)
    # Y1_dummy (100 x 5)
    Y1_dummy <- matrix(0, nrow = 100, ncol = 3)
    Y1_dummy[1:50, 1] <- 1
    Y1_dummy[51:70, 2] <- 1
    Y1_dummy[71:100, 3] <- 1
    # X2 (200 x 150)
    X2 <- matrix(rpois(200 * 150, lambda=98),
            nrow = 200, ncol = 150)
    X2[1:120, 1:20] <- rpois(120*20, lambda=100)
    X2[121:150, 21:50] <- rpois(30*30, lambda=102)
    X2[151:200, 51:70] <- rpois(50*20, lambda=105)
    # Y2 (200 x 50)
    Y2 <- matrix(rpois(200 * 50, lambda=1), nrow = 200, ncol = 50)
    Y2[1:120, 1:10] <- rpois(120*10, lambda=100)
    Y2[121:150, 11:25] <- rpois(30*15, lambda=200)
    Y2[151:200, 26:50] <- rpois(50*25, lambda=150)
    # Y2_dummy (200 x 5)
    Y2_dummy <- matrix(0, nrow = 200, ncol = 3)
    Y2_dummy[1:120, 1] <- 1
    Y2_dummy[121:150, 2] <- 1
    Y2_dummy[151:200, 3] <- 1
    # Color Vector
    col1 <- c(rep("#66C2A5", length=50),
        rep("#FC8D62", length=20),
        rep("#8DA0CB", length=30))
    col2 <- c(rep("#66C2A5", length=120),
        rep("#FC8D62", length=30),
        rep("#8DA0CB", length=50))
    # Output
    list(X1=X1, X2=X2,
        Y1=Y1, Y1_dummy=Y1_dummy,
        Y2=Y2, Y2_dummy=Y2_dummy,
        col1=col1, col2=col2)
}

.flist <- list(
    Easy = .Easy,
    Hard = .Hard
)

toyModel <- function(model="Easy", seeds=123){
    set.seed(seeds)
    out <- .flist[[model]]()
    set.seed(NULL)
    out
}
