## 1. Lalonde data from
## http://users.nber.org/~rdehejia/nswdata2.html <2018-06-25 Mon>

dtr <- read.csv("./nswre74_treated.txt", sep="", header=FALSE)
dct <- read.csv("./psid_controls.txt", sep="", header=FALSE)

cleandata <- function(d) {
    names(d) <- c("treated", "age", "education", "black", "hispanic",
                  "married", "nodegree", "re74", "re75", "re78")

    ## Convert 1-0 variables to logical
    d[, c(1, 4, 5, 6, 7)] <- matrix(as.logical(unname(unlist(
        d[, c(1, 4, 5, 6, 7)]))), nrow=nrow(d))
    ## add unemployment, defined as zero earnings
    d$ue74 <- d$re74==0
    d$ue75 <- d$re75==0
    ## exclude nodegree, which is not in Abadie and Imbens, 2011, JBES
    d$nodegree <- NULL
    d[, c(2, 3)] <- matrix(as.integer(unname(unlist(
        d[, c(2, 3)]))), nrow=nrow(d))

    ## earnings in 1000s of dollars
    d[, 7:9] <- d[, 7:9]/1000
    d <- cbind(d[, -9], d[, 9, drop=FALSE])
}
NSW <- cleandata(rbind(dct, dtr))
## Check we match Abadie and Imbens, 2011, JBES, Table 1
stargazer::stargazer(NSW[NSW$treated==TRUE, ], type="text")
stargazer::stargazer(NSW[NSW$treated==FALSE, ], type="text")
usethis::use_data(NSW, overwrite=TRUE, internal=FALSE)

## 2. Experimental data
## http://users.nber.org/~rdehejia/nswdata2.html <2018-11-30 Fri>
dtr <- read.csv("./nswre74_treated.txt", sep="", header=FALSE)
dct <- read.csv("./nswre74_control.txt", sep="", header=FALSE)
NSWexper <- cleandata(rbind(dct, dtr))
## Check we match Abadie and Imbens, 2011, JBES, Table 1
usethis::use_data(NSWexper, overwrite=TRUE, internal=FALSE)
