

# convert a mothur otu table into a usefule dataframe to work with in R

convert.otu.table <- function (csv = data.frame()) {
    otus <- dim(csv)[1]
    csv$Isolate <- as.character(csv$Isolate)
    pattern <- "(?<=^)Rh-[0-9ab]{1,4}(-((rep[23])|([014]{1,2})))?(?=,)"
    for (i in 1:otus) {
        m <- regexpr(pattern, csv[i,2], perl=T)
        while (m > 0) {
            csv <- rbind(csv, data.frame(otu = csv[i,1], Isolate = regmatches(csv[i,2], m)))
            csv[i, 2] <- sub("(?<=^)Rh-[0-9ab]{1,4}(-((rep[23])|([014]{1,2})))?,", "", csv[i, 2], perl=T)
            m <- regexpr(pattern, csv[i,2], perl=T)
        }
    }
    return (csv)
}

convert.otu.table_ <- function (csv = data.frame()) {
    otus <- dim(csv)[1]
    csv$Isolate <- as.character(csv$Isolate)
    pattern <- "(?<=^)Rh_[0-9]{1,3}(_([123]))?(?=,)"
    for (i in 1:otus) {
        m <- regexpr(pattern, csv[i,2], perl=T)
        while (m > 0) {
            csv <- rbind(csv, data.frame(otu = csv[i,1], Isolate = regmatches(csv[i,2], m)))
            csv[i, 2] <- sub("(?<=^)Rh_[0-9]{1,3}(_([123]))?,", "", csv[i, 2], perl=T)
            m <- regexpr(pattern, csv[i,2], perl=T)
        }
    }
    return (csv)
}