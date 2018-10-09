##########################
##
##   STABILITYSOFT: a new online program to calculate parametric and non-parametric stability statistics
##   
##   Authors:
##   Alireza Pour-Aboughadareh (a.poraboghadareh@gmail.com)
##   Mohsen Yousefian (contact@mohsenyousefian.com)
##
##########################
##
##   Usage:
##
##
##   1. load table data to a dataframe variable named "df"
##   2. results <- Calculate(df)
##   3. print(results$statistics)
##   4. print(results$ranks)
##
##########################

(init <- function()
{
    BPJ <- function(a, b, pureMatBPJ, pureMatAvgCol, pureMatTotalAvg)
    {
        powAvg <- apply((pureMatAvgCol - pureMatTotalAvg) ^ 2, 1, sum)
        powAvgTotal <- sum(powAvg)
        sum <- apply(pureMatBPJ, 1, sum)

        output <- sum / powAvgTotal
        return(output)
    }
    np1 <- function(a, b, mat, rankMid)
    {
        mat = abs(mat - rankMid)
        sum <- apply(mat, 1, sum)
        output <- sum / b
        return(output)
    }
    np2 <- function(a, b, mat, rankMid, rankMidO)
    {
        mat = abs((mat - rankMid) / rankMidO)
        sum <- apply(mat, 1, sum)
        output <- sum / b
        return(output)
    }
    np3 <- function(a, b, mat, rankAvg, rankAvgO)
    {
        mat <- (mat - rankAvg) ^ 2 / b
        sum <- apply(mat, 1, sum)
        output <- sqrt(sum) / rankAvgO
        return(output);
    }
    np4 <- function(a, b, mat, rankAvgO)
    {
        mat <- pmax(mat, 0)
        sum <- vector(length = a)
        for (i in 1:a)
            for (j in 1:b)
                for (k in j:b)
                    sum[i] <- sum[i] + abs(mat[i, j] - mat[i, k])
        output <- (2 / (b * (b - 1))) * (sum / rankAvgO)
        return(output)
    }
    ShuklaEquivalance <- function(a, b, pureMatShu)
    {
        equivalance <- apply(pureMatShu, 1, sum)
        equivalanceTotal <- sum(pureMatShu)
        output <- matrix(nrow = a, ncol = 2)
        shukla <- ((a * equivalance) / ((b - 1) * (a - 2))) - ((equivalanceTotal) / ((a - 1) * (a - 2) * (b - 1)))
        output[, 1] <- equivalance
        output[, 2] <- shukla
        return(output)
    }
    z1 <- function(a, b, mat, eS1, varS1)
    {
        mat <- pmax(mat, 0)
        sum <- vector(length = a)
        for (i in 1:a)
            for (j in 1:b)
                for (k in j:b)
                    sum[i] <- sum[i] + abs(floor(mat[i, j] - mat[i, k]))
        output <- (((2 / (b * (b - 1))) * sum) - eS1) ^ 2 / varS1
        return(output)
    }
    z2 <- function(a, b, mat, matAvg, eS2, varS2)
    {
        mat = pmax(mat, 0)
        sum <- apply((mat - matAvg) ^ 2, 1, sum)
        output <- (((sum / (b - 1)) - eS2) ^ 2) / varS2
        return(output)
    }
    s1 <- function(a, b, mat)
    {
        mat = pmax(mat, 0)
        sum <- vector(length = a)
        for (i in 1:a)
            for (j in 1:b)
                for (k in j:b)
                    sum[i] <- sum[i] + abs(floor(mat[i, j] - mat[i, k]))
        output <- (2 / (b * (b - 1))) * sum
        return(output)
    }
    s2 <- function(a, b, mat, matAvg)
    {
        mat = pmax(mat, 0)
        sum <- apply((mat - matAvg) ^ 2, 1, sum)
        output <- sum / (b - 1)
        return(output)
    }
    s3 <- function(a, b, mat, matAvg)
    {
        mat = pmax(mat, 0)
        sum <- apply((mat - matAvg) ^ 2, 1, sum)
        output <- sum / matAvg
        return(output)
    }
    s6 <- function(a, b, mat, matAvg)
    {
        mat = pmax(mat, 0)
        sum <- apply(abs(mat - matAvg), 1, sum)
        output <- sum / matAvg
        return(output)
    }
    BI <- function(a, b, pureMatBI, pureMatAvgCol, pureMatTotalAvg)
    {
        powAvg <- ((pureMatAvgCol - pureMatTotalAvg) ^ 2)[1,]
        powAvgTotal <- sum(powAvg)
        sum <- apply(pureMatBI, 1, sum)
        output <- sum / powAvgTotal
        return(output);
    }
    SDI <- function(a, b, pureMatSDI, pureMatAvgCol, pureMatTotalAvg, pureMatBI)
    {
        powAvg <- ((pureMatAvgCol - pureMatTotalAvg) ^ 2)[1,]
        powAvgTotal <- sum(powAvg)
        sum <- apply(pureMatSDI, 1, sum)
        bi <- BI(a, b, pureMatBI, pureMatAvgCol, pureMatTotalAvg)
        output <- ((sum) - ((bi * bi) * powAvgTotal)) / 7
        return(output)
    }
    CVR <- function(a, b, pureMatSDI, pureMatCVR, pureMatAvgCol, pureMatAvg, pureMatTotalAvg, pureMatBI, pureMatBPJ, pureMat)
    {
        # calculate bi
        bFWT <- BI(a, b, pureMatBI, pureMatAvgCol, pureMatTotalAvg)
        bPJT <- BPJ(a, b, pureMatBPJ, pureMatAvgCol, pureMatTotalAvg)

        t <- (pureMatAvgCol - pureMatTotalAvg) ^ 2
        d <- (pureMat - pureMatAvg) ^ 2

        # calculate sumCVR & sumSDI & sumXbi
        sumSDI <- apply(pureMatSDI, 1, sum)
        sumCVR <- apply(pureMatCVR, 1, sum)
        sumXbFW <- bFWT ^ 2 * sumCVR
        sumXbPJ <- bPJT ^ 2 * sumCVR
        RFW <- bFWT * t / d
        RPJ <- bPJT * t / d
        CV <- ((sqrt(sumSDI / (b - 1))) / pureMatAvg) * 100
        output <- CV
        return(output)
    }
    Kang <- function(a, b, pureMatShu, PureMatAvg)
    {
        equivalance <- apply(pureMatShu, 1, sum)
        equivalanceTotal <- sum(equivalance)

        shukla <- ((a * equivalance) / ((b - 1) * (a - 2))) - ((equivalanceTotal) / ((a - 1) * (a - 2) * (b - 1)))

        rankAvgR <- vector()
        tmp <- sort(PureMatAvg)
        for (i in 1:a)
            for (j in 1:a)
                if (tmp[j] == PureMatAvg[i])
                {
                    rankAvgR[i] <- a - j
                    break
                }

        rankShu <- vector()
        tmp <- sort(shukla)
        for (i in 1:a)
            for (j in 1:a)
                if (tmp[j] == shukla[i])
                {
                    rankShu[i] <- j + 1
                    break
                }

        kang <- rankAvgR + rankShu

        output <- kang
        return(output)
    }
    P <- function(a, b, pureMatShu)
    {
        equivalance <- apply(pureMatShu, 1, sum)
        equivalanceTotal <- sum(equivalance)

        shukla <- (((-1 * a) * equivalance) / ((b - 1) * (a - 2) * (a - 1))) + ((equivalanceTotal) / ((a - 2) * (b - 1)))

        output <- shukla
        return(output)
    }
    PaP <- function(a, b, pureMatShu)
    {
        equivalance <- apply(pureMatShu, 1, sum)
        equivalanceTotal <- sum(equivalance)

        shukla <- ((a * equivalance) / (2 * (b - 1) * (a - 1))) + ((equivalanceTotal) / (2 * (a - 2) * (b - 1)))

        output <- shukla
        return(output)
    }
    calRankMid <- function(a, b, mat)
    {
        temp <- apply(mat, 1, sort)
        bd2 <- ceiling(b / 2)
        if (b %% 2 == 0)
            rankMid <- ave(temp[bd2:bd2 + 1,])
        else
            rankMid <- temp[bd2 + 1,]
        return(rankMid)
    }
    getranks_df <- function(a, b, df_orig)
    {
        df <- t(df_orig[, c(3, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 20, 19, 18)])
        ranks <- apply(df, 1, rank, ties.method = "min")

        Y <- data.frame(rank(-df_orig[2]))
        colnames(Y) <- c("Y")



        ranks <- cbind(Y, ranks)

        SR <- data.frame(apply(ranks, 1, sum))
        colnames(SR) <- "SR"

        AR <- SR / length(ranks)
        colnames(AR) <- "AR"

        STD <- data.frame(apply(ranks, 1, sd))
        colnames(STD) <- "Std."

        ranks <- cbind(df_orig[1], ranks, SR, AR, STD)


        return(ranks)
    }
    Calculate <<- function(table_original)
    {
        table <- cbind(table_original)[, -1]
        a <- nrow(table)
        b <- ncol(table)
        K <- a
        N <- b
        eS1 <- ((K * K) - 1) / (3 * K)
        varS1 <- (((K * K) - 1) * (((K * K) - 4) * (N + 3) + 30)) / (45 * K * K * N * (N - 1))
        eS2 <- ((K * K) - 1) / 12
        varS2 <- (((K * K) - 1) * (2 * ((K * K) - 4) * (N - 1) + (5 * ((K * K) - 1)))) / (360 * N * (N - 1))

        pureMat <- data.matrix(table)

        pureMatAvg <- rowMeans(pureMat)
        pureMatAvgCol <- t(matrix(colMeans(pureMat), nrow = b, ncol = a))
        pureMatTotalAvg <- ave(pureMat)[1, 1]

        pureMatNP <- pureMat - pureMatAvg + pureMatTotalAvg

        pureMatShu <- (pureMat - pureMatAvg - pureMatAvgCol + pureMatTotalAvg) ^ 2

        pureMatBI <- (pureMat - pureMatAvg) * (pureMatAvgCol - pureMatTotalAvg)

        pureMatBPJ <- (pureMat - pureMatAvg - pureMatAvgCol + pureMatTotalAvg) * (pureMatAvgCol - pureMatTotalAvg)
        BJP <- BPJ(a, b, pureMatBPJ, pureMatAvgCol, pureMatTotalAvg)

        pureMatSDI <- (pureMat - pureMatAvg) ^ 2

        pureMatCVR <- (pureMat - pureMatTotalAvg) ^ 2


        mat <- apply(pureMat, 2, rank, ties.method = "min")
        matNP <- apply(pureMatNP, 2, rank, ties.method = "min")

        matAvg <- apply(mat, 1, ave)[1,]
        matAvgNP <- apply(matNP, 1, ave)[1,]

        rankMid <- calRankMid(a, b, matNP)
        rankMidO <- calRankMid(a, b, mat)

        rankAvg <- apply(matNP, 1, ave)[1,]
        rankAvgO <- apply(mat, 1, ave)[1,]

        Y <- pureMatAvg

        np1 <- np1(a, b, matNP, rankMid)
        np2 <- np2(a, b, mat, rankMid, rankMidO)
        np3 <- np3(a, b, matNP, rankAvg, rankAvgO)
        np4 <- np4(a, b, mat, rankAvgO)

        z1 <- z1(a, b, mat, eS1, varS1)
        z2 <- z2(a, b, mat, matAvg, eS2, varS2)

        s1 <- s1(a, b, mat)
        s2 <- s2(a, b, mat, matAvg)
        s3 <- s3(a, b, mat, matAvg)
        s6 <- s6(a, b, mat, matAvg)

        ShuklaEquivalance <- ShuklaEquivalance(a, b, pureMatShu) # wri Shu

        SDI <- SDI(a, b, pureMatSDI, pureMatAvgCol, pureMatTotalAvg, pureMatBI)

        BI <- BI(a, b, pureMatBI, pureMatAvgCol, pureMatTotalAvg)

        CVR <- CVR(a, b, pureMatSDI, pureMatCVR, pureMatAvgCol, pureMatAvg, pureMatTotalAvg, pureMatBI, pureMatBPJ, pureMat)

        Kang <- Kang(a, b, pureMatShu, pureMatAvg)

        P <- P(a, b, pureMatShu)

        PaP <- PaP(a, b, pureMatShu)

        stats_df <- data.frame(table_original[, 1], Y, s1, z1, s2, z2, s3, s6, np1, np2, np3, np4, ShuklaEquivalance[, 1], ShuklaEquivalance[, 2], SDI, BI, CVR, P, PaP, Kang)
        colnames(stats_df) <- c("Genotype", "Y", "S1", "Z1", "S2", "Z2", "S3", "S6", "NP1", "NP2", "NP3", "NP4", "Wricke’s ecovalence", "Shukla’s stability variance", "Deviation from regression", "Regression coefficient", "Coefficient of variance", "GE variance component", "Mean variance component", "Kang’s rank-sum")

        ranks_df <- getranks_df(a, b, stats_df)

        output <- list(statistics = stats_df, ranks = ranks_df, correlation_matrix = cor(data.matrix(stats_df[-1][, 1:length(stats_df[-1])])))
        return(output)
    }
})()
