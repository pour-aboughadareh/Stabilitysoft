#########################################################################################################
##
##   STABILITYSOFT: a new online program to calculate parametric and non-parametric stability statistics
##
##   Authors:
##   Alireza Pour-Aboughadareh (a.poraboghadareh@gmail.com)
##   Mohsen Yousefian (contact@mohsenyousefian.com)
##
#########################################################################################################
##
##   Usage:
##
##
##   1. load table data to a dataframe variable named "df"
##
##   If you have average matrix:
##
##   2. results <- Calculate(df)
##   3. print(results$statistics)
##   4. print(results$ranks)
##
##   Or, if you have raw data:
##
##   2. results <- CalculateRaw(df)
##   3. print(results$enviroments)
##   4. print(results$average_matrix)
##   5. print(results$statistics)
##   6. print(results$ranks)
##
#########################################################################################################
##
##  np1, np2, np3, np4: Thennarasu’s non-parametric statistics
##
##  ShuklaEquivalance: Shukla’s stability variance
##
##  s1, s2, s3, s6, z1, z2: Huhn’s and Nassar and Huhn’s non-parametric statistics
##
##  BI: Regression coefficient
##
##  SDI: Deviation from regression
##
##  CVI: Coefficient of variance
##
##  Kang: Kang’s rank-sum
##
##  P: Plaisted’s GE variance component
##
##  PaP: Plaisted and Peterson’s mean variance component
##
#########################################################################################################

(init <- function()
{
    np1 <- function(a, b, mat, rankMid)
    {
        mat <- abs(mat - rankMid)
        sum <- apply(mat, 1, sum)
        output <- sum / b
        return(output)
    }
    np2 <- function(a, b, mat, rankMid, rankMidO)
    {
        mat <- abs((mat - rankMid) / rankMidO)
        sum <- apply(mat, 1, sum)
        output <- sum / b
        return(output)
    }
    np3 <- function(a, b, mat, rankAvg, rankAvgO)
    {
        mat <- (mat - rankAvg) ^ 2 / b
        sum <- apply(mat, 1, sum)
        output <- sqrt(sum) / rankAvgO
        return(output)
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
        mat <- pmax(mat, 0)
        sum <- apply((mat - matAvg) ^ 2, 1, sum)
        output <- (((sum / (b - 1)) - eS2) ^ 2) / varS2
        return(output)
    }
    s1 <- function(a, b, mat)
    {
        mat <- pmax(mat, 0)
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
        mat <- pmax(mat, 0)
        sum <- apply((mat - matAvg) ^ 2, 1, sum)
        output <- sum / (b - 1)
        return(output)
    }
    s3 <- function(a, b, mat, matAvg)
    {
        mat <- pmax(mat, 0)
        sum <- apply((mat - matAvg) ^ 2, 1, sum)
        output <- sum / matAvg
        return(output)
    }
    s6 <- function(a, b, mat, matAvg)
    {
        mat <- pmax(mat, 0)
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
        return(output)
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
    CVI <- function(a, b, pureMatSDI, pureMatCVI, pureMatAvgCol, pureMatAvg, pureMatTotalAvg, pureMatBI, pureMatBPJ, pureMat)
    {
        sumSDI <- apply(pureMatSDI, 1, sum)
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

        # Sum Rank
        SR <- data.frame(apply(ranks, 1, sum))
        colnames(SR) <- "SR"

        # Average Rank
        AR <- SR / length(ranks)
        colnames(AR) <- "AR"

        # Standard deviation
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


        eS1 <- ((a * a) - 1) / (3 * a)
        varS1 <- (((a * a) - 1) * (((a * a) - 4) * (b + 3) + 30)) / (45 * a * a * b * (b - 1))
        eS2 <- ((a * a) - 1) / 12
        varS2 <- (((a * a) - 1) * (2 * ((a * a) - 4) * (b - 1) + (5 * ((a * a) - 1)))) / (360 * b * (b - 1))

        pureMat <- data.matrix(table)

        pureMatAvg <- rowMeans(pureMat)
        pureMatAvgCol <- t(matrix(colMeans(pureMat), nrow = b, ncol = a))
        pureMatTotalAvg <- ave(pureMat)[1, 1]

        pureMatNP <- pureMat - pureMatAvg + pureMatTotalAvg

        pureMatShu <- (pureMat - pureMatAvg - pureMatAvgCol + pureMatTotalAvg) ^ 2

        pureMatBI <- (pureMat - pureMatAvg) * (pureMatAvgCol - pureMatTotalAvg)

        pureMatBPJ <- (pureMat - pureMatAvg - pureMatAvgCol + pureMatTotalAvg) * (pureMatAvgCol - pureMatTotalAvg)

        pureMatSDI <- (pureMat - pureMatAvg) ^ 2

        pureMatCVI <- (pureMat - pureMatTotalAvg) ^ 2


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

        CVI <- CVI(a, b, pureMatSDI, pureMatCVI, pureMatAvgCol, pureMatAvg, pureMatTotalAvg, pureMatBI, pureMatBPJ, pureMat)

        Kang <- Kang(a, b, pureMatShu, pureMatAvg)

        P <- P(a, b, pureMatShu)

        PaP <- PaP(a, b, pureMatShu)

        stats_df <- data.frame(table_original[, 1], Y, s1, z1, s2, z2, s3, s6, np1, np2, np3, np4, ShuklaEquivalance[, 1], ShuklaEquivalance[, 2], SDI, BI, CVI, P, PaP, Kang)
        colnames(stats_df) <- c("Genotype", "Y", "S1", "Z1", "S2", "Z2", "S3", "S6", "NP1", "NP2", "NP3", "NP4", "Wricke’s ecovalence", "Shukla’s stability variance", "Deviation from regression", "Regression coefficient", "Coefficient of variance", "GE variance component", "Mean variance component", "Kang’s rank-sum")

        ranks_df <- getranks_df(a, b, stats_df)

        output <- list(statistics = stats_df, ranks = ranks_df, correlation_matrix = cor(data.matrix(stats_df[c(-1, -4, -6)][, 1:(length(stats_df[c(-1)]) - 2)])))
        return(output)
    }

    CalculateRaw <<- function(original_table)
    {
        objkeys <- colnames(original_table)
        specs_keys <- c()
        arr_not_specs <- c("replication", "yield", "genotype")
        for (i in 1:length(objkeys))
        {
            key <- objkeys[i]
            if (!(tolower(key) %in% arr_not_specs))
                specs_keys <- append(specs_keys, key)
        }
        genotypekey <- if ("genotype" %in% objkeys) "genotype" else (if ("Genotype" %in% objkeys) "Genotype" else objkeys[3])
        yieldkey <- if ("yield" %in% objkeys) "yield" else (if ("Yield" %in% objkeys) "Yield" else objkeys[3])
        if (length(specs_keys) < 1)
            stop("At least one property for enviroments (location, year, ...) should be provided in input data")
        specs_columns <- list()
        for (i in 1:length(specs_keys))
            specs_columns <- append(specs_columns, original_table[specs_keys[i]])
        sorted_df <- original_table[do.call(order, specs_columns),]

        genotypes_list <- list()
        enviroments_array <- list()
        for (i in 1:nrow(sorted_df))
        {
            enviroment_case <- list(as.numeric(sorted_df[i, specs_keys]))
            if (length(intersect(enviroment_case, enviroments_array)) < 1)
                enviroments_array <- append(enviroments_array, enviroment_case)

            genotype <- sorted_df[i, genotypekey]
            if (length(intersect(genotype, genotypes_list)) < 1)
                genotypes_list <- append(genotypes_list, genotype)
        }
        enviroments <- data.frame(c(1:length(specs_keys)))
        for (i in 1:length(enviroments_array))
        {
            for (j in 1:length(enviroments_array[[i]]))
            {
                enviroments[i, j] <- enviroments_array[[i]][j]
            }
        }
        colnames(enviroments) <- specs_keys
        rownames(enviroments) <- lapply(list(1:length(enviroments[, 1])), function(x) { return(paste0('E', as.character(x))) })[[1]]
        sum_matrix <- data.frame(matrix(0, ncol = length(enviroments[, 1]), nrow = length(genotypes_list)))
        colnames(sum_matrix) <- rownames(enviroments)
        counter_matrix <- data.frame(sum_matrix)
        specs_keys_count <- length(specs_keys)
        for (i in 1:nrow(original_table))
        {
            row_data <- original_table[i,]
            target_specs <- list()
            for (key in specs_keys)
                target_specs <- append(target_specs, as.numeric(row_data[key]))
            enviroment_index <- 0
            for (i in 1:nrow(enviroments))
            {
                enviroment <- enviroments[i,]
                if (all(as.numeric(enviroment) == as.numeric(target_specs)))
                {
                    enviroment_index <- i
                    break
                }
            }
            genotype_index <- which(as.character(row_data[genotypekey]) == genotypes_list)

            counter_matrix[genotype_index, enviroment_index] <- counter_matrix[genotype_index, enviroment_index] + 1
            sum_matrix[genotype_index, enviroment_index] <- sum_matrix[genotype_index, enviroment_index] + as.numeric(row_data[yieldkey])
        }
        average_matrix <- sum_matrix / counter_matrix
        average_matrix <- cbind(Genotype = as.array(genotypes_list), average_matrix)
        calc_results <- Calculate(average_matrix)
        output <- list(average_matrix = average_matrix, enviroments = enviroments, statistics = calc_results$statistics, ranks = calc_results$ranks, correlation_matrix = calc_results$correlation_matrix)
        return(output)
    }
})()
