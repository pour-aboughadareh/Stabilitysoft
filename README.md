How to use STABILITYSOFT in R
================================================

STABILITYSOFT calculates the several parametric and non-parametric statistics as defined in STABILITYSOFT: a new online program to calculate parametric and non-parametric stability statistics. Applications in Plant Sciences (submitted).

Step 1: Loading the library
--------------------------
To get started, execute the library code ([STABILITYSOFT.R](STABILITYSOFT.R)) in your RStudio console.

```R
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
```
Step 2: Loading the data
--------------------
After execution, function `Calculator` can be used to get the results.
It takes a dataframe as input and returns the result as an object containing dataframes `result$statistics` and `result$ranks`.

First, you need to convert your data to a dataframe (given that the data is already like examples provided [here](Examples.zip)).

### Option 1: Copy data from your favorite app and get it from clipboard

#### 1. Copy data as table in your favorite app
![Copying table data from excel](https://raw.githubusercontent.com/pour-aboughadareh/Stabilitysoft/master/Screenshot%201.jpg)

#### 2. Execute the following command in your RStudio console
```R
df <- read.delim("clipboard")
```

### Option 2: Read from Excel file directly

#### 1. Install `XLConnect` library
```R
install.packages("XLConnect")
```
#### 2. Load `XLConnect` library
```R
library("XLConnect")
```
#### 3. Load data from `xlsx` file
```R
df <- readWorksheetFromFile("<Path to .xlsx file>", sheet=1)
```

Step 3: Getting the results
-------------------
The table data in `df` variable should now look like the following:
```R
> print(df)
   Genotype     E1     E2     E3     E4    E5     E6     E7     E8     E9    E10    E11    E12    E13    E14   E15
1        G1 2835.3 3101.8 3075.5 2483.8 752.3 2960.5 1832.5 1851.3 2667.5 2403.8 2270.7 2006.0 1547.0 2001.6 784.4
2        G2 2296.0 3170.5 3533.4 2830.3 683.0 3644.8 2023.8 2147.5 3392.5 2127.5 2137.5 1892.7 1523.6 1767.8 732.6
3        G3 2252.3 2616.8 2983.4 2405.5 376.3 3033.8 2477.5 2027.5 3086.3 2292.5 2314.6 1687.5 1457.9 1626.6 247.2
4        G4 2326.8 2988.3 3320.8 2297.8 482.0 2899.5 1963.8 2726.3 2526.3 1650.0 1929.3 1686.7 1336.3 1654.8 510.8
5        G5 2403.0 2640.8 3050.9 2431.5 709.0 3055.8 2102.5 2452.5 3060.0 2242.5 2016.0 2079.2 1786.1 2104.7 604.8
6        G6 1736.5 2501.0 2572.1 2375.8 562.3 2741.5 2105.0 2616.3 2718.8 1903.7 2179.4 1598.2  937.5  956.4 407.4
7        G7 2347.5 2852.8 3358.7 2306.5 639.5 3358.3 2183.8 2153.8 2720.0 2241.3 1684.5 1907.3 1495.5 1987.7 663.0
8        G8 2218.0 2621.5 2834.9 2161.0 542.8 2848.8 2172.5 1827.5 2753.8 2198.7 2314.7 1738.7 1753.0 1167.2 501.6
9        G9 2708.0 2969.6 3058.7 2035.3 682.0 3485.0 2955.0 2043.8 3295.0 2418.8 2561.8 1915.5 1767.3 1772.0 915.2
10      G10 1869.5 2587.5 3224.7 1886.8 442.0 2889.3 1742.5 2508.8 2553.8 2395.0 2511.0 1875.2 1640.5 1500.2 528.0
11      G11 2027.5 2795.0 3551.5 2110.0 565.0 3010.5 2142.5 1532.5 2471.3 2182.5 2530.6 2023.1 1640.7 1992.3 643.8
12      G12 2690.5 3041.9 3102.4 1619.8 355.8 2806.3 1927.5 2226.3 2392.5 2137.5 2557.6 1860.9 1532.9 1854.0 458.1
13      G13 2357.8 2847.8 3378.3 2118.0 632.0 2972.8 1693.8 2206.3 2295.0 2350.0 2378.8 1936.4 1532.9 2676.6 755.6
14      G14 2118.5 3099.8 3171.4 2663.8 818.0 3122.0 2057.5 1726.3 2573.8 1985.0 2307.3 2096.7 1157.4 1769.6 557.3
15      G15 2193.8 2829.3 2657.4 2579.8 664.8 3099.5 2113.8 2000.0 2975.0 1514.1 2296.3 1751.9 1322.2 1865.9 621.2
16      G16 2154.5 3048.0 2568.4 2457.3 613.5 3336.8 1991.3 2148.8 2796.3 1946.3 1827.9 1803.2 1584.5 1814.3 481.8
17      G17 2287.0 2911.8 3265.1 2566.5 618.0 3050.3 2213.8 2348.8 2833.8 2200.0 2221.7 1999.5 1284.5 1664.1 713.7
18      G18 2410.3 2553.0 3381.6 2100.3 518.8 2618.0 2040.0 2873.8 3031.3 1795.0 2212.9 1674.7 1884.7 1959.5 511.4
19      G19 2011.0 2956.0 3125.6 2492.5 692.0 3076.8 1630.0 1762.5 2455.0 2060.0 2300.8 1879.2 1584.3 2582.9 663.3
20      G20 2292.0 2621.0 3004.2 2397.0 675.8 3016.3 2191.3 2223.8 2770.0 1737.5 2259.1 1950.8 1589.1 2212.5 700.8
```

#### 1. Use `Calculator` function to get the output
```R
results <- Calculator(df)
```
#### 2. Get the output values available in `results` variable (`results$statistics` and `results$ranks`)
```R
> print(results$statistics)

   Genotype        Y       S1           Z1       S2           Z2       S3        S6      NP1       NP2       NP3       NP4 Wricke’s ecovalence Shukla’s stability variance Deviation from regression Regression coefficient Coefficient of variance GE variance component Mean variance component Kang’s rank-sum
1        G1 2171.600 6.647619 7.502891e-06 32.60000 6.118828e-03 34.57576  5.818182 5.066667 0.3250000 0.4334028 0.5036075                   695993.0                           51760.78                  95686.37              0.9426881                34.27273              63151.13                59194.34                     13
2        G2 2260.233 6.514286 2.437689e-02 31.54286 4.220669e-02 32.95522  5.701493 4.866667 0.4708333 0.4749616 0.4861407                   915588.4                           69188.99                  92365.91              1.1836912                40.64888              62233.86                67449.80                     14
3        G3 2059.047 6.876190 6.771359e-02 34.60000 2.639424e-02 55.04545  8.500000 4.933333 0.6000000 0.6747955 0.7813853                   645775.2                           47775.24                  77638.11              1.1132781                41.90916              63360.90                57306.45                     19
4        G4 2019.967 5.695238 1.206472e+00 27.38095 4.988579e-01 50.00000  8.260870 3.800000 0.7222222 0.6637238 0.7428571                   811231.9                           60906.73                 110962.64              1.0657746                41.48922              62669.76                63526.62                     28
5        G5 2182.620 5.466667 1.853282e+00 22.40952 1.701917e+00 23.29703  4.346535 4.800000 0.2666667 0.3955442 0.4059406                   381286.8                           26784.10                  51588.45              0.9497059                33.65086              64465.69                47363.27                      6
6        G6 1860.793 5.542857 1.622313e+00 26.45714 6.682616e-01 71.23077 11.846154 5.866667 1.5200000 1.2313033 1.0659341                  1477612.4                          113794.07                 210095.13              0.9704829                43.04011              59886.22                88578.52                     40
7        G7 2126.680 4.571429 5.718171e+00 17.26667 3.699783e+00 20.37079  3.741573 3.000000 0.2410256 0.3405254 0.3852327                   516027.4                           37477.80                  72418.25              1.0337834                37.75915              63902.87                52428.71                     11
8        G8 1976.980 5.390476 2.099616e+00 22.49524 1.675109e+00 39.04132  7.355372 4.866667 0.8571429 0.6468493 0.6682409                   685343.2                           50915.56                  94810.62              0.9478676                37.82086              63195.62                58793.96                     26
9        G9 2305.533 6.228571 2.350581e-01 31.42857 4.804695e-02 29.33333  4.666667 6.333333 0.2561404 0.4427969 0.4152381                  1236043.4                           94621.93                 176254.42              1.0168455                35.68733              60895.28                79496.98                     18
10      G10 2010.320 6.952381 1.210141e-01 37.20952 2.270531e-01 63.01613  9.532258 5.666667 0.6750000 0.7750331 0.8410138                   958461.6                           72591.62                 136687.19              0.9856088                39.21739              62054.77                69061.57                     32
11      G11 2081.253 6.990476 1.534266e-01 34.69524 3.024965e-02 44.69939  6.822086 4.600000 0.4333333 0.5203525 0.6432954                   956036.6                           72399.16                 136474.71              0.9905394                38.04627              62064.90                68970.41                     25
12      G12 2037.600 7.200000 4.003618e-01 37.52381 2.645281e-01 60.61538  8.692308 5.533333 0.5866667 0.7432537 0.8307692                  1064631.8                           81017.83                 151396.47              1.0246806                40.28090              61611.29                73052.94                     30
13      G13 2142.140 6.742857 1.141190e-02 33.20952 2.372689e-05 39.62500  5.988636 4.533333 0.3380952 0.4909433 0.5746753                  1264303.9                           96864.82                 173679.95              0.9219704                35.27307              60777.23                80559.41                     25
14      G14 2081.627 7.390476 7.256871e-01 39.17143 5.078022e-01 48.96429  7.071429 5.466667 0.4388889 0.5777184 0.6598639                   709462.2                           52829.77                 101034.86              1.0166797                38.40658              63094.87                59700.70                     21
15      G15 2032.333 5.847619 8.520958e-01 24.54286 1.097976e+00 35.79167  6.500000 4.866667 0.5272727 0.5866309 0.6091270                   695341.4                           51709.07                  98340.03              0.9704518                37.65059              63153.85                59169.84                     24
16      G16 2038.193 5.904762 7.350507e-01 25.06667 9.698444e-01 38.70588  6.735294 4.533333 0.4800000 0.5903455 0.6512605                   635652.9                           46971.89                  90804.27              0.9983007                38.41249              63403.18                56925.91                     19
17      G17 2145.240 4.876190 4.164292e+00 17.92381 3.401810e+00 20.45652  4.054348 3.333333 0.2380952 0.3139028 0.3975155                   234316.3                           15119.78                  32091.03              1.0348424                36.87934              65079.60                41838.07                      6
18      G18 2104.353 8.133333 2.912090e+00 48.66667 3.442092e+00 70.48276  9.517241 6.266667 0.6592593 0.7336323 0.8413793                  1256202.1                           96221.82                 179326.69              0.9892858                38.22442              60811.08                80254.83                     27
19      G19 2084.793 6.552381 1.261236e-02 31.82857 2.926117e-02 42.84615  6.615385 4.733333 0.4111111 0.5276696 0.6300366                  1144622.3                           87366.28                 157867.08              0.9295666                36.24661              61277.16                76060.10                     26
20      G20 2109.413 5.066667 3.317966e+00 18.98095 2.948706e+00 23.04046  4.462428 3.600000 0.2615385 0.3918717 0.4393064                   362776.2                           25315.00                  43392.67              0.9139566                33.43444              64543.01                46667.39                     10


> print(results$ranks)

   Genotype  Y S1 S2 S3 S6 NP1 NP2 NP3 NP4 Wricke’s ecovalence Shukla’s stability variance Deviation from regression Coefficient of variance Kang’s rank-sum Mean variance component GE variance component  SR      AR     Std.
1        G1  4 13 13  7  7  14   6   5   7                          9                                  9                         9                       3                      5                       9                    12 132  8.2500 3.376389
2        G2  2 11 11  6  6  10  11   7   6                         12                                 12                         7                      17                      6                      12                     9 145  9.0625 3.623419
3        G3 13 15 15 16 16  13  15  16  16                          6                                  6                         5                      19                      8                       6                    15 200 12.5000 4.618802
4        G4 17  7  9 15 15   4  18  15  15                         11                                 11                        12                      18                     17                      11                    10 205 12.8125 4.102337
5        G5  3  5  4  4  3   9   5   4   3                          3                                  3                         3                       2                      1                       3                    18  73  4.5625 3.982775
6        G6 20  6  8 20 20  18  20  20  20                         20                                 20                        20                      20                     20                      20                     1 273 17.0625 6.147832
7        G7  7  1  1  1  1   1   2   2   1                          4                                  4                         4                       9                      4                       4                    17  63  3.9375 4.202678
8        G8 19  4  5 10 14  10  19  14  14                          7                                  7                         8                      10                     14                       7                    14 176 11.0000 4.604346
9        G9  1 10 10  5  5  20   3   6   4                         17                                 17                        18                       5                      7                      17                     4 149  9.3125 6.353149
10      G10 18 16 17 18 19  17  17  19  18                         14                                 14                        14                      15                     19                      14                     7 256 16.0000 3.055050
11      G11 12 17 16 13 12   7   9   9  11                         13                                 13                        13                      11                     12                      13                     8 189 11.8125 2.663801
12      G12 15 18 18 17 17  16  14  18  17                         15                                 15                        15                      16                     18                      15                     6 250 15.6250 2.895399
13      G13  6 14 14 11  8   5   7   8   8                         19                                 19                        17                       4                     12                      19                     2 173 10.8125 5.659432
14      G14 11 19 19 14 13  15  10  11  13                         10                                 10                        11                      13                     10                      10                    11 200 12.5000 2.988868
15      G15 16  8  6  8  9  10  13  12   9                          8                                  8                        10                       8                     11                       8                    13 157  9.8125 2.587631
16      G16 14  9  7  9 11   5  12  13  12                          5                                  5                         6                      14                      8                       5                    16 151  9.4375 3.758878
17      G17  5  2  2  2  2   2   1   1   2                          1                                  1                         1                       7                      1                       1                    20  51  3.1875 4.777988
18      G18  9 20 20 19 18  19  16  17  19                         18                                 18                        19                      12                     16                      18                     3 261 16.3125 4.600272
19      G19 10 12 12 12 10   8   8  10  10                         16                                 16                        16                       6                     14                      16                     5 181 11.3125 3.591077
20      G20  8  3  3  3  4   3   4   3   5                          2                                  2                         2                       1                      3                       2                    19  67  4.1875 4.261748
```
Step 3: Plotting the correlation matrix
-------------------
##### 1. Install R library `ggcorrplot`
```R
install.packages("ggcorrplot")
```
##### 2. Load library `ggcorrplot`
```R
library("ggcorrplot")
```
##### 3. Plot the heatmap for pearson's correlation matrix available in `results$correlation_matrix` variable
```R
ggcorrplot(results$correlation_matrix)
```
![Spearman correlation heatmap](https://raw.githubusercontent.com/pour-aboughadareh/STABILITYSOFT/master/Screenshot%201.jpg)
