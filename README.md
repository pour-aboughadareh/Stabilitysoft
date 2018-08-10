How to use STABILITYSOFT in R
================================================

STABILITYSOFT calculates the several parametric and non-parametric statistics as defined in STABILITYSOFT: a new online program to calculate parametric and non-parametric stability statistics. Applications in Plant Sciences (submitted).

Step 1: loading the library
--------------------------
To get started, execute the library code ([STABILITYSOFT.R](STABILITYSOFT.R)) in your RStudio console.

```R
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
    ShuklaEquvalance <- function(a, b, pureMatShu)
    {
        equvalance <- apply(pureMatShu, 1, sum)
        equvalanceTotal <- sum(pureMatShu)
        output <- matrix(nrow = a, ncol = 2)
        shukla <- ((a * equvalance) / ((b - 1) * (a - 2))) - ((equvalanceTotal) / ((a - 1) * (a - 2) * (b - 1)))
        output[, 1] <- equvalance
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
        equvalance <- apply(pureMatShu, 1, sum)
        equvalanceTotal <- sum(equvalance)

        shukla <- ((a * equvalance) / ((b - 1) * (a - 2))) - ((equvalanceTotal) / ((a - 1) * (a - 2) * (b - 1)))

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
        equvalance <- apply(pureMatShu, 1, sum)
        equvalanceTotal <- sum(equvalance)

        shukla <- (((-1 * a) * equvalance) / ((b - 1) * (a - 2) * (a - 1))) + ((equvalanceTotal) / ((a - 2) * (b - 1)))

        output <- shukla
        return(output)
    }
    PaP <- function(a, b, pureMatShu)
    {
        equvalance <- apply(pureMatShu, 1, sum)
        equvalanceTotal <- sum(equvalance)

        shukla <- ((a * equvalance) / (2 * (b - 1) * (a - 1))) + ((equvalanceTotal) / (2 * (a - 2) * (b - 1)))

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

        ShuklaEquvalance <- ShuklaEquvalance(a, b, pureMatShu) # wri Shu

        SDI <- SDI(a, b, pureMatSDI, pureMatAvgCol, pureMatTotalAvg, pureMatBI)

        BI <- BI(a, b, pureMatBI, pureMatAvgCol, pureMatTotalAvg)

        CVR <- CVR(a, b, pureMatSDI, pureMatCVR, pureMatAvgCol, pureMatAvg, pureMatTotalAvg, pureMatBI, pureMatBPJ, pureMat)

        Kang <- Kang(a, b, pureMatShu, pureMatAvg)

        P <- P(a, b, pureMatShu)

        PaP <- PaP(a, b, pureMatShu)

        stats_df <- data.frame(table_original[, 1], Y, s1, z1, s2, z2, s3, s6, np1, np2, np3, np4, ShuklaEquvalance[, 1], ShuklaEquvalance[, 2], SDI, BI, CVR, P, PaP, Kang)
        colnames(stats_df) <- c("Genotype", "Y", "S1", "Z1", "S2", "Z2", "S3", "S6", "NP1", "NP2", "NP3", "NP4", "Wri", "Shu", "Sd", "bFW", "CV", "Pla", "Pla & Pet", "Kang")

        ranks_df <- getranks_df(a, b, stats_df)

        output <- list(statistics = stats_df, ranks = ranks_df)
        return(output)
    }
})()
```
Step 2: loading the data
--------------------
After execution, function `Calculator` can be used to get the results.
It takes a input as dataframe and returns the result as a dataframe containing `result$statistics` and `result$ranks`.

First, you need to convert your data to a dataframe (given that the data is already like examples provided [here](Examples.zip)).

### Option 1: Copy data from your favorite app and get it from clipboard

#### 1. Copy data as table in your favorite app
![Copying table data from excel](SC1.jpg)

#### 2. Execute the followin command in your RStudio console
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

Step 3: getting the results
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

   Genotype         Y       S1         Z1        S2          Z2        S3        S6       NP1       NP2        NP3       NP4 Wricke’s ecovalence Shukla’s stability variance Deviation from regression Regression coefficient Coefficient of variance GE variance component Mean variance component Kang’s rank-sum
1        G1  50952.18 5.029240  4.6360964 18.777778  3.97691278  42.25000  8.250000 0.6315789 0.6140351 0.18232114 0.6286550         17550529389                   750727832                18983990.4              0.8362537                311.8356            6263106144              3673236078              19
2        G2  63643.48 5.111111  4.1795492 21.099415  2.80330261  37.58333  6.739583 0.6842105 0.3301435 0.17368054 0.5057870          7193577069                   111409788               150720497.4              1.0972096                327.7294            6296754462              3370401215              11
3        G3  47057.32 6.081871  0.5696503 28.672515  0.39785934  44.17117  7.954955 1.5263158 0.3298246 0.30464574 0.5205205         33269980625                  1721064328                55330141.5              0.7750115                313.0281            6212035802              4132869155              25
4        G4  58499.32 5.029240  4.6360964 18.982456  3.86521853  27.86266  5.484979 0.2105263 0.3036437 0.04248710 0.4101097           421355496                  -306628581                  479722.2              0.9746333                316.5035            6318756481              3172383040              11
5        G5  30198.76 5.801170  1.2716183 32.257310  0.01871126  88.25600 13.104000 2.4210526 2.8552632 0.85033640 0.8817778        192163246246                 11529290601               209001404.3              0.4581833                289.7705            5695813367              8778871074              35
6        G6 102738.85 7.099415  0.3564597 48.029240  4.14743712  57.03472  7.402778 3.0000000 0.7078947 0.45706896 0.4683642        393948413527                 23985165125               442433698.9              1.7756801                328.5830            5040241023             14679022164              21
7        G7  39183.73 7.391813  0.9711887 48.912281  4.65785138 132.76190 16.952381 1.9473684 2.7500000 0.68050638 1.1146384         94406627627                  5494931427               101931462.4              0.6202213                301.1251            6013411218              5920490412              30
8        G8  50272.94 5.485380  2.3937749 22.315789  2.27012872  39.13846  7.384615 1.1052632 0.4578947 0.25070578 0.5344729         20384561781                   925668103                22029039.9              0.8235269                311.2479            6253898761              3756102522              22
9        G9  86829.03 6.923977  0.1324774 36.578947  0.21042119  63.50254  9.725888 2.1052632 0.6473684 0.46926182 0.6677947        154325052509                  9193599630               177951463.1              1.4854396                325.1364            5818744470              7672491140              18
10      G10  68477.07 5.450292  2.5401873 24.099415  1.58991285  31.21970  5.250000 1.0526316 0.4638158 0.17547325 0.3922559         15143013439                   602115737                21173897.9              1.1519312                319.5992            6270927833              3602840875              12
11      G11  68325.82 3.543860 17.0277421  9.555556 10.66028171  15.63636  4.181818 0.7368421 0.2938596 0.15563116 0.3221691         14571780520                   566854445                20726988.0              1.1490254                319.4982            6272783691              3586138158              12
12      G12  75470.61 4.128655 11.2196542 15.479532  5.99616181  19.32117  3.897810 1.4210526 0.5263158 0.22756872 0.2862936         51163313971                  2825591078                61414983.5              1.2794627                322.1250            6153902815              4656066036              17
13      G13  88341.06 7.368421  0.9109048 39.362573  0.70945320  71.98930 10.449198 2.6315789 0.7224880 0.60444369 0.7486631        176220631746                 10545178595               191380401.2              1.5188575                326.7649            5747608735              8312712755              18
14      G14  52019.76 2.947368 24.1955732  8.093567 12.01638224  11.87983  3.407725 0.1052632 0.1818182 0.02502554 0.2403433         14204829584                   544203153                16530067.9              0.8527285                311.4459            6273975864              3575408598              14
15      G15  51046.86 2.877193 25.1214149  6.257310 13.83464300  11.95531  4.167598 0.5789474 0.3333333 0.15135438 0.3054004         17671139429                   758172896                19709662.4              0.8357115                311.0573            6262714298              3676762687              19
16      G16  47154.69 4.070175 11.7461423 12.152047  8.45193738  23.08889  6.066667 1.4736842 0.5119617 0.36141013 0.4296296         36438649795                  1916661191                40150031.0              0.7640736                307.9343            6201741230              4225520300              25
17      G17  20562.80 7.391813  0.9711887 55.654971  9.53156873 164.08621 19.775862 3.0000000 6.9473684 1.13479191 1.2107280        337714602093                 20513942197               366363986.0              0.2817156                266.6597            5222936967             13034758672              39
18      G18  88165.19 4.198830 10.6038026 18.339181  4.22161477  20.23226  3.541935 2.4736842 0.7397661 0.34454271 0.2573477        170126885628                 10169021427               196064335.3              1.5096880                325.4454            5767406481              8134533044              18
19      G19  80496.58 5.192982  3.7466617 19.766082  3.45229449  35.20833  6.739583 1.7368421 0.5119617 0.39692831 0.5138889         92132800745                  5354571743               100490623.6              1.3751627                324.6294            6020798570              5854004246              17
20      G20  28679.50 5.789474  1.3069026 36.561404  0.20820916 115.77778 15.574074 2.6315789 4.4385965 1.06967490 1.0185185        208669640886                 12548203851               236598327.1              0.4354842                290.3940            5642186354              9261514192              37


> print(results$ranks)

   Genotype  Y S1 S2 S3 S6 NP1 NP2 NP3 NP4 Wricke’s ecovalence Shukla’s stability variance Deviation from regression Coefficient of variance Kang’s rank-sum Mean variance component GE variance component  SR      AR     Std.
1        G1 13  7  7 12 14   4  12   7  14                   6                           6                         3                       9              11                       6                    15 146  9.1250 3.862210
2        G2  9  9 10 10  9   5   5   5  10                   2                           2                        13                      19               1                       2                    19 130  8.1250 5.572253
3        G3 16 15 13 13 13  11   4  10  12                   9                           9                         9                      10              15                       9                    12 180 11.2500 3.022141
4        G4 10  7  8  7  7   2   3   2   7                   1                           1                         1                      11               1                       1                    20  89  5.5625 5.214962
5        G5 18 14 14 17 17  15  18  18  17                  17                          17                        17                       2              18                      17                     4 240 15.0000 4.871687
6        G6  1 17 18 14 12  19  14  14   9                  20                          20                        20                      20              13                      20                     1 232 14.5000 6.292853
7        G7 17 19 19 19 19  13  17  17  19                  13                          13                        12                       4              17                      13                     8 239 14.9375 4.404070
8        G8 14 12 11 11 11   8   7   9  13                   8                           8                         7                       7              14                       8                    13 161 10.0625 2.594064
9        G9  4 16 16 15 15  14  13  15  15                  14                          14                        14                      16               8                      14                     7 210 13.1250 3.556684
10      G10  7 11 12  8  6   7   8   6   6                   5                           5                         6                      13               3                       5                    16 124  7.7500 3.492850
11      G11  8  3  3  3  5   6   2   4   5                   4                           4                         5                      12               3                       4                    17  88  5.5000 3.898718
12      G12  6  5  5  4  3   9  11   8   3                  11                          11                        10                      14               6                      11                    10 127  7.9375 3.395463
13      G13  2 18 17 16 16  17  15  16  16                  16                          16                        15                      18               8                      16                     5 227 14.1875 4.764014
14      G14 11  2  2  1  1   1   1   1   1                   3                           3                         2                       8               5                       3                    18  63  3.9375 4.697074
15      G15 12  1  1  2  4   3   6   3   4                   7                           7                         4                       6              11                       7                    14  92  5.7500 3.855732
16      G16 15  4  4  6  8  10   9  12   8                  10                          10                         8                       5              15                      10                    11 145  9.0625 3.336041
17      G17 20 19 20 20 20  19  20  20  20                  19                          19                        19                       1              20                      19                     2 277 17.3125 6.193747
18      G18  3  6  6  5  2  16  16  11   2                  15                          15                        16                      17               8                      15                     6 159  9.9375 5.697587
19      G19  5 10  9  9  9  12   9  13  11                  12                          12                        11                      15               6                      12                     9 164 10.2500 2.542964
20      G20 19 13 15 18 18  17  19  19  18                  18                          18                        18                       3              19                      18                     3 253 15.8125 5.243647
```
