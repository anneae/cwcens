
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cwcens

<!-- badges: start -->

<!-- badges: end -->

This pacakage implements non-parametric estimation via a kernel
estimator for the probability in state and restricted mean time in state
in an illness-death model under component-wise censoring. Component-wise
censoring arises when illness can only be measured at a finite set of
times, while death is right censored and thus observed continuously up
to the right censoring time. Component-wise censored composite endpoints
arise often in biostatistical practice. For example, in many oncology
studies, progression-free survival is component-wise censored.

## Installation

You can install the development version of cwcens from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("anneae/cwcens")
```

## Example

First, we will simulate data from an illness-death model. Our dataset
will contain 150 patients. Since `scale21` is

``` r
library(cwcens)
simdat(150, scale12=1/.0008, scale13=1/.0002, scale23=1/.0016,
       scale21=1/.0008,
       vital.lfu=c(30.4*36, 30.4*48),
       visit.schedule = 30.4*c(6, 12, 18, 24, 30, 36, 42, 48),
       scatter.sd=10)
#>          dtime dstatus state2obs laststate1       t1       t2       t3       t4
#> 1    230.37275       1  181.9534     0.0000 181.9534       NA       NA       NA
#> 2   1289.19470       1       Inf  1089.0193 169.7163 376.0208 547.8971 744.7745
#> 3    617.88765       1       Inf   575.3473 181.0816 371.7089 575.3473       NA
#> 4   1320.45696       0       Inf  1261.7693 186.4723 358.2408 560.3400 726.7678
#> 5    357.34318       1  185.0992     0.0000 185.0992       NA       NA       NA
#> 6   1268.43230       0 1109.8591   898.9347 185.7930 375.5581 548.0406 746.1186
#> 7   1364.23589       1  719.7921   546.1382 171.2736 371.3974 546.1382 719.7921
#> 8   1304.80987       0       Inf  1284.8238 186.0548 346.4605 535.3857 731.3124
#> 9    639.82729       1  560.8456   352.3708 193.9494 352.3708 560.8456       NA
#> 10  1208.20909       0       Inf  1072.5866 178.3273 366.6908 537.8255 723.0632
#> 11  1304.49986       0       Inf  1290.6352 191.2869 361.0893 550.7671 729.1892
#> 12   752.60089       1  724.4381   549.9709 197.8821 362.2724 549.9709 724.4381
#> 13  1444.15495       0  550.2189   361.2271 187.4117 361.2271 550.2189 734.7715
#> 14  1275.00768       0       Inf  1094.3556 201.3952 364.3516 540.4535 745.8580
#> 15  1146.85401       1  564.1296   366.8275 183.2994 366.8275 564.1296 731.5374
#> 16   349.22745       1       Inf   179.0015 179.0015       NA       NA       NA
#> 17   357.83013       1       Inf   172.3997 172.3997       NA       NA       NA
#> 18   932.98787       1  547.7523   384.5074 187.9287 384.5074 547.7523 733.4108
#> 19   593.84372       1       Inf   555.7929 195.8043 368.3861 555.7929       NA
#> 20  1201.98279       0       Inf  1086.7552 189.7396 351.3935 537.1488 750.8831
#> 21   734.57016       1  551.5007   336.3858 177.8816 336.3858 551.5007 717.8998
#> 22   963.62935       1  532.8584   351.0006 173.5504 351.0006 532.8584 737.5098
#> 23  1265.00271       0  550.7277   359.6269 197.1684 359.6269 550.7277 737.1831
#> 24    72.18823       1       Inf     0.0000       NA       NA       NA       NA
#> 25  1384.29017       0 1088.5820   908.7403 188.2947 365.2388 547.4399 738.1370
#> 26   682.66746       1       Inf   540.8272 184.5972 355.9991 540.8272       NA
#> 27  1308.55077       0  909.3374   734.0507 188.1702 374.3055 536.0885 734.0507
#> 28   734.00399       1  181.1594     0.0000 181.1594 358.9206 546.3320 706.9617
#> 29  1216.00850       0       Inf  1096.7638 186.6330 365.6906 547.7815 742.1339
#> 30    68.76831       1       Inf     0.0000       NA       NA       NA       NA
#> 31  1199.11214       0       Inf  1098.8672 184.0961 375.9146 574.1203 732.3715
#> 32   728.82696       1       Inf   720.8170 184.4729 359.1317 551.0295 720.8170
#> 33  1224.90959       0       Inf  1086.2243 184.7226 363.7551 566.8623 719.6409
#> 34  1449.91812       0       Inf  1281.6340 177.6266 382.1399 551.0936 722.2433
#> 35  1125.27924       1  194.5984     0.0000 194.5984 361.2023 546.4954 723.8799
#> 36  1195.98941       0       Inf  1089.4407 201.6923 361.7416 543.5534 723.3744
#> 37  1425.72623       0       Inf  1272.5793 184.8340 381.5346 538.8925 732.1336
#> 38   605.24553       1       Inf   556.3987 173.9031 366.9060 556.3987       NA
#> 39   609.33447       1       Inf   543.3130 170.2912 362.1860 543.3130       NA
#> 40  1435.45695       0 1097.5720   891.8018 188.9606 368.1223 543.6004 726.9615
#> 41   780.14009       1  727.2887   551.3149 163.4336 349.7689 551.3149 727.2887
#> 42  1114.49596       0       Inf  1114.0161 186.9416 374.1117 533.8913 717.4026
#> 43   497.90931       1       Inf   384.9730 186.7706 384.9730       NA       NA
#> 44  1148.50504       0  731.8311   558.0196 194.5297 359.3352 558.0196 731.8311
#> 45  1248.94527       0  915.7333   714.6602 168.7312 359.6013 560.7346 714.6602
#> 46   708.58379       1  186.2406     0.0000 186.2406 369.9826 542.7325       NA
#> 47   581.87086       1  529.0478   380.2631 172.0568 380.2631 529.0478       NA
#> 48  1282.81576       0  367.5174   181.7288 181.7288 367.5174 554.0655 732.0684
#> 49  1252.82763       0       Inf  1106.0343 179.4193 383.4228 550.7212 730.3968
#> 50  1140.52263       0       Inf  1099.3471 179.0732 341.3408 555.1988 724.2399
#> 51  1186.50059       0  928.2224   725.9363 189.0014 365.5813 565.7747 725.9363
#> 52  1190.62969       1 1098.1888   928.0330 177.4773 363.9107 572.3074 750.9698
#> 53   188.85975       1  180.3638     0.0000 180.3638       NA       NA       NA
#> 54  1198.24033       0       Inf  1095.9312 164.4641 352.0873 562.2915 738.2557
#> 55   658.81122       1       Inf   555.7539 170.3308 380.7771 555.7539       NA
#> 56  1268.94673       0  721.7219   546.1323 169.5688 371.2885 546.1323 721.7219
#> 57   328.41783       1  166.2130     0.0000 166.2130       NA       NA       NA
#> 58  1416.88640       0 1091.1906   914.2970 184.9836 358.7954 541.4164 736.9382
#> 59  1203.93949       0  172.9539     0.0000 172.9539 347.4920 567.7560 727.9091
#> 60    26.96347       1       Inf     0.0000       NA       NA       NA       NA
#> 61  1322.11373       0       Inf  1282.7132 206.4053 355.9732 557.0563 739.6207
#> 62  1409.48989       0  176.5324     0.0000 176.5324 377.8038 539.3376 728.6982
#> 63  1305.00026       0       Inf  1295.8365 196.9166 364.9035 554.9671 719.1214
#> 64   966.74044       1  561.6585   375.3771 187.6066 375.3771 561.6585 726.4543
#> 65   856.70162       1       Inf   723.5379 207.7987 359.4492 534.6696 723.5379
#> 66  1304.28570       0  919.3869   729.2945 160.4370 358.4209 550.3047 729.2945
#> 67   526.29506       1  191.7549     0.0000 191.7549 378.0425       NA       NA
#> 68  1422.01612       0       Inf  1276.8398 175.0122 358.7035 541.5058 728.6205
#> 69  1121.78064       1 1078.2491   906.0268 185.3398 350.9885 556.0308 723.6909
#> 70  1213.89032       0       Inf  1079.4330 180.0305 370.8546 554.9154 727.2297
#> 71  1205.41776       0  901.5867   737.3121 174.8381 363.2928 561.8859 737.3121
#> 72  1427.47696       0       Inf  1284.6588 182.5551 367.2839 536.2012 735.2253
#> 73   594.09972       1       Inf   553.1984 178.1394 362.0821 553.1984       NA
#> 74  1353.69270       0       Inf  1284.4040 193.7740 363.5336 551.4035 735.4045
#> 75  1328.57383       0       Inf  1276.5962 192.5904 363.5950 543.7329 728.2179
#> 76   472.36094       1       Inf   355.5558 164.9214 355.5558       NA       NA
#> 77  1210.45226       0       Inf  1103.7603 191.3497 365.2833 535.1532 728.9695
#> 78  1093.36376       1  906.9110   730.5777 176.1969 369.3318 540.1578 730.5777
#> 79   320.95937       1       Inf   200.1491 200.1491       NA       NA       NA
#> 80   773.90513       1  749.3577   545.7006 182.4299 357.9232 545.7006 749.3577
#> 81   326.94099       1       Inf   185.3955 185.3955       NA       NA       NA
#> 82  1284.03545       0  551.2121   368.7625 177.7893 368.7625 551.2121 722.5551
#> 83  1186.68988       0       Inf  1084.8051 184.6879 363.0210 554.7870 720.9196
#> 84  1364.54065       0       Inf  1257.0563 171.5431 364.8862 566.4428 723.0808
#> 85  1431.33957       0  720.7303   554.2132 184.5566 353.6230 554.2132 720.7303
#> 86  1428.28259       0       Inf  1261.2020 185.9755 365.7002 538.7571 730.5346
#> 87  1202.03744       0       Inf  1110.0540 172.2772 363.1813 549.5054 720.3567
#> 88  1150.43894       0       Inf  1101.0677 171.9504 369.0196 531.8501 711.6122
#> 89  1300.69898       0 1288.6755  1101.4544 186.8457 360.0702 546.8156 731.6850
#> 90    64.22349       1       Inf     0.0000       NA       NA       NA       NA
#> 91  1317.09692       0  187.4812     0.0000 187.4812 363.4027 557.3459 721.9147
#> 92  1164.98798       0  914.1715   727.7720 170.8733 362.5440 568.1645 727.7720
#> 93   460.37350       1       Inf   351.7050 198.2553 351.7050       NA       NA
#> 94  1219.70384       0  716.6958   544.3145 176.1506 353.2872 544.3145 716.6958
#> 95  1407.02866       0       Inf  1279.2403 183.0777 359.6394 556.8921 730.9216
#> 96  1294.48706       1 1085.0423   911.7235 184.8980 351.7779 552.3542 722.7691
#> 97  1189.97458       0       Inf  1112.5939 165.1682 345.3389 527.9844 727.9520
#> 98   257.28745       1       Inf   190.2116 190.2116       NA       NA       NA
#> 99  1359.30940       0 1111.0082   916.0461 181.0725 374.2818 542.2395 720.5138
#> 100  594.19213       1  170.5941     0.0000 170.5941 379.0193 541.6241       NA
#> 101  760.96090       1       Inf   710.1861 179.7058 371.5552 550.9101 710.1861
#> 102 1457.88726       0       Inf  1452.9098 172.9102 368.2738 553.4847 729.1304
#> 103  575.03519       1  351.7058   193.6843 193.6843 351.7058 553.6700       NA
#> 104  159.84169       1       Inf     0.0000       NA       NA       NA       NA
#> 105  448.40098       1  374.2829   174.2343 174.2343 374.2829       NA       NA
#> 106 1059.67786       1  718.6897   540.9073 183.0614 374.0459 540.9073 718.6897
#> 107  329.79307       1       Inf   186.1817 186.1817       NA       NA       NA
#> 108 1127.02382       0       Inf  1085.3565 191.2356 361.1510 537.6755 719.5871
#> 109  732.13218       1       Inf   730.7132 175.1911 360.2933 537.6542 730.7132
#> 110 1096.46642       0       Inf  1091.3434 187.4173 361.8673 536.6312 730.9190
#> 111 1123.84775       1  172.1369     0.0000 172.1369 366.7233 545.9360 723.7665
#> 112 1147.45707       0       Inf  1091.1052 171.7382 372.0127 534.1392 746.0194
#> 113 1169.18023       0  900.2429   722.6827 178.0208 366.9461 535.2433 722.6827
#> 114 1352.28608       0       Inf  1273.6187 184.5140 344.8536 555.6964 738.1582
#> 115 1165.39957       0       Inf  1092.1473 173.8371 348.3159 556.4944 722.0353
#> 116 1402.11639       0  717.1140   545.8597 194.1206 371.0026 545.8597 717.1140
#> 117  587.07195       1       Inf   558.0800 179.5168 346.7216 558.0800       NA
#> 118 1400.73614       0 1265.6785  1083.3528 197.1856 359.3489 534.3480 753.6195
#> 119  812.68753       1  540.0343   374.6195 171.9752 374.6195 540.0343 743.1632
#> 120  866.01556       1  179.4071     0.0000 179.4071 360.2612 537.0749 726.2283
#> 121 1191.14744       0       Inf  1088.5027 182.2688 374.7590 559.5017 713.7534
#> 122 1009.10474       1  549.6578   370.6557 171.7459 370.6557 549.6578 749.4008
#> 123 1173.03727       0 1086.6007   919.0318 177.8049 373.2399 548.4955 732.4602
#> 124 1162.85503       0       Inf  1105.7987 162.5823 368.3321 530.6344 722.4366
#> 125  489.35957       1       Inf   363.7584 198.0690 363.7584       NA       NA
#> 126  979.94547       1  368.7924   171.7390 171.7390 368.7924 546.7469 731.9572
#> 127 1436.43995       0       Inf  1295.9736 177.6063 367.9323 540.7671 730.9133
#> 128 1059.06722       1  910.5796   730.1661 194.9643 361.7521 533.0982 730.1661
#> 129  359.22578       1  193.1507     0.0000 193.1507       NA       NA       NA
#> 130 1207.05333       0       Inf  1098.9841 191.6174 353.2284 545.7324 735.5025
#> 131  382.52073       1       Inf   358.1549 160.4654 358.1549       NA       NA
#> 132 1197.95356       0       Inf  1101.9552 186.8786 373.4752 523.3592 727.8927
#> 133 1109.62795       0  183.4179     0.0000 183.4179 353.0874 566.8423 728.4054
#> 134 1277.97760       0       Inf  1101.2955 195.3581 373.9091 552.0955 724.2882
#> 135  755.71309       1  167.1572     0.0000 167.1572 355.7026 548.5224 729.1991
#> 136  378.98326       1  192.5236     0.0000 192.5236 366.5731       NA       NA
#> 137 1179.18227       0       Inf  1081.0390 177.4421 343.8176 548.0472 731.9474
#> 138 1106.17565       0       Inf  1093.0258 177.9526 351.0235 546.7026 746.8668
#> 139  694.39936       1       Inf   532.2538 193.4197 356.7880 532.2538       NA
#> 140  546.47325       1  164.1011     0.0000 164.1011 378.0456       NA       NA
#> 141 1106.56923       0       Inf  1093.3777 185.7883 360.5954 545.6590 717.9909
#> 142  290.82934       1       Inf   179.2435 179.2435       NA       NA       NA
#> 143 1160.31208       0       Inf  1092.6675 202.9032 367.9244 550.2098 716.5636
#> 144 1119.81810       1 1092.0683   928.4467 175.4285 366.0209 537.6529 744.2866
#> 145 1152.62051       0       Inf  1079.7176 177.7386 356.4475 549.4354 732.6345
#> 146 1112.28156       0       Inf  1109.1743 184.6641 362.7659 540.2386 743.8058
#> 147 1329.96253       0 1277.9193  1102.1593 179.4467 344.4837 542.6736 714.0849
#> 148 1305.30732       0       Inf  1288.2201 186.3509 378.2251 535.3457 719.7863
#> 149 1240.37714       0       Inf  1101.6955 188.4761 357.3266 542.1394 708.8118
#> 150 1446.51836       0       Inf  1259.1700 164.9883 359.3787 537.4179 730.5927
#>           t5       t6       t7      t8 x1 x2 x3 x4 x5 x6 x7 x8 nvisits
#> 1         NA       NA       NA      NA  2 NA NA NA NA NA NA NA       8
#> 2   920.7092 1089.019       NA      NA  1  1  1  1  1  1 NA NA       8
#> 3         NA       NA       NA      NA  1  1  1 NA NA NA NA NA       8
#> 4   882.2011 1090.453 1261.769      NA  1  1  1  1  1  1  1 NA       8
#> 5         NA       NA       NA      NA  2 NA NA NA NA NA NA NA       8
#> 6   898.9347 1109.859       NA      NA  1  1  1  1  1  2 NA NA       8
#> 7   886.9166 1101.962 1291.110      NA  1  1  1  2  2  1  1 NA       8
#> 8   924.5874 1106.283 1284.824      NA  1  1  1  1  1  1  1 NA       8
#> 9         NA       NA       NA      NA  1  1  2 NA NA NA NA NA       8
#> 10  910.0832 1072.587       NA      NA  1  1  1  1  1  1 NA NA       8
#> 11  916.4234 1087.064 1290.635      NA  1  1  1  1  1  1  1 NA       8
#> 12        NA       NA       NA      NA  1  1  1  2 NA NA NA NA       8
#> 13  911.4355 1096.267 1264.885      NA  1  1  2  1  1  1  2 NA       8
#> 14  898.9411 1094.356       NA      NA  1  1  1  1  1  1 NA NA       8
#> 15  914.9536 1112.360       NA      NA  1  1  2  2  1  1 NA NA       8
#> 16        NA       NA       NA      NA  1 NA NA NA NA NA NA NA       8
#> 17        NA       NA       NA      NA  1 NA NA NA NA NA NA NA       8
#> 18  916.9037       NA       NA      NA  1  1  2  2  2 NA NA NA       8
#> 19        NA       NA       NA      NA  1  1  1 NA NA NA NA NA       8
#> 20  924.6115 1086.755       NA      NA  1  1  1  1  1  1 NA NA       8
#> 21        NA       NA       NA      NA  1  1  2  2 NA NA NA NA       8
#> 22  910.9994       NA       NA      NA  1  1  2  2  2 NA NA NA       8
#> 23  921.8130 1102.900       NA      NA  1  1  2  2  2  1 NA NA       8
#> 24        NA       NA       NA      NA NA NA NA NA NA NA NA NA       8
#> 25  908.7403 1088.582 1282.751      NA  1  1  1  1  1  2  1 NA       8
#> 26        NA       NA       NA      NA  1  1  1 NA NA NA NA NA       8
#> 27  909.3374 1104.446 1283.923      NA  1  1  1  1  2  1  1 NA       8
#> 28        NA       NA       NA      NA  2  2  2  2 NA NA NA NA       8
#> 29  898.4296 1096.764       NA      NA  1  1  1  1  1  1 NA NA       8
#> 30        NA       NA       NA      NA NA NA NA NA NA NA NA NA       8
#> 31  899.6830 1098.867       NA      NA  1  1  1  1  1  1 NA NA       8
#> 32        NA       NA       NA      NA  1  1  1  1 NA NA NA NA       8
#> 33  925.2190 1086.224       NA      NA  1  1  1  1  1  1 NA NA       8
#> 34  910.6784 1084.776 1281.634      NA  1  1  1  1  1  1  1 NA       8
#> 35  897.9269 1107.962       NA      NA  2  2  2  1  1  2 NA NA       8
#> 36  919.9559 1089.441       NA      NA  1  1  1  1  1  1 NA NA       8
#> 37  909.8775 1097.665 1272.579      NA  1  1  1  1  1  1  1 NA       8
#> 38        NA       NA       NA      NA  1  1  1 NA NA NA NA NA       8
#> 39        NA       NA       NA      NA  1  1  1 NA NA NA NA NA       8
#> 40  891.8018 1097.572 1279.844      NA  1  1  1  1  1  2  2 NA       8
#> 41        NA       NA       NA      NA  1  1  1  2 NA NA NA NA       8
#> 42  911.2066 1114.016       NA      NA  1  1  1  1  1  1 NA NA       8
#> 43        NA       NA       NA      NA  1  1 NA NA NA NA NA NA       8
#> 44  931.1837 1079.447       NA      NA  1  1  1  2  2  2 NA NA       8
#> 45  915.7333 1099.756       NA      NA  1  1  1  1  2  2 NA NA       8
#> 46        NA       NA       NA      NA  2  2  2 NA NA NA NA NA       8
#> 47        NA       NA       NA      NA  1  1  2 NA NA NA NA NA       8
#> 48  899.4449 1085.043       NA      NA  1  2  2  1  1  1 NA NA       8
#> 49  906.6488 1106.034       NA      NA  1  1  1  1  1  1 NA NA       8
#> 50  941.7558 1099.347       NA      NA  1  1  1  1  1  1 NA NA       8
#> 51  928.2224 1098.546       NA      NA  1  1  1  1  2  2 NA NA       8
#> 52  928.0330 1098.189       NA      NA  1  1  1  1  1  2 NA NA       8
#> 53        NA       NA       NA      NA  2 NA NA NA NA NA NA NA       8
#> 54  911.4289 1095.931       NA      NA  1  1  1  1  1  1 NA NA       8
#> 55        NA       NA       NA      NA  1  1  1 NA NA NA NA NA       8
#> 56  930.4944 1091.480       NA      NA  1  1  1  2  2  1 NA NA       8
#> 57        NA       NA       NA      NA  2 NA NA NA NA NA NA NA       8
#> 58  914.2970 1091.191 1282.727      NA  1  1  1  1  1  2  2 NA       8
#> 59  900.0248 1084.931       NA      NA  2  2  2  2  2  2 NA NA       8
#> 60        NA       NA       NA      NA NA NA NA NA NA NA NA NA       8
#> 61  910.1771 1101.488 1282.713      NA  1  1  1  1  1  1  1 NA       8
#> 62  909.4172 1098.528 1266.960      NA  2  2  2  2  2  2  1 NA       8
#> 63  928.8337 1105.219 1295.836      NA  1  1  1  1  1  1  1 NA       8
#> 64  912.0065       NA       NA      NA  1  1  2  2  2 NA NA NA       8
#> 65        NA       NA       NA      NA  1  1  1  1 NA NA NA NA       8
#> 66  919.3869 1085.224 1277.117      NA  1  1  1  1  2  2  2 NA       8
#> 67        NA       NA       NA      NA  2  2 NA NA NA NA NA NA       8
#> 68  906.7254 1105.121 1276.840      NA  1  1  1  1  1  1  1 NA       8
#> 69  906.0268 1078.249       NA      NA  1  1  1  1  1  2 NA NA       8
#> 70  913.6883 1079.433       NA      NA  1  1  1  1  1  1 NA NA       8
#> 71  901.5867 1092.870       NA      NA  1  1  1  1  2  2 NA NA       8
#> 72  913.9082 1086.303 1284.659      NA  1  1  1  1  1  1  1 NA       8
#> 73        NA       NA       NA      NA  1  1  1 NA NA NA NA NA       8
#> 74  912.0958 1081.389 1284.404      NA  1  1  1  1  1  1  1 NA       8
#> 75  914.9389 1093.738 1276.596      NA  1  1  1  1  1  1  1 NA       8
#> 76        NA       NA       NA      NA  1  1 NA NA NA NA NA NA       8
#> 77  910.7991 1103.760       NA      NA  1  1  1  1  1  1 NA NA       8
#> 78  906.9110       NA       NA      NA  1  1  1  1  2 NA NA NA       8
#> 79        NA       NA       NA      NA  1 NA NA NA NA NA NA NA       8
#> 80        NA       NA       NA      NA  1  1  1  2 NA NA NA NA       8
#> 81        NA       NA       NA      NA  1 NA NA NA NA NA NA NA       8
#> 82  898.7032 1090.231 1266.820      NA  1  1  2  2  2  1  1 NA       8
#> 83  900.0705 1084.805       NA      NA  1  1  1  1  1  1 NA NA       8
#> 84  917.7154 1102.575 1257.056      NA  1  1  1  1  1  1  1 NA       8
#> 85  910.0484 1095.040 1288.523      NA  1  1  1  2  2  2  2 NA       8
#> 86  894.3231 1097.376 1261.202      NA  1  1  1  1  1  1  1 NA       8
#> 87  918.9940 1110.054       NA      NA  1  1  1  1  1  1 NA NA       8
#> 88  912.4428 1101.068       NA      NA  1  1  1  1  1  1 NA NA       8
#> 89  902.1081 1101.454 1288.675      NA  1  1  1  1  1  1  2 NA       8
#> 90        NA       NA       NA      NA NA NA NA NA NA NA NA NA       8
#> 91  902.3358 1096.978 1276.756      NA  2  2  1  1  1  1  1 NA       8
#> 92  914.1715 1089.805       NA      NA  1  1  1  1  2  2 NA NA       8
#> 93        NA       NA       NA      NA  1  1 NA NA NA NA NA NA       8
#> 94  897.9363 1082.441       NA      NA  1  1  1  2  2  2 NA NA       8
#> 95  915.8898 1091.639 1279.240      NA  1  1  1  1  1  1  1 NA       8
#> 96  911.7235 1085.042 1274.815      NA  1  1  1  1  1  2  2 NA       8
#> 97  912.2864 1112.594       NA      NA  1  1  1  1  1  1 NA NA       8
#> 98        NA       NA       NA      NA  1 NA NA NA NA NA NA NA       8
#> 99  916.0461 1111.008 1272.393      NA  1  1  1  1  1  2  2 NA       8
#> 100       NA       NA       NA      NA  2  2  1 NA NA NA NA NA       8
#> 101       NA       NA       NA      NA  1  1  1  1 NA NA NA NA       8
#> 102 907.0000 1091.988 1273.672 1452.91  1  1  1  1  1  1  1  1       8
#> 103       NA       NA       NA      NA  1  2  2 NA NA NA NA NA       8
#> 104       NA       NA       NA      NA NA NA NA NA NA NA NA NA       8
#> 105       NA       NA       NA      NA  1  2 NA NA NA NA NA NA       8
#> 106 918.1152       NA       NA      NA  1  1  1  2  2 NA NA NA       8
#> 107       NA       NA       NA      NA  1 NA NA NA NA NA NA NA       8
#> 108 899.1334 1085.356       NA      NA  1  1  1  1  1  1 NA NA       8
#> 109       NA       NA       NA      NA  1  1  1  1 NA NA NA NA       8
#> 110 912.5367 1091.343       NA      NA  1  1  1  1  1  1 NA NA       8
#> 111 910.2380 1109.817       NA      NA  2  2  1  1  2  2 NA NA       8
#> 112 909.1584 1091.105       NA      NA  1  1  1  1  1  1 NA NA       8
#> 113 900.2429 1071.060       NA      NA  1  1  1  1  2  2 NA NA       8
#> 114 920.4301 1083.221 1273.619      NA  1  1  1  1  1  1  1 NA       8
#> 115 917.1733 1092.147       NA      NA  1  1  1  1  1  1 NA NA       8
#> 116 897.7338 1090.658 1275.908      NA  1  1  1  2  2  2  1 NA       8
#> 117       NA       NA       NA      NA  1  1  1 NA NA NA NA NA       8
#> 118 917.1267 1083.353 1265.678      NA  1  1  1  1  1  1  2 NA       8
#> 119       NA       NA       NA      NA  1  1  2  1 NA NA NA NA       8
#> 120       NA       NA       NA      NA  2  2  2  2 NA NA NA NA       8
#> 121 926.1880 1088.503       NA      NA  1  1  1  1  1  1 NA NA       8
#> 122 913.3279       NA       NA      NA  1  1  2  2  2 NA NA NA       8
#> 123 919.0318 1086.601       NA      NA  1  1  1  1  1  2 NA NA       8
#> 124 910.4165 1105.799       NA      NA  1  1  1  1  1  1 NA NA       8
#> 125       NA       NA       NA      NA  1  1 NA NA NA NA NA NA       8
#> 126 892.3439       NA       NA      NA  1  2  2  2  2 NA NA NA       8
#> 127 908.4775 1096.397 1295.974      NA  1  1  1  1  1  1  1 NA       8
#> 128 910.5796       NA       NA      NA  1  1  1  1  2 NA NA NA       8
#> 129       NA       NA       NA      NA  2 NA NA NA NA NA NA NA       8
#> 130 920.4606 1098.984       NA      NA  1  1  1  1  1  1 NA NA       8
#> 131       NA       NA       NA      NA  1  1 NA NA NA NA NA NA       8
#> 132 902.0414 1101.955       NA      NA  1  1  1  1  1  1 NA NA       8
#> 133 903.7261 1086.551       NA      NA  2  2  2  2  2  2 NA NA       8
#> 134 922.6944 1101.295       NA      NA  1  1  1  1  1  1 NA NA       8
#> 135       NA       NA       NA      NA  2  1  1  2 NA NA NA NA       8
#> 136       NA       NA       NA      NA  2  2 NA NA NA NA NA NA       8
#> 137 910.6981 1081.039       NA      NA  1  1  1  1  1  1 NA NA       8
#> 138 904.2585 1093.026       NA      NA  1  1  1  1  1  1 NA NA       8
#> 139       NA       NA       NA      NA  1  1  1 NA NA NA NA NA       8
#> 140       NA       NA       NA      NA  2  2 NA NA NA NA NA NA       8
#> 141 921.1083 1093.378       NA      NA  1  1  1  1  1  1 NA NA       8
#> 142       NA       NA       NA      NA  1 NA NA NA NA NA NA NA       8
#> 143 912.8681 1092.668       NA      NA  1  1  1  1  1  1 NA NA       8
#> 144 928.4467 1092.068       NA      NA  1  1  1  1  1  2 NA NA       8
#> 145 917.2325 1079.718       NA      NA  1  1  1  1  1  1 NA NA       8
#> 146 921.6367 1109.174       NA      NA  1  1  1  1  1  1 NA NA       8
#> 147 903.2132 1102.159 1277.919      NA  1  1  1  1  1  1  2 NA       8
#> 148 913.5177 1097.123 1288.220      NA  1  1  1  1  1  1  1 NA       8
#> 149 905.4016 1101.695       NA      NA  1  1  1  1  1  1 NA NA       8
#> 150 905.6806 1082.754 1259.170      NA  1  1  1  1  1  1  1 NA       8
```

The standard approach would …

The kernel approach does not …

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub\!
