Calculate management legacies, tree demography and climate data
================

<style type="text/css">
pre {
  font-size: 10px
}
</style>

``` r
library(tidyverse)
library(here)
library(testthat)
library(factoextra)
library(RColorBrewer)
library(ggh4x)
library(sf)
library(patchwork)
library(janitor)
```

``` r
tree23 <- read_csv(
  file = here("01-data", "ifn", "tree23.csv")
  )

Pme2 <- read_csv(
  here("01-data", "ifn", "Pme2_clean.csv")
  )

Pme3 <- read_csv(
  here("01-data", "ifn", "Pme3_clean.csv")
)
```

# Forest type definition

``` r
map(tree23 |> select(Cla3, Subclase3),
    ~levels(as.factor(.x)))
```

    ## $Cla3
    ## [1] "A"
    ## 
    ## $Subclase3
    ## [1] "1"  "3C"

``` r
# only pinus sylvestris trees (Psy)
# I select those plots of pinus sylvestris and
# where the % of pinus sylvestris is >= 50%
Psy23 <- tree23 |> 
  group_by(IDPC3) |> 
  mutate(
    # initial basal area
    AB2m2ha.tot.plot = sum(AB2m2ha, na.rm = T),
  ) |> 
  group_by(IDPC3, sppcompa) |> 
  mutate(
    AB2m2ha.plot = sum(AB2m2ha, na.rm = T), 
    AB2m2ha.rel.plot = AB2m2ha.plot * 100 / AB2m2ha.tot.plot
  ) |> 
  filter(
    # Psy = 21; Ppa = 26; Pni = 25; Pha = 24
    sppcompa == "021" & AB2m2ha.rel.plot >= 50
  )

nrow(Psy23) / nrow(tree23) # 14% Psy
```

    ## [1] 0.1469545

``` r
#  16% Ppa
#  7% Pni
#  9% Pha

# all trees in plots where there is more than one Psy
Psy23 <- Psy23 |> 
  group_by(IDPC3) |> 
  mutate(ntrees = length(IDPC3)) |> 
  filter(ntrees > 1)
```

# Management legacy

``` r
# forest structure, sp richness, regeneration
# join with information from Ceballos
planted <- read_tsv(
  here("01-data", "ifn", "Plotcode_RP.txt"),
  col_types = cols(.default = "c",
                   PLOTCODE = col_integer())
)

names(planted)
```

    ##  [1] "PLOTCODE" "CX"       "CY"       "Aal"      "Api"      "Fsy"     
    ##  [7] "Pca"      "Pha"      "Pni"      "Ppa"      "Ppe"      "Pra"     
    ## [13] "Psy"      "Pun"      "Qca"      "Qfa"      "Qhumilis" "Qil"     
    ## [19] "Qpe"      "Qpy"      "Qro"      "Qsu"      "Pinus_P"

``` r
dat <- planted |>
  dplyr::select(PLOTCODE, Psy) |>
  mutate(Psy = if_else(is.na(Psy), "A", Psy))

test_that("NA is natural", {
  expect_equal(sum(is.na(planted$Psy)) + sum(planted$Psy == "A", na.rm = T),
               sum(dat$Psy == "A", na.rm = T))
})
```

    ## Test passed 🥇

``` r
Psy23_planted <- Psy23 |>
  left_join(dat, by = c("Plotcode3" = "PLOTCODE"))

# A = natural
table(Psy23_planted$Psy, useNA = "always")
```

    ## 
    ##     A     P  <NA> 
    ## 96201 46087     0

``` r
Psy23_planted <- Psy23_planted |> 
  mutate(
    hdbh2 = h2 * 1000 / dbh2,
    hdbh3 = h3 * 1000 / dbh3,
    Psy = factor(Psy, levels = c("P", "A"), 
                 labels = c("Planted", "Natural"))
  )

summary(Psy23_planted$Psy)
```

    ## Planted Natural 
    ##   46087   96201

``` r
# richness: I group by plot and species and 
# select only the alive trees from ifn2
table(tree23$state3)
```

    ## 
    ##     MA  MA-nc     MP  MP-nc      R      V   V-nc 
    ##  97473      3  26031    626 264260 577089   3117

``` r
spp2 <- tree23 |>
  filter(state3 == "V" | state3 == "MA" | state3 == "MP") |>
  group_by(Plotcode2, sppcompa) |>
  summarise(NoTrees2 = length(Plotcode2))

names(spp2)
```

    ## [1] "Plotcode2" "sppcompa"  "NoTrees2"

``` r
sppriqueza2 <- spp2 |>
  group_by(Plotcode2) |>
  summarise(Riqueza2 = length(Plotcode2))

# species richness ifn3
table(tree23$state3)
```

    ## 
    ##     MA  MA-nc     MP  MP-nc      R      V   V-nc 
    ##  97473      3  26031    626 264260 577089   3117

``` r
spp3 <- tree23 |>
  filter(state3 == "V" | state3 == "R") |>
  group_by(Plotcode3, sppcompa) |>
  summarise(NoTrees3 = length(Plotcode3))

names(spp3)
```

    ## [1] "Plotcode3" "sppcompa"  "NoTrees3"

``` r
sppriqueza3 <- spp3 |>
  group_by(Plotcode3) |>
  summarise(Riqueza3 = length(Plotcode3))

Psy2 <- Psy23_planted |>
  group_by(Plotcode2) |>
  summarise(type2 = first(Psy),
            IDPC3 = first(IDPC3),
            ba_ha2 = sum(AB2m2ha, na.rm = TRUE),
            dens2 = sum(dens2, na.rm = TRUE),
            mdbh2 = mean(dbh2, na.rm = TRUE),
            sddbh2 = sd(dbh2, na.rm = TRUE),
            mh2 = mean(h2, na.rm = TRUE),
            sdh2 = sd(h2, na.rm = TRUE),
            mhdbh2 = mean(hdbh2, na.rm = TRUE),
            mqH2 = mean(Calidad2, na.rm = TRUE),
            cvdbh2 = sddbh2 / mdbh2,
            cvh2 = sdh2 / mh2) |>
  dplyr::select(!c(sddbh2, sdh2))

Psy3 <- Psy23_planted |>
  group_by(Plotcode3) |>
  summarise(type3 = first(Psy),
            IDPC3 = first(IDPC3),
            ba_ha3 = sum(AB3m2ha, na.rm = TRUE),
            dens3 = sum(dens3, na.rm = TRUE),
            mdbh3 = mean(dbh3, na.rm = TRUE),
            sddbh3 = sd(dbh3, na.rm = TRUE),
            mh3 = mean(h3, na.rm = TRUE),
            sdh3 = sd(h3, na.rm = TRUE),
            mhdbh3 = mean(hdbh3, na.rm = TRUE),
            mqH3 = mean(Calidad3, na.rm = TRUE),
            cvdbh3 = sddbh3 / mdbh3,
            cvh3 = sdh3 / mh3) |>
  dplyr::select(!c(sddbh3, sdh3))

map(Psy2, ~sum(is.na(.)))
```

    ## $Plotcode2
    ## [1] 1
    ## 
    ## $type2
    ## [1] 0
    ## 
    ## $IDPC3
    ## [1] 0
    ## 
    ## $ba_ha2
    ## [1] 0
    ## 
    ## $dens2
    ## [1] 0
    ## 
    ## $mdbh2
    ## [1] 1
    ## 
    ## $mh2
    ## [1] 1
    ## 
    ## $mhdbh2
    ## [1] 1
    ## 
    ## $mqH2
    ## [1] 1
    ## 
    ## $cvdbh2
    ## [1] 129
    ## 
    ## $cvh2
    ## [1] 129

``` r
map(Psy3, ~sum(is.na(.)))
```

    ## $Plotcode3
    ## [1] 0
    ## 
    ## $type3
    ## [1] 0
    ## 
    ## $IDPC3
    ## [1] 0
    ## 
    ## $ba_ha3
    ## [1] 0
    ## 
    ## $dens3
    ## [1] 0
    ## 
    ## $mdbh3
    ## [1] 57
    ## 
    ## $mh3
    ## [1] 57
    ## 
    ## $mhdbh3
    ## [1] 57
    ## 
    ## $mqH3
    ## [1] 38
    ## 
    ## $cvdbh3
    ## [1] 77
    ## 
    ## $cvh3
    ## [1] 77

``` r
Psy2 <- Psy2 |> 
  mutate(
    type2 = recode(type2, Planted = "1", Natural = "0")
  )
summary(Psy2$type2)
```

    ##    1    0 
    ## 1441 3926

``` r
Psy3 <- Psy3 |> 
  mutate(
    type3 = recode(type3, Planted = "1", Natural = "0")
  )
summary(Psy3$type3)
```

    ##    1    0 
    ## 1440 3926

``` r
Psy2 <- Psy2 |>
  left_join(sppriqueza2, by = "Plotcode2")
glimpse(Psy2)
```

    ## Rows: 5,367
    ## Columns: 12
    ## $ Plotcode2 <dbl> 10067, 10094, 10103, 10116, 10123, 10126, 10135, 10153, 1015…
    ## $ type2     <fct> 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, …
    ## $ IDPC3     <chr> "10067A1", "10094A1", "10103A1", "10116A1", "10123A1", "1012…
    ## $ ba_ha2    <dbl> 4.409798, 14.924834, 22.244708, 13.458388, 19.272495, 21.685…
    ## $ dens2     <dbl> 74.27231, 247.57436, 640.15655, 119.40158, 208.24540, 615.39…
    ## $ mdbh2     <dbl> 290.7500, 284.6000, 223.2759, 401.5909, 350.0938, 239.0652, …
    ## $ mh2       <dbl> 10.125000, 13.700000, 15.034483, 12.454545, 19.906250, 12.10…
    ## $ mhdbh2    <dbl> 36.26612, 49.96737, 70.87539, 31.84902, 57.46009, 53.04158, …
    ## $ mqH2      <dbl> 2.000000, 2.000000, 2.068966, 2.090909, 2.125000, 2.043478, …
    ## $ cvdbh2    <dbl> 0.31343216, 0.21566813, 0.28037706, 0.28235858, 0.16789753, …
    ## $ cvh2      <dbl> 0.28043992, 0.13780523, 0.16071099, 0.21201835, 0.18371327, …
    ## $ Riqueza2  <int> 1, 2, 2, 3, 3, 4, 2, 1, 2, 2, 4, 2, 3, 2, 3, 2, 2, 5, 3, 2, …

``` r
Psy3 <- Psy3 |>
  left_join(sppriqueza3, by = "Plotcode3")
glimpse(Psy3)
```

    ## Rows: 5,366
    ## Columns: 12
    ## $ Plotcode3 <dbl> 10067, 10094, 10103, 10116, 10123, 10126, 10135, 10153, 1015…
    ## $ type3     <fct> 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, …
    ## $ IDPC3     <chr> "10067A1", "10094A1", "10103A1", "10116A1", "10123A1", "1012…
    ## $ ba_ha3    <dbl> 6.123150, 19.390306, 15.466808, 15.311306, 1.903539, 23.6540…
    ## $ dens3     <dbl> 84.88264, 263.27764, 477.46483, 115.44039, 28.29421, 436.579…
    ## $ mdbh3     <dbl> 299.5000, 313.4167, 213.0000, 438.5833, 289.5000, 277.5208, …
    ## $ mh3       <dbl> 12.333333, 15.250000, 14.860000, 17.333333, 13.750000, 15.41…
    ## $ mhdbh3    <dbl> 42.10745, 50.36489, 71.35980, 40.76423, 48.17368, 57.89454, …
    ## $ mqH3      <dbl> 2.000000, 3.368421, 2.100000, 2.000000, 2.000000, 2.720000, …
    ## $ cvdbh3    <dbl> 0.16947663, 0.21478268, 0.19229229, 0.26986840, 0.21005590, …
    ## $ cvh3      <dbl> 0.15943655, 0.12280609, 0.15845923, 0.17146687, 0.07713892, …
    ## $ Riqueza3  <int> 1, 2, 2, 3, 3, 5, 2, 1, 2, 2, 4, 4, 3, 2, 3, 2, 2, 5, 3, 2, …

``` r
# all species and size classes
pPme2 <- Pme2 |> 
  group_by(Plotcode2) |> 
  summarise(Pmedens2 = sum(Pmedens_ha2, na.rm = T)) 

pPme3 <- Pme3 |>
  group_by(Plotcode3) |> 
  summarise(Pmedens3 = sum(Pmedens_ha3, na.rm = T)) 

Psy2 <- Psy2 |> 
  left_join(pPme2, by = "Plotcode2")

Psy3 <- Psy3 |>
  left_join(pPme3, by = "Plotcode3")

summary(Psy2$Pmedens2)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    ##   127.3  1655.2  3437.7  3878.5  5602.3 23682.3     422

``` r
summary(Psy3$Pmedens3)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    ##   127.3  2037.2  5474.9  6171.0  9231.0 35523.4     591

``` r
Psy2 <- Psy2 |> 
  mutate(Pmedens2 = replace_na(Pmedens2, 0))

Psy3 <- Psy3 |> 
  mutate(Pmedens3 = replace_na(Pmedens3, 0))

summary(Psy2$Pmedens2)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##     0.0   763.9  3055.8  3573.6  5347.6 23682.3

``` r
summary(Psy3$Pmedens3)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##     0.0   891.3  4583.7  5491.4  8658.0 35523.4

``` r
Psy23_final <- full_join(Psy2, Psy3, by = "IDPC3")

names(Psy23_final)
```

    ##  [1] "Plotcode2" "type2"     "IDPC3"     "ba_ha2"    "dens2"     "mdbh2"    
    ##  [7] "mh2"       "mhdbh2"    "mqH2"      "cvdbh2"    "cvh2"      "Riqueza2" 
    ## [13] "Pmedens2"  "Plotcode3" "type3"     "ba_ha3"    "dens3"     "mdbh3"    
    ## [19] "mh3"       "mhdbh3"    "mqH3"      "cvdbh3"    "cvh3"      "Riqueza3" 
    ## [25] "Pmedens3"

``` r
nrow(Psy23_final)
```

    ## [1] 5367

``` r
test_that("sum is number of NA", {
  expect_equal(sum(is.na(Psy23_final$Plotcode2)), 1)
  })
```

    ## Test passed 🎉

``` r
test_that("type2 == type3", {
  expect_equal(identical(Psy23_final[["type2"]],
                         Psy23_final[["type3"]]), TRUE)
})
```

    ## Test passed 😀

``` r
Psy23_final <- Psy23_final |> 
  drop_na(Plotcode2)

test_that("Plotcode2 == Plotcode3", {
  expect_equal(identical(Psy23_final[["Plotcode2"]],
                         Psy23_final[["Plotcode3"]]), TRUE)
})
```

    ## Test passed 🥇

``` r
map(Psy23_final, ~sum(is.na(.)))
```

    ## $Plotcode2
    ## [1] 0
    ## 
    ## $type2
    ## [1] 0
    ## 
    ## $IDPC3
    ## [1] 0
    ## 
    ## $ba_ha2
    ## [1] 0
    ## 
    ## $dens2
    ## [1] 0
    ## 
    ## $mdbh2
    ## [1] 0
    ## 
    ## $mh2
    ## [1] 0
    ## 
    ## $mhdbh2
    ## [1] 0
    ## 
    ## $mqH2
    ## [1] 0
    ## 
    ## $cvdbh2
    ## [1] 128
    ## 
    ## $cvh2
    ## [1] 128
    ## 
    ## $Riqueza2
    ## [1] 0
    ## 
    ## $Pmedens2
    ## [1] 0
    ## 
    ## $Plotcode3
    ## [1] 0
    ## 
    ## $type3
    ## [1] 0
    ## 
    ## $ba_ha3
    ## [1] 0
    ## 
    ## $dens3
    ## [1] 0
    ## 
    ## $mdbh3
    ## [1] 57
    ## 
    ## $mh3
    ## [1] 57
    ## 
    ## $mhdbh3
    ## [1] 57
    ## 
    ## $mqH3
    ## [1] 38
    ## 
    ## $cvdbh3
    ## [1] 77
    ## 
    ## $cvh3
    ## [1] 77
    ## 
    ## $Riqueza3
    ## [1] 43
    ## 
    ## $Pmedens3
    ## [1] 0

``` r
nrow(Psy23_final)
```

    ## [1] 5366

``` r
Psy23_final <- Psy23_final |>
  drop_na()
nrow(Psy23_final)
```

    ## [1] 5163

``` r
table(Psy23_final$type2)
```

    ## 
    ##    1    0 
    ## 1376 3787

``` r
# long format
Psy23_final_ba_long <- Psy23_final |>
  filter(mqH3 <= 6 ) |> 
  pivot_longer(
    cols = c(
      "ba_ha2", "ba_ha3"
    ),
    names_to = "ba_ifn",
    values_to = "basal_area"
  ) |> 
  dplyr::select(
    IDPC3, ba_ifn, basal_area, type2
    )
Psy23_final_dens_long <- Psy23_final |>
  filter(mqH3 <= 6 ) |> 
  pivot_longer(
    cols = c(
      "dens2", "dens3"
    ),
    names_to = "dens_ifn",
    values_to = "density"
  ) |> 
  dplyr::select(
    dens_ifn, density
    )
Psy23_final_mdbh_long <- Psy23_final |>
  filter(mqH3 <= 6 ) |> 
  pivot_longer(
    cols = c(
      "mdbh2", "mdbh3"
    ),
    names_to = "mdbh_ifn",
    values_to = "dbh"
  ) |> 
  dplyr::select(
    mdbh_ifn, dbh
    )
Psy23_final_mh_long <- Psy23_final |>
  filter(mqH3 <= 6 ) |> 
  pivot_longer(
    cols = c(
      "mh2", "mh3"
    ),
    names_to = "mh_ifn",
    values_to = "height"
  ) |> 
  dplyr::select(
    mh_ifn, height
  )
Psy23_final_mhdbh_long <- Psy23_final |>
  filter(mqH3 <= 6 ) |> 
  pivot_longer(
    cols = c(
      "mhdbh2", "mhdbh3"
    ),
    names_to = "mhdbh_ifn",
    values_to = "height_dbh"
  ) |> 
  dplyr::select(
    mhdbh_ifn, height_dbh
  )
Psy23_final_cvdbh_long <- Psy23_final |>
  filter(mqH3 <= 6 ) |> 
  pivot_longer(
    cols = c(
      "cvdbh2", "cvdbh3"
    ),
    names_to = "cvdbh_ifn",
    values_to = "cv_dbh"
  ) |> 
  dplyr::select(
    cvdbh_ifn, cv_dbh
  )
Psy23_final_cvh_long <- Psy23_final |>
  filter(mqH3 <= 6 ) |> 
  pivot_longer(
    cols = c(
      "cvh2", "cvh3"
    ),
    names_to = "cvh_ifn",
    values_to = "cv_height"
  ) |> 
  dplyr::select(
    cvh_ifn, cv_height
  )
Psy23_final_Pmedens_long <- Psy23_final |>
  filter(mqH3 <= 6 ) |> 
  pivot_longer(
    cols = c(
      "Pmedens2", "Pmedens3"
    ),
    names_to = "Pmedens_ifn",
    values_to = "saplings"
  ) |> 
  dplyr::select(
    Pmedens_ifn, saplings
  )
Psy23_final_Riqueza_long <- Psy23_final |>
  filter(mqH3 <= 6 ) |> 
  pivot_longer(
    cols = c(
      "Riqueza2", "Riqueza3"
    ),
    names_to = "Riqueza_ifn",
    values_to = "sp_richness"
  ) |> 
  dplyr::select(
    Riqueza_ifn, sp_richness
  )
Psy23_final_mqH_long <- Psy23_final |>
  filter(mqH3 <= 6 ) |> 
  pivot_longer(
    cols = c(
      "mqH2", "mqH3"
    ),
    names_to = "mqH_ifn",
    values_to = "tree_vigour"
  ) |> 
  dplyr::select(
    mqH_ifn, tree_vigour
  )

Psy23_final_long <- bind_cols(
  Psy23_final_ba_long,
  Psy23_final_dens_long,
  Psy23_final_mdbh_long,
  Psy23_final_mh_long,
  Psy23_final_mhdbh_long,
  Psy23_final_cvdbh_long,
  Psy23_final_cvh_long,
  Psy23_final_Pmedens_long,
  Psy23_final_Riqueza_long,
  Psy23_final_mqH_long
)

# for boxplot cluster x origin
Psy23_final_long_str <- Psy23_final_long

test_that("Plotcode2 == Plotcode3", {
  expect_equal(unname(map_dbl(Psy23_final_long, ~sum(is.na(.)))), rep(0, length(names(Psy23_final_long))))
               })
```

    ## Test passed 😀

# PCA and clustering: only planted forests

``` r
names(Psy23_final_long)
```

    ##  [1] "IDPC3"       "ba_ifn"      "basal_area"  "type2"       "dens_ifn"   
    ##  [6] "density"     "mdbh_ifn"    "dbh"         "mh_ifn"      "height"     
    ## [11] "mhdbh_ifn"   "height_dbh"  "cvdbh_ifn"   "cv_dbh"      "cvh_ifn"    
    ## [16] "cv_height"   "Pmedens_ifn" "saplings"    "Riqueza_ifn" "sp_richness"
    ## [21] "mqH_ifn"     "tree_vigour"

``` r
# type = 1 is planted
Psy23_final_long <- Psy23_final_long |> 
  mutate(
    type2 = str_replace(type2, "1", "planted"),
    type2 = str_replace(type2, "0", "natural"),
  )

# ja: clusters para plantado y natural por serapado
table(Psy23_final_long$type2, useNA = "always")
```

    ## 
    ## natural planted    <NA> 
    ##    7574    2750       0

``` r
nrow(Psy23_final_long |> 
       distinct(IDPC3))
```

    ## [1] 5162

``` r
Psy23_final_long_num <- Psy23_final_long |>
  dplyr::select(
    !c(ends_with("ifn"), IDPC3, type2)
    )
glimpse(Psy23_final_long_num)
```

    ## Rows: 10,324
    ## Columns: 10
    ## $ basal_area  <dbl> 4.409798, 6.123150, 14.924834, 19.390306, 22.244708, 15.46…
    ## $ density     <dbl> 74.27231, 84.88264, 247.57436, 263.27764, 640.15655, 477.4…
    ## $ dbh         <dbl> 290.7500, 299.5000, 284.6000, 313.4167, 223.2759, 213.0000…
    ## $ height      <dbl> 10.125000, 12.333333, 13.700000, 15.250000, 15.034483, 14.…
    ## $ height_dbh  <dbl> 36.26612, 42.10745, 49.96737, 50.36489, 70.87539, 71.35980…
    ## $ cv_dbh      <dbl> 0.31343216, 0.16947663, 0.21566813, 0.21478268, 0.28037706…
    ## $ cv_height   <dbl> 0.28043992, 0.15943655, 0.13780523, 0.12280609, 0.16071099…
    ## $ saplings    <dbl> 2291.8312, 6238.8738, 2291.8312, 1273.2395, 4710.9863, 241…
    ## $ sp_richness <int> 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 5, 2, 2, 1, 1, 2, 2, 2, 2…
    ## $ tree_vigour <dbl> 2.000000, 2.000000, 2.000000, 3.368421, 2.068966, 2.100000…

``` r
pca_Psy <- prcomp(Psy23_final_long_num, scale = TRUE)

PCAPsy <- fviz_pca_biplot(pca_Psy, label = "var", 
                          col.ind = "grey",
                          col.var = "black",
                          title = "",
                          labelsize = 3,
                          habillage = Psy23_final_long$type2)

# https://uc-r.github.io/kmeans_clustering
Psy23_final_long_num_sc <- scale(Psy23_final_long_num)

set.seed(123)
k3 <- kmeans(Psy23_final_long_num_sc, centers = 2, nstart = 25)

myColors <- brewer.pal(4, "Paired")[c(2, 4)]

cluster_ifn23 <- fviz_cluster(
  k3, data = Psy23_final_long_num_sc,
  geom = "point"
  ) + 
  ggtitle(label = '') +
  scale_color_manual(values = myColors) +
  scale_fill_manual(values = myColors)
cluster_ifn23
```

![](01-clustering_files/figure-gfm/cluster-1.png)<!-- -->

``` r
# # determining number of clusters
# set.seed(123)
# fviz_nbclust(Psy23_final_long_num_sc, kmeans, method = "wss")
# 
# # quality of the cluster
# set.seed(123)
# fviz_nbclust(Psy23_final_long_num_sc, kmeans, k.max = 6,
#              method = "silhouette") +
#   ggtitle(label = "")

# add the cluster number to the original data
Psy23_final_long <- Psy23_final_long |>
  mutate(
    k3 = k3$cluster
  )
table(Psy23_final_long$k3)
```

    ## 
    ##    1    2 
    ## 6365 3959

``` r
names(Psy23_final_long)
```

    ##  [1] "IDPC3"       "ba_ifn"      "basal_area"  "type2"       "dens_ifn"   
    ##  [6] "density"     "mdbh_ifn"    "dbh"         "mh_ifn"      "height"     
    ## [11] "mhdbh_ifn"   "height_dbh"  "cvdbh_ifn"   "cv_dbh"      "cvh_ifn"    
    ## [16] "cv_height"   "Pmedens_ifn" "saplings"    "Riqueza_ifn" "sp_richness"
    ## [21] "mqH_ifn"     "tree_vigour" "k3"

``` r
# for the contingency table
Psy23_final_long_k <- Psy23_final_long

allPsy23 <- Psy23_final_long |>
   pivot_longer(
    cols = c(
      basal_area,
      density, dbh, height,
      height_dbh, cv_dbh, cv_height,
      saplings, sp_richness,
      tree_vigour
    ),
    names_to = "structure",
    values_to = "str_value"
  ) |> 
  dplyr::select(
    !ends_with("ifn")
  )

gg_str <- function(data, k){
  ggplot(data, aes(x = structure, y = str_value,
                   fill = as.character(.data[[k]]))) +
    geom_boxplot(outlier.size = 0.2, width = 0.5) +
    labs(x = "", title = "", y = "")  +
    scale_fill_manual(name = "cluster", values = myColors) +
    facet_wrap(~structure, scale = "free",
               strip.position = "bottom") +
    ylab("Structure") +
    theme(
      axis.text = element_text(colour = "grey50", size = 10),
      axis.title = element_text(colour = "grey50", size = 12),
      panel.grid.major = element_line(colour = "grey90", size = 0.5), 
      panel.background = element_blank(),
      legend.key = element_blank(),
      strip.text.x = element_blank(),
      axis.ticks = element_blank()
    )
}

plot_str <- gg_str(data = allPsy23,
                   k = "k3")
plot_str
```

![](01-clustering_files/figure-gfm/cluster-2.png)<!-- -->

``` r
plot_str_plan_nat <- gg_str(data = allPsy23,
                   k = "type2") +
  scale_color_manual(values = c("#F8766D", "#619CFF")) +
  scale_fill_manual(values = c("#F8766D", "#619CFF"))
plot_str_plan_nat
```

![](01-clustering_files/figure-gfm/cluster-3.png)<!-- -->

``` r
ggsave( 
  plot = plot_str_plan_nat,
  here("03-results", "others", "ifn23_str_plan_nat.png"),
  width = 10, height = 8
)

plot_PCA <- PCAPsy + theme(
  legend.title = element_blank(),
  text = element_text(colour = "grey50", size = 12),
  panel.grid.minor = element_blank(),
  axis.text = element_text(colour = "grey50", size = 10),
  axis.title = element_text(colour = "grey50", size = 12),
  axis.ticks = element_blank()
)
plot_PCA
```

![](01-clustering_files/figure-gfm/cluster-4.png)<!-- -->

``` r
plot_cluster <- cluster_ifn23 +
  theme(
    axis.text = element_text(colour = "grey50", size = 10),
    axis.title = element_text(colour = "grey50", size = 12),
    panel.grid.major = element_line(colour = "grey90", size = 0.5), 
    panel.background = element_blank(),
    axis.ticks = element_blank()
  )
plot_cluster
```

![](01-clustering_files/figure-gfm/cluster-5.png)<!-- -->

# Demography (2nd & 3rd forest inventory (ifn23)

``` r
names(Psy23_planted)
```

    ##  [1] "ID_Pma3"           "ID_Pma3c"          "ID_Pma2"          
    ##  [4] "IDPC3"             "IDPCc3"            "Plotcode3"        
    ##  [7] "Plotcode2"         "Provincia3"        "Cla3"             
    ## [10] "Subclase3"         "Orden33"           "Orden23"          
    ## [13] "rumbo32"           "distancia32"       "state3"           
    ## [16] "R2"                "R3"                "sppcompa"         
    ## [19] "Especie2"          "Nombre2"           "genero2"          
    ## [22] "Especie3"          "Nombre3"           "genero3"          
    ## [25] "sp_err23"          "sp_mixsp23"        "sp_spmix23"       
    ## [28] "cod_change23"      "intragen_23"       "intergen_23"      
    ## [31] "dbh2"              "dbh3"              "h2"               
    ## [34] "h3"                "dens2"             "dens3"            
    ## [37] "AB2"               "AB3"               "AB2m2ha"          
    ## [40] "AB3m2ha"           "AB32"              "AB32ha"           
    ## [43] "AB2m2ha_muertosc"  "AB2m2ha_presente"  "AB2m2ha_ausente"  
    ## [46] "reclutamiento32"   "reclutamiento32_5" "Agente3"          
    ## [49] "Importancia3"      "Elemento3"         "Calidad3"         
    ## [52] "Calidad2"          "Forma3"            "Forma2"           
    ## [55] "Paramesp3"         "Paramesp2"         "VCCha3"           
    ## [58] "VCCha2"            "VSCha3"            "VSCha2"           
    ## [61] "IAVCha3"           "IAVCha2"           "VLEha3"           
    ## [64] "VLEha2"            "Error23"           "Error3"           
    ## [67] "Error2"            "AB2m2ha.tot.plot"  "AB2m2ha.plot"     
    ## [70] "AB2m2ha.rel.plot"  "ntrees"            "Psy"              
    ## [73] "hdbh2"             "hdbh3"

``` r
# for the dbh of dead treees
Psy23_planted_m <- Psy23_planted

# I obtain the initial basal area and
# the missing dead trees for each plot
dem23_ifn2 <- Psy23_planted |>  
  group_by(IDPC3) |>
  summarise(
    ba_ifn2 = sum(AB2m2ha, na.rm = TRUE),
    Mortality_ausente = sum(AB2m2ha_ausente, na.rm = TRUE)
  ) |> 
  select(IDPC3, ba_ifn2, Mortality_ausente)

table(Psy23_planted$state3)
```

    ## 
    ##    MA    MP MP-nc     R     V  V-nc 
    ## 10924  3493    72 35411 92069   319

``` r
Psy23_planted <- Psy23_planted |> 
  mutate(AB32ha_all = AB3m2ha - AB2m2ha) |> 
  # delete no comparable & missing dead trees 2-3
  filter(!c(
    state3 == "MA-nc" | state3 == "MP-nc" | state3 == "V-nc" 
    | state3 == "MA"
    ))

table(Psy23_planted$state3)
```

    ## 
    ##    MP     R     V 
    ##  3493 35411 92069

``` r
# reclutamiento32_5 recruited trees in
# the subplot of 5 m
dem23 <- Psy23_planted |>  
  group_by(IDPC3) |>
  summarise(
    type2 = first(Psy),
    BAc = sum(AB32ha_all, na.rm = TRUE),
    Ingrowth = sum(reclutamiento32_5, na.rm = TRUE),
    Growth = sum(AB32ha, na.rm = TRUE),
    Mortality = sum(AB2m2ha_muertosc, na.rm = TRUE)
    ) |>
  left_join(dem23_ifn2, by = "IDPC3") |> 
  mutate(
    BAc = BAc / ba_ifn2,
    Ingrowth = Ingrowth / ba_ifn2,
    Growth = Growth / ba_ifn2,
    Mortality_nd = Mortality,
    Mortality = Mortality / ba_ifn2,
    Mortality_ausente = Mortality_ausente / ba_ifn2,
  )

# to get census interval
plot23 <- read_csv(
  here("01-data", "ifn", "plot23.csv"),
  col_types = cols(
    MEROSIVA2 = col_character(),
    MEROSIVA3 = col_character(),  
    FCCTOT = col_character(),
    MEROSIVA4 = col_character(),  
    MejVue13 = col_character(),
    CortaReg3 = col_character(),  
    MejVue23 = col_character(),
    Tratasuelo3 = col_character(),
    MejSue13 = col_character(),
    MejSue23 = col_character(),   
    Nivel33 = col_character(),
    Nivel23 = col_character()
  )
)

plot_dem23 <- plot23 |> 
  dplyr::select(IDPC3, Plotcode3, year32, Cut32) |> 
  right_join(dem23, by = "IDPC3") |> 
  mutate(
    BAc23 = BAc / year32,
    Ingrowth23 = Ingrowth / year32,
    Growth23 = Growth / year32,
    Mortality23 = Mortality / year32,
    Mortality23_nd = Mortality_nd / year32,
    Mortality_ausente23 = Mortality_ausente / year32,
  ) |> 
  dplyr::select(
    !c(BAc, Ingrowth, Growth, Mortality, Mortality_nd, Mortality_ausente)
  )

names(plot_dem23)
```

    ##  [1] "IDPC3"               "Plotcode3"           "year32"             
    ##  [4] "Cut32"               "type2"               "ba_ifn2"            
    ##  [7] "BAc23"               "Ingrowth23"          "Growth23"           
    ## [10] "Mortality23"         "Mortality23_nd"      "Mortality_ausente23"

``` r
Psy_dem23 <- plot_dem23

nrow(Psy_dem23)
```

    ## [1] 5328

``` r
# planted ifn 23
nrow(Psy23_final)
```

    ## [1] 5163

``` r
names(Psy23_final)
```

    ##  [1] "Plotcode2" "type2"     "IDPC3"     "ba_ha2"    "dens2"     "mdbh2"    
    ##  [7] "mh2"       "mhdbh2"    "mqH2"      "cvdbh2"    "cvh2"      "Riqueza2" 
    ## [13] "Pmedens2"  "Plotcode3" "type3"     "ba_ha3"    "dens3"     "mdbh3"    
    ## [19] "mh3"       "mhdbh3"    "mqH3"      "cvdbh3"    "cvh3"      "Riqueza3" 
    ## [25] "Pmedens3"

``` r
table(Psy_dem23$type2, useNA = "always")
```

    ## 
    ## Planted Natural    <NA> 
    ##    1437    3891       0

``` r
# keep the cluster of ifn2
Psy23_final_long <- Psy23_final_long |> 
  distinct(IDPC3, .keep_all = T)

# there are some plots with NA in the planted and natural cluster
# because they only contain one LIVE tree. I deleted them
# kk <- Psy23_planted |>
#   filter(
#   IDPC3 == "10913A1" | IDPC3 == "10634A1"
# )

# demography for natural and planted, all together
Psy_dem23 <- Psy_dem23 |>
  left_join(Psy23_final_long,
            by = "IDPC3")

nrow(Psy_dem23)
```

    ## [1] 5328

``` r
map(Psy_dem23, ~sum(is.na(.)))
```

    ## $IDPC3
    ## [1] 0
    ## 
    ## $Plotcode3
    ## [1] 0
    ## 
    ## $year32
    ## [1] 0
    ## 
    ## $Cut32
    ## [1] 0
    ## 
    ## $type2.x
    ## [1] 0
    ## 
    ## $ba_ifn2
    ## [1] 0
    ## 
    ## $BAc23
    ## [1] 0
    ## 
    ## $Ingrowth23
    ## [1] 0
    ## 
    ## $Growth23
    ## [1] 0
    ## 
    ## $Mortality23
    ## [1] 0
    ## 
    ## $Mortality23_nd
    ## [1] 0
    ## 
    ## $Mortality_ausente23
    ## [1] 0
    ## 
    ## $ba_ifn
    ## [1] 166
    ## 
    ## $basal_area
    ## [1] 166
    ## 
    ## $type2.y
    ## [1] 166
    ## 
    ## $dens_ifn
    ## [1] 166
    ## 
    ## $density
    ## [1] 166
    ## 
    ## $mdbh_ifn
    ## [1] 166
    ## 
    ## $dbh
    ## [1] 166
    ## 
    ## $mh_ifn
    ## [1] 166
    ## 
    ## $height
    ## [1] 166
    ## 
    ## $mhdbh_ifn
    ## [1] 166
    ## 
    ## $height_dbh
    ## [1] 166
    ## 
    ## $cvdbh_ifn
    ## [1] 166
    ## 
    ## $cv_dbh
    ## [1] 166
    ## 
    ## $cvh_ifn
    ## [1] 166
    ## 
    ## $cv_height
    ## [1] 166
    ## 
    ## $Pmedens_ifn
    ## [1] 166
    ## 
    ## $saplings
    ## [1] 166
    ## 
    ## $Riqueza_ifn
    ## [1] 166
    ## 
    ## $sp_richness
    ## [1] 166
    ## 
    ## $mqH_ifn
    ## [1] 166
    ## 
    ## $tree_vigour
    ## [1] 166
    ## 
    ## $k3
    ## [1] 166

``` r
names(Psy_dem23)
```

    ##  [1] "IDPC3"               "Plotcode3"           "year32"             
    ##  [4] "Cut32"               "type2.x"             "ba_ifn2"            
    ##  [7] "BAc23"               "Ingrowth23"          "Growth23"           
    ## [10] "Mortality23"         "Mortality23_nd"      "Mortality_ausente23"
    ## [13] "ba_ifn"              "basal_area"          "type2.y"            
    ## [16] "dens_ifn"            "density"             "mdbh_ifn"           
    ## [19] "dbh"                 "mh_ifn"              "height"             
    ## [22] "mhdbh_ifn"           "height_dbh"          "cvdbh_ifn"          
    ## [25] "cv_dbh"              "cvh_ifn"             "cv_height"          
    ## [28] "Pmedens_ifn"         "saplings"            "Riqueza_ifn"        
    ## [31] "sp_richness"         "mqH_ifn"             "tree_vigour"        
    ## [34] "k3"

``` r
Psy_dem23_sel <- Psy_dem23 |>
  dplyr::select(
    IDPC3,
    Ingrowth23,
    Growth23, 
    Mortality23, 
    Mortality23_nd,
    Mortality_ausente23,
    BAc23,
    ba_ifn2,
    k3,
    type2.x,
    Cut32
  )

Psy_dem23_sel <- Psy_dem23_sel |>
  filter(!is.na(k3))

Psy_dem23_sel_lon <- Psy_dem23_sel |>
  pivot_longer(
    cols = c(
      Ingrowth23, 
      Growth23, 
      Mortality23,
      BAc23
    ),
    names_to = "demography",
    values_to = "dem_value"
  )

myColors <- brewer.pal(4, "Paired")[c(2, 4)]

col_strips <- strip_themed(
  # http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
  background_x = elem_list_rect(fill = c("#F8766D", "#619CFF")),
  by_layer_x = FALSE,
  text_x = elem_list_text(face = "bold")
)

gg_dem <- function(k){
  ggplot(Psy_dem23_sel_lon, aes(x = demography, y = dem_value,
                                fill = as.character(.data[[k]]))) +
    geom_boxplot(outlier.shape = NA) +
    labs(x = "", title = "", y = "Demography")  +
    scale_color_manual(values = myColors) +
    scale_fill_manual(values = myColors) +
    facet_wrap2(
      ~type2.x, strip = col_strips
      ) +
    coord_cartesian(ylim = c(0.0, 0.25)) +
    theme(
      legend.title = element_blank(),
      axis.title.x = element_blank(),
      axis.text = element_text(colour = "grey50", size = 10),
      axis.title = element_text(colour = "grey50", size = 12),
      panel.grid.major = element_line(colour = "grey90", size = 0.5), 
      panel.background = element_blank(),
      legend.key = element_blank()
    )
}

plot_dem <- gg_dem(k = "k3")
plot_dem
```

![](01-clustering_files/figure-gfm/demography-1.png)<!-- -->

``` r
ggsave( 
  plot = plot_dem,
  here("03-results", "others", "ifn23_demo.png"),
  width = 10, height = 8
)

summary(Psy_dem23_sel_lon |> 
  filter(type2.x == "Planted" & demography == "Ingrowth23") |> 
    select(dem_value))
```

    ##    dem_value       
    ##  Min.   :0.000000  
    ##  1st Qu.:0.000000  
    ##  Median :0.000000  
    ##  Mean   :0.012748  
    ##  3rd Qu.:0.004254  
    ##  Max.   :0.859408

``` r
summary(Psy_dem23_sel_lon |> 
  filter(type2.x == "Natural" & demography == "Ingrowth23") |> 
    select(dem_value))
```

    ##    dem_value       
    ##  Min.   :0.000000  
    ##  1st Qu.:0.000000  
    ##  Median :0.000000  
    ##  Mean   :0.005991  
    ##  3rd Qu.:0.003577  
    ##  Max.   :0.446939

# Legacies map

``` r
cPsy <- left_join(Psy_dem23_sel, plot23 |> 
                    select(IDPC3, Plotcode3, CX_teo, CY_teo), by = "IDPC3")

table(cPsy$k3, useNA = "always")
```

    ## 
    ##    1    2 <NA> 
    ## 3430 1732    0

``` r
map_loc <- map_data("world", region = c("Portugal", "Spain", "France",
                                        "Andorra")) |> 
  ggplot(aes(x = long, y = lat)) + 
  geom_polygon(aes(group = group), 
               fill = "white", 
               color = "grey") +
  coord_fixed(xlim = c(-11, 5),  
              ylim = c(35, 45), 
              ratio = 1.3) +
  xlab("Longitude") +
  ylab("Latitude") + 
  theme(legend.title = element_text(colour = "black", face = "bold", size = 12),
        legend.position = c(0.9, 0.2),
        legend.text = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),  
        axis.text.x = element_text(colour = "grey50", size = 8),
        axis.text.y = element_text(colour = "grey50", size = 8),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(colour = "grey90", size = 0.5), 
        axis.line = element_blank(),
        panel.background = element_blank()) 

cPsy_sf <- st_as_sf(x = cPsy, 
                    coords = c("CX_teo", "CY_teo"),
                    crs = "+proj=utm +zone=30 +ellps=intl +units=m +no_defs")
st_crs(cPsy_sf)
```

    ## Coordinate Reference System:
    ##   User input: +proj=utm +zone=30 +ellps=intl +units=m +no_defs 
    ##   wkt:
    ## PROJCRS["unknown",
    ##     BASEGEOGCRS["unknown",
    ##         DATUM["Unknown based on International 1909 (Hayford) ellipsoid",
    ##             ELLIPSOID["International 1909 (Hayford)",6378388,297,
    ##                 LENGTHUNIT["metre",1,
    ##                     ID["EPSG",9001]]]],
    ##         PRIMEM["Greenwich",0,
    ##             ANGLEUNIT["degree",0.0174532925199433],
    ##             ID["EPSG",8901]]],
    ##     CONVERSION["UTM zone 30N",
    ##         METHOD["Transverse Mercator",
    ##             ID["EPSG",9807]],
    ##         PARAMETER["Latitude of natural origin",0,
    ##             ANGLEUNIT["degree",0.0174532925199433],
    ##             ID["EPSG",8801]],
    ##         PARAMETER["Longitude of natural origin",-3,
    ##             ANGLEUNIT["degree",0.0174532925199433],
    ##             ID["EPSG",8802]],
    ##         PARAMETER["Scale factor at natural origin",0.9996,
    ##             SCALEUNIT["unity",1],
    ##             ID["EPSG",8805]],
    ##         PARAMETER["False easting",500000,
    ##             LENGTHUNIT["metre",1],
    ##             ID["EPSG",8806]],
    ##         PARAMETER["False northing",0,
    ##             LENGTHUNIT["metre",1],
    ##             ID["EPSG",8807]],
    ##         ID["EPSG",16030]],
    ##     CS[Cartesian,2],
    ##         AXIS["(E)",east,
    ##             ORDER[1],
    ##             LENGTHUNIT["metre",1,
    ##                 ID["EPSG",9001]]],
    ##         AXIS["(N)",north,
    ##             ORDER[2],
    ##             LENGTHUNIT["metre",1,
    ##                 ID["EPSG",9001]]]]

``` r
ggplot() + 
  geom_sf(data = cPsy_sf)
```

![](01-clustering_files/figure-gfm/map-legacies-1.png)<!-- -->

``` r
cPsy_sf_map <- st_transform(
  cPsy_sf,
  "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
)

st_crs(cPsy_sf_map)
```

    ## Coordinate Reference System:
    ##   User input: +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs 
    ##   wkt:
    ## GEOGCRS["unknown",
    ##     DATUM["World Geodetic System 1984",
    ##         ELLIPSOID["WGS 84",6378137,298.257223563,
    ##             LENGTHUNIT["metre",1]],
    ##         ID["EPSG",6326]],
    ##     PRIMEM["Greenwich",0,
    ##         ANGLEUNIT["degree",0.0174532925199433],
    ##         ID["EPSG",8901]],
    ##     CS[ellipsoidal,2],
    ##         AXIS["longitude",east,
    ##             ORDER[1],
    ##             ANGLEUNIT["degree",0.0174532925199433,
    ##                 ID["EPSG",9122]]],
    ##         AXIS["latitude",north,
    ##             ORDER[2],
    ##             ANGLEUNIT["degree",0.0174532925199433,
    ##                 ID["EPSG",9122]]]]

``` r
cPsy_sp_tb <- cPsy_sf_map |>
  ungroup() |> 
  mutate(CX = unlist(map(cPsy_sf_map$geometry, 1)),
         CY = unlist(map(cPsy_sf_map$geometry, 2))) |> 
  as_tibble() |> 
  dplyr::select(!geometry)

map_plots <- function(col.var){
  map_loc + 
    geom_point(data = cPsy_sp_tb,
               aes(CX, CY, color = as.factor(.data[[col.var]]),
                   shape = as.factor(.data[[col.var]])),
               size = .5) + 
    scale_color_manual(values = myColors, name = "cluster") +
    xlab("") +
    ylab("") + 
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    theme(plot.margin = margin(0, 0, 0, 0),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 9),
          legend.key = element_blank(),
          legend.justification = "left",
          legend.background = element_rect(fill = "white", colour = NA),
          axis.title = element_blank(),
          axis.text.x = element_text(colour = "grey50", size = 7),
          axis.text.y = element_text(colour = "grey50", size = 7),
          panel.grid.major = element_line(colour = "grey90", size = 0.5), 
          axis.line = element_blank(),
          panel.background = element_blank())
}

map_cluster <- map_plots(col.var = "k3") +
  scale_shape_manual(values = c(19, 17), name = "cluster")
map_cluster
```

![](01-clustering_files/figure-gfm/map-legacies-2.png)<!-- -->

``` r
map_planted <- map_plots(col.var = "type2.x") +
  scale_color_manual(
    limits = c("Natural", "Planted"),
    labels = c("natural", "planted"),
    values = c("#F8766D", "#00BFC4"),
    name = "cluster"
    ) +
  scale_shape_manual(
    limits = c("Natural", "Planted"),
    labels = c("natural", "planted"),
    values = c(19, 17),
    name = "cluster"
    )
map_planted
```

![](01-clustering_files/figure-gfm/map-legacies-3.png)<!-- -->

``` r
# save paper plots --------------------------------------------------------
plot_PCA + theme(legend.position = "none") + map_planted +
  plot_annotation(tag_levels = "a")
```

![](01-clustering_files/figure-gfm/map-legacies-4.png)<!-- -->

``` r
ggsave(
  here("03-results", "si_figures", "s1_psy_plan_nat.png"),
  width = 16, height = 8
)

(plot_cluster + plot_str) /
  (map_cluster + plot_spacer()) +
  plot_annotation(tag_levels = "a")
```

![](01-clustering_files/figure-gfm/map-legacies-5.png)<!-- -->

``` r
ggsave(
  here("03-results", "si_figures", "s2_psy_cluster.png"),
  width = 18, height = 10
)
```

# Merge with climate data

``` r
clima <- read_csv(here("01-data", "climate", "completeclimate.csv")) |> 
  dplyr::select(!c(...1, Cut, yearini, yearfin, IFNcode, plotcode))

# remove duplicated plots
clima <- clima |> 
  distinct(IDPC, .keep_all = T)
nrow(clima)
```

    ## [1] 58819

``` r
dem_clim <- left_join(
  Psy_dem23_sel, clima, by = c("IDPC3" = "IDPC")
  )
nrow(dem_clim)
```

    ## [1] 5162

``` r
map(dem_clim, ~sum(is.na(.)))
```

    ## $IDPC3
    ## [1] 0
    ## 
    ## $Ingrowth23
    ## [1] 0
    ## 
    ## $Growth23
    ## [1] 0
    ## 
    ## $Mortality23
    ## [1] 0
    ## 
    ## $Mortality23_nd
    ## [1] 0
    ## 
    ## $Mortality_ausente23
    ## [1] 0
    ## 
    ## $BAc23
    ## [1] 0
    ## 
    ## $ba_ifn2
    ## [1] 0
    ## 
    ## $k3
    ## [1] 0
    ## 
    ## $type2.x
    ## [1] 0
    ## 
    ## $Cut32
    ## [1] 0
    ## 
    ## $CX
    ## [1] 351
    ## 
    ## $CY
    ## [1] 351
    ## 
    ## $sf_nfi
    ## [1] 351
    ## 
    ## $sgdd_nfi
    ## [1] 351
    ## 
    ## $PPplot
    ## [1] 351
    ## 
    ## $PETplot
    ## [1] 351
    ## 
    ## $WAIplot
    ## [1] 351
    ## 
    ## $speimean
    ## [1] 351
    ## 
    ## $speimin
    ## [1] 351

``` r
table(dem_clim$k3)
```

    ## 
    ##    1    2 
    ## 3430 1732

``` r
dem_clim_final <- dem_clim |> 
  dplyr::select(
    IDPC3, CX, CY,
    k3, type2.x, Cut32, everything()
  )

glimpse(dem_clim_final)
```

    ## Rows: 5,162
    ## Columns: 20
    ## $ IDPC3               <chr> "100303A1", "10067A1", "100773A1", "100775A1", "10…
    ## $ CX                  <dbl> 262000, 490000, 263000, 263000, 262000, 258000, 49…
    ## $ CY                  <dbl> 4470000, 4774000, 4468000, 4467000, 4463000, 44590…
    ## $ k3                  <int> 1, 1, 1, 1, 2, 1, 1, 2, 1, 2, 1, 1, 2, 2, 1, 2, 1,…
    ## $ type2.x             <fct> Planted, Planted, Planted, Planted, Planted, Plant…
    ## $ Cut32               <dbl> 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1,…
    ## $ Ingrowth23          <dbl> 0.15401126, 0.00000000, 0.00368975, 0.00000000, 0.…
    ## $ Growth23            <dbl> 0.3532402214, 0.0218017605, 0.1356360226, 0.205318…
    ## $ Mortality23         <dbl> 0.000000000, 0.000000000, 0.008337227, 0.000000000…
    ## $ Mortality23_nd      <dbl> 0.00000000, 0.00000000, 0.12658182, 0.00000000, 0.…
    ## $ Mortality_ausente23 <dbl> 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.…
    ## $ BAc23               <dbl> 0.0201282372, 0.0023934322, 0.0093027502, -0.00958…
    ## $ ba_ifn2             <dbl> 1.381850, 4.409798, 15.182725, 4.547981, 26.897508…
    ## $ sf_nfi              <dbl> 14.5185200, 0.1111111, 14.3333300, 17.7407400, 15.…
    ## $ sgdd_nfi            <dbl> 2031.359, 3091.041, 2044.532, 1888.251, 1993.029, …
    ## $ PPplot              <dbl> 1050, 1028, 1041, 1112, 1067, 978, 1016, 1009, 100…
    ## $ PETplot             <dbl> 957, 798, 958, 937, 956, 994, 804, 807, 816, 814, …
    ## $ WAIplot             <dbl> 9.717868, 28.822055, 8.663883, 18.676628, 11.61087…
    ## $ speimean            <dbl> -0.05407972, -0.38432246, -0.05407972, -0.05407972…
    ## $ speimin             <dbl> -1.713884, -1.885054, -1.713884, -1.713884, -1.713…

``` r
map(dem_clim_final, ~sum(is.na(.)))
```

    ## $IDPC3
    ## [1] 0
    ## 
    ## $CX
    ## [1] 351
    ## 
    ## $CY
    ## [1] 351
    ## 
    ## $k3
    ## [1] 0
    ## 
    ## $type2.x
    ## [1] 0
    ## 
    ## $Cut32
    ## [1] 0
    ## 
    ## $Ingrowth23
    ## [1] 0
    ## 
    ## $Growth23
    ## [1] 0
    ## 
    ## $Mortality23
    ## [1] 0
    ## 
    ## $Mortality23_nd
    ## [1] 0
    ## 
    ## $Mortality_ausente23
    ## [1] 0
    ## 
    ## $BAc23
    ## [1] 0
    ## 
    ## $ba_ifn2
    ## [1] 0
    ## 
    ## $sf_nfi
    ## [1] 351
    ## 
    ## $sgdd_nfi
    ## [1] 351
    ## 
    ## $PPplot
    ## [1] 351
    ## 
    ## $PETplot
    ## [1] 351
    ## 
    ## $WAIplot
    ## [1] 351
    ## 
    ## $speimean
    ## [1] 351
    ## 
    ## $speimin
    ## [1] 351

``` r
# delete 3C plot types
dem_clim_final <- dem_clim_final |> 
  drop_na(
    CX, CY, 
    sf_nfi, sgdd_nfi, 
    PPplot, PETplot,
    WAIplot, speimean, speimin
    )

map(dem_clim_final, ~sum(is.na(.)))
```

    ## $IDPC3
    ## [1] 0
    ## 
    ## $CX
    ## [1] 0
    ## 
    ## $CY
    ## [1] 0
    ## 
    ## $k3
    ## [1] 0
    ## 
    ## $type2.x
    ## [1] 0
    ## 
    ## $Cut32
    ## [1] 0
    ## 
    ## $Ingrowth23
    ## [1] 0
    ## 
    ## $Growth23
    ## [1] 0
    ## 
    ## $Mortality23
    ## [1] 0
    ## 
    ## $Mortality23_nd
    ## [1] 0
    ## 
    ## $Mortality_ausente23
    ## [1] 0
    ## 
    ## $BAc23
    ## [1] 0
    ## 
    ## $ba_ifn2
    ## [1] 0
    ## 
    ## $sf_nfi
    ## [1] 0
    ## 
    ## $sgdd_nfi
    ## [1] 0
    ## 
    ## $PPplot
    ## [1] 0
    ## 
    ## $PETplot
    ## [1] 0
    ## 
    ## $WAIplot
    ## [1] 0
    ## 
    ## $speimean
    ## [1] 0
    ## 
    ## $speimin
    ## [1] 0

``` r
write_csv(dem_clim_final, here("01-data", "legacies", "dem_clim_ps_cut.csv"))
```

# Legacy transitions

``` r
dem_clim_final |> 
  tabyl(k3, type2.x) |> 
  adorn_percentages("row") |>
  adorn_pct_formatting(digits = 2) |>
  adorn_ns()
```

    ##  k3      Planted       Natural
    ##   1 24.33% (772) 75.67% (2401)
    ##   2 31.75% (520) 68.25% (1118)

``` r
# 3C included
k_ifn2 <- Psy23_final_long_k |>
  dplyr::select(IDPC3, ba_ifn, k3) |> 
  filter(ba_ifn == "ba_ha2") |> 
  rename(k3_ifn2 = k3) |> 
  dplyr::select(!ba_ifn)

k_ifn3 <- Psy23_final_long_k |>
  dplyr::select(IDPC3, ba_ifn, k3) |> 
  filter(ba_ifn == "ba_ha3") |> 
  rename(k3_ifn3 = k3) |> 
  dplyr::select(!ba_ifn)

k_ifn23 <- inner_join(k_ifn2, k_ifn3, by = "IDPC3")

cont_table_abs <- Psy23_final_long_k |> 
  tabyl(k3, ba_ifn)
cont_table_abs
```

    ##  k3 ba_ha2 ba_ha3
    ##   1   3430   2935
    ##   2   1732   2227

``` r
write_csv(cont_table_abs, here("01-data", "legacies", "cont_table_abs.csv"))

cont_table_chan <- k_ifn23 |> 
  tabyl(k3_ifn2, k3_ifn3)
cont_table_chan
```

    ##  k3_ifn2    1    2
    ##        1 2787  643
    ##        2  148 1584

``` r
write_csv(cont_table_chan, here("01-data", "legacies", "cont_table_chan.csv"))
```

# Boxplot cluster x sp origen

``` r
dem_clim <- read_csv(
  here("01-data", "legacies", "dem_clim_ps_cut.csv")
)

table(dem_clim$k3)
```

    ## 
    ##    1    2 
    ## 3173 1638

``` r
table(dem_clim$type2.x)
```

    ## 
    ## Natural Planted 
    ##    3519    1292

``` r
summary(dem_clim$Mortality23)
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ## 0.000000 0.000000 0.000000 0.002400 0.002289 0.090909

``` r
dem_clim <- dem_clim |> 
  unite("k_type", k3:type2.x) |> 
  mutate(
    mort_real = if_else(
      Cut32 == 0, Mortality_ausente23 + Mortality23, 
      Mortality23
    )
  )

table(dem_clim$k_type)
```

    ## 
    ## 1_Natural 1_Planted 2_Natural 2_Planted 
    ##      2401       772      1118       520

``` r
myColors <- brewer.pal(4, "Paired")

# growth
gg_growth_clorig <- dem_clim |> 
  mutate(
    k_type = recode_factor(
      k_type, 
      "1_Natural" = "C1 Natural",
      "1_Planted" = "C1 Planted",
      "2_Natural" = "C2 Natural",
      "2_Planted" = "C2 Planted",
      ),
    k_type = fct_relevel(
      k_type, c("C1 Natural", 
              "C1 Planted",
              "C2 Natural",
              "C2 Planted")
    ),
    ) |>
  ggplot(aes(y = Growth23, x = k_type, fill = k_type)) +
  geom_boxplot(outlier.size = 0.2, width = 0.5) +
  labs(x = "", title = "", y = "")  +
  scale_fill_manual(name = "cluster", values = myColors) +
  ylab(expression(paste("Growth (%", ~ha^-1~year^-1,")"))) +
  # coord_flip() +
  theme(
    legend.position = "none",
    axis.text = element_text(colour = "grey50", size = 10),
    axis.title = element_text(colour = "grey50", size = 12),
    panel.grid.major = element_line(colour = "grey90", size = 0.5), 
    panel.background = element_blank(),
    legend.key = element_blank(),
    strip.text.x = element_blank(),
    axis.ticks = element_blank()
  )

# structure
dem_clim_str23 <- left_join(
  dem_clim, Psy23_final_long_str, by = "IDPC3"
  )

test_that("nrow is correct", {
  expect_equal(nrow(dem_clim) * 2,
               nrow(dem_clim_str23))
})
```

    ## Test passed 🌈

``` r
dem_clim_str23_l <- dem_clim_str23 |>
  pivot_longer(
    cols = c(
      basal_area,
      density, dbh, height,
      height_dbh, cv_dbh, cv_height,
      saplings, sp_richness,
      tree_vigour
    ),
    names_to = "structure",
    values_to = "str_value"
  ) |> 
  dplyr::select(
    !ends_with("ifn")
  )

gg_str_clorig <- dem_clim_str23_l |> 
  mutate(
    k_type = recode_factor(
      k_type, 
      "1_Natural" = "C1 Natural",
      "1_Planted" = "C1 Planted",
      "2_Natural" = "C2 Natural",
      "2_Planted" = "C2 Planted",
    )
    ) |> 
  ggplot(aes(x = structure, y = str_value,
             fill = k_type)) +
  geom_boxplot(outlier.size = 0.2, width = 0.5) +
  labs(x = "", title = "", y = "")  +
  scale_fill_manual(name = "", values = myColors) +
  facet_wrap(~structure, scale = "free",
             strip.position = "bottom") +
  ylab("Structure") +
  theme(
    axis.text = element_text(colour = "grey50", size = 10),
    axis.title = element_text(colour = "grey50", size = 12),
    panel.grid.major = element_line(colour = "grey90", size = 0.5), 
    panel.background = element_blank(),
    legend.key = element_blank(),
    strip.text.x = element_blank(),
    axis.ticks = element_blank()
  )

# dbh of plots with dead trees
gg_dbh_mort <- dem_clim_str23_l |> 
  filter(mort_real > 0 & 
           structure == "dbh") |> 
  mutate(
    k_type = recode_factor(
      k_type, 
      "1_Natural" = "C1 Natural",
      "1_Planted" = "C1 Planted",
      "2_Natural" = "C2 Natural",
      "2_Planted" = "C2 Planted",
    )
  ) |> 
  ggplot(aes(x = structure, y = str_value,
             fill = k_type)) +
  geom_boxplot(outlier.size = 0.2, width = 0.5) +
  labs(x = "", title = "", y = "")  +
  scale_fill_manual(name = "", values = myColors) +
  facet_wrap(~structure, scale = "free",
             strip.position = "bottom") +
  ylab("d.b.h. (mm)") +
  theme(
    axis.text.y = element_text(colour = "grey50", size = 10),
    axis.text.x = element_blank(),
    axis.title = element_text(colour = "grey50", size = 12),
    panel.grid.major = element_line(colour = "grey90", size = 0.5), 
    panel.background = element_blank(),
    legend.key = element_blank(),
    strip.text.x = element_blank(),
    axis.ticks = element_blank()
  )
gg_dbh_mort
```

![](01-clustering_files/figure-gfm/cluster-origin-1.png)<!-- -->

``` r
# basal area of plots with dead trees
gg_ba_mort <- dem_clim_str23_l |> 
  filter(mort_real > 0 & 
           structure == "basal_area") |> 
  mutate(
    k_type = recode_factor(
      k_type, 
      "1_Natural" = "C1 Natural",
      "1_Planted" = "C1 Planted",
      "2_Natural" = "C2 Natural",
      "2_Planted" = "C2 Planted",
    )
  ) |> 
  ggplot(aes(x = structure, y = str_value,
             fill = k_type)) +
  geom_boxplot(outlier.size = 0.2, width = 0.5) +
  labs(x = "", title = "", y = "")  +
  scale_fill_manual(name = "", values = myColors) +
  facet_wrap(~structure, scale = "free",
             strip.position = "bottom") +
  ylab(bquote("Basal area"~(m^2~ha^-1))) +
  theme(
    axis.text.y = element_text(colour = "grey50", size = 10),
    axis.text.x = element_blank(),
    axis.title = element_text(colour = "grey50", size = 12),
    panel.grid.major = element_line(colour = "grey90", size = 0.5), 
    panel.background = element_blank(),
    legend.key = element_blank(),
    strip.text.x = element_blank(),
    axis.ticks = element_blank()
  )
gg_ba_mort
```

![](01-clustering_files/figure-gfm/cluster-origin-2.png)<!-- -->

``` r
# tree level
psy_m <- Psy23_planted_m |> 
  select(
    IDPC3, AB2m2ha_muertosc,
    AB2m2ha_ausente, AB2m2ha_presente,
    dbh2, AB2m2ha
  )

psy_m_plots <- left_join(
  psy_m, dem_clim, by = "IDPC3"
  ) |>
  filter(
    (Cut32 == 0 &
    AB2m2ha_muertosc > 0) |
      (Cut32 == 1 &
         AB2m2ha_presente > 0) 
    )

kk <- psy_m_plots |> 
  filter(Cut32 == 1) |> 
  select(AB2m2ha_ausente)

test_that("there are no muerto_ausente when the plot was cut", {
  expect_equal(sum(kk$AB2m2ha_ausente, na.rm = T),
               0)
})
```

    ## Test passed 🎉

``` r
gg_mort_dbh_ind <- psy_m_plots |> 
  mutate(
    k_type = recode_factor(
      k_type, 
      "1_Natural" = "C1 Natural",
      "1_Planted" = "C1 Planted",
      "2_Natural" = "C2 Natural",
      "2_Planted" = "C2 Planted",
    ),
    k_type = fct_relevel(
      k_type, c("C1 Natural", 
                "C1 Planted",
                "C2 Natural",
                "C2 Planted")
    ),
  ) |>
  ggplot(aes(y = dbh2, x = k_type, fill = k_type)) +
  geom_boxplot(outlier.size = 0.2, width = 0.5) +
  labs(x = "", title = "", y = "")  +
  scale_fill_manual(name = "cluster", values = myColors) +
  ylab("d.b.h. (mm)") +
  theme(
    legend.position = "none",
    axis.text = element_text(colour = "grey50", size = 10),
    axis.title = element_text(colour = "grey50", size = 12),
    panel.grid.major = element_line(colour = "grey90", size = 0.5), 
    panel.background = element_blank(),
    legend.key = element_blank(),
    strip.text.x = element_blank(),
    axis.ticks = element_blank()
  )
gg_mort_dbh_ind
```

![](01-clustering_files/figure-gfm/cluster-origin-3.png)<!-- -->

``` r
gg_mort_ba_ind <- psy_m_plots |> 
  mutate(
    k_type = recode_factor(
      k_type, 
      "1_Natural" = "C1 Natural",
      "1_Planted" = "C1 Planted",
      "2_Natural" = "C2 Natural",
      "2_Planted" = "C2 Planted",
    ),
    k_type = fct_relevel(
      k_type, c("C1 Natural", 
                "C1 Planted",
                "C2 Natural",
                "C2 Planted")
    ),
  ) |>
  ggplot(aes(y = AB2m2ha, x = k_type, fill = k_type)) +
  geom_boxplot(outlier.size = 0.2, width = 0.5) +
  labs(x = "", title = "", y = "")  +
  scale_fill_manual(name = "cluster", values = myColors) +
  ylab(bquote("Basal area"~(m^2~ha^-1))) +
  theme(
    legend.position = "none",
    axis.text = element_text(colour = "grey50", size = 10),
    axis.title = element_text(colour = "grey50", size = 12),
    panel.grid.major = element_line(colour = "grey90", size = 0.5), 
    panel.background = element_blank(),
    legend.key = element_blank(),
    strip.text.x = element_blank(),
    axis.ticks = element_blank()
  )
gg_mort_ba_ind
```

![](01-clustering_files/figure-gfm/cluster-origin-4.png)<!-- -->

``` r
ggsave(
  plot = gg_growth_clorig,
  here("03-results", "si_figures",
       "s7_growth.png"),
  width = 6, height = 4
)

ggsave(
  plot = gg_str_clorig,
  here("03-results", "si_figures",
       "s8_str_clorig.png"),
  width = 14, height = 6
)

ggsave(
  plot = gg_dbh_mort,
  here("03-results", "si_figures",
       "s9_mort_dbh.png"),
  width = 6, height = 4
)

ggsave(
  plot = gg_ba_mort,
  here("03-results", "si_figures",
       "s10_mort_ba.png"),
  width = 6, height = 4
)

ggsave(
  plot = gg_mort_dbh_ind,
  here("03-results", "si_figures",
       "s11_mort_dbh_ind.png"),
  width = 6, height = 4
)

ggsave(
  plot = gg_mort_ba_ind,
  here("03-results", "si_figures",
       "s12_mort_ba_ind.png"),
  width = 6, height = 4
)
```
