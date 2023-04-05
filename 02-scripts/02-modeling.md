Modelling management legacies & climate change impacts on tree mortality
================

<style type="text/css">
pre {
  font-size: 10px
}
</style>

``` r
library(tidyverse)
library(here)
library(patchwork)
library(sf)
library(viridis)
library(glmmTMB)
library(DHARMa)
library(MASS, exclude = "select")
library(ggdist)
library(RColorBrewer)
library(testthat)
```

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
# absolute mortality
dem_clim_mort_abs <- dem_clim |> 
  mutate(
    mort_abs = ba_ifn2 * mort_real
  ) |> 
  filter(mort_abs > 0) |> 
  group_by(k_type) |> 
  summarise(
    mort_abs_type = mean(mort_abs, na.rm = T)
  )
```

# Exploratory data analyses

``` r
glimpse(dem_clim)
```

    ## Rows: 4,811
    ## Columns: 20
    ## $ IDPC3               <chr> "100303A1", "10067A1", "100773A1", "100775A1", "10…
    ## $ CX                  <dbl> 262000, 490000, 263000, 263000, 262000, 258000, 49…
    ## $ CY                  <dbl> 4470000, 4774000, 4468000, 4467000, 4463000, 44590…
    ## $ k_type              <chr> "1_Planted", "1_Planted", "1_Planted", "1_Planted"…
    ## $ Cut32               <dbl> 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0,…
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
    ## $ mort_real           <dbl> 0.000000000, 0.000000000, 0.008337227, 0.000000000…

``` r
summary(dem_clim$mort_real)
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ## 0.000000 0.000000 0.000000 0.003459 0.003593 0.090909

``` r
summary(dem_clim$WAIplot)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ## -67.186 -29.780 -11.345  -6.280   9.726 127.395

``` r
summary(dem_clim$speimin)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  -2.565  -1.692  -1.623  -1.619  -1.582  -1.123

``` r
ggplot(dem_clim) +
  geom_histogram(aes(mort_real), bins = 50)
```

![](02-modeling_files/figure-gfm/EDA-1.png)<!-- -->

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
        legend.text = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),  
        axis.text.x = element_text(colour = "grey50", size = 8),
        axis.text.y = element_text(colour = "grey50", size = 8),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(colour = "grey90", size = 0.5), 
        axis.line = element_blank(),
        panel.background = element_blank()) 

# set as sf and assign crs
dem_clim_sf <- st_as_sf(x = dem_clim, 
                        coords = c("CX", "CY"),
                        crs = "+proj=utm +zone=30 +ellps=intl +units=m +no_defs")
st_crs(dem_clim_sf)
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
# ggplot() + 
#   geom_sf(data = dem_clim_sf)

# transform to long lat
dem_clim_sf_map <- st_transform(
  dem_clim_sf,
  "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
)

st_crs(dem_clim_sf_map)
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
# convert to tibble
dem_clim_sp_tb <- dem_clim_sf_map |>
  ungroup() |> 
  mutate(CX = unlist(map(dem_clim_sf_map$geometry, 1)),
         CY = unlist(map(dem_clim_sf_map$geometry, 2))) |> 
  as_tibble() |> 
  dplyr::select(!geometry)

map_plot <- function(col){
  map_loc + 
    geom_point(data = dem_clim_sp_tb,
               aes(CX, CY, color = .data[[col]]), size = .1, alpha = .8) + 
    scale_color_viridis() +
    xlab("") +
    ylab("") + 
    theme(plot.margin = margin(0, 0, 0, 0),
          legend.title = element_text(color = "black", size = 9),
          legend.text = element_text(size = 9),
          legend.key = element_blank(),
          legend.position = "top",
          legend.justification = "left",
          legend.background = element_rect(fill = "white", colour = NA),
          axis.title = element_blank(),
          axis.text.x = element_text(colour = "grey50", size = 7),
          axis.text.y = element_text(colour = "grey50", size = 7),
          panel.grid.major = element_line(colour = "grey90", size = 0.5), 
          axis.line = element_blank(),
          panel.background = element_blank())
}

map_plot("WAIplot") / map_plot("speimin")
```

![](02-modeling_files/figure-gfm/EDA-2.png)<!-- -->

``` r
# full data
gg_mort <- function(x){
  dem_clim |> 
    ggplot(aes(x = .data[[x]], y = Mortality23)) + 
    geom_hex() + 
    ggtitle("Full data = ")
}

var_int_all <- c(
  "WAIplot", "speimin", "k_type"
)

map(var_int_all, ~gg_mort(.)) |> 
  reduce(`/`)
```

![](02-modeling_files/figure-gfm/EDA-3.png)<!-- -->

``` r
ggsave(
  here("03-results", "others", "distri_explanatory.png"),
  width = 12, height = 4
)

# n. of plots
nrow(dem_clim)
```

    ## [1] 4811

``` r
# n. of plots with any dead trees
sum(dem_clim$Mortality23 > 0)
```

    ## [1] 1383

``` r
## correlations ------------------------------------------------------------
corm <- dem_clim |>
  dplyr::select(
    WAIplot, speimin
  ) |>
  corrr::correlate(diagonal = 1) |>
  corrr::shave(upper = FALSE) |>
  pivot_longer(
    cols = -term,
    names_to = "colname",
    values_to = "corr"
  ) |>
  mutate(term = fct_inorder(term),
         colname = fct_inorder(colname))

gg_cor <- ggplot(corm, aes(term, fct_rev(colname),
                 fill = corr)) +
  geom_tile() +
  geom_text(aes(
    label = format(round(corr, 2), nsmall = 2),
    color = abs(corr) < .60)) +
  coord_fixed(expand = FALSE) +
  scale_color_manual(values = c("white", "black"),
                     guide = "none") +
  scale_fill_distiller(
    palette = "PuOr", na.value = "white",
    direction = 1, limits = c(-1, 1)) +
  labs(x = NULL, y = NULL) +
  theme(panel.border = element_rect(color = NA, fill = NA),
        legend.position = c(.85, .8))

gg_cor
```

![](02-modeling_files/figure-gfm/EDA-4.png)<!-- -->

``` r
ggsave(
  plot = gg_cor,
  here("03-results", "si_figures", "s4_correlations.png"),
  width = 3.5, height = 5,
  dpi = 600
)
```

# Modelling

``` r
# https://stackoverflow.com/questions/65745148/is-there-a-difference-between-gamma-hurdle-two-part-models-and-zero-inflated-g

zigamma_model_1 <- glmmTMB(
  mort_real ~ k_type * WAIplot * speimin,
  family = ziGamma(link = "log"),
  # we assume that absences will at least vary by:
  ziformula = ~ k_type * WAIplot * speimin,
  data = dem_clim
  )

# model output interpretation:
# https://journal.r-project.org/archive/2017/RJ-2017-066/RJ-2017-066.pdf
# https://cran.r-project.org/web/packages/glmmTMB/glmmTMB.pdf
# https://cran.r-project.org/web/packages/glmmTMB/vignettes/glmmTMB.pdf
summary(zigamma_model_1)
```

    ##  Family: Gamma  ( log )
    ## Formula:          mort_real ~ k_type * WAIplot * speimin
    ## Zero inflation:             ~k_type * WAIplot * speimin
    ## Data: dem_clim
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##  -6208.9  -5995.1   3137.5  -6274.9     4778 
    ## 
    ## 
    ## Dispersion estimate for Gamma family (sigma^2): 0.68 
    ## 
    ## Conditional model:
    ##                                   Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                     -5.110e+00  3.714e-01 -13.757  < 2e-16 ***
    ## k_type1_Planted                  1.171e+00  7.681e-01   1.524  0.12754    
    ## k_type2_Natural                  1.051e+00  6.187e-01   1.699  0.08932 .  
    ## k_type2_Planted                  1.719e+00  6.641e-01   2.588  0.00965 ** 
    ## WAIplot                          6.748e-03  1.031e-02   0.655  0.51275    
    ## speimin                         -4.107e-01  2.256e-01  -1.820  0.06875 .  
    ## k_type1_Planted:WAIplot         -2.043e-05  2.445e-02  -0.001  0.99933    
    ## k_type2_Natural:WAIplot         -3.427e-03  1.659e-02  -0.207  0.83632    
    ## k_type2_Planted:WAIplot         -4.456e-02  2.233e-02  -1.996  0.04594 *  
    ## k_type1_Planted:speimin          6.550e-01  4.856e-01   1.349  0.17734    
    ## k_type2_Natural:speimin          9.812e-01  3.758e-01   2.611  0.00902 ** 
    ## k_type2_Planted:speimin          1.113e+00  4.200e-01   2.650  0.00805 ** 
    ## WAIplot:speimin                  6.079e-03  6.163e-03   0.986  0.32397    
    ## k_type1_Planted:WAIplot:speimin  3.237e-03  1.581e-02   0.205  0.83779    
    ## k_type2_Natural:WAIplot:speimin -3.705e-03  9.853e-03  -0.376  0.70692    
    ## k_type2_Planted:WAIplot:speimin -2.725e-02  1.443e-02  -1.889  0.05893 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Zero-inflation model:
    ##                                  Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)                      1.239139   0.636473   1.947  0.05155 . 
    ## k_type1_Planted                 -1.112860   1.078765  -1.032  0.30226   
    ## k_type2_Natural                 -3.108057   1.160303  -2.679  0.00739 **
    ## k_type2_Planted                 -0.085724   1.112357  -0.077  0.93857   
    ## WAIplot                         -0.016889   0.017210  -0.981  0.32641   
    ## speimin                          0.152914   0.387175   0.395  0.69288   
    ## k_type1_Planted:WAIplot          0.088537   0.032187   2.751  0.00595 **
    ## k_type2_Natural:WAIplot          0.062538   0.033318   1.877  0.06052 . 
    ## k_type2_Planted:WAIplot          0.014420   0.034406   0.419  0.67513   
    ## k_type1_Planted:speimin         -0.912515   0.677803  -1.346  0.17821   
    ## k_type2_Natural:speimin         -1.182461   0.703455  -1.681  0.09278 . 
    ## k_type2_Planted:speimin          0.548035   0.706782   0.775  0.43811   
    ## WAIplot:speimin                 -0.004183   0.010270  -0.407  0.68378   
    ## k_type1_Planted:WAIplot:speimin  0.054493   0.020585   2.647  0.00812 **
    ## k_type2_Natural:WAIplot:speimin  0.036733   0.019966   1.840  0.06580 . 
    ## k_type2_Planted:WAIplot:speimin  0.006597   0.022230   0.297  0.76666   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# This summary can be broken down into 3 sections

# The top section is a general overview
# containing a description of the model specification
# (Family, Formula, Zero inflation, Dispersion
# Data) and resulting information criteria

# The second section describes the coefficients 
# of the Conditional model
# including Wald Z statistics and p-values. 
# Apart from the intercept, the estimates are all
# contrasts as is standard in regression models.
# This model has a log link as stated in the 
# top line of the summary

# The third section describes the Zero-inflation model 
# similarly to the Conditional model except that
# this model has a logit link. The zero-inflation model 
# estimates the probability of an extra zero such that
# a positive contrast indicates a higher chance of
# absence; this is the opposite of
# the conditional model where a positive contrast 
# indicates a higher abundance

# https://cran.r-project.org/web/packages/glmmTMB/vignettes/model_evaluation.pdf

## model checking and diagnostics ------------------------------------------
# https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
zigamma_model_simres_1 <- simulateResiduals(zigamma_model_1)
# x11()
# par(mfrow = c(2, 2))
# plot(zigamma_model_simres_1)

plotQQunif(zigamma_model_simres_1)
```

![](02-modeling_files/figure-gfm/modelling-1.png)<!-- -->

``` r
plotResiduals(zigamma_model_simres_1)
```

![](02-modeling_files/figure-gfm/modelling-2.png)<!-- -->

``` r
testZeroInflation(zigamma_model_simres_1)
```

![](02-modeling_files/figure-gfm/modelling-3.png)<!-- -->

    ## 
    ##  DHARMa zero-inflation test via comparison to expected zeros with
    ##  simulation under H0 = fitted model
    ## 
    ## data:  simulationOutput
    ## ratioObsSim = 0.99946, p-value = 0.992
    ## alternative hypothesis: two.sided

``` r
# testZeroInflation(zigamma_model_simres_1)
# shows the simulated number
# of zeroes with your actual number shown as red. 
# we have the same number as expected
# so there is no zero-inflation now
# makes sense to fit a zero-inflated model

## predictions: average predictive comparisons -------------------------------------------
summary(zigamma_model_1)
```

    ##  Family: Gamma  ( log )
    ## Formula:          mort_real ~ k_type * WAIplot * speimin
    ## Zero inflation:             ~k_type * WAIplot * speimin
    ## Data: dem_clim
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##  -6208.9  -5995.1   3137.5  -6274.9     4778 
    ## 
    ## 
    ## Dispersion estimate for Gamma family (sigma^2): 0.68 
    ## 
    ## Conditional model:
    ##                                   Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                     -5.110e+00  3.714e-01 -13.757  < 2e-16 ***
    ## k_type1_Planted                  1.171e+00  7.681e-01   1.524  0.12754    
    ## k_type2_Natural                  1.051e+00  6.187e-01   1.699  0.08932 .  
    ## k_type2_Planted                  1.719e+00  6.641e-01   2.588  0.00965 ** 
    ## WAIplot                          6.748e-03  1.031e-02   0.655  0.51275    
    ## speimin                         -4.107e-01  2.256e-01  -1.820  0.06875 .  
    ## k_type1_Planted:WAIplot         -2.043e-05  2.445e-02  -0.001  0.99933    
    ## k_type2_Natural:WAIplot         -3.427e-03  1.659e-02  -0.207  0.83632    
    ## k_type2_Planted:WAIplot         -4.456e-02  2.233e-02  -1.996  0.04594 *  
    ## k_type1_Planted:speimin          6.550e-01  4.856e-01   1.349  0.17734    
    ## k_type2_Natural:speimin          9.812e-01  3.758e-01   2.611  0.00902 ** 
    ## k_type2_Planted:speimin          1.113e+00  4.200e-01   2.650  0.00805 ** 
    ## WAIplot:speimin                  6.079e-03  6.163e-03   0.986  0.32397    
    ## k_type1_Planted:WAIplot:speimin  3.237e-03  1.581e-02   0.205  0.83779    
    ## k_type2_Natural:WAIplot:speimin -3.705e-03  9.853e-03  -0.376  0.70692    
    ## k_type2_Planted:WAIplot:speimin -2.725e-02  1.443e-02  -1.889  0.05893 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Zero-inflation model:
    ##                                  Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)                      1.239139   0.636473   1.947  0.05155 . 
    ## k_type1_Planted                 -1.112860   1.078765  -1.032  0.30226   
    ## k_type2_Natural                 -3.108057   1.160303  -2.679  0.00739 **
    ## k_type2_Planted                 -0.085724   1.112357  -0.077  0.93857   
    ## WAIplot                         -0.016889   0.017210  -0.981  0.32641   
    ## speimin                          0.152914   0.387175   0.395  0.69288   
    ## k_type1_Planted:WAIplot          0.088537   0.032187   2.751  0.00595 **
    ## k_type2_Natural:WAIplot          0.062538   0.033318   1.877  0.06052 . 
    ## k_type2_Planted:WAIplot          0.014420   0.034406   0.419  0.67513   
    ## k_type1_Planted:speimin         -0.912515   0.677803  -1.346  0.17821   
    ## k_type2_Natural:speimin         -1.182461   0.703455  -1.681  0.09278 . 
    ## k_type2_Planted:speimin          0.548035   0.706782   0.775  0.43811   
    ## WAIplot:speimin                 -0.004183   0.010270  -0.407  0.68378   
    ## k_type1_Planted:WAIplot:speimin  0.054493   0.020585   2.647  0.00812 **
    ## k_type2_Natural:WAIplot:speimin  0.036733   0.019966   1.840  0.06580 . 
    ## k_type2_Planted:WAIplot:speimin  0.006597   0.022230   0.297  0.76666   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# reduction in WAI & spei per plot between the median and the 30 percentile
# http://www.ub.edu/irbio/droughts-and-rising-sea-levels-are-the-impacts-of-climate-change-that-will-most-affect-the-mediterranean-basin-n-952-en
summary(dem_clim$WAIplot)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ## -67.186 -29.780 -11.345  -6.280   9.726 127.395

``` r
quantile(dem_clim$WAIplot, probs = c(.25, .3, .5))
```

    ##       25%       30%       50% 
    ## -29.77998 -26.94356 -11.34503

``` r
summary(dem_clim$speimin)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  -2.565  -1.692  -1.623  -1.619  -1.582  -1.123

``` r
quantile(dem_clim$speimin, probs = c(.25, .3, .5))
```

    ##       25%       30%       50% 
    ## -1.692370 -1.692370 -1.622561

``` r
# climate impact
dem_clim_cc <- dem_clim |> 
  mutate(
    WAIplot = WAIplot - 15.59853,
    speimin = speimin - 0.069809809
  )

# without climate impact
dem_clim_no_cc <- dem_clim

dem_clim_new <- list(
  dem_clim_cc,
  dem_clim_no_cc
)

# https://www.rdocumentation.org/packages/glmmTMB/versions/1.1.3/topics/predict.glmmTMB
dem_clim_pred <- map(dem_clim_new, ~predict(
  object = zigamma_model_1, newdata = .x,
  se.fit = TRUE, type = "conditional"
))

# add predictions to the corresponding dataset
add_predictions <- function(dat, pred){
  dat |> 
    mutate(
      predFE = pred$fit
    )
}

arg2_dat_pred <- list(
  dat = dem_clim_new, 
  pred = dem_clim_pred
  )

dem_clim_new_pred <- arg2_dat_pred |>
  pmap(add_predictions)

# length(dem_clim_new_pred)

names(dem_clim_new_pred) <- c(
  "dem_clim_cc",
  "dem_clim_no_cc"
)

# add prediction name to each dataset
mutate_name <- function(x, names_df){
  dem_clim_new_pred[[x]] |> 
    dplyr::mutate(
      df = names_df
    )
}

arg2_x_names_df <- list(
  x = 1:length(dem_clim_new_pred), 
  names_df = names(dem_clim_new_pred)
  )

dem_clim_new_pred_n <- arg2_x_names_df |>
  pmap_df(mutate_name)

table(dem_clim_new_pred_n$df)
```

    ## 
    ##    dem_clim_cc dem_clim_no_cc 
    ##           4811           4811

``` r
table(dem_clim_new_pred_n$k_type)
```

    ## 
    ## 1_Natural 1_Planted 2_Natural 2_Planted 
    ##      4802      1544      2236      1040

``` r
# compare predictions to actual data
dem_clim_c <- dem_clim_new_pred_n |> 
  filter(df == "dem_clim_no_cc") |> 
  filter(predFE > 0 & mort_real > 0)
cor(dem_clim_c$mort_real, dem_clim_c$predFE)
```

    ## [1] 0.2697114

``` r
# mean difference between predicted and 
# predicted under cc
dem_clim_diff <- dem_clim_new_pred_n |> 
  filter(mort_real > 0) |> 
  group_by(df, k_type) |> 
  summarise(mort_df_type = mean(predFE, na.rm = T) * 100)

# c1 natural
dem_clim_diff$mort_df_type[dem_clim_diff$df == "dem_clim_cc" & dem_clim_diff$k_type == "1_Natural"] / dem_clim_diff$mort_df_type[dem_clim_diff$df == "dem_clim_no_cc" & dem_clim_diff$k_type == "1_Natural"]
```

    ## [1] 1.089844

``` r
# c1 planted
dem_clim_diff$mort_df_type[dem_clim_diff$df == "dem_clim_cc" & dem_clim_diff$k_type == "1_Planted"] / dem_clim_diff$mort_df_type[dem_clim_diff$df == "dem_clim_no_cc" & dem_clim_diff$k_type == "1_Planted"]
```

    ## [1] 1.128943

``` r
# c2 natural
dem_clim_diff$mort_df_type[dem_clim_diff$df == "dem_clim_cc" & dem_clim_diff$k_type == "2_Natural"] / dem_clim_diff$mort_df_type[dem_clim_diff$df == "dem_clim_no_cc" & dem_clim_diff$k_type == "2_Natural"]
```

    ## [1] 0.9729331

``` r
# c2 planted
dem_clim_diff$mort_df_type[dem_clim_diff$df == "dem_clim_cc" & dem_clim_diff$k_type == "2_Planted"] / dem_clim_diff$mort_df_type[dem_clim_diff$df == "dem_clim_no_cc" & dem_clim_diff$k_type == "2_Planted"]
```

    ## [1] 0.9924395

``` r
# https://cran.r-project.org/web/packages/ggdist/vignettes/slabinterval.html#multiple-slabs-and-intervals-in-composite-plots
# Quantile dotplots
gg_dem_clim_pred <- dem_clim_new_pred_n |>
  # when there is mortality (conditional part)
  filter(predFE > 0) |> 
  mutate(
    k_type = recode_factor(
      k_type, 
      "1_Natural" = "C1 Natural",
      "1_Planted" = "C1 Planted",
      "2_Natural" = "C2 Natural",
      "2_Planted" = "C2 Planted"
    ),
    df = recode_factor(
      df, 
      "dem_clim_no_cc" = "Predicted",
      "dem_clim_cc" = "Predicted under climate change"
    )
  ) |>
  ggplot(aes(x = k_type, y = predFE * 100,
             color = df, fill = df
             )) +
  stat_dots(
    scale = 0.5,
    quantiles = 50,
    .width = 0,
    alpha = .8,
    color = NA,
    size = 0
    ) +
    stat_halfeye(
    .width = 0,
    slab_fill = NA,
    size = 1,
    point_interval = "mean_qi"
  ) +
  ylab(expression(paste("Mortality (%",~ha^-1~year^-1,")"))) +
  scale_fill_manual(values = c("#E1BE6A", "#40B0A6"),
                    guide = "none") +
  scale_color_manual(values = c("#E1BE6A", "#40B0A6")) +
  coord_cartesian(ylim = c(0, 3)) +
  theme(
    text = element_text(family = "sans"),
    legend.text = element_text(size = 7),
    legend.key = element_blank(),
    legend.title = element_blank(),
    panel.grid.major = element_line(colour = "grey90", size = 0.5),
    panel.background = element_blank(),
    axis.text = element_text(size = 7, color = "black"),
    axis.title.y = element_text(size = 7),
    axis.title.x = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(linetype = 0)))

gg_dem_clim_pred
```

![](02-modeling_files/figure-gfm/modelling-4.png)<!-- -->

``` r
# save plots --------------------------------------------------------------
ggsave(
  plot = map_plot("WAIplot") / map_plot("speimin") +
    plot_annotation(tag_levels = "a"),
  here("03-results", "si_figures", "s3_wai_spei_distribution.png"),
  width = 5, height = 8
)

ggsave(
  plot = gg_dem_clim_pred + 
    theme(
      legend.position = "bottom",
      plot.tag = element_text(face = "bold",
                              size = 9)
      ),
  here("03-results", "main_figure", "main_plot.png"),
  width = 9, height = 9,
  unit = "cm",
  dpi = 600
)

library(svglite)
ggsave(
  plot = gg_dem_clim_pred + 
    theme(
      legend.position = "bottom",
      plot.tag = element_text(face = "bold",
                              size = 9)
      ),
  here("03-results", "main_figure", "main_plot.svg"),
  width = 9, height = 9,
  unit = "cm",
  dpi = 600
)
```

------------------------------------------------------------------------

<details>
<summary>

Session Info

</summary>

``` r
Sys.time()
```

    ## [1] "2023-04-05 12:34:16 CEST"

``` r
git2r::repository()
```

    ## Local:    main C:/Users/julen/OneDrive/Escritorio/GitHub-col/legacies
    ## Remote:   main @ origin (https://github.com/Julenasti/legacies.git)
    ## Head:     [b96bdda] 2023-03-14: remove unnecesary models

``` r
sessionInfo()
```

    ## R version 4.2.2 (2022-10-31 ucrt)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19044)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United Kingdom.utf8 
    ## [2] LC_CTYPE=English_United Kingdom.utf8   
    ## [3] LC_MONETARY=English_United Kingdom.utf8
    ## [4] LC_NUMERIC=C                           
    ## [5] LC_TIME=English_United Kingdom.utf8    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] svglite_2.1.0      testthat_3.1.4     RColorBrewer_1.1-3 ggdist_3.1.1      
    ##  [5] MASS_7.3-58.1      DHARMa_0.4.5       glmmTMB_1.1.3      viridis_0.6.2     
    ##  [9] viridisLite_0.4.1  sf_1.0-7           patchwork_1.1.1    here_1.0.1        
    ## [13] forcats_0.5.1      stringr_1.4.1      dplyr_1.0.9        purrr_0.3.4       
    ## [17] readr_2.1.2        tidyr_1.2.0        tibble_3.1.7       ggplot2_3.3.6     
    ## [21] tidyverse_1.3.2   
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] TH.data_1.1-1        googledrive_2.0.0    minqa_1.2.4         
    ##   [4] colorspace_2.0-3     ellipsis_0.3.2       class_7.3-20        
    ##   [7] rprojroot_2.0.3      estimability_1.4.1   fs_1.5.2            
    ##  [10] rstudioapi_0.13      proxy_0.4-27         hexbin_1.28.2       
    ##  [13] farver_2.1.1         bit64_4.0.5          fansi_1.0.3         
    ##  [16] mvtnorm_1.1-3        lubridate_1.8.0      xml2_1.3.3          
    ##  [19] codetools_0.2-18     splines_4.2.2        knitr_1.40.1        
    ##  [22] jsonlite_1.8.0       nloptr_2.0.3         broom_1.0.0         
    ##  [25] dbplyr_2.2.1         compiler_4.2.2       httr_1.4.3          
    ##  [28] emmeans_1.8.1-1      backports_1.4.1      assertthat_0.2.1    
    ##  [31] Matrix_1.5-1         fastmap_1.1.0        gargle_1.2.0        
    ##  [34] cli_3.3.0            htmltools_0.5.3      tools_4.2.2         
    ##  [37] coda_0.19-4          gtable_0.3.0         glue_1.6.2          
    ##  [40] corrr_0.4.3          maps_3.4.0           Rcpp_1.0.9          
    ##  [43] cellranger_1.1.0     vctrs_0.5.0          nlme_3.1-160        
    ##  [46] xfun_0.32            brio_1.1.3           lme4_1.1-30         
    ##  [49] rvest_1.0.2          lifecycle_1.0.3      gap.datasets_0.0.5  
    ##  [52] googlesheets4_1.0.0  zoo_1.8-10           scales_1.2.1        
    ##  [55] vroom_1.5.7          ragg_1.2.5           hms_1.1.1           
    ##  [58] parallel_4.2.2       sandwich_3.0-2       TMB_1.9.0           
    ##  [61] yaml_2.3.5           gridExtra_2.3        stringi_1.7.8       
    ##  [64] highr_0.9            gap_1.2.3-6          e1071_1.7-11        
    ##  [67] boot_1.3-28          systemfonts_1.0.4    rlang_1.0.6         
    ##  [70] pkgconfig_2.0.3      distributional_0.3.0 evaluate_0.18       
    ##  [73] lattice_0.20-45      labeling_0.4.2       bit_4.0.4           
    ##  [76] tidyselect_1.1.2     magrittr_2.0.3       R6_2.5.1            
    ##  [79] generics_0.1.3       multcomp_1.4-23      DBI_1.1.3           
    ##  [82] pillar_1.8.1         haven_2.5.0          withr_2.5.0         
    ##  [85] units_0.8-0          survival_3.4-0       modelr_0.1.8        
    ##  [88] crayon_1.5.2         KernSmooth_2.23-20   utf8_1.2.2          
    ##  [91] tzdb_0.3.0           rmarkdown_2.16       grid_4.2.2          
    ##  [94] readxl_1.4.0         git2r_0.30.1         reprex_2.0.1        
    ##  [97] digest_0.6.29        classInt_0.4-7       xtable_1.8-4        
    ## [100] numDeriv_2016.8-1.1  textshaping_0.3.6    munsell_0.5.0

</details>
