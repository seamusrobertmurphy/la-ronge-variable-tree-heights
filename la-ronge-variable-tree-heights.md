La Ronge Tree Heights
================
SMurphy
2023-10-15

``` r
dsm = raster::raster("/media/seamus/USB1/la_ronge/dsm_la_ronge.tif", select = 'xyzcr', filter = '-drop_class 19')
plot(dsm)
```

![](la-ronge-variable-tree-heights_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->
