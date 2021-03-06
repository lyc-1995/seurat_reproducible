---
title: "Test: Outputs of Seurat standard pipeline"
output: html_document
---



-   Tested pipeline:

`NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunUMAP(dims = 1:30)`

-   All parameters that weren't mentioned above were kept defaults.

-   Run the pipeline on two independent devices (`DESKTOP-5EA6KQM` and `ubuntu` as the hostnames), using the exactly same input and the same docker image.


```r
library(Seurat)
library(dplyr)
```

## Use original Seurat


```r
obj1 <- readRDS("../output/obj_std_DESKTOP-5EA6KQM_9629.Rds")
obj2 <- readRDS("../output/obj_std_ubuntu_254512.Rds")
```

### UMAP

-   The umap outputs were not reproducible!


```r
umap1 <- obj1[['umap']]@cell.embeddings
umap2 <- obj2[['umap']]@cell.embeddings
print(identical(umap1, umap2))
```

```
## [1] FALSE
```


```r
library(patchwork)
p1 <- DimPlot(obj1, group.by = "Final_anno", label = T) + NoLegend()
```

```
## Rasterizing points since number of points exceeds 100,000.
## To disable this behavior set `raster=FALSE`
```

```r
p2 <- DimPlot(obj2, group.by = "Final_anno", label = T) + NoLegend()
```

```
## Rasterizing points since number of points exceeds 100,000.
## To disable this behavior set `raster=FALSE`
```

```r
p1 + p2
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16-1.png)

### Variable features

-   The variable features were identical!


```r
print(all.equal(VariableFeatures(obj1), VariableFeatures(obj2)))
```

```
## [1] TRUE
```

### Scaled data

-   The scaled data were also identical!


```r
cm1 <- obj1[['RNA']]@scale.data
cm2 <- obj2[['RNA']]@scale.data
cm2 <- cm2[rownames(cm1), ]
print(all.equal(cm1, cm2))
```

```
## [1] TRUE
```


```r
print(identical(cm1, cm2))
```

```
## [1] TRUE
```

### PCA

-   The PCA outputs were the main source of differences!


```r
pc1 <- obj1[['pca']]@cell.embeddings
pc2 <- obj2[['pca']]@cell.embeddings
print(all.equal(pc1, pc2))
```

```
## [1] "Mean relative difference: 0.8369105"
```


```r
gc(reset = TRUE)
```

```
##              used    (Mb)  gc trigger    (Mb)   max used    (Mb)
## Ncells    3458452   184.8     6086921   325.1    3458452   184.8
## Vcells 4890243427 37309.6 11773350978 89823.6 4890243427 37309.6
```

## Use modified Seurat

-   In the original Seurat, `RunPCA()` use `irlba::irlba()` to calculate singular values. In the modified version, I replaced it with `RSpectra::svds()`.

`pca.results <- RSpectra::svds(A = object, k = npcs, ...)`

`pca.results <- RSpectra::svds(A = t(x = object), k = npcs, ...)`


```r
obj1 <- readRDS("../output/obj_std_modify_DESKTOP-5EA6KQM_9702.Rds")
obj2 <- readRDS("../output/obj_std_modify_ubuntu_254641.Rds")
gc(reset = TRUE)
```

```
##              used    (Mb)  gc trigger     (Mb)   max used    (Mb)
## Ncells    3459391   184.8     6086921    325.1    3459391   184.8
## Vcells 9253968999 70602.2 14128101173 107788.9 9253968999 70602.2
```

### UMAP

-   The umap outputs were identical! Hurray!


```r
umap1 <- obj1[['umap']]@cell.embeddings
umap2 <- obj2[['umap']]@cell.embeddings
print(identical(umap1, umap2))
```

```
## [1] TRUE
```


```r
p1 <- DimPlot(obj1, group.by = "Final_anno", label = T) + NoLegend()
```

```
## Rasterizing points since number of points exceeds 100,000.
## To disable this behavior set `raster=FALSE`
```

```r
p2 <- DimPlot(obj2, group.by = "Final_anno", label = T) + NoLegend()
```

```
## Rasterizing points since number of points exceeds 100,000.
## To disable this behavior set `raster=FALSE`
```

```r
p1 + p2
```

![plot of chunk unnamed-chunk-24](figure/unnamed-chunk-24-1.png)


```r
cm1 <- obj1[['RNA']]@scale.data
cm2 <- obj2[['RNA']]@scale.data
cm2 <- cm2[rownames(cm1), ]
print(identical(cm1, cm2))
```

```
## [1] TRUE
```

### PCA

-   The PCA outputs also passed `all.equal()`!


```r
pc1 <- obj1[['pca']]@cell.embeddings
pc2 <- obj2[['pca']]@cell.embeddings
print(all.equal(pc1, pc2))
```

```
## [1] TRUE
```

-   But did not pass `identical()`! Whatever, it dosen't matter anymore.


```r
print(identical(pc1, pc2))
```

```
## [1] FALSE
```
