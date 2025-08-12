# SENNA  
Spatial Expression aNalysis with spliNe Axes


## Description  
SENNA is a user-interactive framework that maps spots from spatially resolved transcriptomics and proteomics datasets onto a ‘curve axis’—a user-drawn path or region boundary on the tissue.  

## Installation  

```
install.packages("devtools")
library(devtools)
install_github("Dalmooree/SENNA")
```
or  
```
install.packages("remotes")
library(remotes)
install_github("Dalmooree/SENNA")
```

## Run  
### Load packages  
```
library(Seurat)
library(SENNA)
```

### Pre-processing  
The sample dataset is available at [10X Genomics repository](https://www.10xgenomics.com/datasets/mouse-brain-coronal-section-1-ffpe-2-standard).


```
# Set data path
dpath <- "../dataset/st/CytAssist_FFPE_Mouse_Brain_Rep1_filtered"

# Load SRT data as `Seurat` object
surt <- Load10X_Spatial(data.dir = dpath,
                        filename = "CytAssist_FFPE_Mouse_Brain_Rep1_filtered_feature_bc_matrix.h5",
                        slice = "mb",
                        filter.matrix = TRUE,
                        image = NULL)

if(min(surt$nCount_Spatial) == 0) surt <- subset(surt, nCount_Spatial > 0)
# `NormalizeData()` is also available
surt <- SCTransform(surt, assay = "Spatial", verbose = FALSE)

# Optional pre-processing
surt <- RunPCA(surt, verbose = FALSE)
surt <- FindNeighbors(surt, verbose = FALSE)
surt <- FindClusters(surt, verbose = FALSE)
```

### `SENNA` object creation and knots picking    
```
# Create SENNA object
sen <- SENNA_Visium(surt,
                    slice_name = "mb",
                    annotation = TRUE)

# Use `SENNA_Xenium()` for Xenium dataset  
# Use `SENNA_CODEX()` for CODEX dataset  
# Use `SENNA_CosMx()` for CosMx dataset  
```

1. Knots picker with default option  
```
AppDat(sen)
knot_picker()
```

![Default kp](images/1_kp1.png)

2. Knots picker with tissue image (Only Visium and VisiumHD are supported currently)
```
# Image path
ipath <- "../dataset/st/CytAssist_FFPE_Mouse_Brain_Rep1_filtered/spatial/"
AppDat(sen,
       image_path = ipath,
       image_resolution = "lowres") # or "hires"
knot_picker()
```

![Image kp](images/1_kp2.png)

3. Knots picker with discrete spot attribute (e.g. Clusters, cell types, ...)  
```
# List of available reference value
names(sen@Gene$Reference)
AppDat(sen,
       reference_value = "Annotation",
       colorset = ggsci::pal_npg("nrc")(9),
       image_path = ipath,
       image_resolution = "lowres")
knot_picker()
```

![Disc kp](images/1_kp31.png)

![Disc kp with image](images/1_kp32.png)


4. Knots picker with continuous spot attribute (e.g. Expression level, CNV score, ...)  
```
# Add reference value
sen <- AddReference(sen,
                    var_name = "Ttr Expression", # Variable name
                    reference = sen@Gene[["Spatial"]][["Ttr"]] # Value
                    )
AppDat(sen,
       reference_value = "Ttr Expression",
       colorset = c("darkblue", "darkred"),
       image_path = ipath,
       image_resolution = "lowres")
knot_picker()
```

![Cont kp](images/1_kp41.png)

![Cont kp with image](images/1_kp42.png)


### Curve axis generation and projection  
```
```

### (Optional) CSD computation
This step is performed only in regionation and islet analysis scenarios.  
```
```

### SVG detection
```
```

### Cell attributes modeling  
```
```

## Citation  


## License  
**To be determined. Please do not use, redistribute, or modify without permission. Patent pending.**  

