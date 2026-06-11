---
title: "Prion Deposition Quantification Pipeline"
author: "Kelsey Martin"
output:
  github_document:
    toc: true
    toc_depth: 3
---

# Overview

This repository contains a CellProfiler-based image analysis workflow and downstream R analysis scripts for quantifying prion protein deposition in brain tissue sections.

The primary goal of this pipeline is to identify, segment, and quantify prion deposits within tissue regions while excluding background signal, edge artifacts, and non-biologically relevant debris.

The workflow generates measurements describing:

- Prion deposition burden
- Deposit count
- Deposit size
- Deposit morphology
- Deposit staining intensity
- Tissue area analyzed

The resulting measurements can be compared across experimental groups, brain regions, and treatments to characterize disease-associated pathology.

# Experimental Goals

Prion disease pathology is characterized by the accumulation of misfolded prion protein deposits throughout the brain, which are darkly stained with DAB

This pipeline was designed to address the following biological questions:

1. Do experimental groups differ in the number of deposits?
2. Are deposits larger in size in one group compared to another?
3. Do deposits differ in staining intensity?
4. Are deposits morphologically distinct between groups?
5. Does overall deposition burden vary across brain regions?

To answer these questions, tissue regions are segmented first and individual deposits are quantified only within the segmented tissue area.

# Repository Structure

Data can be found in Alpine in a folder named CellProfilerData
The different subsets were run separately due to data size, but all run on the same pipeline with same settings, and later joined in R

The pipeline cpproj file can be uploaded to and run by the CellProfiler software

# Input Files
Input Images are the tifs provided by the Telling Lab, and the CellProfiler generated data used in analysis is described below:

## Depositions.csv

Object-level measurements for individual prion deposits.

Each row represents a single segmented deposit.

Measurements include:

- Area
- Shape descriptors
- Intensity metrics
- Spatial coordinates
- Parent tissue assignment

## Image.csv

Image-level summary measurements.

Each row represents one analyzed image.

Measurements include:

- Number of deposits
- Mean deposit area
- Total tissue area
- Area occupied by deposits

## Tissue.csv

Object-level measurements describing segmented tissue regions.

Each row represents a single tissue object.

Measurements include:

- Tissue area
- Tissue morphology
- Parent-child relationships with deposits

## pRions_v3_Control_Experiment.csv

CellProfiler experiment export file.

This file contains:

- Complete pipeline configuration
- Module settings
- Segmentation parameters
- Measurement settings

This file serves as the permanent record of the image analysis workflow and allows complete reproduction of the analysis.

# CellProfiler Pipeline Step by Step Description

**Software Version:** CellProfiler 4.2.8

## 1. Image Import

### Images

Loads grayscale histological images for analysis.

### Metadata

Metadata extraction disabled.

Images are processed independently.

### NamesAndTypes

Assigns imported images as:

```text
DNA
```

Although labeled as "DNA" within CellProfiler, this channel contains the grayscale image used for prion deposition analysis.

# Tissue Segmentation

## 2. Smooth

Gaussian smoothing is applied to reduce image noise and remove small punctate features prior to tissue segmentation.

**Output:**

```text
FilteredImage
```

## 3. EnhanceOrSuppressFeatures

Enhances large tissue structures and improves contrast between tissue and background.

## 4. Threshold

Global thresholding is used to separate tissue from background.

**Output:**

Binary tissue mask.

## 5. ConvertImageToObjects

Converts the binary tissue mask into individual tissue objects.

**Output:**

```text
TissueObjects
```

## 6. ExpandOrShrinkObjects

Tissue objects are shrunk inward.

### Purpose

- Remove edge-associated staining artifacts
- Exclude deposits touching tissue boundaries
- Restrict analysis to internal tissue regions

**Output:**

```text
ShrunkenTissue
```

## 7. MeasureObjectSizeShape

Measures tissue morphology including:

- Area
- Perimeter
- Solidity
- Eccentricity
- Major axis length
- Minor axis length

## 8. FilterObjects

Removes unwanted tissue objects based on size and morphology.

Examples include:

- Small tissue fragments
- Debris
- Sectioning artifacts

# Deposit Identification

## 9. InvertForPrinting

Prion deposits appear dark in the original image.

Images are inverted such that:

```text
Dark deposits → Bright objects
```

allowing standard object identification algorithms to be applied.

## 10. IdentifyPrimaryObjects

Identifies individual prion deposits.

**Output:**

```text
Depositions
```

Each segmented deposit becomes an independent object.

# Deposit Measurements

## 11. MeasureObjectIntensityDistribution

Measures the distribution of pixel intensities within each deposit.

Examples include:

- Lower quartile intensity
- Median intensity
- Upper quartile intensity
- Intensity heterogeneity

## 12. RelateObjects

Assigns deposits to tissue regions.

### Parent Object

```text
ShrunkenTissue
```

### Child Object

```text
Depositions
```

Purpose:

- Restrict analysis to deposits located within tissue
- Generate deposit counts per tissue object

## 13. MeasureObjectSizeShape

Measures deposit morphology.

Examples:

- Area
- Perimeter
- Circularity
- Solidity
- Compactness
- Feret diameters

## 14. MeasureObjectIntensity

Measures deposit staining intensity.

Examples:

- Mean intensity
- Integrated intensity
- Maximum intensity
- Median intensity

# Image-Level Measurements

## 15. MeasureImageAreaOccupied

Measures total tissue area.

Outputs include:

- Tissue area
- Percent image occupied by tissue

## 16. MeasureImageAreaOccupied

Measures total deposit burden.

Outputs include:

- Total deposit area
- Percent image occupied by deposits

# Data Export

## 17. ExportToSpreadsheet

Exports measurements to CSV files.

Generated outputs:

```text
Depositions.csv
Tissue.csv
Image.csv
```

# Downstream R Analysis

CellProfiler outputs are imported into R and metadata are extracted from image filenames.


```{r eval=FALSE}
pRion_Image <- pRions_Image %>%
  separate(
    FileName_DNA,
    into = c(
      "Mouse",
      "Treatment",
      "Region",
      "Magnification",
      "ImageNum"
    ),
    sep = "_"
  )
```

# AI Use Statement
ChatGPT model GPT-5.5 was used minorly for CellProfiler troubleshooting and also wrote the majority of this readme file when provided the experimental pipeline and context.
```

contains the complete CellProfiler workflow and should be retained alongside all raw image data and analysis outputs to ensure reproducibility of all measurements reported from this pipeline.
