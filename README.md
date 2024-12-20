# Non-negative Matrix Factorization

This code was originally developed and written by the Kaspar Podgorski lab at Janelia.

Non-negative matrix factorization (NMF) is a technique used to find correlated regions in a series of images. The goal is to find a limited number basis set of images that can be combined with different weightings to recreate - as best as possible - each of the different images in the dataset. [Here](https://blog.acolyer.org/2019/02/18/the-why-and-how-of-nonnegative-matrix-factorization/) is a good explanation. 

In this repository, the goal is to use non-negative matrix factorization to find correlated clusters of pixels in a time series of fluorescent images during neuronal firing. In other words, it is used to find clusters of pixels that fire similarly throughout the time series. This process involves several stages:

1. The code takes in a flourescent time series as well as an optional sigma to be used in a guassian smoothing of images; the default sigma is 2.
2. A smoothed fluorescent baseline F0 is calculated for each pixel using a leaky cumulative minimum and repeated lowpass filtering to converge on a smooth F0 that obeys the minima.
3. &Delta;F is calculated for each image in the series by subtracting F0 from the image; F0 is then ensured to be positive.
3. A high pass filter is applied to the &Delta;F time series. The resultant images are smoothed and used to create a correlation image.
4. Peaks in the correlation image are used as input to the NMF to help initialize cluster finding.
5. NMF is then applied to the &Delta;F time series, with additional contiguous contraints on clusters such that each correlated cluster of pixels must be contiguous.
6. Small correlated clusters are then merged.
7. Cluster metrics are calculated using the most sensitive pixels per cluster.
8. All relevant information is returned as a struct.

 `processMiniData_edited.m` is the function used to perform the NMF.
