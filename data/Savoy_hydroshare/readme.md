## Supporting data for Savoy et al. (2019): Metabolic rhythms in flowing waters: An approach for classifying river productivity regimes

This document describes several of the derived datasets used in Savoy et al. (2019) as well as Koenig et al. (in review). The analysis presented in Savoy et al. (2019) describes identifying similar characteristic regimes of gross primary productivity (GPP) across 47 U.S. streams and rivers through the use of clustering analysis. This resource contains the following files:

1. **site basic.csv:** A table containing basic site information and clustering assignments for each site
2. **avg_gpp.csv:** Representative time series of GPP for each of the 47 rivers
3. **avg_gpp_filled.csv:** Representative time series of GPP that has been gap-filled for each of the 47 rivers
4. **normalized.csv:** Representative z-normalized time series of GPP for each of the 47 rivers

Below, each of these datasets is described individually including detailed information about the data contained within each file and the derivation of this data.

#### 1. site_basic.csv

This file contains the following columns:

* ***Site_ID:*** The National Water Information System (NWIS) unique identifier for each site. Sites represent USGS gauged sites, and each Site_ID corresponds to a location in the dataset available from Appling et al. (2018b).
* ***Site_name:*** The full site name for each location
* ***Lat:*** Latitude (decimal degrees, NAD83)
* ***Lon:*** Longitude (decimal degrees, NAD83)
* ***WS_area:*** Watershed area (Km<sup>2</sup>)
* ***Width:*** Channel width (m) derived from regional hydraulic geometry coefficients (Gomez-Velez et al. 2015)
* ***two_clus:*** The cluster each site was assigned to for a two cluster solution (summer peak, spring peak)
* ***four_clus:*** The cluster each site was assigned to for a four cluster solution (summer peak, spring peak, summer decline, aseasonal)

The columns ***two_clus*** and ***four_clus***are the clustering assignments for each site based on Savoy et al. (2019). Here, a brief overview of this clustering analysis is provided. Clustering identifies groups based on a measure of dissimilarity, and then this dissimilarity is used to identify a classification. Dynamic time warping (DTW) (Sakoe and Chiba 1978; Berndt and Clifford 1994) was used to define the similarity between time series because of its widespread application in time series analysis. The DTW dissimilarity matrix was then used to perform a hierarchical agglomerative clustering.  We defined between two and ten clusters because no a priori optimal number of clusters exists. A suite of indices that effectively describes a combination of within cluster cohesion and between cluster separation was used to assess the clusters and determine a final set of clusters. The results presented within Savoy et al. (2019) focus primarily on a two clustering solution as the most conservative and robust result across several different clustering methodologies. These clusters were named *Summer Peak Rivers* and *Spring Peak Rivers* based on the temporal patterns of the resulting representative GPP regimes. However, several results indicated the possibility of a less conservative solution of four clusters and these are also presented in the paper. These clusters were named *Summer Peak Rivers*, *Spring Peak Rivers*, *Summer Decline Rivers*, and *Aseasonal* based on the temporal patterns of the resulting representative GPP regimes. Because of this, both the two cluster solution (***two_clus***) and four cluster solution (***four_clus***) are provided for each site. For a full description of the methods used to derive these clusters and the interpretation of these results please refer to Savoy et al. (2019).

#### 2. avg_gpp.csv

This file contains representative time series of GPP for each of the 47 sites used. This data is derived from a subset of daily estimates of stream metabolism described in Appling et al. (2018a) and the original data are freely available to download (Appling et al. 2018b) and full descriptions of the original datasets can be found within these sources. The original set of 356 sites was filtered based on a combination of data quality and coverage to select a subset of 47 rivers that all had data for the time period of 2013-2016. This file consists of the a mean time series of GPP for each site that was calculated by taking the mean GPP for each day of the year across all four years of data. The first column (***DOY***) is the day of year and each subsequent column corresponds to a specific ***Site_ID***.

#### 3. avg_gpp_filled.csv

This file contains representative time series of GPP for each of the 47 sites used; however, the time series of GPP has been gap-filled. To create these series the original daily GPP estimates were gap-filled using a generalized additive model with both seasonal and trend components. These gap-filled series were then used to calculate the mean time series of GPP for each site that was calculated by taking the mean GPP for each day of the year across all four years of data. The first column (***DOY***) is the day of year and each subsequent column corresponds to a specific ***Site_ID***.

#### 4. normalized.csv

To calculate similarity with DTW it is necessary to z-normalize each time series. The representative time series of gap-filled GPP as described above were thus z-normalized. The first column (***DOY***) is the day of year and each subsequent column corresponds to a specific ***Site_ID***.

#### 5. Metadata (LO_letters).pdf

An accompanying set of metadata using the format from Limnology & Oceanography letters.  Note, this metadata largely reiterates the information covered in this readme file but is provided as a separate resource to conform with journal data policy guidelines.

### References

Appling, A. P., and others. 2018a. The metabolic regimes of 356 rivers in the United States. Sci. Data 5: 180292. https://doi.org/10.1038/sdata.2018.292

Appling, A.P., Read, J.S., Winslow, L.A., Arroita, M., Bernhardt, E.S., Griffiths, N.A., Hall, R.O., Jr., Harvey, J.W., Heffernan, J.B., Stanley, E.H., Stets, E.G., and Yackulic, C.B. 2018b, Metabolism estimates for 356 U.S. rivers (2007-2017): U.S. Geological Survey data release. https://doi.org/10.5066/F70864KX

Berndt, D. J., and J. Clifford. 1994. Using dynamic time warping to find patterns in time series. AAAI technical report WS-94-03. Association for the Advancement of Artificial Intelligence.

Gomez-Velez, J. D., J.W. Harvey, M. B. Cardenas, and B. Kiel. 2015. Denitrification in the Mississippi River network controlled by flow through river bedforms. Nature Geoscience 8:941. [10.1038/ngeo2567](http://dx.doi.org/10.1038/ngeo2567)

Sakoe, H., and S. Chiba. 1978. Dynamic programming algorithm optimization for spoken word recognition. IEEE Trans. Acoust. Speech Signal Process. 26:43–49. https://doi.org/10.1109/TASSP.1978.1163055

Savoy, P. , Appling, A. P., Heffernan, J. B., Stets, E. G., Read, J. S., Harvey, J. W. and Bernhardt, E. S. 2019. Metabolic rhythms in flowing waters: An approach for classifying river productivity regimes. Limnol Oceanogr. <https://doi.org/10.1002/lno.11154>

