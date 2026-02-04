# Kansas SSC Geospatial Context Pipeline

This repository contains a fully reproducible, modular Python workflow for building **watershed-scale geospatial context** for suspended sediment concentration (SSC) monitoring locations in Kansas, USA.

The pipeline integrates **water-quality monitoring locations**, **hydrologic unit boundaries**, **geometric watershed attributes**, and **climate covariates** into a single, analysis-ready geospatial dataset. The primary goal is to prepare consistent, watershed-scale geospatial features that support applied geospatial water-quality analysis, modeling, and decision-oriented workflows, while adhering to best practices in geospatial data acquisition, processing, and feature engineering.

---

## Project Overview

Accurate interpretation and modeling of riverine suspended sediment often requires **spatial context** beyond point measurements. This project prepares that context by:

- Discovering SSC monitoring locations from the Water Quality Portal (WQP)
- Assigning each site to a WBD HUC12 watershed
- Computing watershed geometry metrics in a projected CRS
- Extracting mean annual precipitation from PRISM rasters via zonal statistics
- Producing clean, reproducible geospatial outputs suitable for modeling or exploratory analysis

This work is designed as **data preparation and infrastructure**, not hypothesis testing.

---

## Data Sources

All data used in this repository are publicly available and programmatically retrieved:

- **Water Quality Portal (WQP)**  
  Monitoring locations with SSC (USGS pCode `80154`)  
  https://www.waterqualitydata.us/

- **Watershed Boundary Dataset (WBD)**  
  HUC12 hydrologic unit polygons  
  https://www.usgs.gov/national-hydrography/watershed-boundary-dataset

- **PRISM Climate Group**  
  Annual precipitation rasters (4 km resolution)  
  https://prism.oregonstate.edu/

- **U.S. Census TIGER/Cartographic Boundaries**  
  Kansas state boundary  
  https://www.census.gov/geographies/mapping-files/time-series/geo/cartographic-boundary.html

---

## Repository Structure

```text
.
├── README.md
├── requirements.txt
├── data/                  # Raw and intermediate data (not versioned)
├── outputs/               # Figures (optional)
└── src/
    ├── st0_config.py
    ├── st1_get_ks_boundary.py
    ├── st2_get_wqp_ssc_stations.py
    ├── st3_get_wbd_huc12.py
    ├── st4_join_sites_huc12.py
    ├── st5_geometry_metrics.py
    ├── st6_prism_download.py
    ├── st7_zonal_stats.py
    └── st8_generate_figures.py
