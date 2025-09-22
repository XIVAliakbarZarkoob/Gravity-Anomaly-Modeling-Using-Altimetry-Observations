# Gravity Anomaly Modeling Using Altimetry Observations

This repository contains the implementation and results of **gravity
anomaly modeling** using **satellite altimetry data**. The project
leverages altimetry-derived **Sea Surface Height (SSH)** and **Dynamic
Ocean Topography (DOT)** to compute geoid heights and estimate gravity
anomalies using the **Inverse Stokes Formula (ISF)**.

## Overview

Gravity anomaly modeling provides insights into Earth's internal mass
distribution and density variations. Altimetry data, obtained from
satellites like **SARAL**, are a valuable resource for regions with
limited direct measurements. This project:
- Collects and processes SSH and DOT data from **Open Altimeter Database
(OpenADB)**.
- Computes **geoid heights** on a 0.25° grid using SSH and DOT
differences.
- Estimates **gravity anomalies** with ISF.
- Validates results against the **XGM2019** global gravity model.

## Study Area

-   **Longitude:** 326° to 360°
-   **Latitude:** 0° to 33°
-   **Data period:** First two months of 2015

## Methodology

1.  **Data Collection** -- SSH and DOT retrieved and interpolated on a
    regular grid.

2.  **Geoid Height Calculation** -- Geoid heights computed as
    `N = SSH − DOT`.

3.  **Gravity Anomaly Calculation** -- ISF applied in discrete form:

4.  **Validation** -- Comparison with **XGM2019** anomalies.



