# 📄 README Documentation Draft
   
# climate-ch4-hotspots-2030

## **Objective:**
The project aims to develop a machine learning–driven framework to identify and forecast methane (CH₄) emission hotspots toward 2030. By integrating satellite-based methane concentration data with relevant environmental, industrial, and climatic drivers, the study will provide insights into the spatial distribution and future risks of CH₄ emissions.


## **Research Novelty**

**Research Title:** *“2030 Satellite-Based Identification of Global Climate-Driven Methane (CH₄) Emission Hotspots and Their Policy Implications”*

The novelty of this research lies in its integration of multiple dimensions:  
* Utilizing satellite-based data  
* Predicting global methane (CH₄) emission hotspots  
* Producing projections up to the year 2030  
* Linking scientific findings with potential **policy implications** for mitigation and adaptation strategies  

### 📚 Literature Review:

1. **Global Methane Assessment: 2030 Baseline Report**

   * **Summary:** Provides projections of anthropogenic methane emissions through 2030 under various baseline scenarios and assesses the climate benefits of achieving the Global Methane Pledge target.
   * **Key Insight:** While it offers projections, it doesn't focus on satellite-based identification of emission hotspots or detailed policy implications.
   * **Reference:** [Global Methane Assessment: 2030 Baseline Report](https://www.ccacoalition.org/resources/global-methane-assessment-2030-baseline-report)

2. **Automated Detection and Monitoring of Methane Super-Emitters Using Satellite Data**

   * **Summary:** Discusses an automated system for detecting and monitoring methane super-emitters using satellite data.
   * **Key Insight:** Focuses on current emissions rather than future projections or policy implications.
   * **Reference:** [Automated Detection and Monitoring of Methane Super-Emitters](https://www.ghgsat.com/en/scientific-publications/automated-detection-and-monitoring-of-methane-super-emitters-using-satellite-data/)

3. **Deep Learning-Based Quantifications of Methane Emissions with Multispectral Satellite Data**

   * **Summary:** Explores the use of deep learning models to quantify methane emissions using multispectral satellite data.
   * **Key Insight:** Concentrates on current emission quantification without projecting future hotspots or discussing policy implications.
   * **Reference:** [Deep Learning-Based Quantifications of Methane Emissions](https://www.sciencedirect.com/science/article/pii/S1569843224003728)

4. **CH4Net: A Deep Learning Model for Monitoring Methane Super-Emitters**

   * **Summary:** Introduces a deep learning model for automated monitoring of methane super-emitters from Sentinel-2 data.
   * **Key Insight:** Aims at current monitoring rather than future projections or policy analysis.
   * **Reference:** [CH4Net: A Deep Learning Model](https://amt.copernicus.org/articles/17/2583/2024/)

### 🧠 Analysis:

* **Current Research Focus:** Existing studies predominantly concentrate on current methane emissions, detection of super-emitters, and the development of monitoring tools using satellite data.
* **Gap Identified:** There is a lack of research that combines satellite-based data to project global methane emission hotspots up to 2030 with an emphasis on climate-driven factors and policy implications.

## **Scope & Limitations**

### 🔭 Scope
This project focuses on building a machine learning–based framework for forecasting methane (CH₄) emission hotspots toward 2030. The scope includes:
- Processing **satellite-derived CH₄ concentration data** (2010–2024)  
- Integrating **environmental drivers** (temperature, precipitation, soil moisture, permafrost, wetland fraction)  
- Incorporating **industrial drivers** (fuel exploitation CH₄ emissions, land cover/land use changes)  
- Generating **annual, Canada-focused ML-ready datasets** at 0.1° × 0.1° resolution  
- Producing **forecasts up to 2030** using statistical and machine learning models  

### ⚠️ Limitations
While the study provides valuable insights, several constraints must be acknowledged:
- **Geographic focus**: Initial implementation targets Canada, with possible extension to global scale later.  
- **Temporal resolution**: Data are aggregated to **annual means**, which may smooth short-term emission spikes.  
- **Drivers considered**: Focused on key climate and land-use factors; other potential influences (e.g., socio-economic policies, real-time industrial incidents) are not explicitly modeled.  
- **Forecast horizon**: Projections are limited to **2030**, balancing data availability and model uncertainty.  
- **Policy analysis**: The study provides **science-based insights** for policymakers but does not perform in-depth economic or legal policy modeling.  

---

## **Identifying Required Data and Sources**

### **a) Required Data**

To effectively model and forecast methane (CH₄) emission hotspots toward 2030, we require both the target variable (methane concentration) and the causal factors that influence emissions.

### **1. Target Variable**

* **Historical Methane Concentration (2010–2024)**
  *Satellite-based atmospheric CH₄ concentration data.*

  * Provides the dependent variable for hotspot prediction.
  * Harmonized to annual means over Canada.

### **2. Causal Factors Influencing Emissions**

* **Industrial Activities (Anthropogenic Sources):**

  * Agriculture CH₄ emissions (EDGAR v8.0).
  * Fuel exploitation (oil & gas sector, EDGAR v8.0).
  * Waste-related emissions (EDGAR v8.0).

* **Climate Variables (ERA5):**

  * Temperature (K → °C).
  * Precipitation (m → mm).
  * Soil Moisture (m³/m³)

* **Land Cover & Land Use (ESA-CCI):**

  * Annual/static maps of vegetation, croplands, forests, etc.
  * Replicated across years for consistency.

* **Wetlands (GIEMS):**

  * Fractional inundation extent dataset.
  * Extended to 2024 using AR(2) forecasting.

* **Permafrost (NSIDC / ESA CCI):**

  * Frozen ground extent, harmonized and extended to 2024.

* **Elevation (ERA5 Geopotential):**

  * Static topography (m²/s² → meters).
  * Provides terrain context for emissions distribution.

---

## Coverage Check

This list **fully reflects all needed features**:

* CH₄ concentration (target)
* Industrial emissions (agriculture, fuel exploitation, waste)
* Climate (temperature, precipitation, soil moisture)
* Land cover
* Wetlands
* Permafrost
* Elevation


### **b) Data Sources**

All data used in this project were downloaded **manually via official web portals through registration** . 
Every source file is in **NetCDF (`.nc`)** format. 

---

## Feature-by-Feature Download Details

### CH₄ Concentration

* **Portal:** Copernicus Climate Data Store
* **Data Page URL:** [ERA5 hourly data on single levels](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=overview) ([Climate Data Store][1])
* **Steps:**

  1. Log in, search *“ERA5 hourly data on single levels”*.
  2. Choose variable **2 m temperature (incorrect)** — but for methane, use methane dataset.
     (Note: the methane data source is inside the *Satellite Methane* section, not ERA5.)
* **Accurate CH₄ dataset URL:** [Satellite Methane (C3S)](https://cds.climate.copernicus.eu/datasets/satellite-methane)
* Download CH₄ data as NetCDF for 2010–2024.

---

### Industrial CH₄ Emissions (Fuel Exploitation)

* **Portal:** EDGAR v8.0 - JRC GHG Emissions
* **URL:** [EDGAR v8.0](https://edgar.jrc.ec.europa.eu/dataset_ghg80) ([ecmwf-projects.github.io][2], [Nature][3], [edgar.jrc.ec.europa.eu][4])
* **Steps:**

  1. Navigate to CH₄ emissions for Fuel Exploitation.
  2. Select NetCDF format; download data for 2010–2022.

---

### Permafrost Fraction

* **Portal:** ESA CCI Permafrost via CEDA
* **URL:** [ESA CCI Permafrost v4](https://catalogue.ceda.ac.uk/uuid/93444bc1c4364a59869e004bf9bfd94a/) ([Data.europa.eu][5], [ecmwf-projects.github.io][2])
* **Steps:** Download annual NetCDFs (2010–2021), later extended to 2024.

---

### Wetland Fraction (GIEMS-MC)

* **Portal:** ORNL DAAC (GIEMS-D3) via search
* **URL:** [ORNL DAAC GIEMS-D3](https://daac.ornl.gov/) *(navigate to Wetlands GIEMS dataset)*
* Download one comprehensive NetCDF (`GIEMS‑MC_compressed.nc`) which contains wetlands and land cover variables.

---

### Land Cover (Dominant Class)

* **Portal:** Same GIEMS file as wetlands
* **Source:** `dom_land_cover_class` variable extracted from already downloaded `GIEMS‑MC_compressed.nc`.

---

### Temperature & Precipitation (ERA5)

* **Portal:** Copernicus CDS
* **URL:** [ERA5 hourly data on single levels (Temperature, Precipitation)](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=overview) ([Climate Data Store][1])
* Download `t2m` (temperature) and `tp` (precipitation) variables as NetCDF for 2010–2024.

---

### Soil Moisture (ERA5)

* **Portal:** Copernicus CDS (same as above)
* Same dataset page — select `swvl1` (volumetric soil water layer 1) for 2010–2024 as NetCDF.

---

### Elevation (Geopotential)

* **Portal:** Copernicus CDS (same as above)
* Download `geopotential` variable (single timestamp) as NetCDF.

---


## **File Inventory**

### **Input Files:**
- `data_stream-oper_stepType-instant.nc` (Temperature)
- `data_stream-oper_stepType-accum.nc` (Precipitation)  
- `era5_geopotential_canada.nc` (Elevation)
- `_data_stream-oper_stepType-instant.nc'` (Soil_moisture)

###  **Preprocess file:**
- `v8.0_FT2022_GHG_CH4_2000_FUEL_EXPLOITATION_emi.nc` → `CH4_FE_Canada_Annual_2010_2024.nc` (CH₄ Emissions)
- `ESACCI-PERMAFROST-L4-PFR-MODISLST_CRYOGRID-AREA4_PP-2010-fv04.0.nc` → `PERMAFROST_Canada_Annual_2010_2024.nc` (Permafrost)
- `data_sfc_Ch4.nc` → `data_sfc_Ch4_annual.nc` (CH₄_Concentration)
- `GIEMS-MC_compressed.nc` → `WETLAND_Canada_Annual_2010_2024.nc` (Wetland)
- `GIEMS-MC_compressed.nc` → `LandCover_Canada_Static.nc` (LandCover)

---

# **Data Harmonization and Preprocessing**

This step ensured that all raw datasets, obtained from multiple sources (Copernicus, ESA CCI, EDGAR, etc.), were standardized into a **uniform spatiotemporal framework** before being fed into the modeling pipeline.

The preprocessing was carried out in two stages:

---

## **1. Raw Preprocessing (Per Dataset)**

Each dataset required **initial cleaning and subsetting** directly after download:

* **File Loading & Inspection**
  Opened original NetCDF files (`.nc`) and checked dimensions, variables, and metadata.

* **Subsetting to Canada**
  All datasets were clipped to Canada’s geographic extent:

  * Latitude: **40°N – 85°N**
  * Longitude: **-141°W – -52°W**

* **Variable Selection & Renaming**
  Selected only the variables of interest (e.g., `tc_ch4 → ch4_concentration`) and renamed them for consistency across datasets.

* **Temporal Aggregation**
  Where needed, higher-frequency data (e.g., hourly ERA5) were resampled to **annual means (2010–2024)**.

* **Static Features (Land Cover, Elevation)**
  These datasets were treated as static and replicated across all years.

* **Forecasting Missing Years (2023–2024)**
  For datasets available only until 2021/2022 (e.g., Permafrost, Wetlands, Emissions), we applied a **lightweight AR(2) model** to extend time series to 2024:

  $$
  y_t = a \cdot y_{t-1} + b \cdot y_{t-2}
  $$

---

## **2. Harmonization Across All Features**

After raw preprocessing, all datasets were standardized to meet the **benchmark requirements**:

* **Spatial Resolution**
  Interpolated to a uniform **0.1° × 0.1° (\~10 km)** grid.

* **Coordinate System**
  All datasets confirmed in **EPSG:4326 (WGS84)**.

* **Temporal Coverage**
  Restricted to **2010–2024** annual values.

* **Valid Pixel Filtering**
  Excluded NaNs and invalid values (e.g., negative soil moisture, out-of-range concentration).

* **ML-Ready CSV Export**
  Each dataset was flattened into tabular format with the following structure:

  ```
  pixel_id, latitude, longitude, year, <feature_value>
  ```
  Example: `CH4_Concentration_2010-2024.csv`

---

## **3. Final Output Summary**

By the end of harmonization, all features were fully aligned in **space and time**, and stored as **ML-ready CSVs** for direct use in the methane hotspot prediction models.

Features include:

* Methane concentration (`ch4_concentration`)
* Soil moisture (`soil_moisture`)
* Precipitation (`precipitation`)
* Temperature (`temperature`)
* Wetland fraction (`wetland_fraction`)
* Permafrost fraction (`permafrost_fraction`)
* Industrial CH₄ emissions (`ch4_emissions`)
* Land cover (`land_cover_class`)
* Elevation (`elevation`)

  
---

## **Key Code Development**

After successfully preprocessing some of the Downloaded data using the folowing script,
Ch4_raw_process.ipynb
Fuel_Exploitation_emi_raw_ process.ipynb
land_cover_raw_ process.ipynb
permafrost_raw_ process.ipynb
wetland_raw_process.ipynb,


the obtained data was fed to era5_processor.py and climate_processor.py modules.


### **Reusable Modules**

* `era5_processor.py` – 
* `climate_processor.py` – 
---

### **Main Functions**

```python
# Generic (works for ERA5 + non-ERA5 NetCDFs)
process_climate_timeseries(
    file_path: str,
    variable_name: str,
    output_name: str,
    unit_conversion=None,     
    description: str = ""
) -> pd.DataFrame
```


```python
# ERA5 convenience variant, same interface, adds ERA5-specific time/coord fixes
process_era5_timeseries(
    file_path: str,
    variable_name: str,
    output_name: str,
    unit_conversion=None,
    description: str = ""
) -> pd.DataFrame
```

* Loads a NetCDF, converts longitudes (0–360 → -180–180), **clips to Canada**,
* **aggregates to annual** (mean or sum—see `agg` below),
* **resamples to 0.1°**, flattens to ML-ready long table,
* filters years **2010–2024**, saves `*feature_Name_2010-2024.csv`.


---

### **Convenience Functions (wrappers used)**

```python
# ERA5
process_temperature(file_path)          # t2m → temperature (K→°C, annual mean)
process_precipitation(file_path)        # tp  → precipitation (m→mm, annual sum)
process_soil_moisture(file_path)        # swvl1 → soil_moisture (fraction, annual mean)

# Non-ERA5 / other sources
process_ch4_concentration(file_path)    # tc_ch4 → ch4_concentration (ppm or fraction→ppm)
process_permafrost(file_path)           # PFR → permafrost_fraction (annual; AR(2) for 2022–2024 if needed)
process_wetlands(file_path)             # inund_sat_wetland_frac → wetland_fraction (%), monthly→annual + AR(2)
process_land_cover(file_path)           # dom_land_cover_class (static, replicate 2010–2024)
process_elevation(file_path)            # z (geopotential) → elevation (m), static, replicate 2010–2024
process_industrial_emissions(file_path) # fuel_emi/emissions → ch4_emissions (annual; AR(2) for 2023–2024)
```

Each wrapper:

* Passes the correct `variable_name`/`output_name`
* Chooses **aggregation** (mean vs sum)
* Supplies the **unit conversion** (if any)
* Applies the **year filter (2010–2024)**

---

### **Unit Conversions (used in wrappers)**

```python
kelvin_to_celsius(temp_k)          # K → °C
meters_to_millimeters(precip_m)    # m → mm
geopotential_to_elevation(geo)     # m²/s² → m  (geo / 9.80665)
fraction_to_percent(x)             # 0–1 → 0–100 (wetlands)
```

---

### **Common Utilities**

```python
# 1) Coordinate normalization
fix_longitudes_to_west(lon)  # (lon + 180) % 360 - 180

# 2) Spatial subsetting to Canada (benchmark)
CANADA_BOUNDS = dict(lat_min=40.0, lat_max=85.0, lon_min=-145.0, lon_max=-50.0)
subset_canada(ds, bounds=CANADA_BOUNDS)

# 3) Temporal aggregation
annual_aggregate(da, mode="mean")  # "mean" for t2m/swvl1/wetland/conc; "sum" for tp

# 4) Regrid to benchmark grid
TARGET_GRID = (np.arange(40.0, 85.1, 0.1), np.arange(-145.0, -49.9, 0.1))
regrid_to_0p1(da, method="linear")     # continuous fields
regrid_to_0p1_nearest(da)              # categorical (land cover)

# 5) Year filter
filter_years(da_or_df, start=2010, end=2024)

# 6) Flatten to ML long table (consistent schema)
to_long_dataframe(da, var_name)  # -> ['pixel_id','longitude','latitude','year', var_name]

# 7) IO helpers
save_csv(df, out_name)
save_netcdf(da, out_nc)
```

---

### **Usage Examples:**
```python
from era5_processor import process_temperature, process_precipitation

temp_df = process_temperature('data_stream-oper_stepType-instant.nc')
precip_df = process_precipitation('data_stream-oper_stepType-accum.nc')
```
---
 
### **Output Naming (consistent across features)**
* **CSV**: `{OutputName}_2010-2024.csv`
  
  ### **Output Files:**
- `Ch4_Concentration_2010–2024.csv`   
- `Ch4_Emissions_2010–2024.csv`       
- `Elevation_2010–2024.csv`           
- `LandCover_2010–2024.csv`           
- `Permafrost_Fraction_2010–2024.csv`
- `Precipitation_2010–2024.csv`       
- `Soil_Moisture_2010–2024.csv`       
- `Temperature_2010–2024.csv`         
- `Wetland_Fraction_2010–2024.csv`    


---

# **Data Unification**

## **🎯 Goal**

The objective of the unification step was to merge all processed climate features into a single dataset covering **2010–2024**, standardized to:

* **Resolution:** 0.1° × 0.1° (\~10 km)
* **Projection:** EPSG:4326 (WGS84)
* **Temporal:** Annual means
* **Spatial:** Canada boundaries
* **Format:** ML-ready tabular CSV

This unified dataset would serve as the foundation for machine learning models predicting CH₄ emission hotspots.

---

## **🛠 What Was Tried**

To guarantee alignment across all features, I first attempted to compute the **strict spatial and temporal intersection**:

1. Loaded each feature (`ch4_concentration`, `emissions`, `temperature`, `precipitation`, `soil_moisture`, `elevation`, `land_cover`, `permafrost_fraction`, `wetland_fraction`).
2. Identified the set of unique `(pixel_id, year)` pairs in each dataset.
3. Calculated the **intersection** across all features (i.e., only pixels present everywhere).

This approach was designed to ensure a dataset with **no missing values**.

---

## **⚠️ What Happened**

The strict intersection revealed a critical limitation:

* **CH₄ concentration:** \~358k unique pixels
* **ERA5-based features (temperature, precipitation, soil moisture):** \~428k pixels
* **Wetlands:** only \~6,415 unique pixels

When intersected, the dataset shrank to **6,415 pixels**, or less than **2% of ERA5 coverage**.

This meant that the wetlands dataset—being spatially limited—acted as the bottleneck, forcing a drastic reduction in usable data for national-scale modeling.


---

**Decision:** Wetland fraction was excluded from the unified v1 dataset because its limited spatial footprint forced the strict intersection to \~6.4k pixels. We retain national coverage by excluding wetlands and adding a simple wetness proxy from land cover. A wetlands-enhanced variant will be explored as a follow-up model for wetland-dominant regions.

---



































---























































