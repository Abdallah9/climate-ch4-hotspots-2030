"""
Climate Data Processor Module
Processor for non-ERA5 climate datasets (CH4, permafrost, wetlands, etc.)
Author: Your Name
Date: June 2025
"""

#     Process climate time series data with existing time dimensions
#     For datasets that already have annual time series but need spatial harmonization
    
#     Parameters:
#     -----------
#     file_path : str
#         Path to NetCDF file
#     variable_name : str
#         Name of variable in NetCDF (e.g., 'ch4_concentration', 'permafrost_zones')
#     output_name : str
#         Name for output variable in CSV (e.g., 'ch4_concentration', 'permafrost_zones')
#     unit_conversion : function, optional
#         Function to convert units if needed
#     description : str, optional
#         Description for logging
#        
#     Returns:
#     --------
#     pandas.DataFrame
#         Processed data in ML-ready format matching other datasets


import xarray as xr
import numpy as np
import pandas as pd

def process_climate_timeseries(
    file_path, 
    variable_name, 
    output_name, 
    unit_conversion=None,
    description=""
):

    print(f"=== PROCESSING {description.upper()} ===")
    print(f"File: {file_path}")
    print(f"Variable: {variable_name} → {output_name}")
    
    # Step 1: Load and examine data
    print("\n1. Loading data...")
    ds = xr.open_dataset(file_path)
    # Auto-fix: swap axes if lat/lon are mislabeled
    if 'latitude' in ds.coords and 'longitude' in ds.coords:
        lat_vals = ds['latitude'].values
        lon_vals = ds['longitude'].values
        if abs(lat_vals).max() > 90:
            print("⚠️ WARNING: Latitude values exceed ±90° — likely mislabeled. Swapping 'latitude' and 'longitude'.")
            ds = ds.rename({'latitude': 'longitude', 'longitude': 'latitude'})

    print(f"   Dimensions: {dict(ds.dims)}")

    # ✅ Auto-detect time coordinate
    # time_coord = None
    # for coord in ds.coords:
    #     if 'time' in coord.lower():
    #         time_coord = coord
    #         break

    # if not time_coord:
    #     raise ValueError("Time coordinate not found in dataset.")

    # print(f"   Time range: {ds[time_coord].min().values} to {ds[time_coord].max().values}")
    # print(f"   Time steps: {len(ds[time_coord])}")


    # ✅ Auto-detect time-like coordinate
    possible_time_coords = ['time', 'year', 'years']
    time_coord = next((c for c in possible_time_coords if c in ds.coords or c in ds.dims), None)
    
    if not time_coord:
        raise ValueError("No time-like coordinate (e.g., 'time', 'year') found in dataset.")
        
    print(f"   Time coordinate: {time_coord}")
    print(f"   Time range: {ds[time_coord].min().values} to {ds[time_coord].max().values}")
    print(f"   Time steps: {len(ds[time_coord])}")


    # Step 2: Get and clean the target variable
    data_var = ds[variable_name]
    # ✅ Mask absurd fill values (e.g., corrupted CH₄ emissions)
    data_var = data_var.where((data_var >= 0) & (data_var < 1e6))
    print(f"   Variable shape: {data_var.shape}")
    print(f"   Variable dimensions: {data_var.dims}")


    # Step 3: Identify lat/lon coordinates
    print("\n2. Processing coordinates...")
    
    # Try known common names first
    known_lat_names = ['lat', 'latitude']
    known_lon_names = ['lon', 'longitude']
    
    lat_coord = next((c for c in data_var.dims if c.lower() in known_lat_names), None)
    lon_coord = next((c for c in data_var.dims if c.lower() in known_lon_names), None)
    
    # Fallback: look in dataset if not found in variable
    if lat_coord is None:
        lat_coord = next((c for c in ds.coords if c.lower() in known_lat_names), None)
    if lon_coord is None:
        lon_coord = next((c for c in ds.coords if c.lower() in known_lon_names), None)
    
    if lat_coord is None or lon_coord is None:
        raise ValueError(f"Could not detect latitude/longitude coordinates in dataset. Found dims: {list(ds.dims)}")
    
    print(f"   Latitude coordinate: {lat_coord}")
    print(f"   Longitude coordinate: {lon_coord}")

    print("DEBUG lat/lon names:", lat_coord, lon_coord)
    print("DEBUG data_var.coords:", list(data_var.coords))


    # Step 4: Resample to 0.1° grid
    print("\n3. Resampling to 0.1° resolution...")
    target_lats = np.arange(42.0, 82.6, 0.1)
    target_lons = np.arange(-141.0, -52.4, 0.1)
    print(f"   Target grid: {len(target_lats)} lats × {len(target_lons)} lons")

    resampled_data = data_var.interp(
        {lat_coord: target_lats, lon_coord: target_lons},
        method='nearest'
    )
    print(f"   Resampled shape: {resampled_data.shape}")
    print("   Resampled value range:", float(resampled_data.min()), "to", float(resampled_data.max()))

       # Step 5: Extract to CSV
    print("\n4. Extracting to CSV...")
    
    lons = resampled_data[lon_coord].values
    lats = resampled_data[lat_coord].values
    times = resampled_data[time_coord].values
    
    # ✅ Robust year conversion
    if np.issubdtype(times.dtype, np.signedinteger):
        years = times
    else:
        years = pd.to_datetime(times).year
    
    print(f"   Years: {years.min()} to {years.max()}")

    lon_2d, lat_2d = np.meshgrid(lons, lats)
    first_time_data = resampled_data.isel({time_coord: 0}).values
    valid_mask = ~np.isnan(first_time_data.flatten())

    valid_lats = lon_2d.flatten()[valid_mask]
    valid_lons = lat_2d.flatten()[valid_mask]
    n_valid_pixels = len(valid_lats)
    print(f"   Valid pixels: {n_valid_pixels:,}")

    all_records = []
    for i, year in enumerate(years):
        year_data = resampled_data.isel({time_coord: i}).values.flatten()[valid_mask]
        if unit_conversion:
            year_data = unit_conversion(year_data)
        year_df = pd.DataFrame({
            'pixel_id': range(n_valid_pixels),
            'longitude': valid_lons,
            'latitude': valid_lats,
            'year': year,
            output_name: year_data
        })
        all_records.append(year_df)

    final_df = pd.concat(all_records, ignore_index=True)
    final_df = final_df.sort_values(['year', 'pixel_id']).reset_index(drop=True)

    final_df['longitude'] = final_df['longitude'].round(6)
    final_df['latitude'] = final_df['latitude'].round(5)
    final_df[output_name] = final_df[output_name].round(7)

    # Step 6: Save
    output_file = f'{output_name.title()}_2010-2024.csv'
    final_df.to_csv(output_file, index=False)

    print(f"\n✅ {description.upper()} COMPLETE!")
    print(f"   File: {output_file}")
    print(f"   Records: {len(final_df):,}")
    print(f"   Unique pixels: {final_df['pixel_id'].nunique():,}")
    print(f"   Years: {final_df['year'].min()}-{final_df['year'].max()}")
    print(f"   Value range: {final_df[output_name].min():.3f} to {final_df[output_name].max():.3f}")

    print(f"\n   Validation:")
    print(f"   - Missing values: {final_df.isnull().sum().sum()}")
    print(f"   - Records per year: {final_df['year'].value_counts().iloc[0]:,}")
    print(f"   - Consistent pixel count: {final_df['pixel_id'].nunique() == n_valid_pixels}")

    return final_df
    

def process_ch4_concentration(file_path):
    """Process CH4 concentration data (target variable)"""
    return process_climate_timeseries(
        file_path=file_path,
        variable_name='ch4_concentration',
        output_name='ch4_concentration',
        unit_conversion=None,  # ✅ FIXED HERE
        description='CH4 Concentration Data (converted to ppm)'
    )


def process_permafrost(file_path):
    return process_climate_timeseries(
        file_path=file_path,
        variable_name='permafrost_fraction',
        output_name='permafrost_fraction',
        unit_conversion=lambda x: x / 100.0,  # ✅ Convert % to 0–1
        description='Permafrost Fraction'
    )


def process_wetlands(file_path):
    """Process wetland fraction data"""
    return process_climate_timeseries(
        file_path=file_path,
        variable_name='wetland_fraction',
        output_name='wetland_fraction',
        unit_conversion=lambda x: x * 100,  # Convert fraction to percentage
        description='Wetland Fraction Data'
    )

def process_industrial_emissions(file_path):
    """Process industrial CH4 emissions data"""
    return process_climate_timeseries(
        file_path=file_path,
        variable_name='fuel_emi',
        output_name='ch4_emissions',
        unit_conversion=None,  # Already in kg/km²/year
        description='Industrial CH4 Emissions Data',
    )


# Update the batch processing function too:
def process_all_climate_data(ch4_file, permafrost_file, wetlands_file, 
                           emissions_file=None, landcover_file=None, suffix=""):
    """Process all climate datasets at once"""
    
    print("="*80)
    print("BATCH PROCESSING ALL CLIMATE DATASETS")
    print("="*80)
    
    results = {}
    
    # Core datasets
    results['ch4'] = process_ch4_concentration(ch4_file, suffix)
    results['permafrost_zones'] = process_permafrost_zones(permafrost_file, suffix)
    results['permafrost_extent'] = process_permafrost_extent(permafrost_file, suffix)
    results['wetlands'] = process_wetlands(wetlands_file, suffix)
    
    # Optional additional datasets
    if emissions_file:
        results['industrial_emissions'] = process_industrial_emissions(emissions_file, suffix)
    
    if landcover_file:
        results['land_cover'] = process_land_cover(landcover_file, suffix)
    
    print(f"\n{'='*80}")
    print("BATCH PROCESSING COMPLETE!")
    print("="*80)
    
    # Summary
    for name, df in results.items():
        print(f"{name.upper()}: {len(df):,} records, {df['pixel_id'].nunique():,} pixels")
    
    return results

