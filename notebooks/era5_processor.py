"""
ERA5 Data Processor Module
Universal processor for ERA5 time series data (temperature, precipitation, etc.)
Author: Your Name
Date: June 2025
"""

import xarray as xr
import numpy as np
import pandas as pd

def process_era5_timeseries(
    file_path, 
    variable_name, 
    output_name, 
    unit_conversion=None,
    description=""
):
    """
    Universal processor for ERA5 time series data (temperature, precipitation, etc.)
    
    Parameters:
    -----------
    file_path : str
        Path to NetCDF file
    variable_name : str
        Name of variable in NetCDF (e.g., 't2m', 'tp')
    output_name : str
        Name for output variable in CSV (e.g., 'temperature', 'precipitation')
    unit_conversion : function, optional
        Function to convert units (e.g., kelvin_to_celsius)
    description : str, optional
        Description for logging
    suffix : str, optional
        Suffix to add to output filename (e.g., "_NEW", "_TEST")
    
    Returns:
    --------
    pandas.DataFrame
        Processed data in ML-ready format
    """
    
    print(f"=== PROCESSING {description.upper()} ===")
    print(f"File: {file_path}")
    print(f"Variable: {variable_name} → {output_name}")
    
    # Step 1: Load and examine data
    print("\n1. Loading data...")
    ds = xr.open_dataset(file_path)
    
    print(f"   Dimensions: {dict(ds.dims)}")
    print(f"   Time range: {ds.valid_time.min().values} to {ds.valid_time.max().values}")
    
    # Step 2: Coordinate conversion
    print("\n2. Converting coordinates...")
    ds_converted = ds.assign_coords(longitude=(ds.longitude + 180) % 360 - 180)
    ds_converted = ds_converted.sortby('longitude')
    
    # Step 3: Subset to Canada
    print("\n3. Subsetting to Canada...")
    canada_subset = ds_converted.sel(
        latitude=slice(85.0, 40.0),
        longitude=slice(-145.0, -50.0)
    )
    print(f"   Canada dimensions: {dict(canada_subset.dims)}")
    
    # Step 4: Annual aggregation
    print("\n4. Creating annual means...")
    canada_subset['year'] = canada_subset.valid_time.dt.year
    annual_data = canada_subset.groupby('year').mean('valid_time')
    print(f"   Annual data: {dict(annual_data.dims)}")
    print(f"   Years: {annual_data.year.values}")
    
    # Step 5: Resample to 0.1° resolution
    print("\n5. Resampling to 0.1° resolution...")
    target_lats = np.arange(40.0, 85.1, 0.1)
    target_lons = np.arange(-145.0, -49.9, 0.1)
    
    resampled_data = annual_data.interp(
        latitude=target_lats,
        longitude=target_lons,
        method='linear'
    )
    print(f"   Resampled dimensions: {dict(resampled_data.dims)}")
    
    # Step 6: Extract to CSV format
    print("\n6. Extracting to CSV...")
    
    # Create coordinate meshgrid
    lons = resampled_data.longitude.values
    lats = resampled_data.latitude.values
    years = resampled_data.year.values
    years = years[(years >= 2010) & (years <= 2024)]
    lon_2d, lat_2d = np.meshgrid(lons, lats)
    
    # Find valid points
    first_year_data = resampled_data[variable_name].isel(year=0).values
    valid_mask = ~np.isnan(first_year_data.flatten())
    
    valid_lats = lat_2d.flatten()[valid_mask]
    valid_lons = lon_2d.flatten()[valid_mask]
    n_valid_pixels = len(valid_lats)
    
    print(f"   Valid land pixels: {n_valid_pixels:,}")
    
    # Extract data year by year
    all_records = []
    
    for year in years:
        year_data = resampled_data[variable_name].sel(year=year).values.flatten()[valid_mask]
        
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
    
    # Combine all years
    final_df = pd.concat(all_records, ignore_index=True)
    
    # Sort by year first, then pixel_id (CH4 format)
    final_df = final_df.sort_values(['year', 'pixel_id']).reset_index(drop=True)
    
    # Round values
    final_df['longitude'] = final_df['longitude'].round(6)
    final_df['latitude'] = final_df['latitude'].round(5)
    final_df[output_name] = final_df[output_name].round(7)
    
    # Step 7: Save and validate
    output_file = f'{output_name.title()}_2010-2024.csv'
    final_df.to_csv(output_file, index=False)
    
    print(f"\n✅ {description.upper()} COMPLETE!")
    print(f"   File: {output_file}")
    print(f"   Records: {len(final_df):,}")
    print(f"   Pixels: {final_df['pixel_id'].nunique():,}")
    print(f"   Range: {final_df[output_name].min():.3f} to {final_df[output_name].max():.3f}")
    
    # Quick validation
    print(f"\n   Validation:")
    print(f"   - Missing values: {final_df.isnull().sum().sum()}")
    print(f"   - Records per year: {final_df['year'].value_counts().iloc[0]:,}")
    print(f"   - Years: {final_df['year'].min()}-{final_df['year'].max()}")
    
    return final_df


# Add this function to era5_processor.py




# Unit conversion functions
def kelvin_to_celsius(temp_k):
    """Convert temperature from Kelvin to Celsius"""
    return temp_k - 273.15

def meters_to_millimeters(precip_m):
    """Convert precipitation from meters to millimeters"""
    return precip_m * 1000

def geopotential_to_elevation(geopotential):
    """Convert geopotential to elevation (m²/s² to meters)"""
    return geopotential / 9.80665

# Update unit conversion functions section by adding:
def volumetric_to_percentage(vol_fraction):
    """Convert volumetric fraction to percentage"""
    return vol_fraction * 100


# Convenience function for common ERA5 variables
def process_temperature(file_path):
    """Process ERA5 temperature data"""
    return process_era5_timeseries(
        file_path=file_path,
        variable_name='t2m',
        output_name='temperature',
        unit_conversion=kelvin_to_celsius,
        description='ERA5 Temperature Data'
    )


def process_precipitation(file_path):
    """Process ERA5 precipitation data"""
    return process_era5_timeseries(
        file_path=file_path,
        variable_name='tp',
        output_name='precipitation',
        unit_conversion=meters_to_millimeters,
        description='ERA5 Precipitation Data'
    )


def process_soil_moisture(file_path):
    """Process ERA5 soil moisture data"""
    return process_era5_timeseries(
        file_path=file_path,
        variable_name='swvl1',
        output_name='soil_moisture',
        unit_conversion=None,  # Already in volumetric fraction
        description='ERA5 Soil Moisture Layer 1 Data'
    )


def process_elevation(file_path):
    """
    Process ERA5 geopotential data to elevation time series
    Special handling for static elevation data (single time step → replicated across years)
    """
    import xarray as xr
    import numpy as np
    import pandas as pd
    
    print("=== PROCESSING ERA5 ELEVATION DATA ===")
    print(f"File: {file_path}")
    print("Note: Elevation is static - replicating single time step across 2010-2024")
    
    # Load geopotential data
    ds_geo = xr.open_dataset(file_path)
    
    # Find geopotential variable
    geo_var = 'z'  # Standard ERA5 geopotential variable
    
    # Get single time step (elevation doesn't change)
    if 'time' in ds_geo[geo_var].dims:
        geopotential = ds_geo[geo_var].isel(time=0)
    elif 'valid_time' in ds_geo[geo_var].dims:
        geopotential = ds_geo[geo_var].isel(valid_time=0)
    else:
        geopotential = ds_geo[geo_var]
    
    # Convert geopotential to elevation
    elevation_data = geopotential / 9.80665
    
    print(f"Elevation range: {float(elevation_data.min()):.1f}m to {float(elevation_data.max()):.1f}m")
    
    # Resample to 0.1° resolution
    target_lats = np.arange(40.0, 85.1, 0.1)
    target_lons = np.arange(-145.0, -49.9, 0.1)
    
    resampled_elevation = elevation_data.interp(
        latitude=target_lats,
        longitude=target_lons,
        method='linear'
    )
    
    # Extract valid land points
    lons = resampled_elevation.longitude.values
    lats = resampled_elevation.latitude.values
    lon_2d, lat_2d = np.meshgrid(lons, lats)
    
    elev_values = resampled_elevation.values
    valid_mask = ~np.isnan(elev_values.flatten()) & (elev_values.flatten() >= 0)
    
    # Create base elevation DataFrame (single snapshot)
    base_elevation_df = pd.DataFrame({
        'pixel_id': range(valid_mask.sum()),
        'longitude': lon_2d.flatten()[valid_mask],
        'latitude': lat_2d.flatten()[valid_mask], 
        'elevation': elev_values.flatten()[valid_mask]
    })
    
    print(f"Valid land pixels: {len(base_elevation_df):,}")
    
    # Replicate for all years 2010-2024
    all_years_data = []
    for year in range(2010, 2025):
        year_data = base_elevation_df.copy()
        year_data['year'] = year
        all_years_data.append(year_data)
    
    # Combine all years
    final_elevation_df = pd.concat(all_years_data, ignore_index=True)
    
    # Sort by year first, then pixel_id (CH4 format)
    final_elevation_df = final_elevation_df.sort_values(['year', 'pixel_id']).reset_index(drop=True)
    
    # Round values
    final_elevation_df['longitude'] = final_elevation_df['longitude'].round(6)
    final_elevation_df['latitude'] = final_elevation_df['latitude'].round(5)
    final_elevation_df['elevation'] = final_elevation_df['elevation'].round(1)
    
    # Reorder columns
    final_elevation_df = final_elevation_df[['pixel_id', 'longitude', 'latitude', 'year', 'elevation']]
    
    # Save
    output_file = f'Elevation_2010-2024.csv'
    final_elevation_df.to_csv(output_file, index=False)
    
    print(f"\n✅ ELEVATION PROCESSING COMPLETE!")
    print(f"   File: {output_file}")
    print(f"   Records: {len(final_elevation_df):,}")
    print(f"   Unique pixels: {final_elevation_df['pixel_id'].nunique():,}")
    print(f"   Elevation range: {final_elevation_df['elevation'].min():.1f}m to {final_elevation_df['elevation'].max():.1f}m")

    return final_elevation_df


def process_land_cover(file_path):
    import xarray as xr
    import numpy as np
    import pandas as pd

    print("=== PROCESSING LAND COVER DATA ===")
    print(f"File: {file_path}")
    print("Note: Land cover is static - replicating single time step across 2010-2024")

    ds = xr.open_dataset(file_path)
    lc_var = 'dom_land_cover_class'

    land_cover = ds[lc_var]  # shape: (latitude, longitude)

    # Optional: convert to integer categories if needed
    land_cover = land_cover.astype(int)

    # Resample to 0.1° resolution
    target_lats = np.arange(40.0, 85.1, 0.1)
    target_lons = np.arange(-145.0, -49.9, 0.1)

    resampled_lc = land_cover.interp(
        latitude=target_lats,
        longitude=target_lons,
        method='nearest'
    )

    # Extract valid land points
    lons = resampled_lc.longitude.values
    lats = resampled_lc.latitude.values
    lon_2d, lat_2d = np.meshgrid(lons, lats)

    lc_values = resampled_lc.values
    valid_mask = ~np.isnan(lc_values.flatten())

    base_df = pd.DataFrame({
        'pixel_id': range(valid_mask.sum()),
        'longitude': lon_2d.flatten()[valid_mask],
        'latitude': lat_2d.flatten()[valid_mask],
        'land_cover_class': lc_values.flatten()[valid_mask].astype(int)
    })

    print(f"Valid land pixels: {len(base_df):,}")

    # Replicate for all years
    all_years = []
    for year in range(2010, 2025):
        df_year = base_df.copy()
        df_year['year'] = year
        all_years.append(df_year)

    final_df = pd.concat(all_years, ignore_index=True)
    final_df = final_df.sort_values(['year', 'pixel_id']).reset_index(drop=True)

    # Round coords
    final_df['longitude'] = final_df['longitude'].round(6)
    final_df['latitude'] = final_df['latitude'].round(5)

    # Save
    final_df = final_df[['pixel_id', 'longitude', 'latitude', 'year', 'land_cover_class']]
    final_df.to_csv("LandCover_2010-2024.csv", index=False)
    output_file = "LandCover_2010-2024.csv"

    print("\n✅ LAND COVER PROCESSING COMPLETE!")
    print(f"   File: {output_file}")
    print(f"   Records: {len(final_df):,}")
    print(f"   Unique pixels: {final_df['pixel_id'].nunique():,}")
    print(f"   Class range: {final_df['land_cover_class'].min()} to {final_df['land_cover_class'].max()}")
    
    return final_df

