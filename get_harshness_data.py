"""
This script is used to download annual wave height, sea ice concentration, and iceberg density data
from Copernicus Marine Service, and compute the parameters to be used in the computation of the 
Fleming-Drover Harshness Index as follows:
1. The mean annual number of days with a significant wave height greater than a given threshold (1-10 metres)
2. The mean annual number of days with a sea ice concentration greater than a given threshold (0-90 percent in increments of 10)
3. The mean annual open water iceberg areal density (# of icebergs per 100 km^2)

Additionally, relevant data is downloaded from Copernicus Climate Data Service in order to calculate an Icing Predictor Index

Data Acknowledgement: Data downloaded are provided by GHRSST, Met Office and CMEMS, as well as CDS
"""

from osgeo import gdal, osr
import os
import numpy as np
import logging
import shutil
from datetime import datetime
import argparse
from tqdm import tqdm
import time
import tempfile
import geopandas as gpd
import pandas as pd
import subprocess
import json
from cmems_utils import download_from_cmems, download_original_files_from_cmems, get_cmems_product_metadata, cmems_products_to_dict
from cds_utils import download_from_cds
from utils import count_bands_greater_than_thresh


logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

##########################
###Function Definitions###
##########################


def get_waves_daily_averages(input_file=None, 
                             output_file=None):
    """Calculates daily average values from a waves data file from CMEMs.
    Parameters:
        input_file (str): File path to a raster file containing 3-hourly data bands. Bands must include a timestamp field "NETCDF_DIM_time".
        output_file (str): File path to an output Geotiff file where the daily average data wil be saved (one band per day within the file).
    Returns:
        output_file (str)"""
    #Computes daily average data values from a cmems geotiff file {input_file} with 3 hour samples and saves it to {output_file} with one band per day
    CMEMS_WAVES_TIME_EPOCH_OFFSET = np.datetime64("1970-01-01") - np.datetime64("1950-01-01") #Offset between Unix epoch (1970-01-01) and CMEMs waves dataset epoch (1950-01-01)
    with gdal.Open(input_file) as gdal_dataset:
        current_date = None #Current day being averaged
        daily_data = []
        daily_band_num = 1
        num_bands = gdal_dataset.RasterCount
        num_daily_bands = int(num_bands / 8) #8 bands per day

        raster_x_size = gdal_dataset.RasterXSize
        raster_y_size = gdal_dataset.RasterYSize
        raster_data_type = gdal_dataset.GetRasterBand(1).DataType
        raster_geo_transform = gdal_dataset.GetGeoTransform()
        raster_projection = gdal_dataset.GetProjection()

        daily_averages_raster = gdal.GetDriverByName("GTiff").Create(output_file, raster_x_size, raster_y_size, num_daily_bands, raster_data_type)
        daily_averages_raster.SetGeoTransform(raster_geo_transform)
        daily_averages_raster.SetProjection(raster_projection)
        for band_num in range(1, num_bands+1):
            band = gdal_dataset.GetRasterBand(band_num)
            band_timestamp = int(band.GetMetadataItem("NETCDF_DIM_time"))
            band_timestamp = np.datetime64(band_timestamp, 'h') - CMEMS_WAVES_TIME_EPOCH_OFFSET
            band_date = np.datetime64(band_timestamp, 'D')
            band_data = band.ReadAsMaskedArray()
            if band_date != current_date: #New date
                if daily_data: #Get average from previous day
                    daily_average_data = np.ma.mean(daily_data, axis=0)
                    daily_averages_raster.GetRasterBand(daily_band_num).WriteArray(daily_average_data)
                    logger.debug(f"Wrote band {daily_band_num}")
                    daily_band_num += 1
                daily_data = [band_data]
                current_date = band_date
            else: #Same date
                daily_data.append(band_data)
        if daily_data: #Get average from final day
            daily_average_data = np.ma.mean(daily_data, axis=0)
            daily_average_data = np.ma.masked_array(daily_average_data, fill_value=-1)
            daily_average_data = daily_average_data.filled()
            daily_averages_raster.GetRasterBand(daily_band_num).WriteArray(daily_average_data)
            daily_averages_raster.GetRasterBand(daily_band_num).SetNoDataValue(-1)
            logger.debug(f"Wrote band {daily_band_num}")   
    del daily_averages_raster #Close the file after writing
    return(output_file)



def process_iceberg_data(input_dir=None, 
                         output_file=None):
    """Generates a Geotiff file containging mean annual iceberg density from a collection of individual NetCDF files containing iceberg concentration data.
    Parameters:
        input_dir (str): Path to a directory of NetCDF files of the same size, geo_transform, and projection, each containing a dataset "ibc" which contains iceberg concentration data.
        output_file (str): File path to a new Geotiff file that will be created containging the average iceberg density over all files in <input_dir>.
    Returns:
        output_file (str)"""
    iceberg_data = []
    raster_x_size = None
    raster_y_size = None
    raster_geo_transform = None
    raster_projection = None
    raster_data_type = gdal.GDT_Float32
    logger.info(f"Processing iceberg data in {input_dir}")
    for iceberg_file in os.listdir(input_dir):
        dataset_name = f"NETCDF:{os.path.join(input_dir, iceberg_file)}:ibc" #NETCDF dataset of the iceberg concentration data
        with gdal.Open(dataset_name) as gdal_dataset:
            band = gdal_dataset.GetRasterBand(1)
            band_data = band.ReadAsMaskedArray()
            iceberg_data.append(band_data)
            if not(raster_x_size): #Just need to do this once
                raster_x_size = gdal_dataset.RasterXSize
                raster_y_size = gdal_dataset.RasterYSize
                raster_geo_transform = gdal_dataset.GetGeoTransform()
                raster_projection = gdal_dataset.GetProjection()
            else:
                assert raster_x_size == gdal_dataset.RasterXSize, "Input files are not the same size."
                assert raster_y_size == gdal_dataset.RasterYSize, "Input files are not the same size."
                assert raster_geo_transform == gdal_dataset.GetGeoTransform(), "Input files are not of the same geotransform."
                assert raster_projection == gdal_dataset.GetProjection(), "Input files are not in the same projection." 
    mean_iceberg_concentration = np.ma.mean(iceberg_data, axis=0)
    mean_iceberg_concentration = np.ma.masked_array(mean_iceberg_concentration, fill_value=-1)
    mean_iceberg_concentration = mean_iceberg_concentration.filled()
    final_iceberg_concentration_raster = gdal.GetDriverByName("GTiff").Create(output_file, raster_x_size, raster_y_size, 1, raster_data_type)
    final_iceberg_concentration_raster.SetGeoTransform(raster_geo_transform)
    final_iceberg_concentration_raster.SetProjection(raster_projection)
    final_iceberg_concentration_raster.GetRasterBand(1).SetNoDataValue(-1)
    final_iceberg_concentration_raster.GetRasterBand(1).WriteArray(mean_iceberg_concentration)
    del final_iceberg_concentration_raster
    logger.info(f"Finished processing iceberg data. Output file is {output_file}")
    return(output_file)


def process_excel_iceberg_data(input_dir=None, 
                               output_file=None, 
                               data_year=None,
                               reference_tif="oilco_iceberg_data/fallback_cmems_raster.tif"):
    """Generates geoTIFF and shapefile of the mean annual iceberg density for oilco data.

    Args:
        input_iceberg_file (str, optional): _description_. Defaults to "Iceberg_Oilco_Insight_CellAll_V1_Final_2022.xlsx".
        shapefile_path (str, optional): _description_. Defaults to "Oilco_Insight_Grid_Cells/Metocean_Nalcor_Mar31.shp".
        data_year (_type_, optional): _description_. Defaults to None.
    """

    input_iceberg_file = os.path.join(input_dir, "Iceberg_Oilco_Insight_CellAll_V1_Final_2022.xlsx")
    ice_density_df = pd.read_excel(input_iceberg_file)
    logger.info("Done.")

    # Drop unrelated columns
    ice_density_df = ice_density_df.drop(columns=["AD_ESAT", "AD_Aerial", "AD_Sentinel"])
    # If any values are missing fill them with nan
    ice_density_df["AD_All"] = ice_density_df["AD_All"].fillna(np.nan)
    
    # Group by cell id and year and average the ice concentration per year per cell
    mean_ice_concentration_df = (
            ice_density_df.groupby(["Cell_ID", "Year"])
            .agg(Avg_AD_All = ("AD_All", lambda x: x.mean(skipna=True)))
            .reset_index()
    )
    # Replace nan values with zero
    mean_ice_concentration_df["Avg_AD_All"] = mean_ice_concentration_df["Avg_AD_All"].fillna(0) #TODO: not sure we want to do this. Could we keep them nans?
    # In preperation of merging with the shapefile that uses the term cell name instead of cell id
    mean_ice_concentration_df.rename(columns={"Cell_ID": "Cell_Name"}, inplace=True)

    # Filter the data for a specific year only
    df_selected_year = mean_ice_concentration_df[mean_ice_concentration_df["Year"] == data_year]
    # icebergs/1km^2 -> icebergs/100km^2
    df_selected_year["Avg_AD_All"] = df_selected_year["Avg_AD_All"] * 100 #TODO: Should this be * 10000?
    logger.info("Reading oilco shapefile...")

    shapefile_path = os.path.join(input_dir, "oilco_shapefile_dir/Metocean_Nalcor_Mar31.shp")
    shapefile_gdf = gpd.read_file(filename=shapefile_path)
    logger.info("Done.")
   
    # If the shapefile does not have a crs set it to EPSG:4326
    if shapefile_gdf.crs is None or shapefile_gdf.crs == "":
        shapefile_gdf = shapefile_gdf.set_crs("EPSG:4326")

    shapefile_gdf = shapefile_gdf.merge(df_selected_year, on="Cell_Name", how="left")
    
    temp_dir = tempfile.mkdtemp()
    updated_shapefile = os.path.join(temp_dir, f"oilco_mean_iceberg_density_{data_year}.shp")
    shapefile_gdf.to_file(updated_shapefile)

    #final_iceberg_concentration_file = f"oilco_mean_iceberg_density_raster_{data_year}.tif"
    resolution = 0.1 
    # Apply shapefile geometry to dataframe and generate raster geotiff
    gdal.Rasterize(output_file, 
                updated_shapefile,
                format="GTiff",
                attribute="Avg_AD_All",
                xRes=resolution, yRes=resolution,
                outputType=gdal.GDT_Float32)
    
    shutil.rmtree(temp_dir)    
    
    ds_oilco = gdal.Open(output_file, gdal.GA_Update)
    if ds_oilco is None:
        raise RuntimeError(f"Could not open {output_file} for update to fix metadata.")

    # Enforce crs to be EPSG:4326 since the data is lat/long. Currently unknown/lambert conic
    srs = osr.SpatialReference()
    srs.SetFromUserInput("EPSG:4326")
    ds_oilco.SetProjection(srs.ExportToWkt())
    ds_oilco.FlushCache()
    del ds_oilco
    logger.info(f"Forced {output_file} to EPSG:4326 (lat/lon).")

    logger.info(f"Warping {output_file} to match {reference_tif} projection.")
    warp_oilco_to_cmems(
        oilco_in=output_file,
        cmems_in=reference_tif,
        oilco_out=output_file,
        dst_nodata=-1
    )
    logger.info(f"Warp complete => final Lambert Conic in {output_file}")


def get_cmems_params(cmems_path="oilco_iceberg_data/fallback_cmems_raster.tif"):
    ds = gdal.Open(cmems_path, gdal.GA_ReadOnly)
    if ds is None:
        raise FileNotFoundError(f"Could not open {cmems_path}")
    
    projection_wkt = ds.GetProjection()

    gt = ds.GetGeoTransform()
    
    xres = abs(gt[1])  
    yres = abs(gt[5])
    
    del ds
    return projection_wkt, xres, yres


def warp_oilco_to_cmems(oilco_in, cmems_in, oilco_out, dst_nodata=-1):

    projection_wkt, x_resolution_cmems, y_resolution_cmems = get_cmems_params(cmems_in)

    warp_options = gdal.WarpOptions(dstSRS=projection_wkt, xRes=x_resolution_cmems, yRes=y_resolution_cmems, dstNodata=dst_nodata, resampleAlg="near")

    out_ds = gdal.Warp(oilco_out, oilco_in, options=warp_options)
    del out_ds
    
    return oilco_out


def mosaic_with_cmems(warped_oilco_file, cmems_file, output_mosaic_file, nodata=-1):
    
    os.environ["CHECK_DISK_FREE_SPACE"] = "FALSE"
    
    with tempfile.NamedTemporaryFile(suffix=".tif", delete=False) as temp_output:
        temp_output_file = temp_output.name
    
    cmd = [
        "gdal_merge.py",
        "-o", temp_output_file,
        "-n", str(nodata),
        "-a_nodata", str(nodata),
        warped_oilco_file,   
        cmems_file        
    ]
    subprocess.run(cmd, check=True)
    
    shutil.move(temp_output_file, output_mosaic_file)
    
    return output_mosaic_file

def create_cmems_metadata_json(product_ids,
                               json_filename):
    """Creates a json file containing metadata from the given list of CMEMs product IDs
    Parameters:
        product_ids (list(str)): A list of CMEMs product IDs to gather metadats for.
        json_filename (str): A path to the file to store the metadata. Must end in ".json"
    Returns:
        json_filename (str)"""
    products = []
    for product_id in product_ids:
        product = get_cmems_product_metadata(product_id)
        products.append(product)
    metadata_dict = cmems_products_to_dict(products)
    with open(json_filename, 'w') as file:
        json.dump(metadata_dict, file)


def download_and_preprocess_icing_predictor_data(): #TODO: Finish this function
    """Downloads sea surface temperature, wind speed, air temperature, sea ice cover, and sea water salinity data for a given year 
    from Copernicus Climate Data Store and Copernicus Marine Service and uses this data to calculate an annual Icing Predictor Index stored in <data_dir>/icing/.
        Data Downloaded from CDS are:
        derived-era5-single-levels-daily-statistics:
            sea_surface_temperature
            10m_u_component_of_wind
            110m_v_component_of_wind
            2m_temperature
            sea_ice_cover
        *cmems_mod_glo_phy_myint_0.083deg_P1D-m (so) from CMEMs for salinity (Not currently implemented)
        
    Output data calculated are:
        Number of days with Icing Predictor Index <= 0 (no icing), < 0 (light icing), > 22 (moderate icing), > 53 (heavy icing), > 83 (extreme icing)
    """
    SEA_ICE_THRESH = 0.15 #The threshold of sea ice cover below which icing can occur (expressed as a fraction of 1)
    FREEZING_POINT = -1.8 #The freezing point of sea water (in degrees Celcius)
    LIGHT_ICING_THRESH = 0 #The threshold of the icing predictor for light icing
    MODERATE_ICING_THRESH = 22 #The threshold of the icing predictor for moderate icing
    HEAVY_ICING_THRESH = 53 #The threshold of the icing predictor for heavy icing
    EXTREME_ICING_THRESH = 83 #The threshold of the icing predictor for extreme icing
    OLDEST_CDS_DATA_YEAR = 1940

    # Note on Salinity data: This data is used to calculate Tf, the freezing point of seawater, but after a brief time searching online, 
    # I couldn't find a universally agreed upon formula for salinity dependant freezing point.
    # At some point we can look into this further, but for now, maybe easiest to just use a constant value of -1.8 for Tf and not worry about salinity
    
    #Name output files
    no_icing_geotiff_name = os.path.join(data_dir, "icing_predictor_annual_data", f"icing_predictor_annual_data_{data_year}_no_icing.tif")
    light_icing_geotiff_name = os.path.join(data_dir, "icing_predictor_annual_data", f"icing_predictor_annual_data_{data_year}_light_icing.tif")
    moderate_icing_geotiff_name = os.path.join(data_dir, "icing_predictor_annual_data", f"icing_predictor_annual_data_{data_year}_moderate_icing.tif")
    heavy_icing_geotiff_name = os.path.join(data_dir, "icing_predictor_annual_data", f"icing_predictor_annual_data_{data_year}_heavy_icing.tif")
    extreme_icing_geotiff_name = os.path.join(data_dir, "icing_predictor_annual_data", f"icing_predictor_annual_data_{data_year}_extreme_icing.tif")
    
    if os.path.exists(no_icing_geotiff_name) \
    and os.path.exists(light_icing_geotiff_name) \
    and os.path.exists(moderate_icing_geotiff_name) \
    and os.path.exists(heavy_icing_geotiff_name) \
    and os.path.exists(extreme_icing_geotiff_name):
        logger.info(f"Icing predictor maps already exist for {data_year}. Skipping Icing Predictor processing.")
    
    else:
        #Download datasets#

        sea_surface_temp_netcdf_name = os.path.join(data_dir, "raw_data", f"sea_surface_temperature_{data_year}.nc")
        wind_speed_u_netcdf_name = os.path.join(data_dir, "raw_data", f"wind_speed_u_{data_year}.nc")
        wind_speed_v_netcdf_name = os.path.join(data_dir, "raw_data", f"wind_speed_v_{data_year}.nc")
        air_temp_netcdf_name = os.path.join(data_dir, "raw_data", f"air_temperature_{data_year}.nc")
        sea_ice_cover_netcdf_name = os.path.join(data_dir, "raw_data", f"sea_ice_cover_{data_year}.nc")

        #Sea Surface Temperature Data
        if os.path.exists(sea_surface_temp_netcdf_name):
            logger.info(f"Sea surface temperature data already exists for {data_year}. Skipping download.")
        else:
            logger.info(f"Downloading sea surface temperature data for {data_year}")
            start = time.time()
            download_from_cds(output_file =   sea_surface_temp_netcdf_name, 
                            data_year =     data_year,
                            dataset =       "derived-era5-single-levels-daily-statistics", 
                            variables =     ["sea_surface_temperature"])
            end = time.time()
            logger.info(f"Sea surface temperature download took {end-start} seconds")
        
        #Wind Speed Data
        if os.path.exists(wind_speed_u_netcdf_name):
            logger.info(f"U component wind speed data already exists for {data_year}. Skipping download.")
        else:
            logger.info(f"Downloading U component wind speed data for {data_year}")
            start = time.time()
            download_from_cds(output_file =   wind_speed_u_netcdf_name, 
                            data_year =     data_year,
                            dataset =       "derived-era5-single-levels-daily-statistics", 
                            variables =     ["10m_u_component_of_wind"])
            end = time.time()
            logger.info(f"U component wind speed download took {end-start} seconds")
        if os.path.exists(wind_speed_v_netcdf_name):
            logger.info(f"V component wind speed data already exists for {data_year}. Skipping download.")
        else:
            logger.info(f"Downloading V component wind speed data for {data_year}")
            start = time.time()
            download_from_cds(output_file =   wind_speed_v_netcdf_name, 
                            data_year =     data_year,
                            dataset =       "derived-era5-single-levels-daily-statistics", 
                            variables =     ["10m_v_component_of_wind"])
            end = time.time()
            logger.info(f"V component wind speed download took {end-start} seconds")
        
        #Air Temperature Data
        if os.path.exists(air_temp_netcdf_name):
            logger.info(f"Air temperature data already exists for {data_year}. Skipping download.")
        else:
            logger.info(f"Downloading air temperature data for {data_year}")
            start = time.time()
            download_from_cds(output_file =   air_temp_netcdf_name, 
                            data_year =     data_year,
                            dataset =       "derived-era5-single-levels-daily-statistics", 
                            variables =     ["2m_temperature"])
            end = time.time()
            logger.info(f"Air temperature download took {end-start} seconds")
        
        #Sea Ice Cover Data
        if os.path.exists(sea_ice_cover_netcdf_name):
            logger.info(f"Sea ice cover data already exists for {data_year}. Skipping download.")
        else:
            logger.info(f"Downloading sea ice cover data for {data_year}")
            start = time.time()
            download_from_cds(output_file =   sea_ice_cover_netcdf_name, 
                            data_year =     data_year,
                            dataset =       "derived-era5-single-levels-daily-statistics", 
                            variables =     ["sea_ice_cover"])
            end = time.time()
            logger.info(f"Sea ice cover download took {end-start} seconds")
        

        #Open Datasets#
        
        
        try:
            sea_surface_temp_netcdf_dataset_name = f"NETCDF:{sea_surface_temp_netcdf_name}:sst"
            wind_speed_U_netcdf_dataset_name = f"NETCDF:{wind_speed_u_netcdf_name}:u10"
            wind_speed_V_netcdf_dataset_name = f"NETCDF:{wind_speed_v_netcdf_name}:v10"
            air_temp_netcdf_dataset_name = f"NETCDF:{air_temp_netcdf_name}:t2m"
            sea_ice_cover_netcdf_dataset_name = f"NETCDF:{sea_ice_cover_netcdf_name}:siconc"

            sea_surface_temp_dataset = gdal.Open(sea_surface_temp_netcdf_dataset_name)
            wind_speed_U_dataset = gdal.Open(wind_speed_U_netcdf_dataset_name)
            wind_speed_V_dataset = gdal.Open(wind_speed_V_netcdf_dataset_name)
            air_temp_dataset = gdal.Open(air_temp_netcdf_dataset_name)
            sea_ice_cover_dataset = gdal.Open(sea_ice_cover_netcdf_dataset_name)
            
            raster_x_size = sea_surface_temp_dataset.RasterXSize
            raster_y_size = sea_surface_temp_dataset.RasterYSize
            ref_projection = sea_surface_temp_dataset.GetProjection()    
            ref_geotransform  = sea_surface_temp_dataset.GetGeoTransform()
            num_bands = sea_surface_temp_dataset.RasterCount
        except Exception as e:
            logger.error(f"Failed to read raw datasets: {e}")
            return None
            
        try:
            assert ref_projection == wind_speed_U_dataset.GetProjection(), "Projections of raw data files do not match"
            assert ref_projection == wind_speed_V_dataset.GetProjection(), "Projections of raw data files do not match"
            assert ref_projection == air_temp_dataset.GetProjection(), "Projections of raw data files do not match"
            assert ref_projection == sea_ice_cover_dataset.GetProjection(), "Projections of raw data files do not match"
            assert ref_geotransform == wind_speed_U_dataset.GetGeoTransform(), "Geotransforms of raw data files do not match"
            assert ref_geotransform == wind_speed_V_dataset.GetGeoTransform(), "Geotransforms of raw data files do not match"
            assert ref_geotransform == air_temp_dataset.GetGeoTransform(), "Geotransforms of raw data files do not match"
            assert ref_geotransform == sea_ice_cover_dataset.GetGeoTransform(), "Geotransforms of raw data files do not match"
            assert num_bands == wind_speed_U_dataset.RasterCount, "Number of bands in raw data files do not match"
            assert num_bands == wind_speed_V_dataset.RasterCount, "Number of bands in raw data files do not match"
            assert num_bands == air_temp_dataset.RasterCount, "Number of bands in raw data files do not match"
            assert num_bands == sea_ice_cover_dataset.RasterCount, "Number of bands in raw data files do not match"
        except AssertionError as a:
            logger.error(f"Raw data files metadata does not match: {a}")
            return None
        
        
        #Calculate daily icing predictor values#
        
        
        if os.path.exists(no_icing_geotiff_name) \
        and os.path.exists(light_icing_geotiff_name) \
        and os.path.exists(moderate_icing_geotiff_name) \
        and os.path.exists(heavy_icing_geotiff_name) \
        and os.path.exists(extreme_icing_geotiff_name):
            logger.info(f"Icing predictor maps already exist for {data_year}. Skipping calculation.")
        
        else:
            logger.info(f"Calculating daily icing predictor values for {data_year}")     
            no_icing_days = np.zeros([raster_y_size, raster_x_size])
            light_icing_days = np.zeros([raster_y_size, raster_x_size])
            moderate_icing_days = np.zeros([raster_y_size, raster_x_size])
            heavy_icing_days = np.zeros([raster_y_size, raster_x_size])
            extreme_icing_days = np.zeros([raster_y_size, raster_x_size])
            
            start = time.time()
            for band_num in tqdm(range(1, num_bands+1)):
                sea_surface_temp = sea_surface_temp_dataset.GetRasterBand(band_num).ReadAsArray()
                wind_speed_U = wind_speed_U_dataset.GetRasterBand(band_num).ReadAsArray()
                wind_speed_V = wind_speed_V_dataset.GetRasterBand(band_num).ReadAsArray()
                air_temp = air_temp_dataset.GetRasterBand(band_num).ReadAsArray()
                sea_ice_cover = sea_ice_cover_dataset.GetRasterBand(band_num).ReadAsArray()
                
                sea_surface_temp = sea_surface_temp - 273.15 #Convert from Kelvin to Celcius
                wind_speed = np.sqrt(wind_speed_U**2 + wind_speed_V**2) #Calculate total wind speed
                air_temp = air_temp - 273.15 #Convert from Kelvin to Celcius
                sea_ice_mask = np.where(np.isnan(sea_ice_cover), np.nan, sea_ice_cover < SEA_ICE_THRESH) #Calculate whether sea_ice_cover is below threshold, preserving nan values
                icing_predictor = sea_ice_mask * (wind_speed * (FREEZING_POINT - air_temp) / (1 + 0.3 * (sea_surface_temp - FREEZING_POINT))) #Calculate icing predictor
                
                no_icing_days += np.where(np.isnan(icing_predictor), np.nan, icing_predictor <= LIGHT_ICING_THRESH) #Count days with no icing
                light_icing_days += np.where(np.isnan(icing_predictor), np.nan, (icing_predictor > LIGHT_ICING_THRESH)) #Count days with light or greater icing
                moderate_icing_days += np.where(np.isnan(icing_predictor), np.nan, (icing_predictor > MODERATE_ICING_THRESH)) #Count days with moderate or greater icing
                heavy_icing_days += np.where(np.isnan(icing_predictor), np.nan, (icing_predictor > HEAVY_ICING_THRESH)) #Count days with heavy or greater icing
                extreme_icing_days += np.where(np.isnan(icing_predictor), np.nan, icing_predictor > EXTREME_ICING_THRESH) #Count days with extreme icing
            end = time.time()
            logger.info(f"Icing predictor calculation took {end-start} seconds")

            
            #Create output files
            logger.info(f"Writing icing predictor maps to {os.path.join(data_dir, 'icing_predictor_annual_data')}")
            no_icing_raster = gdal.GetDriverByName("GTiff").Create(no_icing_geotiff_name, raster_x_size, raster_y_size, bands=1, eType=gdal.GDT_Float32)
            light_icing_raster = gdal.GetDriverByName("GTiff").Create(light_icing_geotiff_name, raster_x_size, raster_y_size, bands=1, eType=gdal.GDT_Float32)
            moderate_icing_raster = gdal.GetDriverByName("GTiff").Create(moderate_icing_geotiff_name, raster_x_size, raster_y_size, bands=1, eType=gdal.GDT_Float32)
            heavy_icing_raster = gdal.GetDriverByName("GTiff").Create(heavy_icing_geotiff_name, raster_x_size, raster_y_size, bands=1, eType=gdal.GDT_Float32)
            extreme_icing_raster = gdal.GetDriverByName("GTiff").Create(extreme_icing_geotiff_name, raster_x_size, raster_y_size, bands=1, eType=gdal.GDT_Float32)
            
            no_icing_raster.SetGeoTransform(ref_geotransform)
            light_icing_raster.SetGeoTransform(ref_geotransform)
            moderate_icing_raster.SetGeoTransform(ref_geotransform)
            heavy_icing_raster.SetGeoTransform(ref_geotransform)
            extreme_icing_raster.SetGeoTransform(ref_geotransform)
            
            no_icing_raster.SetProjection(ref_projection)
            light_icing_raster.SetProjection(ref_projection)
            moderate_icing_raster.SetProjection(ref_projection)
            heavy_icing_raster.SetProjection(ref_projection)
            extreme_icing_raster.SetProjection(ref_projection)
            
            no_icing_raster.GetRasterBand(1).WriteArray(no_icing_days)
            light_icing_raster.GetRasterBand(1).WriteArray(light_icing_days)
            moderate_icing_raster.GetRasterBand(1).WriteArray(moderate_icing_days)
            heavy_icing_raster.GetRasterBand(1).WriteArray(heavy_icing_days)
            extreme_icing_raster.GetRasterBand(1).WriteArray(extreme_icing_days)
            
            #CDS datasets sometimes use x coordinates of 0 to 360 instead of -180 to 180. Warp the rasters to use 0 as center longitude
            default_center_long = gdal.GetConfigOption('CENTER_LONG')
            gdal.SetConfigOption('CENTER_LONG', '0')
            gdal.Warp(no_icing_geotiff_name, no_icing_geotiff_name, srcSRS="+proj=longlat +ellps=WGS84", dstSRS="WGS84")
            gdal.Warp(light_icing_geotiff_name, light_icing_geotiff_name, srcSRS="+proj=longlat +ellps=WGS84", dstSRS="WGS84")
            gdal.Warp(moderate_icing_geotiff_name, moderate_icing_geotiff_name, srcSRS="+proj=longlat +ellps=WGS84", dstSRS="WGS84")
            gdal.Warp(heavy_icing_geotiff_name, heavy_icing_geotiff_name, srcSRS="+proj=longlat +ellps=WGS84", dstSRS="WGS84")
            gdal.Warp(extreme_icing_geotiff_name, extreme_icing_geotiff_name, srcSRS="+proj=longlat +ellps=WGS84", dstSRS="WGS84")
            gdal.SetConfigOption('CENTER_LONG', default_center_long)
            
            #Clean up datasets
            del sea_surface_temp_dataset
            del wind_speed_U_dataset
            del wind_speed_V_dataset
            del air_temp_dataset
            del sea_ice_cover_dataset
            del no_icing_raster
            del light_icing_raster
            del moderate_icing_raster
            del heavy_icing_raster
            del extreme_icing_raster
            
            logger.info(f"Finished calculating daily icing predictor values for {data_year}")






def download_and_preprocess_iceberg_data(): #TODO: Finish this function
    """Downloads iceberg for a given year from Copernicus Marine Service.
    Averages the iceberg density over the year.
    """
    OLDEST_ICEBERG_DATA_YEAR = 2020
    OLDEST_OILCO_ICEBERG_DATA_YEAR = 1998
    NEWEST_OILCO_ICEBERG_DATA_YEAR = 2021

    iceberg_mean_density_name = f"iceberg_annual_data_{data_year}.tif"
    iceberg_mean_density_name = os.path.join(data_dir, "iceberg_annual_data", iceberg_mean_density_name)
    if os.path.exists(iceberg_mean_density_name):
        logger.info(f"Iceberg data already exists for {data_year}. Skipping Iceberg processing.")
    else:
        if (int(data_year) >= OLDEST_ICEBERG_DATA_YEAR):
            #Download Current Year Data
            logger.info(f"Downloading Iceberg Data")
            iceberg_dir_name = f"iceberg_raw_data_{data_year}"
            iceberg_dir_name = os.path.join(data_dir, "raw_data", iceberg_dir_name)
            if not(os.path.exists(iceberg_dir_name)):
                start = time.time()
                download_original_files_from_cmems(output_dir = iceberg_dir_name,
                                                data_year =  data_year,
                                                dataset =    "DMI-ARC-SEAICE_BERG-L4-NRT-OBS",
                                                credentials_file =  CMEMS_CREDENTIALS_FILE)
                end = time.time()
                logger.info(f"Iceberg download took {end-start} seconds")

            #Get average annual iceberg density
            logger.info("Processing iceberg data")
            if not(os.path.exists(iceberg_mean_density_name)):
                start = time.time()
                process_iceberg_data(input_dir =    iceberg_dir_name,
                                    output_file =  iceberg_mean_density_name)
                end = time.time()
                logger.info(f"Iceberg processing took {end-start} seconds")


            logger.info("Finished processing iceberg data")

        else:
            logger.info(f"Iceberg data is only available from the year {OLDEST_ICEBERG_DATA_YEAR} onward. Not available for {data_year}. Skipping Iceberg data download")

    ####################
    #OILCO Iceberg Data#

    if (int(data_year) >= OLDEST_OILCO_ICEBERG_DATA_YEAR and int(data_year) <= NEWEST_OILCO_ICEBERG_DATA_YEAR):
        logger.info("Reading oilco ice density excel data...")
        
        oilco_data_dir = "oilco_iceberg_data"
        oilco_iceberg_mean_density_name = f"oilco_iceberg_annual_data_{data_year}.tif"
        oilco_iceberg_mean_density_name = os.path.join(oilco_data_dir, oilco_iceberg_mean_density_name)
        
        if not(os.path.exists(oilco_iceberg_mean_density_name)):
            process_excel_iceberg_data(input_dir = oilco_data_dir,
                                    output_file = oilco_iceberg_mean_density_name,
                                    data_year = int(data_year),
                                    reference_tif=iceberg_mean_density_name)
        else:
            logger.info(f"{oilco_iceberg_mean_density_name} already exists.")
    else:
        logger.info(f"No Oilco data for year {data_year} (valid range is 1998-2021)")
        
           
    cmems_tif = os.path.join(data_dir, "iceberg_annual_data", f"iceberg_annual_data_{data_year}.tif")
    oilco_tif = os.path.join("oilco_iceberg_data", f"oilco_iceberg_annual_data_{data_year}.tif")     
    
    # Only overlap over 2020-2021
    if os.path.exists(cmems_tif) and os.path.exists(oilco_tif):
        logger.info(f"Merging Oilco & CMEMS iceberg rasters found for {data_year}...")
        cmems_tif = os.path.join(data_dir, "iceberg_annual_data", f"iceberg_annual_data_{data_year}.tif")
        mosaic_with_cmems(oilco_tif, cmems_tif, cmems_tif, nodata=-1)
        logger.info("Mosaic created!")
    else:
        logger.info("No overlapping Oilco & CMEMS data to merge or only one dataset is available.")



def download_and_preprocess_cmems_data(data_year = datetime.today().year-1,
                                 data_dir =  os.path.join(os.getcwd(), "data"),
                                 variables = ["siconc", "VHM0"],
                                 metadata_file = os.path.join(os.getcwd(), "data", "CMEMs_metadata.json"),
                                 clean = True): #TODO: Finish this function
    """Downloads data defined in [variables] for a given year from Copernicus Marine Service.
    Calculates the number of days in the year that they exceed various thresholds.
    Averages the iceberg density over the year.
    Output data is stored in <data_dir>/<variable_name>/
    
    Default data downloaded from CMEMs are:
        cmems_mod_glo_wav_my_0.2deg_PT3H-i (VHM0) for significant wave height data 
        GLOBAL_MULTIYEAR_PHY_001_030 (siconc) for sea ice concentration data
            Note: used to use (make sure they are similar): METOFFICE-GLO-SST-L4-REP-OBS-SST (sea_ice_fraction) for sea ice concentration data,
        DMI-ARC-SEAICE_BERG-L4-NRT-OBS (ibc) for iceberg concentration data
        
    Output data calculated are:
        Number of days with significant wave heights > 1, 2, 3, 4, 5, 6, 7, 8, 9, and 10 metres
        Number of days with sea ice concentration > 0, 10, 20, 30, 40, 50, 60, 70, 80, and 90 percent
    Parameters:
        data_year (int): The year over which to download data and make calculations.
        data_dir (str): Path to a parent directory in which output data will be stored.
        variables (list(str)): a list of variables to be downloaded and preprocessed.
        metadata_file (str): The path to a json file containing CMEMs metadata for the variables provided. Created by create_cmems_metadata_json.
        clean (bool): If true, intermediate files will be removed.
    Returns:
        data_dir (str)"""
    
    NUM_THRESHOLDS = 9
    DEFAULT_TIME_UNITS = "milliseconds since 1970-01-01 00:00:00Z (no leap seconds)"
    MILLISECONDS_IN_A_DAY = 86400000
    NUM_BANDS_TO_ESTIMATE_MIN_MAX = 4 #If min/max values cannot be read from metadata, use a sample from this many bands to estimate min/max of the dataset

    CMEMS_CREDENTIALS_FILE = os.path.join(data_dir, ".copernicusmarine-credentials")
    
    logger.info(f"Beginning Download and Preprocessing of the following variables for year {data_year}: {variables}")

    if not(os.path.exists(data_dir)):
        os.mkdir(data_dir)

    for variable in variables:
        variable_dir = os.path.join(data_dir, variable)
        if not(os.path.exists(variable_dir)):
            os.mkdir(variable_dir)

        existing_file_count = sum(f"{data_year}" in file_name for file_name in os.listdir(variable_dir))
        if existing_file_count == NUM_THRESHOLDS:
            logger.info(f"{variable} files already exist for {data_year}. Skipping.")
            continue

        with open(metadata_file, 'r') as file:
            metadata = json.load(file)

        if variable not in metadata:
            error = f"No entry for {variable} found in metadata file, {metadata_file}. Skipping."
            continue    

        #Select dataset to use
        logger.info(f"Selecting best dataset to use for {variable}")
        for dataset_name in metadata[variable]["datasets"]:

            try:
                dataset = metadata[variable]["datasets"][dataset_name]
                
                time_units = dataset["time_units"]
                if time_units != DEFAULT_TIME_UNITS:
                    logger.debug(f"{dataset_name} uses unknown time units. Dataset will not be used.")
                    continue
                
                time_step = dataset["time_step"]
                if time_step != MILLISECONDS_IN_A_DAY:
                    logger.debug(f"{dataset_name} uses a time step other than 1 day. Dataset will not be used.")
                    continue

                #Check that dataset contains the entire year #TODO: Deal with the case when data_year is split between two datasets
                start_time = dataset["start_time"]
                start_time = datetime.fromtimestamp(start_time / 1000)
                if start_time > datetime(data_year, 1, 1):
                    logger.debug(f"{dataset_name} does not contain all of {data_year}. Dataset will not be used.")
                    continue
                end_time = dataset["end_time"]
                end_time = datetime.fromtimestamp(end_time / 1000)
                if end_time < datetime(data_year + 1, 1, 1):
                    logger.debug(f"{dataset_name} does not contain all of {data_year}. Dataset will not be used.")
                    continue
                
                #All checks passed, dataset is suitable, break from loop:
                break

            except KeyError as e:
                error = f"Error while reading from metadata file {metadata_file}: {e}"
                logger.error(error)
                raise(error)

        #Download data
        logger.info(f"Downloading {variable} Data")
        raw_data_netcdf_name = f"{variable}_{data_year}_raw_data.nc"
        raw_data_netcdf_name = os.path.join(variable_dir, raw_data_netcdf_name)
        if os.path.exists(raw_data_netcdf_name):
            logger.info(f"{raw_data_netcdf_name} already exists. Skipping download")
        else:
            start = time.time()
            download_from_cmems(dataset =           dataset_name, 
                                variables =         [variable],
                                start_datetime =    datetime(data_year, 1, 1),
                                end_datetime =      datetime(data_year + 1, 1, 1),
                                output_file =       raw_data_netcdf_name, 
                                credentials_file =  CMEMS_CREDENTIALS_FILE)
            end = time.time()
            logger.info(f"{variable} download took {end-start} seconds")

        #Compute thresholds
        #Note: Was working on a method to calculate thresholds using mean snd standard dev, but won't work well for non normal distributed data, so will just use uniform thresholds between min and max for now
        # with gdal.Open(raw_data_netcdf_name) as file:
        #     scale_factor = file.GetMetadata()[f"{variable}#scale_factor"]
        #     add_offset = file.GetMetadata()[f"{variable}#add_offset"]
        #     num_bands = file.RasterCount
        #     num_samples = file.RasterXSize * file.RasterYSize
        #     means = []
        #     variances = []
        #     for band_num in range(1, num_bands+1):
        #         band = file.GetRasterBand(band_num)
        #         (minimum, maximum, mean, standard_deviation) = band.GetStatistics(True, True) #approx_ok=True, force=True
        #         means.append(mean)
        #         variances.append(standard_deviation**2)
        # total_samples = num_samples * num_bands
        # overall_mean = np.mean(means)
        # s1 = np.sum([v * num_samples for v in variances])
        # s2 = np.sum([num_samples * (m - overall_mean)**2 for m in means])
        # overall_variance = (s1 + s2) / total_samples
        # overall_standard_deviation = np.sqrt(overall_variance)

        logger.info(f"Processing {variable} data")
        #Compute uniform thresholds
        with gdal.Open(raw_data_netcdf_name) as file:
            try:
                minimum = float(file.GetMetadata()[f"{variable}#valid_min"])
                maximum = float(file.GetMetadata()[f"{variable}#valid_max"])
            except KeyError as e:
                num_bands = min(NUM_BANDS_TO_ESTIMATE_MIN_MAX, file.RasterCount)
                minimum = float('inf')
                maximum = float('-inf')
                logger.warning(f"Could not obtain min/max value from {raw_data_netcdf_name} metadata. Using sample of {num_bands} raster bands to estimate min/max.")
                for i in range(1,num_bands+1):
                    band_num = i * int(file.RasterCount / num_bands) #Pick evenly distributed bands
                    band = file.GetRasterBand(band_num)
                    (band_minimum, band_maximum, band_mean, band_standard_deviation) = band.GetStatistics(True, True) #approx_ok=True, force=True
                    minimum = min(band_minimum, minimum)
                    maximum = max(band_maximum, maximum)
            try:
                scale_factor = float(file.GetMetadata()[f"{variable}#scale_factor"])
                add_offset = float(file.GetMetadata()[f"{variable}#add_offset"])            
            except KeyError as e:
                scale_factor = 1
                add_offset = 0
                logger.warning(f"Scale information not found for {variable}. This will not affect results, but threshold values in file names may not represent expected units.")
        step = (maximum - minimum) / (NUM_THRESHOLDS + 1)
        thresholds = [minimum + step * t for t in range(1,NUM_THRESHOLDS+1)]
        
        #Count days where value is > thresholds
        for threshold in thresholds:
            # Adjust threshold using scale and offset from metadata.
            # Note, unadjusted values and thresholds are used for calculation since we are only counting bands > threshold, it makes no difference.
            # Adjusted threshold is just used for naming output files since it properly reflects the real world units.
            scaled_threshold = threshold * scale_factor + add_offset

            final_count_geotiff_name = f"{variable}_{data_year}_{scaled_threshold:.2g}.tif"
            final_count_geotiff_name = os.path.join(variable_dir, final_count_geotiff_name)
            if os.path.exists(final_count_geotiff_name):
                logger.info(f"{final_count_geotiff_name} already exists. Skipping Processing.")
            else:
                logger.info(f"Counting number of days with {variable} > {scaled_threshold:.2g}")
                start = time.time()
                count_bands_greater_than_thresh(input_file =    raw_data_netcdf_name,
                                                output_file =   final_count_geotiff_name,
                                                thresh =        threshold)
                end = time.time()
                logger.info(f"Count took {end-start:.1f} seconds")

        if clean:
            try:
                os.remove(raw_data_netcdf_name)
            except UnboundLocalError:
                logger.debug(f"raw_data_netcdf_name not set. Skipping deletion.")

        logger.info(f"Finished processing {variable} data")
    
    return(data_dir)
        
        
        
##########
###MAIN###
##########

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=download_and_preprocess_data.__doc__, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("--data_year", help="The year over which to download data and make calculations. (int)")
    parser.add_argument("--data_dir", help="Path to a parent directory in which output data will be stored. (str)")
    parser.add_argument("--clean", help="If true, intermediate files will be removed. (bool)")
         
    #Default values
    data_year=datetime.today().year-1
    data_dir = os.path.join(os.getcwd(), "data")
    clean = True

    args = parser.parse_args()
    if args.data_year is not None:
        data_year = args.data_year
    if args.data_dir is not None:
        data_dir = args.data_dir
    if args.clean is not None:
        clean = args.clean

    download_and_preprocess_data(data_year=data_year,
                                 data_dir = data_dir,
                                 clean=clean)