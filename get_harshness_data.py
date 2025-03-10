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

import copernicusmarine
from osgeo import gdal
import os
import numpy as np
import logging
from osgeo_utils import gdal_calc
import calendar
import shutil
from datetime import datetime
import argparse
import cdsapi
from tqdm import tqdm
import time

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

##########################
###Function Definitions###
##########################

def num_days_in_year(year=None):
    """Returns the number of days in a given year
    Parameters:
        year (int)
    Returns:
        number_of_days (int)"""
    return 365 + calendar.isleap(year)

def download_from_cmems(output_file=None, 
                        data_year=None, 
                        dataset=None, 
                        variables=None):
    """Downloads data from the Copernicus Marine Service using copernicusmarine.subset().
    Parameters:
        output_file (str): File path where the downloaded dataset will be stored. Must be a NetCDF (.nc) file.
        data_year (int): The dataset will contain data from Jan 1 to Dec 31 of this year.
        dataset (str): The name of the CMEMs dataset to download.
        variables (list[str]): A list of variables to download from the dataset.
    Returns:
        output_file (str)"""
    if not os.path.exists(output_file):
        data_year_start = f"{data_year}-01-01"
        data_year_end = f"{data_year}-12-31T21:00:00"
        logger.info(f"Downloading {data_year} {dataset} from Copernicus Marine")
        try:
            copernicusmarine.subset(dataset_id = dataset, 
                                    variables = variables, 
                                    start_datetime = data_year_start, 
                                    end_datetime = data_year_end, 
                                    output_directory = os.path.dirname(output_file), 
                                    output_filename = os.path.basename(output_file), 
                                    force_download = True,
                                    credentials_file=os.path.join(data_dir, ".copernicusmarine-credentials"))
            logger.info(f"Download complete: {output_file}")
        except Exception as e:
            logger.error(f"Download failed: {e}")
            os.remove(output_file)
            raise e

    else: #File already downloaded
        logger.info(f"{output_file} already exists. Skipping download")
    return(output_file)

def download_original_files_from_cmems(output_dir=None, 
                                       data_year=None, 
                                       dataset=None):
    """Downloads individual data files from the Copernicus Marine Service using copernicusmarine.get().
    This function should be used inplace of download_from_cmems for datasets that do not yet support the copernicusmarine.subset function.
    Parameters:
        output_dir (str): Path to the directory where the downloaded data files will be stored.
        data_year (int): The dataset will contain data from Jan 1 to Dec 31 of this year.
        dataset (str): The name of the CMEMs dataset to download.
    Returns:
        output_dir (str)"""
    if not os.path.exists(output_dir):
        logger.info(f"Downloading {data_year} {dataset} from Copernicus Marine")
        copernicusmarine.get(dataset_id = dataset, 
                             filter = f"*_{data_year}*",
                             output_directory = output_dir, 
                             no_directories = True,
                             force_download = True,
                             credentials_file=os.path.join(data_dir, ".copernicusmarine-credentials"))
        logger.info(f"Download complete: {output_dir}")
    else: #File already downloaded
        logger.info(f"{output_dir} already exists. Skipping download")    
    return(output_dir)

def download_from_cds(output_file=None,
                      data_year=None, 
                      dataset=None, 
                      variables=None):
    """Downloads data from the Copernicus Climate Data Service using the cdsapi.
    Parameters:
        output_file (str): File path where the downloaded dataset will be stored. Must be a NetCDF (.nc) file.
        data_year (int): The dataset will contain data from Jan 1 to Dec 31 of this year.
        dataset (str): The name of the CDS dataset to download.
        variables (list[str]): A list of variables to download from the dataset.
    Returns:
        output_file (str)"""
    
    if not os.path.exists(output_file):
        
        logger.info(f"Downloading {data_year} {dataset} {variables} from Copernicues Climate Data Store")

        try:
            client = cdsapi.Client()
            request = {
                "product_type": "reanalysis",
                "variable": variables,
                "year": data_year,
                "month": [
                    "01", "02", "03",
                    "04", "05", "06",
                    "07", "08", "09",
                    "10", "11", "12"
                ],
                "day": [
                    "01", "02", "03",
                    "04", "05", "06",
                    "07", "08", "09",
                    "10", "11", "12",
                    "13", "14", "15",
                    "16", "17", "18",
                    "19", "20", "21",
                    "22", "23", "24",
                    "25", "26", "27",
                    "28", "29", "30",
                    "31"
                ],
                "daily_statistic": "daily_mean",
                "time_zone": "utc+00:00",
                "frequency": "1_hourly"
            }
            client.retrieve(dataset, request, output_file)
            logger.info(f"Download complete: {output_file}")
            
        except Exception as e:
            logger.error(f"Download failed: {e}")
            os.remove(output_file)
            raise e
    
    else: #File already downloaded
        logger.info(f"{output_file} already exists. Skipping download")
        
    return(output_file)

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

def get_average_raster(input_files = [],
                       output_file = None):
    """Calculates the mean of all rasters in the list <input_files> and writes the result to <output_file>.
    Parameters:
        input_files (list[str]): A list of paths to input raster files to be averaged. All rasters must be the same size, geo_transform, and projection.
        output_file (str): The path to the file in which the averaged raster will be saved. Must be a Geotiff.
    Returns: 
        output_file (str)"""
    data = []
    raster_x_size = None
    raster_y_size = None
    raster_geo_transform = None
    raster_projection = None
    raster_data_type = gdal.GDT_Float32
    logger.info(f"Calculating average of {len(input_files)} rasters.")
    for file in input_files:
        with gdal.Open(file) as gdal_dataset:
            band = gdal_dataset.GetRasterBand(1)
            band_data = band.ReadAsMaskedArray()
            data.append(band_data)
            if not(raster_x_size): #Just need to do this once
                raster_x_size = gdal_dataset.RasterXSize
                raster_y_size = gdal_dataset.RasterYSize
                raster_geo_transform = gdal_dataset.GetGeoTransform()
                raster_projection = gdal_dataset.GetProjection()
            else:
                assert raster_x_size == gdal_dataset.RasterXSize, "Raster files are not the same size."
                assert raster_y_size == gdal_dataset.RasterYSize, "Raster files are not the same size."
                assert raster_geo_transform == gdal_dataset.GetGeoTransform(), "Raster files are not of the same geotransform."
                assert raster_projection == gdal_dataset.GetProjection(), "Raster files are not in the same projection." 
    mean_data = np.ma.mean(data, axis=0)
    mean_data = np.ma.masked_array(mean_data, fill_value=-1)
    mean_data = mean_data.filled()
    output_raster = gdal.GetDriverByName("GTiff").Create(output_file, raster_x_size, raster_y_size, 1, raster_data_type)
    output_raster.SetGeoTransform(raster_geo_transform)
    output_raster.SetProjection(raster_projection)
    output_raster.GetRasterBand(1).SetNoDataValue(-1)
    output_raster.GetRasterBand(1).WriteArray(mean_data)
    del output_raster
    logger.info(f"Finished averageing rasters. Output file is {output_file}")
    return output_file


def count_bands_greater_than_thresh(input_file=None, 
                                    output_file=None, 
                                    thresh=None):
    """Counts the number of bands with pxel values greater than a given threshold, given an input raster file with multiple bands and saves the output to <output_file>.
    Parameteres:
        input_file (str): File path to a raster file containging multiple bands of data.
        output_file (str): File path to an output Geotiff file containging a single band with the counts of input bands > thresh.
        thresh (float): The threshold to compare input raster bands to.
    Returns:
        output_file (str)"""
    with gdal.Open(input_file) as input_raster:
        input_data_array = input_raster.ReadAsArray()
        raster_x_size = input_raster.RasterXSize
        raster_y_size = input_raster.RasterYSize
        raster_geo_transform = input_raster.GetGeoTransform()
        raster_projection = input_raster.GetProjection()
        raster_no_data_value = input_raster.GetRasterBand(1).GetNoDataValue()
    input_data_array = np.ma.masked_values(input_data_array, raster_no_data_value) #Mask out nodata value
    output_data_array = np.ma.sum((input_data_array > thresh), axis=0) #Count bands with pixel value > thresh
    del input_data_array
    output_data_array = np.ma.masked_array(output_data_array, fill_value=-1) #Set new nodata value to -1
    output_data_array = output_data_array.filled() #Fill nodata values
    raster_data_type = gdal.GDT_UInt16
    output_raster = gdal.GetDriverByName("GTiff").Create(output_file, raster_x_size, raster_y_size, 1, raster_data_type)
    output_raster.SetGeoTransform(raster_geo_transform)
    output_raster.SetProjection(raster_projection)
    output_raster.GetRasterBand(1).SetNoDataValue(-1)
    output_raster.GetRasterBand(1).WriteArray(output_data_array)
    del output_data_array
    print(f"Finished counting bands > {thresh} with numpy. Output file is {output_file}")
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

def download_and_preprocess_data(data_year = datetime.today().year-1,
                                 data_dir =  os.path.join(os.getcwd(), "data"),
                                 clean = True):
    """Downloads wave height, sea ice, and iceberg data for a given year from Copernicus Marine Service.
    Calculates the number of days in the year that the wave height and sea ice concentration data exceed various thresholds.
    Averages the iceberg density over the year.
    Output data is stored in <data_dir>/waves_annual_data/, <data_dir>/sea_ice_annual_data/, and <data_dir>/iceberg_annual_data/
    
    Also downloads sea surface temperature, wind speed, air temperature, sea ice cover, and sea water salinity data for a given year 
    from Copernicus Climate Data Store and Copernicus Marine Service and uses this data to calculate an annual Icing Predictor Index stored in <data_dir>/icing_predictor_annual_data/.
    Data downloaded from CMEMs are: 
        cmems_mod_glo_wav_my_0.2deg_PT3H-i (VHM0) for significant wave height data 
        METOFFICE-GLO-SST-L4-REP-OBS-SST (sea_ice_fraction) for sea ice concentration data,
        DMI-ARC-SEAICE_BERG-L4-NRT-OBS (ibc) for iceberg concentration data
        cmems_mod_glo_phy_myint_0.083deg_P1D-m (so) for salinity (Not currently implemented)
    Data Downloaded from CDS are:
        derived-era5-single-levels-daily-statistics:
            sea_surface_temperature
            10m_u_component_of_wind
            110m_v_component_of_wind
            2m_temperature
            sea_ice_cover
    Output data calculated are:
        Number of days with significant wave heights > 1, 2, 3, 4, 5, 6, 7, 8, 9, and 10 metres
        Number of days with sea ice concentration > 0, 10, 20, 30, 40, 50, 60, 70, 80, and 90 percent
        Average annual iceberg density
        Number of days with Icing Predictor Index <= 0 (no icing), < 0 (light icing), > 22 (moderate icing), > 53 (heavy icing), > 83 (extreme icing)
    Parameters:
        data_year (int): The year over which to download data and make calculations.
        data_dir (str): Path to a parent directory in which output data will be stored.
        clean (bool): If true, intermediate files will be removed.
    Returns:
        data_dir (str)"""
    
    SEA_ICE_THRESH = 0.15 #The threshold of sea ice cover below which icing can occur (expressed as a fraction of 1)
    FREEZING_POINT = -1.8 #The freezing point of sea water (in degrees Celcius)
    LIGHT_ICING_THRESH = 0 #The threshold of the icing predictor for light icing
    MODERATE_ICING_THRESH = 22 #The threshold of the icing predictor for moderate icing
    HEAVY_ICING_THRESH = 53 #The threshold of the icing predictor for heavy icing
    EXTREME_ICING_THRESH = 83 #The threshold of the icing predictor for extreme icing
    OLDEST_WAVE_HEIGHT_DATA_YEAR = 1993
    OLDEST_SEA_ICE_DATA_YEAR = 2007
    OLDEST_ICEBERG_DATA_YEAR = 2020
    OLDEST_CDS_DATA_YEAR = 1940
    
    # Note on Salinity data: This data is used to calculate Tf, the freezing point of seawater, but after a brief time searching online, 
    # I couldn't find a universally agreed upon formula for salinity dependant freezing point.
    # At some point we can look into this further, but for now, maybe easiest to just use a constant value of -1.8 for Tf and not worry about salinity
    
    logger.info(f"Beginning Download and Preprocessing for Year {data_year}")

    if not(os.path.exists(data_dir)):
        os.mkdir(data_dir)
    if not(os.path.exists(os.path.join(data_dir, "raw_data"))):
        os.mkdir(os.path.join(data_dir, "raw_data"))
    if not(os.path.exists(os.path.join(data_dir, "waves_annual_data"))):
        os.mkdir(os.path.join(data_dir, "waves_annual_data"))
    if not(os.path.exists(os.path.join(data_dir, "sea_ice_annual_data"))):
        os.mkdir(os.path.join(data_dir, "sea_ice_annual_data"))
    if not(os.path.exists(os.path.join(data_dir, "iceberg_annual_data"))):
        os.mkdir(os.path.join(data_dir, "iceberg_annual_data"))
    if not os.path.exists(os.path.join(data_dir, "icing_predicor_annual_data")):
        os.mkdir(os.path.join(data_dir, "icing_predictor_annual_data"))
    

    ###################################
    ##Surface Wave Significant Height##

    if (int(data_year) >= OLDEST_WAVE_HEIGHT_DATA_YEAR):
        #Download current year dataset
        logger.info(f"Downloading Wave Height Data")
        waves_raw_data_netcdf_name = f"waves_raw_data_{data_year}.nc"
        waves_raw_data_netcdf_name = os.path.join(data_dir, "raw_data", waves_raw_data_netcdf_name)
        if not(os.path.exists(waves_raw_data_netcdf_name)):
            start = time.time()
            download_from_cmems(output_file =   waves_raw_data_netcdf_name, 
                                data_year =     data_year, 
                                dataset =       "cmems_mod_glo_wav_my_0.2deg_PT3H-i", 
                                variables =     ["VHM0"])
            end = time.time()
            logger.debug(f"Waves download took {end-start} seconds")

        #Read raw data and generate daily averages
        logger.info("Processing waves data")
        waves_daily_averages_geotiff_name = f"waves_daily_averages_{data_year}.tif"
        waves_daily_averages_geotiff_name = os.path.join(data_dir, "raw_data", waves_daily_averages_geotiff_name)
        if not(os.path.exists(waves_daily_averages_geotiff_name)):
            start = time.time()
            get_waves_daily_averages(input_file =    waves_raw_data_netcdf_name,
                            output_file =   waves_daily_averages_geotiff_name)
            end = time.time()
            logger.debug(f"Waves daily averages took {end-start} seconds")

        #Count days where average is > thresholds (1m - 10m)
        for wave_height_thresh in range(1,11): #Height threshold in metres
            waves_final_count_geotiff_name = f"waves_annual_data_{data_year}_{wave_height_thresh}.tif"
            waves_final_count_geotiff_name = os.path.join(data_dir, "waves_annual_data", waves_final_count_geotiff_name)
            if not(os.path.exists(waves_final_count_geotiff_name)):
                logger.info(f"Counting number of days with wave height > {wave_height_thresh}m")
                start = time.time()
                count_bands_greater_than_thresh(input_file =    waves_daily_averages_geotiff_name,
                                                output_file =   waves_final_count_geotiff_name,
                                                thresh =        wave_height_thresh * 100) #Height threshold in centimetres
                end = time.time()
                logger.debug(f"Wave height count took {end-start} seconds")

        #Cleanup files
        if clean:
            os.remove(waves_daily_averages_geotiff_name)
            os.remove(waves_raw_data_netcdf_name)

        logger.info("Finished processing wave data")
    else:
        logger.info(f"Wave Height data is only available from the year {OLDEST_WAVE_HEIGHT_DATA_YEAR} onward. Not available for {data_year}. Skipping Wave Height data download")

    #########################
    ##Sea ice conecntration##

    if (int(data_year) >= OLDEST_SEA_ICE_DATA_YEAR):
        #Download current year dataset
        logger.info(f"Downloading Sea Ice Concentration Data")
        sea_ice_raw_data_netcdf_name = f"sea_ice_raw_data_{data_year}.nc"
        sea_ice_raw_data_netcdf_name = os.path.join(data_dir, "raw_data", sea_ice_raw_data_netcdf_name)
        if not(os.path.exists(sea_ice_raw_data_netcdf_name)):
            start = time.time()
            download_from_cmems(output_file =   sea_ice_raw_data_netcdf_name,
                                data_year =     data_year, 
                                dataset =       "METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2", 
                                variables =     ["sea_ice_fraction"])
            end = time.time()
            logger.debug(f"Sea ice download took {end-start} seconds")

        #Count days where average is > threshold (concentration 0% - 90%)
        logger.info("Processing sea ice data")
        for sea_ice_concentration_threshold in range(0,100,10):
            sea_ice_final_count_geotiff_name = f"sea_ice_annual_data_{data_year}_{sea_ice_concentration_threshold}.tif"
            sea_ice_final_count_geotiff_name = os.path.join(data_dir, "sea_ice_annual_data", sea_ice_final_count_geotiff_name)
            if not(os.path.exists(sea_ice_final_count_geotiff_name)):
                logger.info(f"Counting number of days with sea ice concentration > {sea_ice_concentration_threshold}%")
                start = time.time()
                count_bands_greater_than_thresh(input_file =    sea_ice_raw_data_netcdf_name,
                                                output_file =   sea_ice_final_count_geotiff_name,
                                                thresh =        sea_ice_concentration_threshold)
                end = time.time()
                logger.debug(f"Sea ice count took {end-start} seconds")

        #Cleanup files
        if clean:
            os.remove(sea_ice_raw_data_netcdf_name)

        logger.info("Finished processing sea ice data")
    else:
        logger.info(f"Sea Ice data is only available from the year {OLDEST_SEA_ICE_DATA_YEAR} onward. Not available for {data_year}. Skipping Sea Ice data download")

    ################
    ##Iceberg data##

    if (int(data_year) >= OLDEST_ICEBERG_DATA_YEAR):
        #Download Current Year Data
        logger.info(f"Downloading Iceberg Data")
        iceberg_dir_name = f"iceberg_raw_data_{data_year}"
        iceberg_dir_name = os.path.join(data_dir, "raw_data", iceberg_dir_name)
        if not(os.path.exists(iceberg_dir_name)):
            start = time.time()
            download_original_files_from_cmems(output_dir = iceberg_dir_name,
                                            data_year =  data_year,
                                            dataset =    "DMI-ARC-SEAICE_BERG-L4-NRT-OBS")
            end = time.time()
            logger.debug(f"Iceberg download took {end-start} seconds")

        #Get average annual iceberg density
        logger.info("Processing iceberg data")
        iceberg_mean_density_name = f"iceberg_annual_data_{data_year}.tif"
        iceberg_mean_density_name = os.path.join(data_dir, "iceberg_annual_data", iceberg_mean_density_name)
        if not(os.path.exists(iceberg_mean_density_name)):
            start = time.time()
            process_iceberg_data(input_dir =    iceberg_dir_name,
                                output_file =  iceberg_mean_density_name)
            end = time.time()
            logger.debug(f"Iceberg processing took {end-start} seconds")

        #Cleanup files
        if clean:
            shutil.rmtree(iceberg_dir_name)
            shutil.rmtree(os.path.join(data_dir, "raw_data"))

        logger.info("Finished processing iceberg data")
    else:
        logger.info(f"Iceberg data is only available from the year {OLDEST_ICEBERG_DATA_YEAR} onward. Not available for {data_year}. Skipping Iceberg data download")

   
    #######################
    #Icing Predictor Index#
    
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
        logger.debug(f"Sea surface temperature download took {end-start} seconds")
    
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
        logger.debug(f"U component wind speed download took {end-start} seconds")
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
        logger.debug(f"V component wind speed download took {end-start} seconds")
    
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
        logger.debug(f"Air temperature download took {end-start} seconds")
    
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
        logger.debug(f"Sea ice cover download took {end-start} seconds")
    

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
    
    
    #Name output files
    no_icing_geotiff_name = os.path.join(data_dir, "icing_predictor_annual_data", f"icing_predictor_{data_year}_no_icing.tif")
    light_icing_geotiff_name = os.path.join(data_dir, "icing_predictor_annual_data", f"icing_predictor_{data_year}_light_or_greater.tif")
    moderate_icing_geotiff_name = os.path.join(data_dir, "icing_predictor_annual_data", f"icing_predictor_{data_year}_moderate_or_greater.tif")
    heavy_icing_geotiff_name = os.path.join(data_dir, "icing_predictor_annual_data", f"icing_predictor_{data_year}_heavy_or_greater.tif")
    extreme_icing_geotiff_name = os.path.join(data_dir, "icing_predictor_annual_data", f"icing_predictor_{data_year}_extreme_icing.tif")
    
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
        logger.debug(f"Icing predictor calculation took {end-start} seconds")

        
        #Create output files
        logger.info(f"Writing icing predictor maps to {os.path.join(data_dir, 'icing_predictor_annual_data')}")
        no_icing_raster = gdal.GetDriverByName("GTiff").Create(no_icing_geotiff_name, raster_x_size, raster_y_size, bands=1, eType=gdal.GDT_Float32)
        light_icing_raster = gdal.GetDriverByName("GTiff").Create(light_icing_geotiff_name, raster_x_size, raster_y_size, bands=1, eType=gdal.GDT_Float32)
        moderate_icing_raster = gdal.GetDriverByName("GTiff").Create(moderate_icing_geotiff_name, raster_x_size, raster_y_size, bands=1, eType=gdal.GDT_Float32)
        heavy_icing_raster = gdal.GetDriverByName("GTiff").Create(heavy_icing_geotiff_name, raster_x_size, raster_y_size, bands=1, eType=gdal.GDT_Float32)
        extreme_icing_raster = gdal.GetDriverByName("GTiff").Create(extreme_icing_geotiff_name, raster_x_size, raster_y_size, bands=1, eType=gdal.GDT_Float32)
        
        #CDS datasets sometimes use x coordinates of 0 to 360 instead of -180 to 180. We will explicitly set the upper left corner to -180,90 for the output rasters
        ref_geotransform[0] = -180
        ref_geotransform[3] = 90
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