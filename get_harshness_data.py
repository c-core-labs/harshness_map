"""
This script is used to download annual wave height, sea ice concentration, and iceberg density data
from Copernicus Marine Service, and compute the parameters to be used in the computation of the 
Fleming-Drover Harshness Index as follows:
1. The mean annual number of days with a significant wave height greater than a given threshold (1-10 metres)
2. The mean annual number of days with a sea ice concentration greater than a given threshold (0-90 percent in increments of 10)
3. The mean annual open water iceberg areal density (# of icebergs per 100 km^2)

Data Acknowledgement: Data downloaded are provided by GHRSST, Met Office and CMEMS
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
        copernicusmarine.subset(dataset_id = dataset, 
                                variables = variables, 
                                start_datetime = data_year_start, 
                                end_datetime = data_year_end, 
                                output_directory = os.path.dirname(output_file), 
                                output_filename = os.path.basename(output_file), 
                                force_download = True,
                                credentials_file=os.path.join(data_dir, ".copernicusmarine-credentials"))
        logger.info(f"Download complete: {output_file}")
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
                    del daily_data
                    del daily_average_data
                daily_data = [band_data]
                current_date = band_date
            else: #Same date
                daily_data.append(band_data)
            del band
            del band_timestamp
            del band_date
            del band_data
            
        if daily_data: #Get average from final day
            daily_average_data = np.ma.mean(daily_data, axis=0)
            daily_average_data = np.ma.masked_array(daily_average_data, fill_value=-1)
            daily_average_data = daily_average_data.filled()
            daily_averages_raster.GetRasterBand(daily_band_num).WriteArray(daily_average_data)
            daily_averages_raster.GetRasterBand(daily_band_num).SetNoDataValue(-1)
            logger.debug(f"Wrote band {daily_band_num}")   
    del daily_averages_raster #Close the file after writing
    return(output_file)

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
        raster_x_size = input_raster.RasterXSize
        raster_y_size = input_raster.RasterYSize
        num_bands = input_raster.RasterCount
        raster_geo_transform = input_raster.GetGeoTransform()
        raster_projection = input_raster.GetProjection()
        raster_no_data_value = input_raster.GetRasterBand(1).GetNoDataValue()
    
        rows_per_chunk = int(104857600 / raster_x_size / num_bands) #TODO: Make this a parameter or optimize. For now, read in 100MB at a time (assuming 8bit input data)
        output_data_array = np.empty((0, raster_x_size), dtype=np.uint16)

        for y_off in range(0, raster_y_size, rows_per_chunk):
            y_size = min(rows_per_chunk, raster_y_size - y_off)
            input_data_array = input_raster.ReadAsArray(xoff=0, yoff=y_off, xsize=raster_x_size, ysize=y_size)
            input_data_array = np.ma.masked_values(input_data_array, raster_no_data_value) #Mask out nodata value
            chunk_output = np.ma.sum((input_data_array > thresh), axis=0) #Count bands with value greater than threshold
            chunk_output.fill_value = -1
            chunk_output = chunk_output.filled() #Fill with new nodata value of -1
            output_data_array = np.vstack((output_data_array, chunk_output))
            del input_data_array
            del chunk_output

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
            if not(mean_iceberg_concentration): #Create running count files
                raster_x_size = gdal_dataset.RasterXSize
                raster_y_size = gdal_dataset.RasterYSize
                raster_geo_transform = gdal_dataset.GetGeoTransform()
                raster_projection = gdal_dataset.GetProjection()
                mean_iceberg_concentration = np.zeros(band_data.shape)
                valid_data_count = np.zeros(band_data.shape)
            else: #Check that data is consistent
                assert raster_x_size == gdal_dataset.RasterXSize, "Input files are not the same size."
                assert raster_y_size == gdal_dataset.RasterYSize, "Input files are not the same size."
                assert raster_geo_transform == gdal_dataset.GetGeoTransform(), "Input files are not of the same geotransform."
                assert raster_projection == gdal_dataset.GetProjection(), "Input files are not in the same projection."

            #Add to running counts
            mask = np.logical_not(band_data.mask)
            mean_iceberg_concentration[mask] += band_data[mask]
            valid_data_count += mask

            del band
            del band_data
            del mask
    #Get average
    mean_iceberg_concentration /= valid_data_count
    mean_iceberg_concentration = np.ma.masked_invalid(mean_iceberg_concentration)
    mean_iceberg_concentration.fill_value= -1
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
    Data downloaded from CMEMs are: 
        cmems_mod_glo_wav_my_0.2deg_PT3H-i (VHM0) for significant wave height data 
        METOFFICE-GLO-SST-L4-REP-OBS-SST (sea_ice_fraction) for sea ice concentration data,
        DMI-ARC-SEAICE_BERG-L4-NRT-OBS (ibc) for iceberg concentration data
    Output data calculated are:
        Number of days with significant wave heights > 1, 2, 3, 4, 5, 6, 7, 8, 9, and 10 metres
        Number of days with sea ice concentration > 0, 10, 20, 30, 40, 50, 60, 70, 80, and 90 percent
        Average annual iceberg density
    Parameters:
        data_year (int): The year over which to download data and make calculations.
        data_dir (str): Path to a parent directory in which output data will be stored.
        clean (bool): If true, intermediate files will be removed.
    Returns:
        data_dir (str)"""
    
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
    
    oldest_wave_height_data_year = 1993
    oldest_sea_ice_data_year = 2007
    oldest_iceberg_data_year = 2020
    ###################################
    ##Surface Wave Significant Height##

    if (int(data_year) >= oldest_wave_height_data_year):
        #Download current year dataset
        logger.info(f"Downloading Wave Height Data")
        waves_raw_data_netcdf_name = f"waves_raw_data_{data_year}.nc"
        waves_raw_data_netcdf_name = os.path.join(data_dir, "raw_data", waves_raw_data_netcdf_name)
        if not(os.path.exists(waves_raw_data_netcdf_name)):
            download_from_cmems(output_file =   waves_raw_data_netcdf_name, 
                                data_year =     data_year, 
                                dataset =       "cmems_mod_glo_wav_my_0.2deg_PT3H-i", 
                                variables =     ["VHM0"])

        #Read raw data and generate daily averages
        logger.info("Processing waves data")
        waves_daily_averages_geotiff_name = f"waves_daily_averages_{data_year}.tif"
        waves_daily_averages_geotiff_name = os.path.join(data_dir, "raw_data", waves_daily_averages_geotiff_name)
        if not(os.path.exists(waves_daily_averages_geotiff_name)):
            get_waves_daily_averages(input_file =    waves_raw_data_netcdf_name,
                            output_file =   waves_daily_averages_geotiff_name)

        #Count days where average is > thresholds (1m - 10m)
        for wave_height_thresh in range(1,11): #Height threshold in metres
            waves_final_count_geotiff_name = f"waves_annual_data_{data_year}_{wave_height_thresh}.tif"
            waves_final_count_geotiff_name = os.path.join(data_dir, "waves_annual_data", waves_final_count_geotiff_name)
            if not(os.path.exists(waves_final_count_geotiff_name)):
                logger.info(f"Counting number of days with wave height > {wave_height_thresh}m")
                count_bands_greater_than_thresh(input_file =    waves_daily_averages_geotiff_name,
                                                output_file =   waves_final_count_geotiff_name,
                                                thresh =        wave_height_thresh * 100) #Height threshold in centimetres

        #Cleanup files
        if clean:
            os.remove(waves_daily_averages_geotiff_name)
            os.remove(waves_raw_data_netcdf_name)

        logger.info("Finished processing wave data")
    else:
        logger.info(f"Wave Height data is only available from the year {oldest_wave_height_data_year} onward. Not available for {data_year}. Skipping Wave Height data download")

    #########################
    ##Sea ice conecntration##

    if (int(data_year) >= oldest_sea_ice_data_year):
        #Download current year dataset
        logger.info(f"Downloading Sea Ice Concentration Data")
        sea_ice_raw_data_netcdf_name = f"sea_ice_raw_data_{data_year}.nc"
        sea_ice_raw_data_netcdf_name = os.path.join(data_dir, "raw_data", sea_ice_raw_data_netcdf_name)
        if not(os.path.exists(sea_ice_raw_data_netcdf_name)):
            download_from_cmems(output_file =   sea_ice_raw_data_netcdf_name,
                                data_year =     data_year, 
                                dataset =       "METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2", 
                                variables =     ["sea_ice_fraction"])

        #Count days where average is > threshold (concentration 0% - 90%)
        logger.info("Processing sea ice data")
        for sea_ice_concentration_threshold in range(0,100,10):
            sea_ice_final_count_geotiff_name = f"sea_ice_annual_data_{data_year}_{sea_ice_concentration_threshold}.tif"
            sea_ice_final_count_geotiff_name = os.path.join(data_dir, "sea_ice_annual_data", sea_ice_final_count_geotiff_name)
            if not(os.path.exists(sea_ice_final_count_geotiff_name)):
                logger.info(f"Counting number of days with sea ice concentration > {sea_ice_concentration_threshold}%")
                count_bands_greater_than_thresh(input_file =    sea_ice_raw_data_netcdf_name,
                                                output_file =   sea_ice_final_count_geotiff_name,
                                                thresh =        sea_ice_concentration_threshold)

        #Cleanup files
        if clean:
            os.remove(sea_ice_raw_data_netcdf_name)

        logger.info("Finished processing sea ice data")
    else:
        logger.info(f"Sea Ice data is only available from the year {oldest_sea_ice_data_year} onward. Not available for {data_year}. Skipping Sea Ice data download")

    ################
    ##Iceberg data##

    if (int(data_year) >= oldest_iceberg_data_year):
        #Download Current Year Data
        logger.info(f"Downloading Iceberg Data")
        iceberg_dir_name = f"iceberg_raw_data_{data_year}"
        iceberg_dir_name = os.path.join(data_dir, "raw_data", iceberg_dir_name)
        if not(os.path.exists(iceberg_dir_name)):
            download_original_files_from_cmems(output_dir = iceberg_dir_name,
                                            data_year =  data_year,
                                            dataset =    "DMI-ARC-SEAICE_BERG-L4-NRT-OBS")

        #Get average annual iceberg density
        logger.info("Processing iceberg data")
        iceberg_mean_density_name = f"iceberg_annual_data_{data_year}.tif"
        iceberg_mean_density_name = os.path.join(data_dir, "iceberg_annual_data", iceberg_mean_density_name)
        if not(os.path.exists(iceberg_mean_density_name)):
            process_iceberg_data(input_dir =    iceberg_dir_name,
                                output_file =  iceberg_mean_density_name)

        #Cleanup files
        if clean:
            shutil.rmtree(iceberg_dir_name)
            shutil.rmtree(os.path.join(data_dir, "raw_data"))

        logger.info("Finished processing iceberg data")
        return(data_dir)
    else:
        logger.info(f"Iceberg data is only available from the year {oldest_iceberg_data_year} onward. Not available for {data_year}. Skipping Iceberg data download")

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