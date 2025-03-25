
from osgeo import gdal
import numpy as np
import logging
import calendar

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

def num_days_in_year(year=None):
    """Returns the number of days in a given year
    Parameters:
        year (int)
    Returns:
        number_of_days (int)"""
    return 365 + calendar.isleap(year)

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
        raster_x_size = input_raster.RasterXSize
        raster_y_size = input_raster.RasterYSize
        raster_geo_transform = input_raster.GetGeoTransform()
        raster_projection = input_raster.GetProjection()
        raster_no_data_value = input_raster.GetRasterBand(1).GetNoDataValue()
        
        rows_per_chunk = 10 #TODO: Make this a parameter or optimize
        output_data_array = np.empty((0, raster_x_size), dtype=np.uint16)
        
        for y_off in range(0, raster_y_size, rows_per_chunk):
            y_size = min(rows_per_chunk, raster_y_size - y_off)
            input_data_array = input_raster.ReadAsArray(xoff=0, yoff=y_off, xsize=raster_x_size, ysize=y_size)
            input_data_array = np.ma.masked_values(input_data_array, raster_no_data_value) #Mask out nodata value
            chunk_output = np.ma.sum((input_data_array > thresh), axis=0) #Count bands with value greater than threshold
            chunk_output = np.ma.masked_array(chunk_output, fill_value=-1).filled() #Fill with new nodata value of -1
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