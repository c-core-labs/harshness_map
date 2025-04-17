"""
This script is used to generate harshness maps from Wave Height, Sea Ice Concentration, and Iceberg Density data.
By default the Fleming-Drover Harshness Index is used, but users may specify custom formulas for harshness calculations.

Data Acknowledgement: Data utilized are provided by GHRSST, Met Office and CMEMS
"""

from datetime import datetime
import logging
import os
import numpy as np
from osgeo import gdal
from osgeo_utils import gdal_calc
import argparse
import base64
from string import ascii_uppercase as alphabet

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

##########################
###Function Definitions###
##########################

def encode_string(original_string):
    """Encodes any string into a unique URL (or filename) safe string
    Parameters:
        original_string (str): Any UTF-8 string.
    Returns:
        encoded_string (str): A unique URL (or filename) safe base64 encoded UTF-8 string."""
    encoded_string = base64.urlsafe_b64encode(original_string.encode('utf-8'))
    return encoded_string.decode('utf-8')

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

    if output_file is None:
        raise ValueError("Output file is required")

    logger.info(f"Calculating average of {len(input_files)} rasters.")

    if len(input_files == 1):
        shutil.copy(input_files[0], output_file)
        logger.info(f"Only one input file given. Input file simply copied to output. Output file is {output_file}")
        return output_file

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

def calculate_harshness_orig(start_year=datetime.today().year-1,
                        end_year=datetime.today().year-1,
                        data_dir = os.path.join(os.getcwd(), "data"),
                        wave_height_thresh=4,
                        sea_ice_concentration_thresh=60,
                        icing_thresh="light",
                        formula="6*S/350 + 2.5*W/110 + 1.5*(I>0.01)*(12 + 2*log10((I/10000)+1e-40))",
                        x_resolution=0.2,
                        y_resolution=0.2,
                        bounds=(-71, 41, 25, 82),
                        clean = True):
    """Calculates a harshness index given by <formula> using the input data defined by input parameters.
    Input annual data files are found in <data_dir>.
    Files containging data between <start_year> and <end_year> for the given <wave_height_thresh>, <sea_ice_concentration_thresh>, and <icing_thresh> are averaged together to generate the input parameters for <formula>
    The resulting harshness map is generated in EPSG:4326 at a resolution <x_resolution> x <y_resolution> degrees for the region defined in <bounds>
    The harshness map is saved as a Geotiff in "<data_dir>/harshness_maps/" with the name:
    "harshness_<start_year>_<end_year>_<wave_height_thresh>_<sea_ice_concentration_thresh>_<icing_thresh>_<x_resolution>_<y_resolution>_<bounds>.tif"
    Parameters:
        start_year (int): The first year of data files to include (inclusive).
        end_year (int): The last year of data to include (inclusive).
        data_dir (str): The parent directory containing the directories "iceberg_annual_data", "waves_annual_data", and "sea_ice_annual_data".
        wave_height_thresh (int 1-10): Used to define wave input files. Input files contain # of days with wave height > <wave_height_thresh> metres.
        sea_ice_concentration_thresh (int 0-90 by 10s): Used to define sea ice input files. Input files contain # of days with sea ice concentration > <sea_ice_concentration_thresh> percent.
        icing_thresh (str): The threshold for the icing predictor index. Options are "light", "moderate", "heavy", or "extreme".
        formula (str): The formula used for calculating the harshness index using gdalCalc where S = Sea Ice Data, W = Wave Data, I = Iceberg Data, P = Icing Predictor Data.
        x_resolution (float): Resolution of the output file in degrees longitude.
        y_resolution (float): Resolution of the output file in degrees latitude.
        bounds (tuple(float)): The bounds of the output file in degrees (lon_min, lat_min, lon_max, lat_max).
        clean (bool): If true, intermediate files will be removed.
    Returns:
        harshness_file_name (str): File path of the output harshness map file"""

    formula = formula.replace(" ", "")
    logger.info("Running Harshness Map Calculator with the following parameters:")
    logger.info(f"Start Year: {start_year}")
    logger.info(f"End Year: {end_year}")
    logger.info(f"Data Directory: {data_dir}")
    logger.info(f"Wave Height Threshold: {wave_height_thresh}")
    logger.info(f"Sea Ice Concentration Threshold: {sea_ice_concentration_thresh}")
    logger.info(f"Icing Threshold: {icing_thresh}")
    logger.info(f"Harshness Formula: {formula}")
    logger.info(f"X Resolution: {x_resolution}")
    logger.info(f"Y Resolution: {y_resolution}")
    logger.info(f"Bounds: {bounds}")    
    logger.info(f"Clean: {clean}")

    encoded_formula = encode_string(formula) #Create a unique identifier to represent the formula in the filename
    intermediate_files = []
    harshness_file_name = f"harshness_{start_year}_{end_year}_{wave_height_thresh}_{sea_ice_concentration_thresh}_{icing_thresh}_{x_resolution}_{y_resolution}_{bounds[0]}_{bounds[1]}_{bounds[2]}_{bounds[3]}_{encoded_formula}.tif".replace(" ", "")
    harshness_file_name = os.path.join(data_dir, 'harshness_maps', harshness_file_name)

    #If the requested harshness map already exists, return it
    if os.path.exists(harshness_file_name):
        logger.info(f"Harshness map with these parameters already ecists: {harshness_file_name}")
        return harshness_file_name
    
    #Ensure files exist for all years
    data_years = list(range(start_year, end_year + 1))
    waves_files = []
    sea_ice_files = []
    iceberg_files = []
    icing_files = []
    for data_year in data_years:
        try:
            waves_file = os.path.join(data_dir, "waves_annual_data", f"waves_annual_data_{data_year}_{wave_height_thresh}.tif")
            assert os.path.exists(waves_file), f"Annual data file in given date range does not exist: {waves_file}."
            waves_files.append(waves_file)
            sea_ice_file = os.path.join(data_dir, "sea_ice_annual_data", f"sea_ice_annual_data_{data_year}_{sea_ice_concentration_thresh}.tif")
            assert os.path.exists(sea_ice_file), f"Annual data file in given date range does not exist: {sea_ice_file}."
            sea_ice_files.append(sea_ice_file)
            iceberg_file = os.path.join(data_dir, "iceberg_annual_data", f"iceberg_annual_data_{data_year}.tif")
            assert os.path.exists(iceberg_file), f"Annual data file in given date range does not exist: {iceberg_file}."
            iceberg_files.append(iceberg_file)
            icing_file = os.path.join(data_dir, "icing_predictor_annual_data", f"icing_predictor_annual_data_{data_year}_{icing_thresh}_icing.tif")
            assert os.path.exists(icing_file), f"Annual data file in given date range does not exist: {icing_file}."
            icing_files.append(icing_file)
        except AssertionError as a:
            logger.error(a)
            return None

    if len(data_years) > 1:
        #Average waves files
        average_waves_file = os.path.join(data_dir, "waves_annual_data", f"waves_average_data_{start_year}-{end_year}_{wave_height_thresh}.tif")
        get_average_raster(waves_files, average_waves_file)
        intermediate_files.append(average_waves_file)

        #Average sea_ice files
        average_sea_ice_file = os.path.join(data_dir, "sea_ice_annual_data", f"sea_ice_average_data_{start_year}-{end_year}_{sea_ice_concentration_thresh}.tif")
        get_average_raster(sea_ice_files, average_sea_ice_file)
        intermediate_files.append(average_sea_ice_file)

        #Average iceberg files
        average_iceberg_file = os.path.join(data_dir, "iceberg_annual_data", f"iceberg_average_data_{start_year}-{end_year}.tif")
        get_average_raster(iceberg_files, average_iceberg_file)
        intermediate_files.append(average_iceberg_file)
        
        #Average icing files
        average_icing_file = os.path.join(data_dir, "icing_predictor_annual_data", f"icing_predictor_average_{start_year}-{end_year}_{icing_thresh}_icing.tif")
        get_average_raster(icing_files, average_icing_file)
        intermediate_files.append(average_icing_file)

    elif len(data_years) == 1:
        average_waves_file = waves_files[0]
        average_sea_ice_file = sea_ice_files[0]
        average_iceberg_file = iceberg_files[0]
        average_icing_file = icing_files[0]
    else:
        logger.error(f"End year must be later than or equal to start year. Got start year: {start_year}, end year: {end_year}")
        return None

    #Warp files to the desired resolution and bounds
    #Waves
    warped_waves_file = average_waves_file.replace(".tif", f"_{x_resolution}_{y_resolution}_{bounds}.tif".replace(" ", ""))
    warped = gdal.Warp(warped_waves_file, average_waves_file, xRes=x_resolution, yRes=y_resolution, srcSRS="EPSG:4326", dstSRS="EPSG:4326", outputBounds=bounds)
    del warped
    intermediate_files.append(warped_waves_file)

    #Sea Ice
    warped_sea_ice_file = average_sea_ice_file.replace(".tif", f"_{x_resolution}_{y_resolution}_{bounds}.tif".replace(" ", ""))
    warped = gdal.Warp(warped_sea_ice_file, average_sea_ice_file, xRes=x_resolution, yRes=y_resolution, srcSRS="EPSG:4326", dstSRS="EPSG:4326", outputBounds=bounds)
    del warped
    intermediate_files.append(warped_sea_ice_file)

    #Icebergs
    warped_iceberg_file = average_iceberg_file.replace(".tif", f"_{x_resolution}_{y_resolution}_{bounds}.tif".replace(" ", ""))
    warped = gdal.Warp(warped_iceberg_file, average_iceberg_file, xRes=x_resolution, yRes=y_resolution, dstSRS="EPSG:4326", outputBounds=bounds) #Iceberg files contain projection information, so no need to specify srcSRS
    del warped
    intermediate_files.append(warped_iceberg_file)
    
    #Icing Predictor
    warped_icing_file = average_icing_file.replace(".tif", f"_{x_resolution}_{y_resolution}_{bounds}.tif".replace(" ", ""))
    warped = gdal.Warp(warped_icing_file, average_icing_file, xRes=x_resolution, yRes=y_resolution, srcSRS="EPSG:4326", dstSRS="EPSG:4326", outputBounds=bounds)
    del warped
    intermediate_files.append(warped_icing_file)
   
    #Calculate Harshness
    if not(os.path.exists(os.path.join(data_dir, 'harshness_maps'))):
           os.mkdir(os.path.join(data_dir, 'harshness_maps'))         
    gdal_calc.Calc(formula, 
                S=warped_sea_ice_file, 
                W=warped_waves_file,
                I=warped_iceberg_file,
                P=warped_icing_file,
                outfile=harshness_file_name)
    
    if clean:
        for intermediate_file in intermediate_files:
            os.remove(intermediate_file)

    logger.info(f"Finished createing harshness map: {harshness_file_name}")
    return harshness_file_name


def calculate_harshness(start_year=datetime.today().year-1,
                        end_year=datetime.today().year-1,
                        data_dir = os.path.join(os.getcwd(), "data"),
                        variables = ["siconc", "VHM0", "ibc", "icing"],
                        thresholds = [60, 4, None, "light"],
                        formula="6*A/350 + 2.5*B/110 + 1.5*(C>0.01)*(12 + 2*log10((C/10000)+1e-40))",
                        x_resolution=0.2,
                        y_resolution=0.2,
                        bounds=(-71, 41, 25, 82),
                        clean = True):
    """Calculates a harshness index given by <formula> using the input data defined by input parameters.
    Input annual data files are found in <data_dir>.
    Files containging data between <start_year> and <end_year> for the given thresholds are averaged together to generate the input parameters for <formula>
    The resulting harshness map is generated in EPSG:4326 at a resolution <x_resolution> x <y_resolution> degrees for the region defined in <bounds>
    The harshness map is saved as a Geotiff in "<data_dir>/harshness_maps/" with the name:
    "harshness_<start_year>_<end_year>_<variables>_<thresholds>_<x_resolution>_<y_resolution>_<bounds>_<formula_id>.tif"
    Parameters:
        start_year (int): The first year of data files to include (inclusive).
        end_year (int): The last year of data to include (inclusive).
        data_dir (str): The parent directory containing the directories "iceberg_annual_data", "waves_annual_data", and "sea_ice_annual_data".
        variables (list(str)): A list of variables to be used in the formula. Also corresponds to the input files that will be searched for in data_dir.
        thresholds (list): A list of thresholds to use corrsponding to the given variables. For more information, see README.md
        formula (str): The formula used for calculating the harshness index using gdalCalc where A, B, C, etc. correspond to the ordered list of variables provided.
        x_resolution (float): Resolution of the output file in degrees longitude.
        y_resolution (float): Resolution of the output file in degrees latitude.
        bounds (tuple(float)): The bounds of the output file in degrees (lon_min, lat_min, lon_max, lat_max).
        clean (bool): If true, intermediate files will be removed.
    Returns:
        harshness_file_name (str): File path of the output harshness map file"""

    if len(variables) != len(thresholds):
        error = f"Lengths of 'variables' and 'thresholds' lists must be equal. Got lengths {len(variables)} and {len(thresholds)}. Harshness map generation cancelled."
        logger.error(error)
        raise ValueError(error)

    if start_year > end_year:
        error = f"End year must be later than or equal to start year. Got start year: {start_year}, end year: {end_year}. Harshness map generation cancelled."
        logger.error(error)
        raise ValueError(error)

    formula = formula.replace(" ", "")
    logger.info("Running Harshness Map Calculator with the following parameters:")
    logger.info(f"Start Year: {start_year}")
    logger.info(f"End Year: {end_year}")
    logger.info(f"Data Directory: {data_dir}")
    logger.info("Variables and corresponding thresholds:")
    for (variable, threshold) in zip(variables, thresholds):
        logger.info(f"Var: {variable}")
        logger.info(f"Thresh: {threshold}")
    logger.info(f"Harshness Formula: {formula}")
    logger.info(f"X Resolution: {x_resolution}")
    logger.info(f"Y Resolution: {y_resolution}")
    logger.info(f"Bounds: {bounds}")    
    logger.info(f"Clean: {clean}")

    encoded_formula = encode_string(formula) #Create a unique identifier to represent the formula in the filename
    variables_string = '_'.join(variables)
    thresholds_string = '_'.join(map(str, thresholds))
    intermediate_files = []
    harshness_file_name = f"harshness_{start_year}_{end_year}_{variables_string}_{thresholds_string}_{x_resolution}_{y_resolution}_{bounds[0]}_{bounds[1]}_{bounds[2]}_{bounds[3]}_{encoded_formula}.tif".replace(" ", "")
    harshness_file_name = os.path.join(data_dir, 'harshness_maps', harshness_file_name)

    #If the requested harshness map already exists, return it
    if os.path.exists(harshness_file_name):
        logger.info(f"Harshness map with these parameters already exists: {harshness_file_name}")
        return harshness_file_name
    
    #Gather all input files
    data_years = list(range(start_year, end_year + 1))
    input_files = [] #This will be a 2 dimensional list of input files per variable
    variables_missing_files = []
    for (variable, threshold) in zip(variables, thresholds):
        variable_files = []
        for data_year in data_years:
            try:
                input_file = os.path.join(data_dir, variable, f"{variable}_{data_year}_{threshold}.tif")
                assert os.path.exists(input_file), f"Annual data file in given date range does not exist: {input_file}."
                variable_files.append(input_file)
            except AssertionError as e:
                variables_missing_files.append(variable)
                continue
        input_files.append(variable_files)
    if variables_missing_files:
        error = f"The following variables are missing input files: {variables_missing_files}. Harshness map generation cancelled."
        logger.error(error)
        raise AssertionError(error)

    #Average and warp input files and create arg dict
    args = {}
    for index, (variable, treshold) in enumerate(zip(variables, thresholds)):
        #Average
        average_file = os.path.join(data_dir, variable, f"average_{variable}_{start_year}-{end_year}_{threshold}.tif")
        get_average_raster(input_files[index], average_file)
        intermediate_files.append(average_file)

        #Warp
        warped_file = average_file.replace(".tif", f"_{x_resolution}_{y_resolution}_{bounds}.tif".replace(" ", ""))
        warped = gdal.Warp(warped_file, average_file, xRes=x_resolution, yRes=y_resolution, srcSRS="EPSG:4326", dstSRS="EPSG:4326", outputBounds=bounds)
        del warped
        intermediate_files.append(warped_file)

        #Add to arg dict for gdal_calc
        args[alphabet[index]] = warped_file
   
    #Calculate Harshness
    if not(os.path.exists(os.path.join(data_dir, 'harshness_maps'))):
           os.mkdir(os.path.join(data_dir, 'harshness_maps'))
    
    gdal_calc.Calc(formula, 
                outfile=harshness_file_name,
                *args)
    
    if clean:
        for intermediate_file in intermediate_files:
            os.remove(intermediate_file)

    logger.info(f"Finished createing harshness map: {harshness_file_name}")
    return harshness_file_name


##########
###Main###
##########

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=calculate_harshness.__doc__, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("--start_year", help="The first year of data files to include (inclusive). (int)")
    parser.add_argument("--end_year", help="The last year of data to include (inclusive). (int)")
    parser.add_argument("--data_dir", help="The parent directory containing the directories 'iceberg_annual_data', 'waves_annual_data', and 'sea_ice_annual_data'. (str)")
    parser.add_argument("--variables", help ="A list of variables to be used in the formula. Also corresponds to the input files that will be searched for in data_dir. (Comma Separated Strings)")
    parser.add_argument("--thresholds", help ="A list of thresholds to use corrsponding to the given variables. For more information, see README.md. Comma Separated Values")
    parser.add_argument("--formula", help="The formula used for calculating the harshness index using gdalCalc where A, B, C, etc. correspond to the ordered list of variables provided. (str)")
    parser.add_argument("--x_resolution", help="Resolution of the output file in degrees longitude. (float)")
    parser.add_argument("--y_resolution", help="Resolution of the output file in degrees latitude. (float)")
    parser.add_argument("--lon_min", help="The minimum longitude bound of the output file in degrees. (float)")
    parser.add_argument("--lat_min", help="The minimum latitude bound of the output file in degrees. (float)")
    parser.add_argument("--lon_max", help="The maximum longitude bound of the output file in degrees. (float)")
    parser.add_argument("--lat_max", help="The maximum latitude bound of the output file in degrees. (float)")
    parser.add_argument("--clean", help="If true, intermediate files will be removed. (bool)") 

    #Default values
    start_year=datetime.today().year-1
    end_year=start_year
    data_dir = os.path.join(os.getcwd(), "data")
    variables = ["siconc", "VHM0", "ibc", "icing"]
    thresholds = [60, 4, None, "light"]
    formula="6*A/350 + 2.5*B/110 + 1.5*(C>0.01)*(12 + 2*log10((C/10000)+1e-40))"
    x_resolution=0.2
    y_resolution=0.2
    lon_min=-71
    lat_min=41
    lon_max=25
    lat_max=82
    clean = True

    args = parser.parse_args()
    
    if args.thresholds is not None:
        #Parse through thresholds list since it can contain strings, ints, floats, or None values
        thresholds = []
        threshold_strings = args.thresholds.split(",")
        for t in threshold_strings:
            if t.isdigit():
                thresholds.append(int(t))
            elif re.match('^[+-]?\d*\.\d+$', t):
                threholds.append(float(t))
            elif t == "None":
                thresholds.append(None)
            else:
                thresholds.append(t)

    if args.start_year is not None:
        start_year = int(args.start_year)
    if args.end_year is not None:
        end_year = int(args.end_year)
    else:
        end_year = start_year
    if args.data_dir is not None:
        data_dir = args.data_dir
    if args.variables is not None:
        variables = args.variables.split(",")
    if args.formula is not None:
        formula = args.formula
    if args.x_resolution is not None:
        x_resolution = float(args.x_resolution)
    if args.y_resolution is not None:
        y_resolution = float(args.y_resolution)
    if args.lon_min is not None:
        lon_min = float(args.lon_min)
    if args.lat_min is not None:
        lat_min = float(args.lat_min)
    if args.lon_max is not None:
        lon_max = float(args.lon_max)
    if args.lat_max is not None:
        lat_max = float(args.lat_max)   
    if args.clean is not None:
        clean = bool(args.clean)

    try:
        calculate_harshness(start_year=start_year,
                            end_year=end_year,
                            data_dir = data_dir,
                            variables=variables,
                            thresholds=thresholds,
                            formula=formula,
                            x_resolution=x_resolution,
                            y_resolution=y_resolution,
                            bounds=(lon_min,lat_min,lon_max,lat_max),
                            clean=clean)
    except ValueError as e:
        logger.error(f"There was an error with the input parameters to calculate_harshness: {e}")
    except AssertionError as e:
        logger.error(e)

