import copernicusmarine
import os
import logging
from datetime import datetime, timedelta

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

#Documentation on data structures returned from CMEMs:
#https://toolbox-docs.marine.copernicus.eu/en/pre-releases-2.0.0a4/response-types.html#copernicusmarine.CopernicusMarineDataset

def cmems_products_to_dict(products):
    """Converts a copernicusmarine.CopernicusMarineProduct into python dict.
    Parameters:
        products (list(copernicusmarine.CopernicusMarineProduct)): A list of the products to be converted.
    Returns:
        product_dict (dict)"""
    product_dict = {} #Dict to be populated and returned
    for product in products:
        for dataset in product.datasets:
            logger.info(f"Extracting data from {dataset.dataset_id}")
            #Extract verion
            version = dataset.versions[0]
            if len(dataset.versions) > 1:
                logger.warning(f"Multiple versions of {dataset.dataset_id} exist; using {version.label}")
            #Extract part
            part = None
            for p in version.parts:
                if p.name == 'default' or p.name == 'mask': #Want "mask" for static datasets. "default" for others
                    part = p
            if part:
                if len(version.parts) > 1:
                    logger.warning(f"Multiple parts of {dataset.dataset_id} v{version.label} exist; using {part.name}")
                service = None
                for s in part.services:
                    if 'arco' in s.service_name:
                        service = s
                        break #Assumption that variables in arco-geo and arco-time series will be the same
                if service:
                    for variable in service.variables:
                        #Add new entry for each variable
                        if variable.short_name not in product_dict:
                            product_dict[variable.short_name] = {
                                "product": product.title,
                                "long_name": variable.standard_name,
                                "units": variable.units,
                                "datasets": {}
                            }
                        #Get time and depth coordinates specific to dataset
                        start_time = None
                        end_time = None
                        time_step = None
                        time_units = None
                        depth_values = None
                        depth_units = None
                        for coordinate in variable.coordinates:
                            if coordinate.coordinate_id == 'time':
                                start_time = coordinate.minimum_value
                                end_time = coordinate.maximum_value
                                time_step = coordinate.step
                                time_units = coordinate.coordinate_unit
                            if coordinate.coordinate_id == 'depth':
                                depth_values = coordinate.values
                                depth_units = coordinate.coordinate_unit
                        #Add an entry for each dataset
                        product_dict[variable.short_name]["datasets"][dataset.dataset_id] = {
                            "start_time": start_time,
                            "end_time": end_time,
                            "time_step": time_step,
                            "time_units": time_units,
                            "depth_values": depth_values,
                            "depth_units": depth_units
                        }
                else: #No arco-geo-series or arco-time-series in list of services
                    logger.warning(f"No arco series found in {dataset.dataset_id} v{version.label} part {part.name}. Skipping")
            else: #No default or mask part found in list of parts
                logger.warning(f"No default or mask part found in {dataset.dataset_id} v{version.label}. Skipping")
    return product_dict



def get_cmems_product_metadata(product_id):
    """Downloads metadata for the given product id from CMEMs
    Parameters:
        product_id (str): A CMEMs product id.
    Returns:
        product (copernicusmarine.CopernicusMarineProduct)"""
    logger.info(f"Querying CMEMs catalog for products with product_id: {product_id}")
    catalog = copernicusmarine.describe(product_id = product_id)
    products = catalog.products
    if len(products) == 0:
        logger.warning(f"CMEMs returned no products with product_id: {product_id}")
        return None
    elif len(products) == 1:
        product = products[0]
    else:
        try:
            index = [p.product_id for p in products].index(product_id)
        except ValueError as e:
            logger.warning(f"CMEMs returned no products with product_id: {product_id}.")
            logger.debug("CMEMs found these similar product_ids:")
            for product in products:
                logger.debug(f"{product.product_id}")
            return None
        product = products[index]    
    return product

def query_cmems():
    catalog = copernicusmarine.describe()
    for product in catalog.products:
        print(product.title)
        for dataset in product.datasets:
            print("\t" + dataset.dataset_name)
            for version in dataset.versions:
                print("\t\t" + version.label)
                for part in version.parts:
                    print("\t\t\t" + part.name)
                    for service in part.services:
                        print("\t\t\t\t" + service.service_name)
                        for variable in service.variables:
                            print("\t\t\t\t\t" + variable.short_name)

def download_from_cmems(dataset,
                        variables,
                        start_datetime=datetime.today() - timedelta(days=1),
                        end_datetime=datetime.today(),
                        minimum_longitude=None,
                        maximum_longitude=None,
                        minimum_latitude=None,
                        maximum_latitude=None,
                        output_file=None,
                        credentials_file=None):
    """Downloads data from the Copernicus Marine Service using copernicusmarine.subset().
    Parameters:
        dataset (str): The name of the CMEMs dataset to download.
        variables (list[str]): A list of variables to download from the dataset.
        start_datetime (datetime)
        end_datetime (datetime)
        minimum_longitude (float)
        maximum_longitude (float)
        minimum_latitude (float)
        maximum_latitude (float)
        output_file (str): File path where the downloaded dataset will be stored. Must be a NetCDF (.nc) file.
        credentials_file (str): Path to the file containing CMEMs credentials
    Returns:
        output_path (str) or None"""
    if output_file and not output_file.endswith(".nc"):
        logger.error(f"output_file must be a path to a NetCDF file ending in '.nc'. Got output_file: {output_file} Skipping download")
        return None
    if start_datetime > end_datetime:
        logger.error(f"end_datetime must be later than start_datetime. Got start_datetime: {start_datetime}, end_datetime: {end_datetime}. Skipping download")
        return None
    if os.path.exists(output_file):
        logger.info(f"{output_file} already exists. Skipping download")
        return None
    
    #Proceed with download
    logger.info(f"Downloading data from from Copernicus Marine")
    logger.info(f"Dataset: {dataset}")
    logger.info(f"variables: {variables}")
    logger.info(f"start_datetime: {start_datetime}")
    logger.info(f"end_datetime: {end_datetime}")
    logger.info(f"minimum_longitude: {minimum_longitude}")
    logger.info(f"maximum_longitude: {maximum_longitude}")
    logger.info(f"minimum_latitude: {minimum_latitude}")
    logger.info(f"maximum_latitude: {maximum_latitude}")
    logger.info(f"output_file: {output_file}")
    try:
        output_path = copernicusmarine.subset(dataset_id = dataset, 
                                variables = variables, 
                                minimum_longitude = minimum_longitude,
                                maximum_longitude = maximum_longitude,
                                minimum_latitude = minimum_latitude,
                                maximum_latitude = maximum_latitude,
                                minimum_depth = None,
                                maximum_depth = None,
                                start_datetime = start_datetime, 
                                end_datetime = end_datetime, 
                                output_directory = os.path.dirname(output_file), 
                                output_filename = os.path.basename(output_file), 
                                credentials_file=credentials_file)
        logger.info(f"Download complete: {output_path}")
        return(str(output_path))
    except Exception as e:
        logger.error(f"Download failed: {e}")
        os.remove(output_file)
        raise e
        
    
    
    # def download_from_cmems(output_file=None, 
    #                     data_year=None, 
    #                     dataset=None, 
    #                     variables=None,
    #                     credentials_file=None):
    # """Downloads data from the Copernicus Marine Service using copernicusmarine.subset().
    # Parameters:
    #     output_file (str): File path where the downloaded dataset will be stored. Must be a NetCDF (.nc) file.
    #     data_year (int): The dataset will contain data from Jan 1 to Dec 31 of this year.
    #     dataset (str): The name of the CMEMs dataset to download.
    #     variables (list[str]): A list of variables to download from the dataset.
    # Returns:
    #     output_file (str)"""
    # if not os.path.exists(output_file):
    #     data_year_start = f"{data_year}-01-01"
    #     data_year_end = f"{data_year}-12-31T21:00:00"
    #     logger.info(f"Downloading {data_year} {dataset} from Copernicus Marine")
    #     try:
    #         copernicusmarine.subset(dataset_id = dataset, 
    #                                 variables = variables, 
    #                                 start_datetime = data_year_start, 
    #                                 end_datetime = data_year_end, 
    #                                 output_directory = os.path.dirname(output_file), 
    #                                 output_filename = os.path.basename(output_file), 
    #                                 credentials_file=credentials_file)
    #         logger.info(f"Download complete: {output_file}")
    #     except Exception as e:
    #         logger.error(f"Download failed: {e}")
    #         os.remove(output_file)
    #         raise e

    # else: #File already downloaded
    #     logger.info(f"{output_file} already exists. Skipping download")
    # return(output_file)

def download_original_files_from_cmems(output_dir=None, 
                                       data_year=None, 
                                       dataset=None,
                                       credentials_file=None):
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
                             force_download = True, #TODO: Replace/remove, this is deprecated
                             credentials_file=credentials_file)
        logger.info(f"Download complete: {output_dir}")
    else: #File already downloaded
        logger.info(f"{output_dir} already exists. Skipping download")    
    return(output_dir)
