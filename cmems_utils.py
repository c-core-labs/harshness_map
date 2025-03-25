
import copernicusmarine
import os
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def download_from_cmems(output_file=None, 
                        data_year=None, 
                        dataset=None, 
                        variables=None,
                        credentials_file=None):
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
                                    credentials_file=credentials_file)
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
                             force_download = True,
                             credentials_file=credentials_file)
        logger.info(f"Download complete: {output_dir}")
    else: #File already downloaded
        logger.info(f"{output_dir} already exists. Skipping download")    
    return(output_dir)
