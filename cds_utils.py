
import os
import logging
import cdsapi

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

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