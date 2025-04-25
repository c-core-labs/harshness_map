import os
import logging
import argparse
import json
import re
from cmems_utils import create_cmems_metadata_json
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--data_dir", help="Path to a parent directory in which CMEMs metadata file is/will be stored. (str)")

    data_dir = os.path.join(os.getcwd(), "data")

    args = parser.parse_args()
    if args.data_dir is not None:
        data_dir = args.data_dir

    #Ensure metadata file exists
    metadata_filename = os.path.join(data_dir, "CMEMs_metadata.json")
    cmems_products_to_use = ['GLOBAL_MULTIYEAR_BGC_001_02', 'GLOBAL_MULTIYEAR_PHY_001_030', 'GLOBAL_MULTIYEAR_WAV_001_032']
    if not(os.path.exists(metadata_filename)):
        logger.info(f"CMEMs metadata file, {metadata_filename}, not found. Getting metadata and creating file.")
        create_cmems_metadata_json(cmems_products_to_use,
                                   metadata_filename)
        
    with open(metadata_filename) as metadata_file:
        metadata = json.load(metadata_file)

    print("The following variables are downloaded and available for the indicated years, and thresholds.")
    for variable in os.listdir(data_dir):
        if variable in metadata or variable in ["ibc", "icing"]:
            if variable == "ibc":
                long_name = "Iceberg Concentration"
            elif variable == "icing":
                long_name = "Icing Predictor Index"
            else:
                long_name = metadata[variable]["long_name"]
            print(f"\nVariable: {variable} ({long_name:})")
            year = None
            for file in sorted(os.listdir(os.path.join(data_dir, variable))):
                regex_pattern = fr"^{variable}_(.+?)_(.+).tif"
                match = re.search(regex_pattern, file)
                if match:
                    if match.group(1) != year:
                        year = match.group(1)
                        print(f"Year: {year}")
                        print("Thresholds:")
                    if match.group(2) != "raw_data":
                        print(f"{match.group(2)}") #Threshold


