import os
import logging
import argparse
from cmems_utils import create_cmems_metadata_json, get_available_cmems_variables

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
        
    (short_names, long_names) = get_available_cmems_variables(metadata_filename)
    print("The following variables are available for download from CMEMs:")
    print(f"{'Short Name':<15}{'Long Name':<10}")
    for (short_name, long_name) in zip(short_names, long_names):
        print(f"{short_name:<15}{long_name:<10}")
    print("\nIn addition, the following variables are available from other data sources:")
    print(f"{'Short Name':<15}{'Long Name':<10}")
    print(f"{'ibc':<15}{'Iceberg Concentration':<10}")
    print(f"{'icing':<15}{'Icing Predictor Index':<10}")