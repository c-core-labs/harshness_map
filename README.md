# C-CORE Harshness Map

## Installation

1. If not already installed, [download and install Docker](https://docs.docker.com/engine/install/)
2. From the directory you would like to install the Harshness Map package `git clone https://github.com/c-core-labs/harshness_map.git`
3. In the *harshness_map* directory, create a file called .cdsapirc and populate it with your Copernicus Clmate Data Service API URL and Key as per Step 1 here: https://cds.climate.copernicus.eu/how-to-api
4. From the *harshness_map* directory `docker build -t harshness_map .`
5. To download data from Copernicus Marine, you will need to set your login credentials. This only needs to be done once by running the command `docker run --rm -it -v $PWD/data:/app/data harshness_map copernicusmarine login --configuration-file-directory=/app/data` and entering your Copernicus Marine credentials

## Determine Data Available for Download and Processing
The Harshness Map software is capable of downloading various datasets from Copernicus Marine Service (as well some additinal datasets from other sources).
The `print_variables_available_for_download` script can be used to print a list of available variables.
Note that this script lists all variables, some variables may not be available for all time periods.

### To Run
1. From the *harshness_map* directory `docker run --rm -it -v $PWD/data:/app/data harshness_map python -m print_variables_available_for_download`

### Parameters
`--data_dir` Path to a parent directory in which CMEMs metadata file is/will be stored. (str)

### Default Values
```
data_dir = os.path.join(os.getcwd(), "data")  
```

### Example
```
docker run --rm -it -v $PWD/data:/app/data harshness_map python -m print_variables_available_for_download \
--data_dir="./data" \
```

## Downloading and Preprocessing Annual Data
The `get_harshness_data` script is used to download and preprocess data. The script performs three main functions:
1. Download data for given variables and given year
2. Compute nine uniform thresholds based on the range of the data
3. Count the number of days in the year that the variable exceeds these thresholds

Note that at this time, there are two datasets that are treated slightly differently than others:
1. Icing Predictor Index (icing): The Icing Predictor Index (described below) uses predefined thresholds of 'none', 'light', 'moderate', 'heavy', and 'extreme'
rather than using uniformly calculated thresholds.
2. Iceberg Concentration (ibc): Instead of computing thresholds and counting number of days which exceed the thresholds, average annual iceberg concentration is computed instead.
This aligns with how iceberg concentration data is utilized in the computation of the Fleming-Drover Harshness Index (described below).

### To Run
1. From the *harshness_map* directory `docker run --rm -it -v $PWD/data:/app/data harshness_map python -m get_harshness_data`

### Parameters
`--data_year` The year over which to download data and make calculations. (int)  
`--data_dir` Path to a parent directory in which output data will be stored. (str)  
`--variables` Comma separated variables to be downloaded and processed (Variable short names from print_variables_available_for_download). (Comma Separated Strings)  
`--no_cleanup` If present, intermediate files will NOT be removed.  

### Default Values
```
data_year=datetime.today().year-1  
data_dir = os.path.join(os.getcwd(), "data")  
variables = None
```

### Example
```
docker run --rm -it -v $PWD/data:/app/data harshness_map python -m get_harshness_data \
--data_year=2020 \
--data_dir="./data" \
--variables="siconc,VHM0,icing,ibc"
```

## Determine Downloaded Data Available for Harshness Map Creation
The `print_downloaded_variables` script can be used to print a list of variables that have been downloaded and preprocessed along with their corresponding available years and thresholds.

### To Run
1. From the *harshness_map* directory `docker run --rm -it -v $PWD/data:/app/data harshness_map python -m print_downloaded_variables`

### Parameters
`--data_dir` Path to a parent directory in which CMEMs metadata file is/will be stored. (str)

### Default Values
```
data_dir = os.path.join(os.getcwd(), "data")  
```

### Example
```
docker run --rm -it -v $PWD/data:/app/data harshness_map python -m print_downloaded_variables \
--data_dir="./data" \
```

## Creating a Harshness Map
The `create_harshness_map` script is used to generate harshness maps from given data.
By default the Fleming-Drover Harshness Index is used, but users may specify custom variables and formulas for harshness calculations.

### To Run
1. From the *harshness_map* directory `docker run --rm -it -v $PWD/data:/app/data harshness_map python -m create_harshness_map`

### Parameters
`--start_year` The first year of data files to include (inclusive). (int)  
`--end_year` The last year of data to include (inclusive). (int)  
`--data_dir` The parent directory containing the given variable data. (str)  
`--variables` A list of variables to be used in the formula. Also corresponds to the input files that will be searched for in data_dir. (Comma Separated Strings)"  
`--thresholds` A list of thresholds to use corrsponding to the given variables. (Comma Separated Values)"  
`--formula` The formula used for calculating the harshness index using gdalCalc where A, B, C, etc. correspond to the ordered list of variables provided. (str)"  
`--x_resolution` Resolution of the output file in degrees longitude. (float)  
`--y_resolution` Resolution of the output file in degrees latitude. (float)  
`--lon_min` The minimum longitude bound of the output file in degrees. (float)  
`--lat_min` The minimum latitude bound of the output file in degrees. (float)  
`--lon_max` The maximum longitude bound of the output file in degrees. (float)  
`--lat_max` The maximum latitude bound of the output file in degrees. (float)  
`--no_cleanup` If true, intermediate files will NOT be deleted.  

**Default Values:**  
```
start_year = datetime.today().year-1  
end_year = start_year  
data_dir = os.path.join(os.getcwd(), "data")  
variables = ["siconc", "VHM0", "ibc", "icing"]
thresholds = [0.64, 4, None, "light"]
formula = "6*S/350 + 2.5*W/110 + 1.5*(I>0.01)*(12 + 2*log10((I/10000)+1e-40))"  
x_resolution = 0.2  
y_resolution = 0.2  
lon_min = -71  
lat_min = 41  
lon_max = 25  
lat_max = 82  
```

### Example
```
docker run --rm -it -v $PWD/data:/app/data harshness_map python -m create_harshness_map \
--start_year=2020 \
--end_year=2022 \
--data_dir="./data" \
--variables="siconc,VHM0,ibc"
--threholds=0.64,4,None
--formula="4*A/350 + 4.5*B/110 + 1.5*(C>0.01)*(12 + 2*log10((C/10000)+1e-40))" \
--x_resolution=0.4 \
--y_resolution=0.4 \
--lon_min=-60 \
--lat_min=35 \
--lon_max=-20 \
--lat_max=75 \
```

## Harshness Index
A harshness index is a single parameter which takes various data into account to provide a measure of environmental harshness in different ocean regions.  
  
The default harshness index used in this project is the Fleming-Drover Harshness Index which takes into account sea ice concentration, significant wave height, and iceberg density, and is given by the following formula:  
  
`6*A/350 + 2.5*B/110 + 1.5*(C>0.01)*(12 + 2*log10((C/10000)+1e-40))`  
  
Where,  
`A` is the average number of days per year with a sea ice concentration > 60%  
`B` is the average number of days per year with a significant wave height > 4 metres  
`C` is the average annual icebergdensity (# of icebergs per 100 km<sup>2</sup>)  
  
When calling the `create_harshness_map` script, a custom harshness index formula (along with custom input variables) may be passed as an argument.  
A good place to start when trying out custom formulas is changing certain parameters of the default formula.  
The default formula can be parameterized as follows:  
  
`A*S/D + B*W/E + C*(I>T)*(F + G*log10((I/H)+J))`  

Where,  
`A`, `B`, and `C` are weights for `S`, `W`, and `I` respectively, and should sum to 10.  
`D`, `E`, `F`, `G`, and 'H' are normalization factors, and should be adjusted based on the data being used.  
`J` is a small value added to the log10 argument simply to avoid taking log10 of zero.  
`T` is a threshold for iceberg density. Densities < T will be considered 0.  

These parameters can be modified to create custom formulas.  
For example, to put more weight on wave height, and less on sea ice and icebergs, you may try:  
`1*A/350 + 8*B/110 + 1*(C>0.01)*(12 + 2*log10((C/10000)+1e-40))`  
Or, you may only want to consider two variables:  
`5*A/200 + 5*B/10`  
Of course, the formula is completely customizable, so feel free to experiment, here is an example of an acceptable (though probably non-sensical) formula:  
`A*B + (A>C/10) * D/10 + B/100 + 15` 



## Data Acknowledgement
Data utilized in this project are provided by GHRSST, Met Office and CMEMS


