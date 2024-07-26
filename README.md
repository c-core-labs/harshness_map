# C-CORE Harshness Map

## Installation

1. If not already installed, [download and install Docker](https://docs.docker.com/engine/install/)
2. From the directory you would like to install the Harshness Map package `git clone https://github.com/c-core-labs/harshness_map.git`
3. From the newly created *harshness_map* directory `docker build -t harshness_map .`
4. To download data from Copernicus Marine, you will need to set your login credentials. This only needs to be done once by running the command `docker run --rm -it -v $PWD/data:/app/data harshness_map copernicusmarine login --configuration-file-directory=/app/data` and entering your Copernicus Marine credentials

## Downloading and Preprocessing Annual Data
The `get_harshness_data` script is used to download annual wave height, sea ice concentration, and iceberg density data
from Copernicus Marine Service, and compute the parameters to be used in the computation of the 
Fleming-Drover Harshness Index as follows:
1. The mean annual number of days with a significant wave height greater than a given threshold (1-10 metres)
2. The mean annual number of days with a sea ice concentration greater than a given threshold (0-90 percent in increments of 10)
3. The mean annual open water iceberg areal density (# of icebergs per 100 km^2)

### To Run
1. From the *harshness_map* directory `docker run --rm -it -v $PWD/data:/app/data harshness_map python -m get_harshness_data`

### Parameters
`--data_year` The year over which to download data and make calculations. (int)  
`--data_dir` Path to a parent directory in which output data will be stored. (str)  
`--clean` If true, intermediate files will be deleted. (bool)  

### Default Values
```
data_year=datetime.today().year-1  
data_dir = os.path.join(os.getcwd(), "data")  
clean = True
```

### Example
```
docker run --rm -it -v $PWD/data:/app/data harshness_map python -m get_harshness_data \  
--data_year=2020 \
--data_dir="./data" \  
--clean=True
```
## Creating a Harshness Map
The `create_harshness_map` script is used to generate harshness maps from Wave Height, Sea Ice Concentration, and Iceberg Density data.
By default the Fleming-Drover Harshness Index is used, but users may specify custom formulas for harshness calculations.

### To Run
1. From the *harshness_map* directory `docker run --rm -it -v $PWD/data:/app/data harshness_map python -m create_harshness_map`

### Parameters
`--start_year` The first year of data files to include (inclusive). (int)  
`--end_year` The last year of data to include (inclusive). (int)  
`--data_dir` The parent directory containing the directories 'iceberg_annual_data', 'waves_annual_data', and 'sea_ice_annual_data'. (str)  
`--wave_height_thresh` Used to define wave input files. Input files contain # of days with wave height > <wave_height_thresh> metres. (int 1-10)  
`--sea_ice_concentration_thresh` Used to define sea ice input files. Input files contain # of days with sea ice concentration > <sea_ice_concentration_thresh> percent. (int 0-90 by 10s)  
`--formula` The formula used for calculating the harshness index using gdalCalc where S = Sea Ice Data, W = Wave Data, I = Iceberg Data. (str)  
`--x_resolution` Resolution of the output file in degrees longitude. (float)  
`--y_resolution` Resolution of the output file in degrees latitude. (float)  
`--lon_min` The minimum longitude bound of the output file in degrees. (float)  
`--lat_min` The minimum latitude bound of the output file in degrees. (float)  
`--lon_max` The maximum longitude bound of the output file in degrees. (float)  
`--lat_max` The maximum latitude bound of the output file in degrees. (float)  
`--clean` If true, intermediate files will be deleted. (bool)  

**Default Values:**  
```
start_year = datetime.today().year-1  
end_year = start_year  
data_dir = os.path.join(os.getcwd(), "data")  
wave_height_thresh = 4  
sea_ice_concentration_thresh = 60  
formula = "6*S/350 + 2.5*W/110 + 1.5*(I>0.01)*(12 + 2*log10((I/10000)+1e-40))"  
x_resolution = 0.2  
y_resolution = 0.2  
lon_min = -71  
lat_min = 41  
lon_max = 25  
lat_max = 82  
clean = True  
```

### Example
```
docker run --rm -it -v $PWD/data:/app/data harshness_map python -m create_harshness_map \
--start_year=2020 \
--end_year=2022 \
--data_dir="./data" \  
--wave_height_thresh=2 \
--sea_ice_concentration_thresh=0.2 \
--formula="4*S/350 + 4.5*W/110 + 1.5*(I>0.01)*(12 + 2*log10((I/10000)+1e-40))" \
--x_resolution=0.4 \
--y_resolution=0.4 \
--lon_min=-60 \
--lat_min=35 \
--lon_max=-20 \
--lat_max=75 \
--clean=True
```

## Harshness Index
A harshness index is a single parameter which takes various data into account to provide a measure of environmental harshness in different ocean regions.  
  
The default harshness index used in this project is the Fleming-Drover Harshness Index which takes into account sea ice concentration, significant wave height, and iceberg density, and is given by the following formula:  
  
`6*S/350 + 2.5*W/110 + 1.5*(I>0.01)*(12 + 2*log10((I/10000)+1e-40))`  
  
Where,  
`S` is the average number of days per year with a sea ice concentration > 60%  
`W` is the average number of days per year with a significant wave height > 4 metres  
`I` is the average annual icebergdensity (# of icebergs per 100 km<sup>2</sup>)  
  
When calling the `create_harshness_map` script, a custom harshness index formula may be passed as an argument.  
A good place to start when trying out custom formulas is changing certain parameters of the default formula.  
The default formula can be parameterized as follows:  
  
`A*S/D + B*W/E + C*(I>T)*(F + G*log10((I/H)+J))`  

Where,  
`A`, `B`, and `C` are weights for `S`, `W`, and `I` respectively, and should sum to 10.  
`D`, `E`, `F`, `G`, and 'H' are normalization factors, and should not generally be changed without a knowledge of the data being used.  
`J` is a small value added to the log10 argument simply to avoid taking log10 of zero.  
`T` is a threshold for iceberg density. Densities < T will be considered 0.  

These parameters can be modified to create custom formulas.  
For example, to put more weight on wave height, and less on sea ice and icebergs, you may try:  
`1*S/350 + 8*W/110 + 1*(I>0.01)*(12 + 2*log10((I/10000)+1e-40))`  
Or, you may not want to consider iceberg data at all:  
`5*S/350 + 5*W/110`  
Of course, the formula is completely customizable, so feel free to experiment, here is an example of an acceptable (though probably non-sensical) formula:  
`S*W + (S>W/10) * I/10 + I/100 + 15` 



## Data Acknowledgement
Data utilized in this project are provided by GHRSST, Met Office and CMEMS


