FROM python:3.12-slim as builder

#Install requirements
RUN apt-get update
RUN apt-get -y install g++
RUN apt-get install -y gdal-bin libgdal-dev
WORKDIR /app
COPY ./requirements.txt .
RUN pip install --upgrade pip
RUN pip install -r requirements.txt
RUN pip install --no-cache gdal[numpy]==3.6.4

#Copy source code
COPY ./create_harshness_map.py .
COPY ./get_harshness_data.py .

#Create data dir
RUN mkdir /app/data
RUN mkdir /app/data/harshness_maps
RUN mkdir /app/data/waves_annual_data
RUN mkdir /app/data/sea_ice_annual_data
RUN mkdir /app/data/iceberg_annual_data

CMD ["/bin/bash"]
