FROM python:3.12-slim as builder

#Install requirements
RUN apt-get update
RUN apt-get -y install g++

#Install gdal with conda
RUN apt-get -y install wget
ENV CONDA_DIR /opt/conda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
RUN /bin/bash ~/miniconda.sh -b -p /opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH
RUN conda install -c conda-forge -y libgdal=3.9.1

WORKDIR /app
COPY ./requirements.txt .
RUN pip install --upgrade pip
RUN pip install -r requirements.txt
RUN pip install --no-cache gdal[numpy]==3.9.1

#Set environment variables
ENV PROJ_LIB "/opt/conda/pkgs/proj-9.4.1-hb784bbd_0/share/proj/"

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
