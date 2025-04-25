FROM python:3.12-slim AS builder


RUN apt-get update && apt-get install -y \
    g++ wget build-essential libsqlite3-dev && \
    apt-get clean


ENV CONDA_DIR=/opt/conda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh
ENV PATH=$CONDA_DIR/bin:$PATH


RUN wget https://www.sqlite.org/2023/sqlite-autoconf-3410200.tar.gz && \
    tar -xzf sqlite-autoconf-3410200.tar.gz && \
    cd sqlite-autoconf-3410200 && \
    ./configure --prefix=/usr && make && make install && \
    cd .. && rm -rf sqlite-autoconf-3410200 sqlite-autoconf-3410200.tar.gz


RUN conda install -c conda-forge -y libgdal=3.9.1 gdal sqlite


ENV PROJ_LIB="/opt/conda/pkgs/proj-9.4.1-h54d7996_1/share/proj/"
ENV GDAL_DATA="/opt/conda/pkgs/libgdal-3.9.1-*/share/gdal"


WORKDIR /app
COPY ./requirements.txt .
RUN pip install --upgrade pip
RUN pip install -r requirements.txt
RUN pip install --no-cache gdal[numpy]==3.9.1


COPY ./*.py .


RUN mkdir /app/data
RUN mkdir /app/data/harshness_maps
RUN mkdir /app/data/oilco

COPY ./data/oilco /app/data/oilco

WORKDIR /app

CMD ["/bin/bash"]
