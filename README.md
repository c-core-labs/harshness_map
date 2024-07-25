# C-CORE Harshness Map

## Installation

1. If not already installed, [download and install Docker](https://docs.docker.com/get-docker/)
2. 'git clone https://github.com/c-core-labs/harshness-map.git'
3. From the **harshness** directory 'docker build -t harshness_map .'

## Downloading and Preprocessing Annual Data

1. From the **harshness** directory 'docker run --rm -it -v ./data:/app/data harshness_map python -m get_harshness_data'

## Creating a Harshness Map

1. From the **harshness** directory 'docker run --rm -it -v ./data:/app/data harshness_map python -m create_harshness_map'
