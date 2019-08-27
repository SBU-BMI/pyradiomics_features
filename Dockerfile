FROM python:3.7.4
MAINTAINER Feiqiao Wang "feiqiao.wang@stonybrook.edu" 

RUN apt-get upgrade && apt-get update && apt-get install -y  \    
  git \
  vim  \
  openslide-tools
  
# copy project source code to app folder
COPY ./app/* /app/
COPY ./radiomics_features_selected.txt /app/radiomics_features_selected.txt
COPY ./run.sh /root

# install python other package such as geojson,pymongo,opencv-python,shapely,openslide-python
RUN pip3.7 install --upgrade pip && pip3.7 install -r /app/requirements.txt

# Setup necessary python packages
RUN pip3.7 install wheel>=0.29.0 && \
    pip3.7 install setuptools>=28.0.0 && \
    pip3.7 install scipy && \
    pip3.7 install trimesh && \
    pip3.7 install scikit-image && \
    pip3.7 install --trusted-host www.itk.org -f https://itk.org/SimpleITKDoxygen/html/PyDownloadPage.html SimpleITK>=0.9.1 && \
    python3.7 -c "import SimpleITK; print('SimpleITK Version:' + SimpleITK.Version_VersionString())"

# install pyradiomics package   
RUN python3.7 -m pip install pyradiomics

RUN chmod 755 /app/computing_patch_level_radiomics_features.sh

WORKDIR /app

CMD ["sh", "/root/run.sh"]
