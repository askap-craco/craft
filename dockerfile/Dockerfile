# Based on the dockerfile from https://github.com/AA-ALERT/AMBER_docker/blob/master/Dockerfile
# Use NVIDIA Docker image
FROM nvidia/cuda

# Install all necessary system packages
WORKDIR /

RUN apt-get -qq -y update && apt-get -qq -y install \
    build-essential \
    git \
    cmake \
    python-dev \
    python3-dev \
    python-pip \
    unzip \
    wget \
    automake \
    autoconf \
    autotools-dev \
    csh \
    libtool \
    m4


    
RUN apt-get clean

# Install CRAFT

ENV CRAFT=/opt/craft
WORKDIR $CRAFT
RUN git clone https://bitbucket.csiro.au/scm/craf/craft.git -b dadain src
WORKDIR $CRAFT/src
RUN pip install -r requirements.txt
ENV PATH=$CRAFT/src/python:$PATH
ENV PYTHONPATH=$CRAFT/src/python:$PYTHONPATH


# install CUB - FREDDA needs this to compile
WORKDIR $CRAFT/include/
RUN wget https://github.com/NVlabs/cub/archive/1.8.0.zip
RUN unzip 1.8.0.zip
RUN rm 1.8.0.zip
ENV CUB_PATH=$CRAFT/include/cub-1.8.0/cub

# install psrdada
WORKDIR $CRAFT/include
RUN git clone git://git.code.sf.net/p/psrdada/code psrdada 
WORKDIR $CRAFT/include/psrdada

COPY PsrdadaMakefile.am Makefile.am

RUN pwd && mkdir build/ && \
    ./bootstrap && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    make clean

ENV CUDA=/usr/local/cuda

# install FREDDA
ENV FREDDA=$CRAFT/src/cuda-fdmt/cudafdmt/src/
WORKDIR $FREDDA
ENV PATH=$FREDDA:$PATH
RUN cp fredda.mk.template fredda.mk
RUN make -j 24
