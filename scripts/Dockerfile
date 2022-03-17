FROM ubuntu:latest
  
RUN apt update
RUN DEBIAN_FRONTEND=noninteractive apt install -y g++ cmake python3 python3-pip python3-venv pkg-config libboost-dev libopenmpi3 libsuitesparse-dev gmsh git git-lfs

RUN CXXFLAGS='-Wno-narrowing' pip install dune-mmesh matplotlib

RUN git clone https://github.com/samuelburbulla/dune-mmesh.git

WORKDIR dune-mmesh/doc/examples
RUN python3 coupling.py

RUN echo 'alias python=python3' >> ~/.bashrc
