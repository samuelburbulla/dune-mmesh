FROM ubuntu:latest

RUN apt update
ENV DEBIAN_FRONTEND=noninteractive
RUN apt install -yq g++ cmake python3 python3-pip python3-venv pkg-config libboost-dev libopenmpi-dev libsuitesparse-dev gmsh git git-lfs > /dev/null

RUN python3 -m venv venv
RUN alias python="venv/bin/python3"
RUN alias pip="venv/bin/python3 -m pip"

ENV CXXFLAGS='-Wno-narrowing'
RUN pip install dune-common
RUN pip install dune-geometry
RUN pip install dune-grid
RUN pip install dune-istl
RUN pip install dune-localfunctions
RUN pip install dune-alugrid
RUN pip install dune-fem
RUN pip install dune-mmesh
RUN pip install matplotlib
RUN pip install mpi4py

RUN git clone https://github.com/samuelburbulla/dune-mmesh.git
WORKDIR dune-mmesh/doc/examples
RUN python coupling.py

RUN echo 'alias python=python3' >> ~/.bashrc
