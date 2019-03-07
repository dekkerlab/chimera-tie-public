FROM continuumio/miniconda
MAINTAINER Mihir Metkar (mihir.metkar@gmail.com)

# Install chimeraTie dependencies using conda
COPY chimeraTie_py2.yml /apps/chimeraTie_py2.yml

COPY scripts chimeratie
COPY test test_data

# Make port 80 available to the world outside this container
EXPOSE 80

RUN conda install python=2.7 && \
    conda env update -n root --file /apps/chimeraTie_py2.yml && \
    rm -rf /opt/conda/pkgs/* &&\
    chmod u+x chimeratie
