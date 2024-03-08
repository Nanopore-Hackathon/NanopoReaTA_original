FROM docker.io/mambaorg/micromamba:1.5-jammy

USER root

COPY environment.yml /opt/conda.yml

RUN micromamba install -y -n base -f /opt/conda.yml && micromamba clean -afy

ENV PATH="${MAMBA_ROOT_PREFIX}/bin:${PATH}" PORT=42420 SHINYDATA="/data"

RUN R -e "install.packages(c('dashboardthemes'), repos='http://cran.rstudio.com/'); BiocManager::install('rnaseqDTU')"

COPY server /opt/server
COPY ui /opt/ui
COPY app.R NanopoReaTA_Rpackages.RDS /opt/

WORKDIR /opt

#ENTRYPOINT ["Rscript", "/opt/app.R"]

USER $MAMBA_USER