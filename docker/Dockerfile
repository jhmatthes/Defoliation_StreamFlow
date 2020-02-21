FROM rocker/geospatial:3.6.2 

RUN R -e "install.packages('dataRetrieval',dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
        R -e "install.packages('daymetr',dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
        R -e "install.packages('cowplot',dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
        R -e "install.packages('lme4',dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
        R -e "install.packages('lmerTest',dependencies=TRUE, repos='http://cran.rstudio.com/')" 

