FROM rocker/geospatial:3.6.2 

RUN wget http://cran.r-project.org/src/contrib/dataRetrieval_2.7.5.tar.gz && \
    R CMD INSTALL dataRetrieval_2.7.5.tar.gz && \
    wget http://cran.r-project.org/src/contrib/daymetr_1.4.tar.gz && \
    R CMD INSTALL daymetr_1.4.tar.gz && \
    wget http://cran.r-project.org/src/contrib/cowplot_1.0.0.tar.gz && \
    R CMD INSTALL cowplot_1.0.0.tar.gz && \
    wget https://cran.r-project.org/src/contrib/lme4_1.1-21.tar.gz && \
    R CMD INSTALL lme4_1.1-21.tar.gz && \
    wget http://cran.r-project.org/src/contrib/lmerTest_3.1-1.tar.gz && \
    R CMD INSTALL lmerTest_3.1-1.tar.gz && \
    wget http://cran.r-project.org/src/contrib/lmerTest_3.1-1.tar.gz && \
    R CMD INSTALL lmerTest_3.1-1.tar.gz && \
    rm -rf *.tar.gz
