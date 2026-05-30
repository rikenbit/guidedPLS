# Base Image
FROM bioconductor/bioconductor_docker:devel

# Install R Packages
RUN R -e "pak::pak('rikenbit/guidedPLS', \
    lib = .libPaths()[1]); \
    tools::testInstalledPackage('guidedPLS')"