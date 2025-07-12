# Use a specific R version from the Rocker project
FROM rocker/r-ver:4.4.1

# Install system dependencies that might be required by R packages
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
gdal-bin \
libgdal-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install renv
RUN R -e "install.packages('renv', repos = 'https://packagemanager.posit.co/cran/latest')"

# Copy the renv.lock file
COPY renv.lock renv.lock

# Restore the R environment
RUN R -e "renv::restore()"

# Copy the rest of the project files
COPY . .

# Expose the port if you are running a Shiny app, for example
# EXPOSE 8080

# Specify the command to run when the container starts
CMD ["R", "-e", "rmarkdown::render('workflow.Rmd')"]
