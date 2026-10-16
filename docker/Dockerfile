FROM rocker/rstudio:4.5.2

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    cmake \
    gcc \
    g++ \
    gfortran \
    make \
    libblas-dev \
    liblapack-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libcairo2-dev \
    libxt-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libharfbuzz-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff5-dev \
    libgit2-dev \
    libglpk-dev \
    libgmp-dev \
    libmagick++-dev \
    libnlopt-dev \
    libuv1-dev \
    libzmq3-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/hcocena

RUN R -q -e "options(timeout = 1200); install.packages(c('remotes', 'BiocManager', 'graphlayouts', 'ellmer'), repos = 'https://cloud.r-project.org')" && \
    R -q -e "options(timeout = 1200); repos <- BiocManager::repositories(); repos['CRAN'] <- 'https://cloud.r-project.org'; options(repos = repos); pkgs <- c('ape', 'treeio', 'ggtree', 'GO.db', 'org.Hs.eg.db', 'org.Mm.eg.db', 'GOSemSim', 'DOSE', 'enrichplot', 'clusterProfiler', 'DESeq2', 'limma', 'sva', 'edgeR', 'tximport', 'apeglm', 'ashr', 'EnhancedVolcano', 'pheatmap', 'CALIBERrfimpute', 'decoupleR', 'dorothea', 'OmnipathR', 'clipr', 'mice', 'progeny', 'RCy3', 'SpatialExperiment', 'GSVA'); for (attempt in seq_len(3)) { missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; if (!length(missing)) break; message('Docker dependency install attempt ', attempt, ': ', paste(missing, collapse = ', ')); BiocManager::install(missing, ask = FALSE, update = FALSE) }; status <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE); print(status); if (!all(status)) stop('Missing Docker dependency packages: ', paste(names(status)[!status], collapse = ', '))"

# Separate layer so that adding a vignette dependency does not invalidate the
# expensive Bioconductor install above.
RUN R -q -e "options(timeout = 1200); repos <- BiocManager::repositories(); repos['CRAN'] <- 'https://cloud.r-project.org'; options(repos = repos); pkgs <- c('BiocStyle', 'knitr', 'rmarkdown'); missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; if (length(missing)) BiocManager::install(missing, ask = FALSE, update = FALSE); status <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE); print(status); if (!all(status)) stop('Missing vignette packages: ', paste(names(status)[!status], collapse = ', '))"

COPY . /opt/hcocena

# build_vignettes = TRUE: the vignette is the package's worked example, so it
# belongs in the image people are told to start from. Building it here also
# proves the whole analysis runs in this container before the image is tagged.
RUN R -q -e "options(timeout = 1200); repos <- BiocManager::repositories(); repos['CRAN'] <- 'https://cloud.r-project.org'; options(repos = repos); remotes::install_local('/opt/hcocena', dependencies = NA, upgrade = 'never', build_vignettes = TRUE)" && \
    R -q -e "status <- vapply(c('graphlayouts', 'ellmer', 'clusterProfiler', 'DESeq2', 'limma', 'sva', 'edgeR', 'tximport', 'apeglm', 'ashr', 'EnhancedVolcano', 'pheatmap', 'org.Hs.eg.db', 'org.Mm.eg.db', 'CALIBERrfimpute', 'decoupleR', 'dorothea', 'OmnipathR', 'clipr', 'mice', 'progeny', 'RCy3', 'SpatialExperiment', 'GSVA', 'hcocena'), requireNamespace, logical(1), quietly = TRUE); print(status); if (!all(status)) stop('Missing Docker packages: ', paste(names(status)[!status], collapse = ', '))" && \
    R -q -e "packageVersion('hcocena')" && \
    R -q -e "v <- vignette(package = 'hcocena')\$results; print(v[, 'Item', drop = FALSE]); if (nrow(v) == 0) stop('The vignette was not installed.')"

# Workspace. The templates live in /usr/local/share/hcocena and the source
# tree is removed, so users only see their own home. The home is filled at
# container start and at the first R session (hcocena-setup-workspace), which
# also works when a platform mounts its own /home/rstudio:
#   ~/README_hcocena.md, ~/reference_files/, ~/projects/example_project/
# Further projects: copy example_project (paths are relative to the project).
RUN T=/usr/local/share/hcocena && \
    mkdir -p $T/project_template/count_data $T/project_template/annotation_data \
             $T/project_template/project_data $T/project_template/output && \
    cp -r /opt/hcocena/docker/reference_files $T/reference_files && \
    cp /opt/hcocena/docker/workspace/README_hcocena.md $T/README_hcocena.md && \
    cp /opt/hcocena/docker/project_workflows/hcocena_main.Rmd $T/project_template/01_hcocena_main.Rmd && \
    cp /opt/hcocena/docker/project_workflows/hcocena_satellite.Rmd $T/project_template/02_hcocena_satellite.Rmd && \
    cp /opt/hcocena/docker/project_data/README.md $T/project_template/project_data/README.md && \
    install -m 0755 /opt/hcocena/docker/workspace/hcocena-setup-workspace /usr/local/bin/ && \
    mkdir -p /etc/cont-init.d && \
    install -m 0755 /opt/hcocena/docker/workspace/99_hcocena_workspace /etc/cont-init.d/ && \
    sed -i 's/\r$//' /usr/local/bin/hcocena-setup-workspace \
                     /etc/cont-init.d/99_hcocena_workspace && \
    tr -d '\r' < /opt/hcocena/docker/workspace/Rprofile.site.append >> "$(R RHOME)/etc/Rprofile.site" && \
    chmod -R a+rX $T && find $T -type f -exec chmod 0644 {} + && \
    rm -rf /opt/hcocena && \
    rm -rf /home/rstudio/* && \
    su rstudio -s /bin/sh -c "/usr/local/bin/hcocena-setup-workspace /home/rstudio" && \
    test -f /home/rstudio/projects/example_project/01_hcocena_main.Rmd

WORKDIR /home/rstudio

EXPOSE 8787
CMD ["/init"]
