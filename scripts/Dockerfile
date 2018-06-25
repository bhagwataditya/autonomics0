# Run like so:
#  docker build --squash -t autonomics:latest .

# 'Squash' (experimental as of June 25th 2018), delivers result as a single layer
# without overlays. This is preferrable for later saving/exchanging.
# To enable '--squash' configure docker daemon like so (ArchLinux):
#  nano /usr/lib/systemd/system/docker.service
#    add --experimental=true to the end of the ExecStart=/usr/bin/dockerd -H fd:// line
#  systemctl daemon-reload
#  systemctl restart docker.service

FROM rocker/rstudio:latest
MAINTAINER Johannes Graumann (johannes.graumann@mpi-bn.mpg.de)

RUN \
    # Install additional tools
    ## git-lfs in current base image requires backports
    echo 'deb http://http.debian.net/debian stretch-backports main' > /etc/apt/sources.list.d/stretch-backports-main.list && \
    ## Assure up-to-date package lists
    apt-get update && \
    ## Assure up-to-date environment
    apt-get dist-upgrade --assume-yes && \
    apt-get install --assume-yes --no-install-recommends \
      ## Fetch the additional tools
      git-lfs ssh-client \
      # Install libraries needed for compilation of the autonomics toolkit
      ## Dependencies of rgl (which is a dpendency of mixOmics)
      libx11-dev libglu1-mesa-dev \
      ## Dependencies of roxygen2
      libxml2-dev zlib1g-dev &&\
    # Clean the cache(s)
    apt-get clean && \
    rm -r /var/lib/apt/lists/* && \
    # Install the autonomics infrastructure
    Rscript -e 'source("https://raw.github.com/bhagwataditya/autonomics/master/scripts/install_autonomics.R")' && \
    # Clean /tmp (downloaded packages)
    rm -r /tmp/*

# To save environment/backup do like so:
#  # Projectify:
#  docker tag autonomics <PROJECT_NAME>
#  # Possibly (if the image is not needed for any other project(s)):
#  docker rmi autonomics
#  # Start the container (& do your work via http://localhost:8787)
#  CONTAINER_ID=$(docker run -d -p 8787:8787 -v /tmp:/tmp -e ROOT=TRUE <PROJECT_NAME>)
#  # After you are done: stop the container
#  docker stop $CONTAINER_ID
#  # Commit the changes to the image:
#  docker commit --message="Work-preserving snapshot" $CONTAINER_ID <PROJECT_NAME>:latest
#  # Save an image that may later be loaded (xz-compressed)
#  docker save --output <PROJECT_NAME>.tar <PROJECT_NAME> | xz
