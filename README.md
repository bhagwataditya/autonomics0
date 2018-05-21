# Install autonomics

1. Make sure you have R 3.5.0. Install if required.

2. Close all active R sessions. Open a fresh R session

3. Install the latest version of BioConductor
     ```
     remove.packages('BiocInstaller')               # We want the latest version of BiocInstaller
     source("https://bioconductor.org/biocLite.R")  # Install latest version of BiocInstaller
     biocLite(suppressUpdates = TRUE)               # Install latest version of Bioconductor
     ```
    
4. Upgrade packages
     ```
     biocLite("BiocUpgrade")
     ```

5. Install autonomics:
    ```
     source("https://raw.github.com/bhagwataditya/autonomics/master/scripts/install_autonomics.R")
    ```
