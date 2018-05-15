
download_human_reference <- function(organism, version){

   # Satisfy CHECK
   . <- NULL

   if(version >= 76){
         create_dir <- dir.create(sprintf("~/.autonomics/gtf/%s","Homo_sapiens"), recursive=TRUE, showWarnings = FALSE)
         #remote <- 'ftp://ftp.ensembl.org/pub/release-91/gtf/homo_sapiens/Homo_sapiens.GRCh38.91.gtf.gz'
         remote <- sprintf('ftp://ftp.ensembl.org/pub/release-%s/gtf/homo_sapiens/Homo_sapiens.GRCh38.%s.gtf.gz',version,version)
         local <- sprintf("~/.autonomics/gtf/Homo_sapiens/%s", basename(remote))
         utils::download.file(url = remote, destfile = local )
         R.utils::gunzip(local,  remove = FALSE, overwrite = TRUE)
         message(sprintf("GTF file version %s downloaded under ~/.autonomics/gtf/Homo_sapiens", version))
   }
   else {
         create_dir <- dir.create(sprintf("~/.autonomics/gtf/%s","Homo_sapiens"), recursive=TRUE, showWarnings = FALSE)
         remote <- sprintf('ftp://ftp.ensembl.org/pub/release-%s/gtf/homo_sapiens/Homo_sapiens.GRCh37.%s.gtf.gz',version,version)
         local <- sprintf("~/.autonomics/gtf/Homo_sapiens/%s", basename(remote))
         utils::download.file(url = remote, destfile = local )
         R.utils::gunzip(local,  remove = FALSE, overwrite = TRUE)
         message(sprintf("GTF file version %s downloaded under ~/.autonomics/gtf/Homo_sapiens", version))
        }

}
download_mouse_reference <- function(organism, version){

   # Satisfy CHECK
   . <- NULL

   if(version >= 68){
         create_dir <- dir.create(sprintf("~/.autonomics/gtf/%s","Mus_musculus"), recursive=TRUE, showWarnings = FALSE)
         remote <- sprintf('ftp://ftp.ensembl.org/pub/release-%s/gtf/mus_musculus/Mus_musculus.GRCm38.%s.gtf.gz',version,version)
         local <- sprintf("~/.autonomics/gtf/Mus_musculus/%s", basename(remote))
         utils::download.file(url = remote, destfile = local )
         R.utils::gunzip(local,  remove = FALSE, overwrite = TRUE)
         message(sprintf("GTF file version %s downloaded under ~/.autonomics/gtf/Mus_musculus",version))
   }
   else {
         create_dir <- dir.create(sprintf("~/.autonomics/gtf/%s","Mus_musculus"), recursive=TRUE, showWarnings = FALSE)
         remote <- sprintf('ftp://ftp.ensembl.org/pub/release-%s/gtf/mus_musculus/Mus_musculus.NCBIM37.%s.gtf.gz',version,version)
         local <- sprintf("~/.autonomics/gtf/Mus_musculus/%s", basename(remote))
         utils::download.file(url = remote, destfile = local )
         R.utils::gunzip(local,  remove = FALSE, overwrite = TRUE)
         message(sprintf("GTF file version %s downloaded under ~/.autonomics/gtf/Mus_musculus",version))
        }
}


download_rat_reference <- function(organism, version){

   # Satisfy CHECK
   . <- NULL
   if(version >= 80){
         create_dir <- dir.create(sprintf("~/.autonomics/gtf/%s","Rattus_norvegicus"), recursive=TRUE, showWarnings = FALSE)
         remote <- sprintf('ftp://ftp.ensembl.org/pub/release-%s/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.%s.gtf.gz',version,version)
         local <- sprintf("~/.autonomics/gtf/Rattus_norvegicus/%s", basename(remote))
         utils::download.file(url = remote, destfile = local )
         R.utils::gunzip(local,  remove = FALSE, overwrite = TRUE)
         message(sprintf("GTF file version %s downloaded under ~/.autonomics/gtf/Rattus_norvegicus",version))
   }
   else{
        create_dir <- dir.create(sprintf("~/.autonomics/gtf/%s","Rattus_norvegicus"), recursive=TRUE, showWarnings = FALSE)
        remote <-sprintf('ftp://ftp.ensembl.org/pub/release-%s/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_5.0.%s.gtf.gz',version,version)
        local <- sprintf("~/.autonomics/gtf/Rattus_norvegicus/%s", basename(remote))
        utils::download.file(url = remote, destfile = local )
        R.utils::gunzip(local,  remove = FALSE, overwrite = TRUE)
        message(sprintf("GTF file version %s downloaded under ~/.autonomics/gtf/Rattus_norvegicus",version))
       }


}

#' Download gene annotations
#' Download gene annotations in GTF file format
#' @param organism    select orgnaism ('Homo_sapiens' or 'Mus_musculus' or 'Rattus_norvegicus')
#' @param version     select gtf version (By default version = 91)
#' @export
download_gtf <- function(organism, version=91){

   # Assert validity
   assertive.sets::assert_is_subset(organism, c('Homo_sapiens', 'Mus_musculus', 'Rattus_norvegicus'))

   # Organism specific parts
   reference_genomes <-  switch(organism,
                                Homo_sapiens  = download_human_reference(organism, version),
                                Mus_musculus = download_mouse_reference(organism,version),
                                Rattus_norvegicus = download_rat_reference(organism,version))
}




