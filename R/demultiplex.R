#' Demultiplex double indexed Illumina sequencing reads.
#' 
#' demultiplex 
#' 
#' Function deploys the read_preprocessing.sh script from the screenforbio pipeline of Alex Crampton-Platt. The original bash script has been edited to eliminate bugs. The bash script uses [bcl2fastq](https://emea.support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html), [cutadapt](https://doi.org/10.14806/ej.17.1.200), [AdapterRemoval](https://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-016-1900-2) and [usearch](https://www.drive5.com/usearch/). All four programmes must be installed on the computer running the analysis and their superior directory where their executable files are saved must be added to the PATH variable of the environment (e.g. use base::Sys.setenv() to add them in R).  
#' 
#' @param seq_dir Run folder of an Illumina sequencing run.
#' @param seq_run Sequencing run name (e.g. seq_run = "LeechSeq001").
#' @param out_dir Custom name of the created output folder that is created in the run folder. Default is "DerepOut_" followed by current date and time (e.g. "DerepOut_20221219_110204")
#' @param loci String decoding the genetic markers that were sequenced as a binary string with each locus coded as either present (1) or absent (0) in the order: 12S, 16S, CytB, COI. Default is set to 16S only loci= "0100", 12S and CytB would be 1010.
#' @param cycles Integer giving the number of sequencing cycles. Default is cycles = 300.
#' @param index Integer giving length of the indices. Default is index = 5.
#' @param tag_dir Directory containing a text file with a list of sample names and the sample tags for each plate tag. File name format: PlateLabel.txt (e.g. "p5001.txt"). PlateLabel must match Sample_ID in SampleSheet.csv extactly. File contents format: sample_name fwdTag revTag
#' @param seq_machine Name of the sequencing machine, needed for the bcl2fastq script that is deployed to convert the base calls files to fastq files. default is set to "M01108", the Illumina MiSeq machine of the Leibniz-IZW. 
#' 

#' @examples
#' Simple query with requesting the results for BiDoup NP: 
#' rep_tab(sample_id = "VBDC")
#' 
#' @export

demultiplex = function(seq_dir = NA,
                       seq_run = NA,
                       out_dir = NA,
                       loci = "0100",
                       cycles = 300,
                       index = 5,
                       tag_dir = "sample_tags",
                       seq_machine = "M01108"){
  
  # test for input variables
  if(dir.exists(seq_dir) == F){
    stop(paste0("seq_dir = \"", seq_dir, "\" does not exist. Please check."))
  }
  
  in_dir = paste0(seq_dir, "/Data/Intensities/BaseCalls/")
  if(dir.exists(in_dir) == F){
    stop(paste0("The folder \"", in_dir, "\" that should contain the basecall files does not exist. Please check the run folder."))
  }
  
  # set output folder name, if out_dir not specified
  if(is.na(out_dir) == T){
    out_dir = paste0("DerepOut_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  }
  
  
  # setting the  bin folders containing bcl2fastq, AdapterRemoval and cutadapt to PATH 
  Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/usr/local/bcl2fastq2-v2.20.0.422/bin/", sep=":"))
  Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/usr/bin/", sep=":"))

  sh_script = system.file("read_preprocessing_v2.sh", package = "RHelperSFB")
  
  # set working directory 
  setwd(paste0(seq_dir, "/"))
  
  # run bash script
  system(paste0("bash ", sh_script," ", in_dir, " ", seq_machine, " ", out_dir, " ", cycles, " ", index, " ", "SampleSheet.csv", " ", loci, " ", tag_dir))
  
  
}