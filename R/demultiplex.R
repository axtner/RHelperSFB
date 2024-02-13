#' Demultiplex double indexed Illumina sequencing reads.
#' 
#' demultiplex 
#' 
#' Function deploys the read_preprocessing.sh script from the screenforbio pipeline of Alex Crampton-Platt. The original bash script has been edited to eliminate bugs. The bash script uses \code{bcl2fastq} (1), \code{cudadapt} (2), \code{adapterremoval} (3) and \code{vsearch} (4). All four programmes must be installed on the computer running the analysis and their superior directory where their executable files are saved must be added to the PATH variable of the environment (e.g. use base::Sys.setenv() to add them in R).  
#' 
#' @param seq_dir Run folder of an Illumina sequencing run.
#' @param seq_run Sequencing run name (e.g. seq_run = "LeechSeq001").
#' @param out_dir Custom name of the created output folder that is created in the run folder. Default is "DerepOut_" followed by current date and time (e.g. "DerepOut_20221219_110204")
#' @param loci String decoding the genetic markers that were sequenced as a binary string with each locus coded as either present (1) or absent (0) in the order: 12S, 16S, CytB, COI. Default is set to 16S only loci= "010", 12S and CytB would be 101.
#' @param cycles Integer giving the number of sequencing cycles. Default is cycles = 300.
#' @param index Integer giving length of the indices. Default is index = 5.
#' @param tag_dir Directory containing a text file with a list of sample names and the sample tags for each plate tag. File name format: PlateLabel.txt (e.g. "p5001.txt"). Plate label must match Sample_ID in SampleSheet.csv exactly. File contents format: sample_name fwdTag revTag
#' @param seq_machine Name of the sequencing machine, needed for the bcl2fastq script that is deployed to convert the base calls files to fastq files. default is set to "M01108", the Illumina MiSeq machine of the Leibniz-IZW. 
#' 
#' @examples
#' Simple example to analyse 16S sequencing reads with using th default settings:
#' demultiplex(seq_dir = "/SequencingData/221114_M01108_0115/", seq_run = "SeqRun001")
#' 
#' @source
#' (1) \code{bcl2fastq}: 
#' \url{https://emea.support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html}
#' 
#' (2) \code{cutadapt}:
#' \url{https://doi.org/10.14806/ej.17.1.200}
#' 
#' (3) \code{adapterremoval}:  
#' \url{https://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-016-1900-2}
#' 
#' (4) \code{vsearch}: 
#' \url{https://github.com/torognes/vsearch}
#' 
#' @export
demultiplex = function(seq_dir = NA,
                       seq_run = NA,
                       out_dir = NA,
                       loci = "010",
                       cores = 24,
                       cycles = 300,
                       index = 5,
                       tag_dir = "sample_tags",
                       seq_machine = "M01108"){
  
  # test for input variables
  if(is.na(seq_dir) == T){
    stop("\'seq_dir\' is not set, please provide Illumina run folder.\n")
  }
  if(dir.exists(seq_dir) == F){
    stop(paste0("seq_dir = \"", seq_dir, "\" does not exist. Please check.\n"))
  }
  
  in_dir = paste0(seq_dir, "/Data/Intensities/BaseCalls/")
  if(dir.exists(in_dir) == F){
    stop(paste0("The folder \"", in_dir, "\" that should contain the basecall files does not exist. \n
                Please check the Illumina run folder!\n"))
  }
  
  if(is.na(seq_dir) == T){
    stop("\'seq_run\' is not defined, please provide sequencing run name. \n
         The run name is written to every read header and allows unique read identification.\n")
  }
  
  # set output folder name, if out_dir not specified
  if(is.na(out_dir) == T){
    out_dir = paste0("DerepOut_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  }
  
  
  # setting the  bin folders containing bcl2fastq, AdapterRemoval and cutadapt to PATH 
  Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/home/bioadmin/mambaforge/bin/", sep=":"))
  Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/usr/bin/", sep=":"))

  sh_script = system.file("read_preprocessing_v2.sh", package = "RHelperSFB")
  
  # set working directory 
  setwd(paste0(seq_dir, "/"))
  
  # run bash script
  system(paste0("bash ", sh_script," ", in_dir, " ", seq_run, " ", seq_machine, " ", out_dir, " ", cycles, " ", index, " ", "SampleSheet.csv", " ", loci, " ", tag_dir, " ", cores))
  
  writeLines(paste0("Demultiplexed and dereplicated reads can be found in ", seq_dir,"/", out_dir, "/derep/"))
}
