#' Taxonomic assignment of metabarcoding sequencing reads.
#' 
#' assign 
#' 
#' Function deploys the bash scripts for PROTAX-assignments of metabarcoding sequencing reads from the screenforbio pipeline of Alex Crampton-Platt. It can be used to call either weighted or unweighted PROTAX models. PROTAX models must have been pre-trained using \code{train_protax}. 
#' 
#' @param in_dir Run folder of an Illumina sequencing run.
#' @param out_dir Custom name of the created output folder that is created in the run folder. Default is "AssignOut_" followed by current date and time (e.g. "AssignOut_20221219_110204")
#' 
#' @examples
#' Simple example to analyse 16S sequencing reads with using th default settings:
#' demultiplex(seq_dir = "/SequencingData/221114_M01108_0115/", seq_run = "SeqRun001")
#' 
#' @source
#' (1) \code{bcl2fastq}: 
#' \url{https://emea.support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html}
#' 
#' @export

assign = function(in_dir = NA,
                  out_dir = NA,
                  marker = "16S",
                  model_dir = NA,
                  sfb_dir = NA
                  ){


                    
  
}