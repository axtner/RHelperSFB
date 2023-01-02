#' Taxonomic assignment of metabarcoding sequencing reads.
#' 
#' assign 
#' 
#' Function deploys the bash scripts for PROTAX-assignments of metabarcoding sequencing reads from the screenforbio pipeline of Alex Crampton-Platt. It can be used to call either weighted or unweighted PROTAX models ((1) Somervuo et al. 2017). PROTAX models must have been pre-trained using \code{train_protax}. The bash script uses \code{last} (2) and \code{usearch} (3), which must be installed on the computer running the analysis and their superior directory containing its their executable files must be added to the PATH variable of the environment (e.g. use base::Sys.setenv() to add them in R).
#' 
#' @param in_dir Folder containing demultiplexed and dereplicated sequencing reads.
#' @param out_dir Custom name of the created output folder that is created. Default is "AssignOut_" followed by marker and current date and time (e.g. "AssignOut_16S_20221219_110213")
#' @param marker String defining the genetic marker used for metabarcoding. Must be can be one or several of "16S", "12S" or "CytB" (e.g. marker = c("12S", "16S"). Default is marker = "16S".
#' @param model_dir Path to a directory containing PROTAX models and clean databases for all loci. Default uses the latest model from '/home/bioadmin/Protax/models'.
#' @param sfb_dir Path to the screenforbio-mbc directory and it's subdirectory '/protaxscripts' (screenforbio pipeline of Alex Crampton-Platt).
#' 
#' @examples
#' Simple example to analyse 16S sequencing reads with using th default settings:
#' assign(in_dir = "/SequencingData/221114_M01108_0115/seqRun001/data/derep")
#' 
#' @source
#' (1) Somervuo et al. Methods in Ecology & Evolution 2017: 
#' \url{https://doi.org/10.1111/2041-210X.12721}
#' 
#' #' @source
#' (2) \code{last}: 
#' \url{https://gitlab.com/mcfrith/last}
#' 
#' #' (3) \code{usearch}: 
#' \url{https://www.drive5.com/usearch/}
#' 
#' @export

assign = function(in_dir = NA,
                  out_dir = NA,
                  marker = "16S",
                  model_dir = NA
                  ){
  # loop throug markers, each loop calls the assignment shell script 
  for(m in marker){
    
    # test for input variables
    if(is.na(in_dir) == T){
    stop("\'in_dir\' not defined. Please provide folder containing dereplicated sequencing reads.")
    }
    in_dir = paste0(in_dir,"/data/derep")
    if(dir.exists(in_dir) == F){
      stop(paste0("in_dir = \"", in_dir, "\" does not exist. Please check."))
    }
  
    # use latest model if model_dir is not defined
    if(is.na(model_dir) == T){
      sort(dir("/home/bioadmin/Protax/models", pattern = "clean_dbs"), decreasing =T)[1]
      model_dir = paste0("/home/bioadmin/Protax/models/", sort(dir("/home/bioadmin/Protax/models", pattern = "clean_dbs"), decreasing =T)[1])
    }
    # check if model_dir exists  
    if(dir.exists(model_dir) == F){
      stop(paste0("The folder \"", model_dir, "\" that should contain the pre-trained PROTAX models and the reference taxonomy does not exist."))
    }
    
    # shell script
    sh_script = system.file("weighted_protax_classify.sh", package = "RHelperSFB", mustWork = T)
    
    # set output folder name, if out_dir not specified
    if(is.na(out_dir) == T){
      out_dir = paste0("AssignOut_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    }

    # set working directory 
    setwd(dirname(in_dir))
    
    
    Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/opt/anaconda/anaconda3/bin/lastal", sep=":"))
    
    # run bash script
    system(paste0("bash ", sh_script," ", in_dir, " ", m, " ", model_dir, " ", dirname(sh_script), " ", out_dir))
                    
  
  }
}
