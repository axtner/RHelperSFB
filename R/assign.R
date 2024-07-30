#' Taxonomic assignment of metabarcoding sequencing reads.
#' 
#' assign 
#' 
#' Function deploys the bash scripts for PROTAX-assignments of metabarcoding sequencing reads from the screenforbio pipeline of Alex Crampton-Platt. It can be used to call either weighted or unweighted PROTAX models ((1) Somervuo et al. 2017). PROTAX models must have been pre-trained using \code{train_protax}. The bash script uses \code{last} (2) and \code{vsearch} (3), which must be installed on the computer running the analysis and their superior directory containing its their executable files must be added to the PATH variable of the environment (e.g. use base::Sys.setenv() to add them in R).
#' 
#' @param in_dir Folder containing de-multiplexed and de-replicated sequencing reads.
#' @param out_dir Custom name of the created output folder that is created. Default is "AssignOut_" followed by marker and current date and time (e.g. "AssignOut_16S_20221219_110213")
#' @param marker String defining the genetic marker used for metabarcoding. Must be can be one or several of "16S", "12S" or "CytB" (e.g. marker = c("12S", "16S"). Default is marker = "16S".
#' @param model_dir Path to a directory containing PROTAX models and clean databases for all loci. Default uses the latest model from '/home/bioadmin/Protax/models'.
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
#' #' (3) \code{vsearch}: 
#' \url{https://github.com/torognes/vsearch}
#' 
#' @export

assign = function(in_dir = NA,
                  out_dir = NA,
                  marker = "16S",
                  model_dir = NA,
                  num_cores = 20
                  ){
  
  # test for input variables
  if(is.na(in_dir) == T){
    stop("\'in_dir\' not defined. Please provide folder containing dereplicated sequencing reads.")
    }
  if(length(list.files(in_dir, pattern="*filter.derep.fq")) == 0){
    in_dir = list.dirs(in_dir, full.names = T)[grepl("*derep", list.dirs(in_dir, full.names = T))]
    if(length(list.files(in_dir, pattern="*filter.derep.fq")) == 0){
      stop(paste0("in_dir = \"", in_dir, "\" contains no result files matching the filter '*filter.derep.fq'. Please check."))
      }
    }
    
  # use latest model if model_dir is not defined
  if(is.na(model_dir) == T){
    sort(dir("/home/bioadmin/Protax/models", pattern = "clean_dbs"), decreasing =T)[1]
    model_dir = paste0("/home/bioadmin/Protax/models/", sort(dir("/home/bioadmin/Protax/models", pattern = "clean_dbs"), decreasing =T)[1])
    }
  # check if model_dir exists  
  if(dir.exists(model_dir) == F){
    stop(paste0("The folder \"", model_dir, "\" that should contain the pre-trained PROTAX models and the reference taxonomy does not exist."))
  } else {
      model_date = strsplit(model_dir, "_")[[1]][3]
    }
  
  # set path to lastal
  Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/opt/anaconda/anaconda3/bin/lastal", sep=":"))
    
  # set output folder name, if out_dir not specified
  if(is.na(out_dir) == T){
    out_dir = paste0("AssignOut_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  }
  
  # set working directory 
  setwd(dirname(in_dir))
  
  # find input files
  files = list.files(in_dir, pattern="*filter.derep.fq", full.names = T)
  
    
  # loop throug markers, each loop calls the assignment shell script 
  for(m in marker){
    #create output dir for each marker
    if(dir.exists(paste0(out_dir, "_", m)) == F){
      dir.create(paste0(out_dir, "_", m))
    }
    
    # loop through the files
    for(file in files[grepl(m, files)]){
      
      # sequencing sample name
      seq_sample = gsub(".filter.derep.fq", "", basename(file))
      message(paste0("\nProcessing ", m, " reads of ", seq_sample,":"))
      
      #create output dir for each marker
      if(dir.exists(paste0(out_dir, "_", m, "/", seq_sample)) == F){
        dir.create(paste0(out_dir, "_", m, "/", seq_sample))
      }
      
      # convert fastq to fasta file
      convert = paste0("usearch -fastq_filter ", file, " -fastaout ", out_dir, "_", m, "/", seq_sample, "/", seq_sample, ".fa")
      writeLines("\nConverting fastq file to fasta file...")
      system(convert)
      if(file.size(paste0(out_dir, "_", m, "/", seq_sample, "/", seq_sample, ".fa")) == 0){
        message(paste0("No results for ", seq_sample, " skipping to next sample."))
        next
      }
      
      
      # run LAST search
      last = paste0("lastal -T 1 -a 1 -f 0 -m 1000 -P20 ", model_dir, "/w_model_", model_date, "_", m, "/lastref_", m, " ", out_dir, "_", m, "/", seq_sample, "/", seq_sample, ".fa >", out_dir, "_", m, "/", seq_sample, "/", seq_sample, ".last")
      writeLines("\nRunning LAST search...")
      system(last)
      
      # convert last result to sequence similarity
      conv_last = paste0("perl ", system.file("last2sim.pl", package = "RHelperSFB", mustWork = T), " ", out_dir, "_", m, "/", seq_sample, "/", seq_sample, ".last >",  out_dir, "_", m, "/", seq_sample, "/", seq_sample, ".lastsim")
      writeLines("\nConverting last result to sequence similarity...")
      system(conv_last)
      if(file.size(paste0(out_dir, "_", m, "/", seq_sample, "/", seq_sample, ".lastsim")) == 0){
        message(paste0("No results for ", seq_sample, " skipping to next sample."))
        next
      }
      
      # get read IDs
      fa_file <- readLines(paste0(out_dir, "_", m, "/", seq_sample, "/", seq_sample, ".fa"))
      ids <- substring(grep("^>", fa_file, value = TRUE), 2)
      writeLines(ids, paste0(out_dir, "_", m, "/", seq_sample, "/", seq_sample, ".ids"))
      
      # Generate base log probability file
      base_log = paste0("perl ", system.file("testsample2init.pl", package = "RHelperSFB", mustWork = T), " ",  out_dir, "_", m, "/", seq_sample, "/", seq_sample, ".ids >", out_dir, "_", m, "/", seq_sample, "/", seq_sample, ".0.logprob")
      system(base_log)
      
      # classify reads at each taxonomic level 1 to 4 
      for(level in c(1:4)){
        if(level == 0){
          l_name = "class"
        }
        if(level == 1){
          l_name = "order"
        }
        if(level == 2){
          l_name = "family"
        }
        if(level == 3){
          l_name = "genus"
        } 
        if(level == 4){
          l_name = "species"
        }
        writeLines(paste0("\nClassify at tax level ", l_name, " ..."))
        prelevel = level - 1
        ifile = paste0(out_dir, "_", m, "/", seq_sample, "/", seq_sample, ".", prelevel, ".logprob")
        ofile = paste0(out_dir, "_", m, "/", seq_sample, "/", seq_sample, ".", level, ".logprob")
        
        protax = paste0("perl ", system.file("classify4.pl", package = "RHelperSFB", mustWork = T), " ", ifile, " ",  model_dir, "/w_model_", model_date, "_", m, "/wtax", level, " ", model_dir, "/w_model_", model_date, "_", m, "/ref.wtax", level, " ", model_dir, "/w_model_", model_date, "_", m, "/rseqs", level, " ", model_dir, "/w_model_", model_date, "_", m, "/w_mcmc", level, " map ", out_dir, "_", m, "/", seq_sample, "/", seq_sample, ".lastsim 0 .05 ", ofile, " 1")
        system(protax)
      }
      
      # Adding taxonomy
      for(level in c(0:4)){
        if(level == 0){
          l_name = "class"
        }
        if(level == 1){
          l_name = "order"
        }
        if(level == 2){
          l_name = "family"
        }
        if(level == 3){
          l_name = "genus"
        } 
        if(level == 4){
          l_name = "species"
        }
        taxonomy = paste0("perl ", system.file("add_taxonomy_info.pl", package = "RHelperSFB", mustWork = T), " ", model_dir, "/w_model_", model_date, "_", m, "/taxonomy ", out_dir, "_", m, "/", seq_sample, "/", seq_sample, ".", level, ".logprob > ", out_dir, "_", m, "/", seq_sample, "/", seq_sample, ".", l_name, "_probs")
        system(taxonomy)
      }
      
      # adding sequence simlarity of assigned /species/genus (for unassigned reads takes best matching sequence)
      writeLines(paste0("\nParallelizing to add sequence simlarity of assigned /species/genus to result tables. This might take a while..."))
      
      # Set the number of cores to use for parallel processing
      while(parallel::detectCores() < num_cores){
        num_cores <- num_cores/2
      }
      if(parallel::detectCores() == 2){
        num_cores == 1
      }
      
      # Read data
      if(file.size(paste0(out_dir, "_", m, "/", seq_sample, "/", seq_sample, ".species_probs")) != 0){
        species_probs <- read.table(paste0(out_dir, "_", m, "/", seq_sample, "/", seq_sample, ".species_probs"), header = FALSE)
        lastsim <- read.table(paste0(out_dir, "_", m, "/", seq_sample, "/", seq_sample, ".lastsim"), header = FALSE)
      
        # Extract unique IDs
        unique_ids <- unique(species_probs$V1)
      
        # Initialize a parallel backend using doParallel
        cl <- parallel::makeCluster(num_cores)
        doParallel::registerDoParallel(cl)
      
        # Initialize an empty list to store results
        bestsim <- foreach::foreach(id = unique_ids, .combine = 'c') %dopar% {
          # Extract species and probability
          sp <- subset(species_probs, V1 == id)[, 5]
          sp <- sub(".*_,", "", sp)
          sp <- paste(strsplit(sp, ",")[[1]][3], strsplit(sp, ",")[[1]][4], sep = "_")
        
        
          # Check if species ends with "_unk"
          if (grepl("_unk$", sp)) {
            gen <- sub("_unk", "", sp)
            # Find gen
            subset_data <- subset(lastsim, V1 == id & V3 == gen)
          } else {
            # Find sp (for unassigned reads will pick up best matching sequence)
            subset_data <- lastsim[grepl(id, lastsim$V1),]
            subset_data <- subset_data[grepl(sp, subset_data$V2),]
          }
        
          # Find the row with the maximum value in column 3
          want <- subset_data[which.max(subset_data$V3), ]
          # Return the result
          want
        }
      
        # Stop the parallel backend
        parallel::stopCluster(cl)
      
        # Combine the list of data frames into a single data frame
        bestsim_df = as.data.frame(matrix(unlist(bestsim), ncol = 3, byrow = TRUE))
      
        # Write results to file
        write.table(bestsim_df, file = paste0(out_dir, "_", m, "/", seq_sample, "/", seq_sample, ".bestsim"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      
        # merge species_probs with bestsim results in single table
        species_probs_sim = (merge(species_probs, bestsim_df, by="V1", all.x = T))
        write.table(species_probs_sim, file = paste0(out_dir, "_", m, "/", seq_sample, "/", seq_sample, ".species_probs_sim"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      } else {next}
    }
  }
  #rm("in_dir", "out_dir", "marker", "model_dir")
  gc()
}



