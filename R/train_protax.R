#' train Protax models for taxonomic asignment of metabarcoding reads.
#' 
#' train_protax 
#' 
#' Function deploys the train_weighted_protax.sh script from the screenforbio pipeline of Alex Crampton-Platt. It will take the taxon from the protax taxonomy file name and assumes curated database FASTA files are in current directory and labelled with format '*.final_database.marker.fa'. The bash script uses \code{last}, which must be installed on the computer running the analysis and their superior directory containing its their executable files must be added to the PATH variable of the environment (e.g. use base::Sys.setenv() to add them in R).
#' 
#' @param in_dir Input directory containing the reference taxonomy and reference databases containing the reference sequences. The latter must be fasta files and follow the name schema 'taxon.database.date.marker.fa' (e.g. 'Tetrapoda.database.20220315.12S.fa').
#' @param out_dir Optional output directory. If not provided an output folder including date and time is created (e.g. 'TrainProtax_20240722_170211')
#' @param splist splist is a list of expected species to use in weighing in the format 'Genus,species' (e.g. 'Homo,sapiens').
#' @param taxon Reference taxonomy used for Protax training. File must follow the name schema '*.protax_taxonomy.date.txt' (e.g. 'Tetrapoda.protax_taxonomy.20220318.txt').
#' @param marker Optinal character giving the sequenced genetic markers. By default it is only "16S", but "12S" and "cytB" are also possible, e.g. 'marker = c("16S", "12S", "cytB")' would use all three options.
#' @param weight Optional numeric variable giving the weight for the expected species. Must be between 0 and 1.
#' @param cores Optional parameter that can be used to parallelize 'lastal' in step 3, which speeds-up processing. Default is set 'cores = 1', which means a single thread is running.
#' 
#' @examples
#' Simple example to train a Protax model for 126S reads:
#' trainProtax(splist = "/training_data/Species_Vietnam.txt/")
#' 
#' @export

train_protax = function(in_dir = NA,
                        splist = NA,
                        out_dir = NA,
                        marker = "16S",
                        weight = 0.9,
                        taxon= "Tetrapoda",
                        cores = 1
                        ){
  # set path to lastal ----
  Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/opt/anaconda/anaconda3/bin/lastal", sep=":"))
  model_date = format(Sys.Date(),"%Y%m%d")
  
  
  # setting path to directory containing the perl scripts and testing for it ----
  scripts = system.file("", package ="RHelperSFB")
  if(dir.exists(scripts) == F){
    stop(paste0("Tried to find folder of the RHelperSFB package at '", scripts, "' and could not find it. Please make sure the package is installed."))
  }

  
  # setting input directory ----
  if(is.na(in_dir) == T){
    message("Please provide input directory containing the relevant data")
    in_dir = readline("Provide path to the needed reference data:")
  }
  setwd(in_dir)
  writeLines(paste0("Input directory is :", in_dir))
  
  
  # files in in_dir'
  files = list.files(in_dir, full.names = T)
  
  # check for species list ----
  if(is.na(splist) == T){
    writeLines("Please provide species list used for weighting of assignments.")
    splist = utils::select.list(basename(files), 
                                title = "\nSelect species list:", graphics = FALSE)
    splist <<- files[grepl(splist, files)]
  }
  writeLines(paste0("List of weighted species is :", splist))
  
  # check for taxonomy file ----
  taxonomy = files[grepl("*.protax_taxonomy.*.txt", files)]
  if(length(taxonomy) == 0){
    stop(paste0("No reference taxonomy in ", in_dir))
  }
  writeLines(paste0("Reference taxonomy file is :", taxonomy))
  
  
  # check for reference databases ----
  refs = files[grepl("*.database.*.fa", files)]
  if(length(refs) == 0){
    stop(paste0("No reference database files in ", in_dir))
  }
  
  
  # setting output directory ----
  if(is.na(out_dir) == T){
    out_dir = paste0("TrainProtax_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  }
  dir.create(out_dir)
  writeLines(paste0("Main output directory is :", in_dir, "/", out_dir))
  
  for(m in marker){
    # create subdir for each marker
    dir.create(paste0(out_dir, "/", m))
    }
  writeLines(paste0("Subdirectories for each marker created inside ", out_dir))
  
  
  # start message ----
  writeLines(paste0("\nYou started to train a new PROTAX model for the taxonomic group '", taxon, "' and for the marker(s) '", paste0(marker, collapse = "', '"),"'."))
  
  # Step 1, making taxonomy priors ----
  writeLines("\nStep 1, making taxonomy priors...")
  
  # building taxonomy
  step1.1 = paste0("perl ", system.file("maketaxonomy.pl", package = "RHelperSFB", mustWork = T), " ", taxonomy, " > ", paste0(out_dir,"/taxonomy"))
  system(step1.1)
  if(file.size(paste0(out_dir,"/taxonomy"))==0 | file.exists(paste0(out_dir,"/taxonomy"))==F){
    stop("Step 1.1 failed!")
  } else {writeLines("Step 1.1 done")}
  
  # setting taxonomy priors
  step1.2 = paste0("perl ", system.file("taxonomy_priors.pl", package = "RHelperSFB", mustWork = T), " ", paste0(out_dir,"/taxonomy"), " > ", paste0(out_dir,"/tax4"))
  system(step1.2)
  if(file.size(paste0(out_dir,"/tax4"))==0 | file.exists(paste0(out_dir,"/tax4"))==F){
    stop("Step 1.2 failed!")
  } else {writeLines("Step 1.2 done")}
    
  # filter tax4 for expected species
  #step1.3 = paste0("grep -f ", splist, " ", paste0(out_dir,"/tax4"), " > ", paste0(out_dir,"/expected_sp"))
  #system(step1.3)
  sp_list = read.table(splist, col.names="species")
  tax4 <- readLines(paste0(out_dir,"/tax4"))
  matching_lines <- c()
  for (species in sp_list$species) {
    matches <- grep(species, tax4, value = TRUE)
    matching_lines <- c(matching_lines, matches)
  }
  writeLines(matching_lines, con= paste0(out_dir,"/expected_sp"))
  #rm(sp_list, tax4, expected_sp)
  if(file.size(paste0(out_dir,"/expected_sp"))==0 | file.exists(paste0(out_dir,"/expected_sp"))==F){
    stop("Step 1.3 failed!")
  } else {writeLines("Step 1.3 done")}
    
  # set weighted taxonomy priors for expected species
  file.exists(paste0(out_dir,"/tax4"))
  step1.4 = paste0("perl ",  system.file("setpriors.pl", package = "RHelperSFB", mustWork = T), " ", weight, " ",  paste0(out_dir,"/expected_sp"), " < ", paste0(out_dir,"/tax4"), " > ", paste0(out_dir,"/wtax4"))
  system(step1.4)
  if(file.size(paste0(out_dir,"/wtax4"))==0 | file.exists(paste0(out_dir,"/wtax4"))==F){
    stop("Step 1.4 failed!")
  } else {writeLines("Step 1.4 done")}
    
  for(n in 1:3){
    step1.5 = paste0("perl ", system.file("thintaxonomy.pl", package = "RHelperSFB", mustWork = T)," ", n, " ", paste0(out_dir,"/wtax4") , " > ",  paste0(out_dir,"/wtax", n))
    system(step1.5)
    if(file.size(paste0(out_dir,"/wtax", n))==0 | file.exists(paste0(out_dir,"/wtax", n))==F){
      stop(paste0("Step 1.5 failed at taxonomic level ", n, "!"))
      } 
  }
  writeLines("Step 1.5 done")
    
  # matching references with taxonomy
  for(m in marker){
    # select reference database
    ref_data = refs[grepl(m, refs)]
      
    reftax = paste0("perl ", system.file("initialseqid2tax.pl", package = "RHelperSFB", mustWork = T), " ", paste0(out_dir,"/taxonomy"), " ",  ref_data, " > ", paste0(out_dir, "/", m, "/", m, ".seq2tax4"))
      system(reftax)
    }
    
  file.remove(paste0(out_dir,"/expected_sp"))
  file.remove(paste0(out_dir,"/tax4"))
    
    
  # Step 2, generating training ----
  writeLines("\nStep 2, generating training data for the different taxonomic levels...")
    
  # loop through markers
  for(m in marker){
    writeLines(paste0("\nAnalysing marker ", m, "..."))
    
    # select reference database
    ref_data = refs[grepl(m, refs)]
    
    # loop through tax levels
    for(level in 1:4){
      writeLines(paste0("\nProcessing taxonomic level ", level, "..."))
      
      step2.1 = paste0("perl ", system.file("makeseqid2tax.pl", package = "RHelperSFB", mustWork = T), " ", level, " ", paste0(out_dir,"/", m, "/", m, ".seq2tax4"), " > ", paste0(out_dir,"/", m, "/", m, ".ref.wtax", level))
      system(step2.1)
      writeLines("Step 2.1 done")
                     
      step2.2 = paste0("perl ", system.file("get1layer_reference_sequences_all.pl", package = "RHelperSFB", mustWork = T), " ", level, " ", paste0(out_dir,"/wtax", level), " ", paste0(out_dir, "/", m, "/", m, ".ref.wtax", level), " ", paste0(out_dir, "/", m, "/", m, ".rseqs", level))
      system(step2.2)
      writeLines("Step 2.2 done")
      
      step2.3 = paste0("perl ", system.file("generate_training_data.pl", package = "RHelperSFB", mustWork = T), " ",  paste0(out_dir, "/wtax", level), " ",  paste0(out_dir,"/", m, "/", m, ".ref.wtax", level), " ", paste0(out_dir, "/", m, "/", m, ".rseqs", level), " 4500 1 no ", paste0(out_dir,"/", m, "/", m, ".train.w.level", level))
      system(step2.3)
      writeLines("Step 2.3 done")
      
      step2.4 = paste0("perl ", system.file("generate_unk_training_data.pl", package = "RHelperSFB", mustWork = T), " ", level, " ", paste0(out_dir,"/wtax", level), " ", paste0(out_dir, "/", m, "/", m, ".ref.wtax", level), " ", paste0(out_dir, "/", m, "/", m, ".rseqs", level), " 500 1 no ", paste0(out_dir, "/", m, "/", m, ".train.w.unk", level))
      system(step2.4)
      writeLines("Step 2.4 done")
      
      step2.5 = paste0("cat ", paste0(out_dir,"/", m, "/", m, ".train.w.level", level), " ", paste0(out_dir,"/", m, "/", m, ".train.w.unk", level), " > ",  paste0(out_dir,"/", m, "/", m, ".train.w", level))
      system(step2.5)
      writeLines("Step 2.5 done")
      
      step2.6 = paste0("cut -f6 -d\" \" ",  paste0(out_dir,"/", m, "/", m, ".train.w", level), " | sort | uniq > ", paste0(out_dir, "/", m, "/", m, ".train.w", level, ".id"))
      system(step2.6)
      writeLines("Step 2.6 done")
      
      step2.7 = paste0("cat ", paste0(out_dir,"/", m, "/", m,".train.w", level, ".id | sort | uniq > "), paste0(out_dir, "/", m, "/", m, ".train.w.ids"))
      system(step2.7)
      writeLines("Step 2.7 done")
      
      step2.8 = paste0("perl ", system.file("fastagrep.pl", package = "RHelperSFB", mustWork = T), " ", paste0(out_dir, "/", m, "/", m, ".train.w.ids"), " ", ref_data, " > ", paste0(out_dir, "/", m, "/", m, ".training.fa"))
      system(step2.8)
      writeLines("Step 2.8 done")
      
      }
    }
      
  # Step 3, run LAST searches ----
  writeLines("\nStep 3, running similarity searches for each marker with LAST. This might take a while...")
    if(cores == 1){
      writeLines("\nCurrently you are running only a single thread. To speed-up processing you could parallelize lastal by setting 'cores' to a higher number (depending on available cores)")
    }
    
  # loop through markers
  for(m in marker){
    # select reference database
    ref_data = refs[grepl(m, refs)]
      
    step3.1 = paste0("lastdb ", paste0(out_dir, "/", m, "/lastref_", m), " ", ref_data)
    system(step3.1)
    writeLines(paste0("Step 3.1 done for marker ", m))
      
    step3.2 = paste0("lastal -T 1 -a 1 -f 0 -m 1000 -P ", cores, " ", paste0(out_dir, "/", m, "/lastref_", m), " ", paste0(out_dir,"/", m, "/", m, ".training.fa"), " >  ", paste0(out_dir, "/", m, "/", m, ".training.last"))
    system(step3.2)
    writeLines(paste0("Step 3.2 done for marker ", m))
      
    step3.3 = paste0("perl ", system.file("last2sim.pl", package = "RHelperSFB", mustWork = T), " ", paste0(out_dir, "/", m, "/", m, ".training.last"), " >  ",paste0(out_dir, "/", m, "/", m, ".train.w.lastsim"))
    system(step3.3)
    writeLines(paste0("Step 3.3 done for marker ", m))
    
    file.remove(paste0(out_dir, "/", m, "/", m, ".training.last"))
    }
    
    
  # Step 4, making x-matrices for R input ----
  writeLines("\nStep 4, making x-matrices for R input...")
    
  # loop through reference files
  for(m in marker){
      # select reference database
      ref_data = refs[grepl(m, refs)]
    # loop through tax levels
    for(level in 1:4){
      step4.1 = paste0("perl ", system.file("create1layer_xdata4.pl", package = "RHelperSFB", mustWork = T), " ", paste0(out_dir,"/", m, "/", m, ".train.w", level), " ", paste0(out_dir,"/wtax", level), " ", paste0(out_dir,"/", m,"/", m, ".ref.wtax", level), " ", paste0(out_dir,"/", m,"/", m, ".rseqs", level), " ", paste0(out_dir,"/", m, "/", m, ".train.w.lastsim"), " ", paste0(out_dir,"/", m, "/", m, ".train.w", level, ".xdat 1"))
      system(step4.1)
      writeLines(paste0("Step 4.1 done for marker ", m," and taxonomic level", level))
    }
  }
  # Step 5, parameterise model in R ----
  writeLines("\nStep 5, parameterise model in R...\nThis will take a while, you better get a coffee or two.")
  
  # load the script that contains necessary functions like logprior, loglikelihood or adaptiveMCMC    
  source(system.file("amcmc.rcode.txt", package = "RHelperSFB", mustWork = T))
  logprior = compiler::cmpfun(logprior)
  loglikelihood = compiler::cmpfun(loglikelihood)
  adaptiveMCMC = compiler::cmpfun(adaptiveMCMC)
  
  # setting parameters
  num.params=1+4
  ind=1001:2000  
  
  for(m in marker){ 
    writeLines(paste0("\nTraining model for marker ", m))
    
    for(level in c(1,2,3,4)){
      writeLines(paste0("\nWorking on taxonomic level ",level))
      dat=read.xdata(paste0(out_dir, "/", m,"/",m,".train.w",level,".xdat"))
      ppa=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1)
      initstate=initialize.adaptation(ppa$params[2000,])
      ppb=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1,prev.state=initstate)
      initstate=initialize.adaptation(ppb$params[2000,])
      ppc=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1,prev.state=initstate)
      initstate=initialize.adaptation(ppc$params[2000,])
      ppd=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1,prev.state=initstate)
      pdf(paste0(out_dir, "/", m,"/",m,".weighted_training_plot_",m,"_level",level,"a_MCMC.pdf"))
      traceplot.all(ppa,ind,num.levels=1, title="iter1")
      amcmc.diagnostic.plot(ppa)
      dev.off()
      pdf(paste0(out_dir, "/", m,"/",m,".weighted_training_plot_",m,"_level",level,"b_MCMC.pdf"))
      traceplot.all(ppb,ind,num.levels=1, title="iter2")
      amcmc.diagnostic.plot(ppb)
      dev.off()
      pdf(paste0(out_dir, "/", m,"/",m,".weighted_training_plot_",m,"_level",level,"c_MCMC.pdf"))
      traceplot.all(ppc,ind,num.levels=1, title="iter3")
      amcmc.diagnostic.plot(ppc)
      dev.off()
      pdf(paste0(out_dir, "/", m,"/",m,".weighted_training_plot_",m,"_level",level,"d_MCMC.pdf"))
      traceplot.all(ppd,ind,num.levels=1, title="iter4")
      amcmc.diagnostic.plot(ppd)
      dev.off()
      k=which.max(ppa$postli[ind])
      write.postparams(ppa,paste0(out_dir, "/", m,"/",m,".w_mcmc",level,"a"),ind[k])
      k=which.max(ppb$postli[ind])
      write.postparams(ppb,paste0(out_dir, "/", m,"/",m,".w_mcmc",level,"b"),ind[k])
      k=which.max(ppc$postli[ind])
      write.postparams(ppc,paste0(out_dir, "/", m,"/",m,".w_mcmc",level,"c"),ind[k])
      k=which.max(ppd$postli[ind])
      write.postparams(ppd,paste0(out_dir, "/", m,"/",m,".w_mcmc",level,"d"),ind[k])
    }
    
  }
   
  writeLines(paste0("\nTraining protax model has finished.\nYou will find the results in teh respective subfolder for each marker under ", out_dir,".\nHave a nice day!"))   
    
    
    

  
}

