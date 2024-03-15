# using UMI setting
  
demultiplex2 = function(seq_dir = NA,
                        seq_run = NA,
                        out_dir = NA,
                        marker = "16S",
                        cores = 24,
                        cycles = 210,
                        paired_end = T,
                        merge = NA,
                        demulti_b = T,
                        demulti_s = T
                        ){
  
  # test for input variables
  if(is.na(seq_dir) == T){
    stop("\'seq_dir\' is not set, please provide Illumina run folder.\n")
  }
  if(dir.exists(seq_dir) == F){
    stop(paste0("seq_dir = \"", seq_dir, "\" does not exist. Please check.\n"))
  } else {
    setwd(seq_dir)
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



  if(dir.exists(out_dir) == F){
    dir.create(out_dir)
  }

  
# demultiplex plates  
dmBatch <- function(){
  w <- paste0("bcl2fastq --input-dir ", in_dir, " --output-dir ", seq_dir, "/", out_dir, "/batches/ --barcode-mismatches 1 --with-failed-reads --minimum-trimmed-read-length ", cycles, " --sample-sheet SampleSheet.csv --min-log-level ERROR --loading-threads ", cores, " --processing-threads ", cores, " --writing-threads ", cores)
  message("\nStep 1\nDeploying 'bcl2fastq' to demultiplex pcr batches.\nThis might take a while...\n")
  system(w)
  # simplifying file names
  message("\nSimplifying file names...\n")
  files = list.files(full.names = T, path = paste0(seq_dir, "/", out_dir, "/batches/"), pattern = ".fastq.gz")
  for(file in files){
    new_name <- paste0(lapply(strsplit(basename(file),"_"), "[[", 1), "_",
                       lapply(strsplit(basename(file),"_"), "[[", 4),
                       lapply(strsplit(basename(file),"_"), "[[", 5))
    new_name = gsub("001.", ".", new_name)
    paste(new_name)
    file.rename(file, paste0(seq_dir, "/", out_dir, "/batches/", new_name))
  }
  message(paste0("\nStep 1 finished. PCR batches are demultiplexed in '", out_dir, "/batches/'."))
}

if(demulti_b == T){
  dir.create(paste0(seq_dir,"/",out_dir, "/batches"))
  dmBatch()
  }


# demultiplex samples of pcr batches  
  batches <- gsub(".txt", "", list.files(paste0(seq_dir, "/sample_tags")))
  dmSample <- function(){
    if(paired_end == T){
      x <- paste0("AdapterRemoval --file1 ", seq_dir,"/", out_dir, "/batches/", batch, "_R1.fastq.gz --file2 ", seq_dir,"/", out_dir, "/batches/", batch, "_R2.fastq.gz --basename ", seq_dir,"/", out_dir, "/samples/", batch, " --barcode-list ", seq_dir,"/sample_tags/", batch, ".txt --barcode-mm-r1 1 --barcode-mm-r2 1 --threads ", cores," --maxn 50")
      } else {
        x <- paste0("AdapterRemoval --file1 ", seq_dir,"/", out_dir, "/batches/", batch, "_R1.fastq.gz --basename ", seq_dir,"/", out_dir, "/samples/", batch, " --barcode-list ", seq_dir,"/sample_tags/", batch, ".txt --barcode-mm-r1 1 --barcode-mm-r2 1 --threads ", cores," --maxn 50")
      }
    system(x)
    }

  if(demulti_s == T){
    dir.create(paste0(out_dir, "/samples"))
    message(paste0("\nStep 2\nDeploying AdapterRemoval to demultiplex to sample level. This might take a while, maybe go and get a coffee..."))
    
    for(batch in batches){
      message(paste0("\nProcessing pcr batch '", batch, "':"))
      dmSample()
    }
    
    
    # simplifying file names
    message("\nSimplifying file names...\n")
    files <- list.files(full.names = T, path = paste0(seq_dir, "/", out_dir, "/samples"), pattern = "truncated")
    samples <<- list()
    for(file in files){
      # processing paired end reads
      if(paired_end == T){
        if(grepl("singleton", file) == F){
          new_name = gsub("pair", "R", basename(file))
          new_name = gsub("truncated", "fq", new_name)
          sample = paste0(lapply(strsplit(basename(file),"\\."), "[[", 1), ".", lapply(strsplit(basename(file),"\\."), "[[", 2))
          samples = unique(c(samples, sample))
          writeLines(paste0("\nRenaming '", basename(file), "' to '", new_name, "'"))
          file.rename(file, paste0(seq_dir, "/", out_dir, "/samples/", new_name))
          
        } else {next}
      }
      # processing single reads
      if(demulti_s == F){
        new_name = gsub("truncated", "R1.fq", basename(file))
        sample = paste0(lapply(strsplit(basename(file),"\\."), "[[", 1), ".",
                        lapply(strsplit(basename(file),"\\."), "[[", 2))
        samples = c(samples, sample)
        writeLines(paste0("\nRenaming '", basename(file), "' to '", new_name, "'"))
        file.rename(file, paste0(seq_dir, "/", out_dir, "/samples/", new_name))
      }
    }
    
    
    # merging read pairs R1/R2 per sample if paired-end reads and primer clipping; looping per sample
    if(paired_end == T){    
        if(exists(paste0(out_dir,"/merge")) == F){
          dir.create(paste0(out_dir, "/merge"))
        }
        for(sample in samples){
          R1 <- paste0(seq_dir, "/", out_dir, "/samples/", sample, ".R1.fq")
          R2 <- paste0(seq_dir, "/", out_dir, "/samples/", sample, ".R2.fq")
          if(sum(grepl("@M",readLines(R1))) == 0){
            message(paste0("\nSample '", sample, "' has no sequence pairs that could be merged, skipping to next sample..."))
            next
          } else {
            # merging reads
            writeLines(paste0("\nSample '", sample, "' has ", sum(grepl("@M",readLines(R1))), " sequence pairs that will be merged."))
            if(!is.na(merge) & merge == "vsearch"){
              y <- paste0("vsearch -fastq_mergepairs ", R1, " -reverse ", R2, " -fastqout ", out_dir, "/merge", sample, ".merge.fq -fastq_minovlen 50 -fastq_maxdiffpct 20 -fastq_maxdiffs 20 -threads", cores)
            }
            if(!is.na(merge) & merge == "flash"){
              y <- paste0("flash ", R1, R2, " --output-prefix=", sample, " --output-directory=", out_dir, "/merge -M 65 --threads ", cores)
            } else {
              y <- paste0("usearch -fastq_mergepairs ", R1, " -reverse ", R2, " -fastqout ", out_dir, "/merge/", sample, ".merge.fq -fastq_minovlen 50 -fastq_maxdiffpct 20 -fastq_maxdiffs 20 -threads ", cores)
            }
            system(y)
            if(length(readLines(paste0(out_dir, "/merge/", sample, ".merge.fq")))==0){
              message("WARNING: Failed to merge existing read pairs.")
              unlink(paste0(out_dir, "/merge/", sample, ".merge.fq"))
            }
          }
          
          # primer clipping
          if(exists(paste0(out_dir,"/primerclip")) == F){
            dir.create(paste0(out_dir, "/primerclip"))
            }
          if("12s" %in% tolower(marker)){
            writeLines("\nClipping of 12S primers")
            z1 <- paste0("cutadapt -j ", cores, " -a AAAAAGCTTCAAACTGGGATTAGATACCCCACTAT...ACACACCGCCCGTCACCCTCTGCAGTCA$ --minimum-length 350 -o ", out_dir, "/primerclip/", sample, ".12S.fq --discard-untrimmed ", out_dir, "/merge/", sample, ".merge.fq")
            writeLines("\nPrimer clipping of 12S...")
            system(z1)
          }
          if("cytb" %in% tolower(marker)){
            writeLines("\nClipping of cytB primers")
            z2 <- paste0("cutadapt -j ", cores, " -a AAAAAGCTTCCATCCAACATCTCAGCATGATGAAA...TGAGGACAAATATCATTCTGAGGGGCTGCAGTTT$ --minimum-length 270 -o ", out_dir, "/primerclip/", sample, ".CytB.fq --discard-untrimmed ", out_dir, "/merge/", sample, ".merge.fq")
            writeLines("\nPrimer clipping of CytB...")
            system(z2)
          }
          if("16s" %in% tolower(marker)){
            writeLines("\nClipping of 16S primers")
            z3 <- paste0("cutadapt -j ", cores, " -a CGGTTGGGGTGACCTCGGA...AGTTACCCTAGGGATAACAGC$ --minimum-length 80 -o ", out_dir, "/primerclip/", sample, ".16S.fq --discard-untrimmed ", out_dir, "/merge/", sample, ".merge.fq")
            writeLines("\nPrimer clipping of 16S...")
            system(z3)
          }
        }
      } #end of paired end merging and primer clipping
  
  # processing single reads
  if(paired_end == F){
    if(exists(paste0(out_dir,"/primerclip")) == F){
      dir.create(paste0(out_dir, "/primerclip"))
    }
    for(sample in samples){
        # primer clipping
        writeLines("\nStep 3\nPrimer clipping")
        if("12s" %in% tolower(marker)){
          writeLines("\nClipping of 12S primers")
          z1 <- paste0("cutadapt -j ", cores, " -a AAAAAGCTTCAAACTGGGATTAGATACCCCACTAT...ACACACCGCCCGTCACCCTCTGCAGTCA$ --minimum-length 350 -o ", out_dir, "/primerclip/", sample, ".12S.fq --discard-untrimmed ", out_dir, "/samples/", sample, ".R1.fq")
          system(z1)
        }
        if("cytb" %in% tolower(marker)){
          writeLines("\nClipping of cytB primers")
          z2 <- paste0("cutadapt -j ", cores, " -a AAAAAGCTTCCATCCAACATCTCAGCATGATGAAA...TGAGGACAAATATCATTCTGAGGGGCTGCAGTTT$ --minimum-length 270 -o ", out_dir, "/primerclip/", sample, ".CytB.fq --discard-untrimmed ", out_dir, "/samples/", sample, ".R1.fq")
          system(z2)
        }
        if("16s" %in% tolower(marker)){
          writeLines("\nClipping of 16S primers")
          z3 <- paste0("cutadapt -j ", cores, " -a CGGTTGGGGTGACCTCGGA...AGTTACCCTAGGGATAACAGC$ --minimum-length 80 -o ", out_dir, "/primerclip/", sample, ".16S.fq --discard-untrimmed ", out_dir, "/samples/", sample, ".R1.fq")
          system(z3)
        }
      }
    } # end of single read primer clipping
  } # end of batch demultiplexing
  
  # quality filtering
  if(exists(paste0(out_dir,"/filter")) == F){
      dir.create(paste0(out_dir, "/filter"))
    }
  for(sample in samples){
    for(m in marker){
      x1 <- paste0("vsearch -fastx_filter ", out_dir, "/primerclip/", sample,".", m, ".fq -fastqout ", out_dir, "/filter/", sample, ".", m, ".filter.fq -fastq_maxee 0.5")
    
      system(x1)
      }
  }
  
  # dereplication
  if(exists(paste0(out_dir,"/derep")) == F){
    dir.create(paste0(out_dir, "/derep"))
  }
  for(sample in samples){
    for(m in marker){
      x2<- paste0("vsearch -fastx_uniques ", out_dir, "/filter/", sample, ".", m, ".filter.fq -fastqout ", out_dir, "/derep/", sample, ".filter.derep.fq -sizeout -strand both -minuniquesize 2 -relabel ", seq_run, ".", m, ".", sample, "_")
      system(x2)
    }
  }
} 
#end of demultiplex2 function

