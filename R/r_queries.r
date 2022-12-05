#' Query sequencing reads from database and convert them to fasta file
#' 
#' The get_reads function queries read sequences assigned to a certain taxon and writes them to fasta file. The query is based in sample_id, assigned taxon and taxonomic level. To connect to the PostgreSQL database from outside the IZW domain it would need ssh connection to the database host and port forwarding to 5433. The db_host must be set to "localhost" then.
#' @param sample_id Sample IDs as they are documented in the database (mandatory). String or vector of strings of complete or incomplete IDs as wildcards are added (e.g. c("VBDC", "VBF") or "VBDC_12_1_").
#' @param tax_level Taxonomic level at which you query (mandatory). If you want to query at other taxonomic levels than the default "species" you can set it to either "genus" or "family".
#' @param tax_query String or vector of strings of Latin taxon names that you wish to query (mandatory). Strings can be complete or incomplete taxonomic names (e.g. "Canis lupus" or c("Canis", "Felis")). Keep in mind that tax_query must fit the tax_level your perform your query on.
#' @param read_freq Numeric threshold for the minimum frequency of a read (optional). Default is 100.
#' @param dir Custom output directory (optional), default is the current working directory (getwd()).
#' @param file_out Custom name of output file (optional). Output file is in fasta format ("*.fas") providing sequence ID, taxon and read frequency in the sequence header. Default creates a file name starting with "ReadsOut", followed by tax_query and date (e.g. "ReadsOut_Canis_20221202.fas"). 
#' @param db_name String giving the database name (optional). Default is the biodiv database.
#' @param db_user String giving the name of the database user (optional). Default user has only SELECT privileges.
#' @param db_pwd String giving the password of the database user (optional). Default user has only SELECT privileges.
#' @param db_host String giving the internal IP address of the host server of the biodiv database. If you wish to connect to another database you have to adjust db_name and db_host accordingly. To connect to the PostgreSQL database from outside the IZW domain it would need ssh connection to the database host and port forwarding to 5433. The db_host must be set to "localhost" then.
#' @param db_port String giving the port the Postgres server listens to. Default is set to "5433", the port bidov database listens to.
#' @examples 
#' Simple use with default settings:
#' get_reads(sample_id = c("VDBC144_2_1", "VDBC_185"), tax_query = c("Canis", "Felis"))
#' 
#' Query at family level from different database with different user from the same database host:  
#' get_reads(sample_id = "VDBC", tax_query = "Canidae", tax_level = "family", db_name = "my_db", db_user = "my_db_user", db_pwd = "my_db_password")
#' @export
get_reads = function(sample_id = NA,
                     tax_level = "species", 
                     tax_query = NA,
                     read_freq = 100,
                     dir = getwd(),
                     file_out = NA,
                     db_name = "biodiv", 
                     db_user = "r_select",
                     db_pwd = "rs#izw22",
                     db_host = "192.168.2.71",
                     db_port = "5433"){
  
  ## needed packages
  packages = c("DBI", "RPostgres")
  
  ## Now load or install & load all
  invisible(lapply(packages, FUN = function(x) {
    if (!require(x, character.only = TRUE) == F) {
      message(paste(x, "package is already installed and loaded\n"))} else {
        install.packages(x, dependencies = TRUE)
        library(x, character.only = TRUE)
      }
  }
  ))
  
  # set output file name, if file_out not specified
  if(is.na(file_out) == T){
    file_out = paste0("ReadsOut_", tax_query[1], "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".fas")
  }
  
  
  # test tax_level and define name of the assignment results table
  if(tax_level %in% c("species", "genus", "family") == F) {
    stop("tax_level must be either \"species\", \"genus\" or \"family\"")
  }
  if(tax_level == "species"){
    level_tab = "spec"
  }
  if(tax_level == "genus"){
    level_tab = "gen"
  }
  if(tax_level == "family"){
    level_tab = "fam"
  }
  
  
  # test tax_query
  if(is.na(tax_query)[1] == T) {
    stop("tax_query must be specified as character (e.g. \"Canis lupus\" or c(\"Canis luspus\", \"Felis catus\"))")
  }
  
  # create tax_query string q_1 for more flexible LIKE ANY query in combination with ARRAY
  for(i in 1 : length(tax_query)){
    q_i = paste0("'", paste0("%", tax_query[i], "%"), "'")
    if(i == 1){
      q_1 = q_i
      }
    if(i > 1){
      q_1 = paste(q_1, q_i, sep = ", ")
    }
  }
    
  # create tax_query string q_1 for more flexible LIKE ANY query in combination with ARRAY
  for(i in 1 : length(sample_id)){
    q_i = paste0("'", paste0("%", sample_id[i], "%"), "'")
    if(i == 1){
      q_2 = q_i
      }
    if(i > 1){
      q_2 = paste(q_2, q_i, sep = ", ")
    }
  }
  
  # establish db connection
  db = DBI::dbConnect(RPostgres::Postgres(), host = db_host, user = db_user, password = db_pwd, dbname = db_name, port = db_port)
  
  
  # query database
  q_tab = DBI::dbGetQuery(db, paste0(
    "
    SELECT 
      reads.sample_id, 
      reads.read_id, 
      reads.frequency, 
      reads.seq_read, 
      taxonomy.taxon 
    FROM 
      reads 
      LEFT JOIN ", level_tab, "_assign ON reads.read_id = ", level_tab, "_assign.read_id 
      LEFT JOIN taxonomy ON ", tax_level, "_id = taxonomy.tax_id  
    WHERE 
    reads.frequency >= ", read_freq, 
    " AND taxonomy.taxon LIKE ANY (ARRAY[", q_1, "])
      AND reads.sample_id LIKE ANY (ARRAY[", q_2, "])
    ")
    )
  
  if(nrow(q_tab) == 0){
    message(paste0("No results for \"", sample_id, "\" and \"", tax_query, "\" at the ", tax_level, " level."))
  } else {
    # create fasta file out of q_tab
    fa = character(2 * nrow(q_tab))
    fa[c(TRUE, FALSE)] = sprintf("> %s", paste(q_tab$read_id, q_tab$taxon, q_tab$frequency))
    fa[c(FALSE, TRUE)] = q_tab$seq_read
    writeLines(fa, paste0(dir, "/", file_out))
  
    message(paste0("\nYou queried ", nrow(q_tab), " sequences from the \"", db_name, "\" database. \nResults are written to the ", paste0("\"", dir, "/", file_out, "\""), " file.", "\n"))
  }

  # disconnect from database
  DBI::dbDisconnect(db)
  
  rm("q_1", "q_2", "q_tab", "db", "fa" )
  gc()
  }

#' Query reference sequences from database and convert them to fasta file
#'
#' get_refs
#'
#' The get_refs function, queries reference sequences and writes them to fasta file.
#' The query is based in taxon and genetic marker. To connect to the PostgreSQL database from outside the IZW domain it would need ssh connection to the database host and port forwarding to 5433. The db_host must be set to "localhost" then.
#' @param tax_query String or vector of strings of Latin taxon names that you wish to query (mandatory). Strings can be complete or incomplete taxonomic names (e.g. "Canis lupus" or c("Canis", "Felis")).
#' @param gen_marker String defining the genetic marker, could be either "16S", "12S" or "CytB". Default is set to "16S".
#' @param dir Custom output directory (optional), default is the current working directory (getwd()).
#' @param file_out Custom name of output file (optional). Output file is in fasta format ("*.fas") providing taxon, its reference ID in the sequence header. Default creates a file name starting with "RefsOut", followed by taxon, genetic marker and date (e.g. "RefsOut_Canis_16S_20221202.fas"). 
#' @param db_name String giving the database name (optional). Default is the biodiv database.
#' @param db_user String giving the name of the database user (optional). Default user has only SELECT privileges.
#' @param db_pwd String giving the password of the database user (optional). Default user has only SELECT privileges.
#' @param db_host String giving the internal IP address of the host server of the biodiv database. If you wish to connect to another database you have to adjust db_name and db_host accordingly. To connect to the PostgreSQL database from outside the IZW domain it would need ssh connection to the database host and port forwarding to 5433. The db_host must be set to "localhost" then.
#' @param db_port String giving the port the Postgres server listens to. Default is set to "5433", the port bidov database listens to.
#' @examples
#' Simple query with default settings: 
#' get_refs(tax_query = c("Canis", "Felis"))
#' @export
get_refs = function(tax_query = NA,
                    gen_marker = "16S",
                    dir = getwd(),
                    file_out = NA,
                    db_name = "biodiv", 
                    db_user = 'R_select',
                    db_pwd = "rs#izw22",
                    db_host = "192.168.2.71",
                    db_port = "5433"){
  
  ## needed packages
  packages = c("DBI", "RPostgres")
  
  ## Now load or install & load all
  invisible(lapply(packages, FUN = function(x) {
    if (!require(x, character.only = TRUE) == F) {
      message(paste(x, "package is already installed and loaded"))} else {
        install.packages(x, dependencies = TRUE)
        library(x, character.only = TRUE)
      }
  }
  ))
  
  # set output file name, if file_out not specified
  if(is.na(file_out) == T){
    file_out = paste0("RefsOut_", gen_marker, "_", tax_query[1], "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".fas")
  }
  
  # test tax_level and define name of the assignment results table
  if(gen_marker %in% c("16S", "12S", "CytB") == F) {
    stop("gen_marker must be either \"16S\", \"12S\" or \"CytB\"")
  }
  
  # test tax_query
  if(is.na(tax_query)[1] == T) {
    stop("tax_query must be specified as character (e.g. \"Canis lupus\" or c(\"Canis luspus\", \"Felis catus\"))")
  }
  
  # create tax_query string q_1 for more flexible LIKE ANY query in combination with ARRAY
  for(i in 1 : length(tax_query)){
    q_i = paste0("'", paste0("%", tax_query[i], "%"), "'")
    if(i == 1){
      q_1 = q_i
    }
    if(i > 1){
      q_1 = paste(q_1, q_i, sep = ", ")
    }
  }
  
  # establish db connection
  db = DBI::dbConnect(RPostgres::Postgres(), host = db_host, user = db_user, password = db_pwd, dbname = db_name, port = db_port)

  
  # query database
  q_tab = DBI::dbGetQuery(db, paste0(
    "
    SELECT 
      ref_id, 
      taxon, 
      ref_seq 
    FROM 
      public.references 
    WHERE 
      taxon LIKE ANY (ARRAY[", q_1, "])
      AND 
      gen_marker = ", paste0("'", gen_marker, "'"))
  )
  
  # test if query result is empty
  if(nrow(q_tab) == 0){
    message(paste0("No results for \"", sample_id, "\" and \"", tax_query, "\" at the ", tax_level, " level."))
  } else {
    
    # create fasta file out of q_tab
    fa = character(2 * nrow(q_tab))
    fa[c(TRUE, FALSE)] = sprintf("> %s", paste(q_tab$taxon, q_tab$ref_id))
    fa[c(FALSE, TRUE)] = q_tab$ref_seq
    writeLines(fa, paste0(dir, "/", file_out))
  
    message(paste0("\nYou queried ", nrow(q_tab), " reference sequences from the \"", db_name, "\" reference database. \nResults are written to the ", paste0("\"", dir, "/", file_out, "\""), " file.", "\n"))
  }
  
  # disconnect from database
  DBI::dbDisconnect(db)
  
  rm("q_1", "q_tab", "db", "fa" )
  gc()
}
