#' Aggregate assignment results per sample and replicate ID.
#' 
#' res_tab 
#' 
#' Lists the assignments for the different taxonomic levels and calculated the mean log probability for those assignments made by PROTAX. Also lists the reference with the highest similarity and the average similarity score. To connect to the PostgreSQL database from outside the IZW domain it would need ssh connection to the database host and port forwarding to 5433. The db_host must be set to "localhost" then. 
#' 
#' @param sample_id String defining a sample ID (optional). Strings can be complete or incomplete sample IDs (e.g. "VDBC" or c("VBDC068_1_1", "VBDC_062_1_2")). Either sample_id or rep_id (or both) must be provided.
#' @param rep_id String defining a replicate ID (optional). Strings can be complete or incomplete replicate IDs (e.g. "rep1" or c("rep1", "rep2")). Either sample_id or rep_id (or both) must be provided.
#' @param dir Custom output directory (optional), default is the current working directory (getwd()).
#' @param file_out Custom name of output file (optional). Output file is a tab delimited text file ("*.txt"). Default creates a file name starting with "AssignResOut", followed by sample_id and rep_id (if provided) and date (e.g. "AssignResOut_20221207_105429.txt" or "AssignResOut_VBDC_068_1_1_rep1_20221207_105429.txt"). 
#' @param db_name String giving the database name (optional). Default is the biodiv database.
#' @param db_user String giving the name of the database user (optional). Default user has only SELECT privileges.
#' @param db_pwd String giving the password of the database user (optional). Default user has only SELECT privileges.
#' @param db_host String giving the internal IP address of the host server of the biodiv database. If you wish to connect to another database you have to adjust db_name and db_host accordingly. To connect to the PostgreSQL database from outside the IZW domain it would need ssh connection to the database host and port forwarding to 5433. The db_host must be set to "localhost" then.
#' @param db_port String giving the port the Postgres server listens to. Default is set to "5433", the port biodiv database listens to.
#' 
#' @export



res_tab = function(sample_id = NA,
                   rep_id = NA,
                   dir = getwd(),
                   file_out = NA,
                   db_name = "biodiv", 
                   db_user = "r_select",
                   db_pwd = "rs#izw22", 
                   db_host = "192.168.2.71",
                   db_port = "5433"
                   ){
  
  ## needed packages
  packages = c("DBI", "RPostgres", "BiocManager", "Biostrings")
  
  ## Now load or install & load all
  invisible(lapply(packages, FUN = function(x) {
    if (!require(x, character.only = TRUE) == F) {
      message(paste(x, "package is already installed and loaded\n"))} else {
        if(x == "Biostrings"){
          BiocManager::install("Biostrings", version= "3.16")} else {
            install.packages(x, dependencies = TRUE)}
        library(x, character.only = TRUE)
      }
  }
  ))
  
  # set output file name, if file_out not specified
  if(is.na(file_out) == T){
    file_out = paste0("AssignResOut", if(is.na(sample_id[1]) == F){paste0("_", sample_id[1])}, if(is.na(rep_id[1]) == F){paste0("_", rep_id[1])}, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
  }
  
  # check if either sample_id or rep_id is provided
  if(is.na(sample_id) == T & is.na(sample_id) == T){
    stop("\nPlease provide either sample_id or rep_id, or both.")
  }
  
  # create sample_id string q_1 for more flexible LIKE ANY query in combination with ARRAY
  for(i in 1 : length(sample_id)){
    q_i = paste0("'", paste0("%", sample_id[i], "%"), "'")
    if(i == 1){
      q_1 = q_i
    }
    if(i > 1){
      q_1 = paste(q_1, q_i, sep = ", ")
    }
  }
  
  # create rep_id string q_2 for more flexible LIKE ANY query in combination with ARRAY
  for(i in 1 : length(rep_id)){
    q_i = paste0("'", paste0("%", rep_id[i], "%"), "'")
    if(i == 1){
      q_2 = q_i
    }
    if(i > 1){
      q_2 = paste(q_2, q_i, sep = ", ")
    }
  }
  
  # establish db connection
  db = DBI::dbConnect(RPostgres::Postgres(), host = db_host, user = db_user, password = db_pwd, dbname = db_name, port = db_port)
  
  query = dbSendQuery(db, paste0(
    "
    SELECT
      tab1.sample_id, 
      tab1.rep_id, 
      tab1.family, 
      to_char((sum(tab1.log_prob_family * tab1.frequency)/sum(tab1.frequency)), '9D99') AS \"mean log prob family\",
      tab1.genus, 
      to_char((sum(tab1.log_prob_genus * tab1.frequency)/sum(tab1.frequency)), '9D99') AS \"mean log prob genus\",
      tab1.species, 
      to_char((sum(tab1.log_prob_species * tab1.frequency)/sum(tab1.frequency)), '9D99') AS \"mean log prob species\",
      sum(tab1.frequency) AS \"assigned reads\", 
      tab1.\"best reference\", 
      to_char((sum(tab1.sim_score * tab1.frequency)/sum(tab1.frequency)), '9D99') AS \"mean sim score\" 
    FROM 
      (SELECT
         reads.sample_id,  
         reads.rep_id, 
         reads.read_id AS \"read ID\", 
         reads.frequency, 
         fam.taxon AS family, 
         fam_assign.log_prob_family, 
         gen.taxon AS genus, 
         gen_assign.log_prob_genus,  
         spec.taxon AS species, 
         spec_assign.log_prob_species, 
         best_sim_spec AS \"best reference\", 
         sim_score 
       FROM  
         reads 
         LEFT JOIN similarity ON reads.read_id = similarity.read_id 
         LEFT JOIN spec_assign ON reads.read_id=spec_assign.read_id 
         LEFT JOIN taxonomy spec ON species_id = spec.tax_id 
         LEFT JOIN gen_assign ON reads.read_id = gen_assign.read_id 
         LEFT JOIN taxonomy gen ON genus_id = gen.tax_id 
         LEFT JOIN fam_assign ON reads.read_id = fam_assign.read_id 
         LEFT JOIN taxonomy fam ON family_id = fam.tax_id
         ",
     if(is.na(sample_id[1]) == F | is.na(rep_id[1]) == F){
       paste0(
      "WHERE ",
         if(is.na(sample_id[1]) == F){paste0(
         "reads.sample_id LIKE ANY (ARRAY[", q_1, "])")},
         if(is.na(rep_id[1]) == F){paste0( 
      " AND 
         reads.read_id LIKE ANY (ARRAY[", q_2, "])")}
       )},
         
       "GROUP BY
         reads.sample_id,  
         reads.rep_id,
         reads.read_id,
         reads.frequency, 
         family, 
         log_prob_family,
         genus, 
         log_prob_genus, 
         species, 
         log_prob_species, 
         \"best reference\",
         sim_score
       
       ) 
       AS 
         tab1
    GROUP BY 
      tab1.sample_id,
      tab1.rep_id, 
      tab1.family,
      tab1.genus, 
      tab1.species, 
      tab1.\"best reference\"
    
    ORDER BY 
      sample_id, 
      rep_id, 
      \"assigned reads\" DESC, 
      family,
      genus, 
      species, 
      \"best reference\"
    ")
    )
  
  q_tab = DBI::dbFetch(query)
  
  DBI::dbClearResult(query)
  
  if(nrow(q_tab) == 0){
    message(paste0("No results for found. Please check for correct sample_id or rep_id if provided"))
  } else {
    # save q_tab as comma tab-separated text file
    write.table(q_tab, file = file_out, sep = "\t")
    
    message(paste0("\nYou queried ", nrow(q_tab), " assignments from the \"", db_name, "\" database. \nResults are written to the ", paste0("\"", dir, "/", file_out, "\""), " file.", "\n"))
  } 
  # disconnect from database
  DBI::dbDisconnect(db)
  
  rm("q_1", "q_2", "q_tab", "db")
  gc()
}
  