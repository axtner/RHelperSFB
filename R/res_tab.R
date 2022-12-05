#' Aggregate assignment results per sample and replicate ID. 
#' 
#' Lists the assignments for the different taxonomic levels and calculated the mean log probability for those assignments made by PROTAX. Also lists the reference with the highest similarity and the average similarity score.
#' 
#' 


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
  
  
  q_tab = DBI::dbGetQuery(db, paste0(
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
  