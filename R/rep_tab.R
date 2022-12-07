#' Aggregate assignment results per sample and replicate ID.
#' 
#' rep_tab 
#' 
#' Lists the assignments for the different taxonomic levels and calculated the mean log probability for those assignments made by PROTAX. Also lists the reference with the highest similarity and the average similarity score. To connect to the PostgreSQL database from outside the IZW domain it would need ssh connection to the database host and port forwarding to 5433. The db_host must be set to "localhost" then. 
#' 
#' @param sample_id String defining a sample ID or study site (e.g. VBDC = BiDoup National Park). For a list of current study sites and their abbreviations use study_sites(). Strings can be complete or incomplete sample IDs (e.g. "VBDC" or c("VBDC_068_1_1", "VBDC_062_1_2")).
#' @param dir Custom output directory (optional), default is the current working directory (getwd()).
#' @param file_out Custom name of output file (optional). Output file is a tab delimited text file ("*.txt"). Default creates a file name starting with "RepAssignResOut", followed by sample_id and date (e.g. "RepAssignResOut_VDBC_068_1_1_rep1_20221207_105429.txt"). 
#' @param db_name String giving the database name (optional). Default is the biodiv database.
#' @param db_user String giving the name of the database user (optional). Default user has only SELECT privileges.
#' @param db_pwd String giving the password of the database user (optional). Default user has only SELECT privileges.
#' @param db_host String giving the internal IP address of the host server of the biodiv database. If you wish to connect to another database you have to adjust db_name and db_host accordingly. To connect to the PostgreSQL database from outside the IZW domain it would need ssh connection to the database host and port forwarding to 5433. The db_host must be set to "localhost" then.
#' @param db_port String giving the port the Postgres server listens to. Default is set to "5433", the port biodiv database listens to.
#' @examples
#' Simple query with requesting the results for BiDoup NP: 
#' rep_tab(sample_id = "VBDC")
#' 
#' @export

rep_tab = function(sample_id = NA,
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
    file_out = paste0("RepAssignResOut", if(is.na(sample_id[1]) == F){paste0("_", sample_id[1])}, if(is.na(rep_id[1]) == F){paste0("_", rep_id[1])}, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
  }
  
  # check if either sample_id or rep_id is provided
  if(is.na(sample_id) == T & is.na(sample_id) == T){
    stop("\nPlease provide at least partial sample_id to define sample or study area.")
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
  
  # establish db connection
  db = DBI::dbConnect(RPostgres::Postgres(), host = db_host, user = db_user, password = db_pwd, dbname = db_name, port = db_port)
  
  query = dbSendQuery(db,paste0( 
                        "
                        SELECT
                         DISTINCT
                          tab2.sample_id,
                          tab2.species,
                          \"prob 1\",
                          \"reads 1\",
                          \"prob 2\",
                          \"reads 2\",
                          \"prob 3\",
                          \"reads 3\",
                          \"prob 4\",
                          \"reads 4\",
                          \"prob 5\",
                          \"reads 5\",
                          \"prob 6\",
                          \"reads 6\"
                         FROM
                          (SELECT 
                             tab1.sample_id, 
                             tab1.rep_id, 
                             tab1.species, 
                             to_char((sum(exp(tab1.log_prob_species) * tab1.frequency)/sum(tab1.frequency)),
                             '9D99') AS \"mean probability\",
                             sum(tab1.frequency) AS \"assigned reads\" 
                           FROM 
                             (SELECT  
                               reads.sample_id,  
                               reads.rep_id, 
                               reads.read_id AS \"read ID\",
                               reads.frequency, 
                               spec.taxon AS species, 
                               spec_assign.log_prob_species 
                              FROM  
                               reads 
                              LEFT JOIN 
                               spec_assign ON reads.read_id=spec_assign.read_id 
                              LEFT JOIN 
                               taxonomy spec ON species_id = spec.tax_id 
                              GROUP BY
                               reads.sample_id,  
                               reads.rep_id,
                               reads.read_id, 
                               reads.frequency, 
                               species, 
                               log_prob_species) 
                              AS 
                               tab1
                            WHERE 
                              tab1.species IS NOT NULL
                            AND
                              tab1.sample_id LIKE ANY (ARRAY[", q_1, "])
                            GROUP BY 
                              tab1.sample_id, 
                              tab1.rep_id, 
                              tab1.species)
                            AS
                              tab2
                         LEFT JOIN
                         (SELECT
                          DISTINCT
                           sample_id,
                           species,
                           \"mean probability\" AS \"prob 1\", 
                           \"assigned reads\" AS \"reads 1\" 
                          FROM
                           (SELECT 
                             tab1.sample_id, 
                             tab1.rep_id, 
                             tab1.species, 
                             to_char((sum(exp(tab1.log_prob_species) * tab1.frequency)/sum(tab1.frequency)),
                             '9D99') AS \"mean probability\",
                             sum(tab1.frequency) AS \"assigned reads\" 
                           FROM 
                             (SELECT  
                               reads.sample_id,  
                               reads.rep_id, 
                               reads.read_id AS \"read ID\",
                               reads.frequency, 
                               spec.taxon AS species, 
                               spec_assign.log_prob_species 
                              FROM  
                               reads 
                              LEFT JOIN 
                               spec_assign ON reads.read_id=spec_assign.read_id 
                              LEFT JOIN 
                               taxonomy spec ON species_id = spec.tax_id 
                              GROUP BY
                               reads.sample_id,  
                               reads.rep_id,
                               reads.read_id, 
                               reads.frequency, 
                               species, 
                               log_prob_species) 
                              AS 
                               tab1
                            WHERE 
                              tab1.species IS NOT NULL
                            GROUP BY 
                              tab1.sample_id, 
                              tab1.rep_id, 
                              tab1.species)
                            AS
                              tab2 
                          WHERE 
                           rep_id LIKE 'rep1%' AND species IS NOT NULL) 
                          AS rep1
                         ON tab2.sample_id=rep1.sample_id AND tab2.species=rep1.species
                         LEFT JOIN
                         (SELECT 
                          DISTINCt
                           sample_id, 
                           species, 
                           \"mean probability\" AS \"prob 2\", 
                           \"assigned reads\" AS \"reads 2\" 
                          FROM 
                           (SELECT 
                             tab1.sample_id, 
                             tab1.rep_id, 
                             tab1.species, 
                             to_char((sum(exp(tab1.log_prob_species) * tab1.frequency)/sum(tab1.frequency)),
                             '9D99') AS \"mean probability\",
                             sum(tab1.frequency) AS \"assigned reads\" 
                           FROM 
                             (SELECT  
                               reads.sample_id,  
                               reads.rep_id, 
                               reads.read_id AS \"read ID\",
                               reads.frequency, 
                               spec.taxon AS species, 
                               spec_assign.log_prob_species 
                              FROM  
                               reads 
                              LEFT JOIN 
                               spec_assign ON reads.read_id=spec_assign.read_id 
                              LEFT JOIN 
                               taxonomy spec ON species_id = spec.tax_id 
                              GROUP BY
                               reads.sample_id,  
                               reads.rep_id,
                               reads.read_id, 
                               reads.frequency, 
                               species, 
                               log_prob_species) 
                              AS 
                               tab1
                            WHERE 
                              tab1.species IS NOT NULL
                            GROUP BY 
                              tab1.sample_id, 
                              tab1.rep_id, 
                              tab1.species)
                            AS
                              tab2 
                          WHERE 
                           rep_id LIKE 'rep2%' AND species IS NOT NULL) 
                         AS rep2
                         ON tab2.sample_id=rep2.sample_id AND tab2.species=rep2.species
                         
                         LEFT JOIN
                         (SELECT 
                          DISTINCt
                           sample_id, 
                           species, 
                           \"mean probability\" AS \"prob 3\", 
                           \"assigned reads\" AS \"reads 3\" 
                          FROM 
                           (SELECT 
                             tab1.sample_id, 
                             tab1.rep_id, 
                             tab1.species, 
                             to_char((sum(exp(tab1.log_prob_species) * tab1.frequency)/sum(tab1.frequency)),
                             '9D99') AS \"mean probability\",
                             sum(tab1.frequency) AS \"assigned reads\" 
                           FROM 
                             (SELECT  
                               reads.sample_id,  
                               reads.rep_id, 
                               reads.read_id AS \"read ID\",
                               reads.frequency, 
                               spec.taxon AS species, 
                               spec_assign.log_prob_species 
                              FROM  
                               reads 
                              LEFT JOIN 
                               spec_assign ON reads.read_id=spec_assign.read_id 
                              LEFT JOIN 
                               taxonomy spec ON species_id = spec.tax_id 
                              GROUP BY
                               reads.sample_id,  
                               reads.rep_id,
                               reads.read_id, 
                               reads.frequency, 
                               species, 
                               log_prob_species) 
                              AS 
                               tab1
                            WHERE 
                              tab1.species IS NOT NULL
                            GROUP BY 
                              tab1.sample_id, 
                              tab1.rep_id, 
                              tab1.species)
                            AS
                              tab2 
                          WHERE 
                           rep_id LIKE 'rep3%' AND species IS NOT NULL) 
                         AS rep3
                         ON tab2.sample_id=rep3.sample_id AND tab2.species=rep3.species
                         
                         LEFT JOIN
                         (SELECT 
                          DISTINCt
                           sample_id, 
                           species, 
                           \"mean probability\" AS \"prob 4\", 
                           \"assigned reads\" AS \"reads 4\" 
                          FROM 
                           (SELECT 
                             tab1.sample_id, 
                             tab1.rep_id, 
                             tab1.species, 
                             to_char((sum(exp(tab1.log_prob_species) * tab1.frequency)/sum(tab1.frequency)),
                             '9D99') AS \"mean probability\",
                             sum(tab1.frequency) AS \"assigned reads\" 
                           FROM 
                             (SELECT  
                               reads.sample_id,  
                               reads.rep_id, 
                               reads.read_id AS \"read ID\",
                               reads.frequency, 
                               spec.taxon AS species, 
                               spec_assign.log_prob_species 
                              FROM  
                               reads 
                              LEFT JOIN 
                               spec_assign ON reads.read_id=spec_assign.read_id 
                              LEFT JOIN 
                               taxonomy spec ON species_id = spec.tax_id 
                              GROUP BY
                               reads.sample_id,  
                               reads.rep_id,
                               reads.read_id, 
                               reads.frequency, 
                               species, 
                               log_prob_species) 
                              AS 
                               tab1
                            WHERE 
                              tab1.species IS NOT NULL
                            GROUP BY 
                              tab1.sample_id, 
                              tab1.rep_id, 
                              tab1.species)
                            AS
                              tab2 
                          WHERE 
                           rep_id LIKE 'rep4%' AND species IS NOT NULL) 
                         AS rep4
                         ON tab2.sample_id=rep4.sample_id AND tab2.species=rep4.species
                         
                         LEFT JOIN
                         (SELECT 
                          DISTINCt
                           sample_id, 
                           species, 
                           \"mean probability\" AS \"prob 5\", 
                           \"assigned reads\" AS \"reads 5\" 
                          FROM 
                           (SELECT 
                             tab1.sample_id, 
                             tab1.rep_id, 
                             tab1.species, 
                             to_char((sum(exp(tab1.log_prob_species) * tab1.frequency)/sum(tab1.frequency)),
                             '9D99') AS \"mean probability\",
                             sum(tab1.frequency) AS \"assigned reads\" 
                           FROM 
                             (SELECT  
                               reads.sample_id,  
                               reads.rep_id, 
                               reads.read_id AS \"read ID\",
                               reads.frequency, 
                               spec.taxon AS species, 
                               spec_assign.log_prob_species 
                              FROM  
                               reads 
                              LEFT JOIN 
                               spec_assign ON reads.read_id=spec_assign.read_id 
                              LEFT JOIN 
                               taxonomy spec ON species_id = spec.tax_id 
                              GROUP BY
                               reads.sample_id,  
                               reads.rep_id,
                               reads.read_id, 
                               reads.frequency, 
                               species, 
                               log_prob_species) 
                              AS 
                               tab1
                            WHERE 
                              tab1.species IS NOT NULL
                            GROUP BY 
                              tab1.sample_id, 
                              tab1.rep_id, 
                              tab1.species)
                            AS
                              tab2 
                          WHERE 
                           rep_id LIKE 'rep5%' AND species IS NOT NULL)
                         AS rep5
                         ON tab2.sample_id=rep5.sample_id AND tab2.species=rep5.species
                      
                         LEFT JOIN
                         (SELECT 
                          DISTINCt
                           sample_id, 
                           species, 
                           \"mean probability\" AS \"prob 6\", 
                           \"assigned reads\" AS \"reads 6\" 
                          FROM 
                           (SELECT 
                             tab1.sample_id, 
                             tab1.rep_id, 
                             tab1.species, 
                             to_char((sum(exp(tab1.log_prob_species) * tab1.frequency)/sum(tab1.frequency)),
                             '9D99') AS \"mean probability\",
                             sum(tab1.frequency) AS \"assigned reads\" 
                           FROM 
                             (SELECT  
                               reads.sample_id,  
                               reads.rep_id, 
                               reads.read_id AS \"read ID\",
                               reads.frequency, 
                               spec.taxon AS species, 
                               spec_assign.log_prob_species 
                              FROM  
                               reads 
                              LEFT JOIN 
                               spec_assign ON reads.read_id=spec_assign.read_id 
                              LEFT JOIN 
                               taxonomy spec ON species_id = spec.tax_id 
                              GROUP BY
                               reads.sample_id,  
                               reads.rep_id,
                               reads.read_id, 
                               reads.frequency, 
                               species, 
                               log_prob_species) 
                              AS 
                               tab1
                            WHERE 
                              tab1.species IS NOT NULL
                            GROUP BY 
                              tab1.sample_id, 
                              tab1.rep_id, 
                              tab1.species)
                            AS
                              tab2 
                          WHERE 
                           rep_id LIKE 'rep6%' AND species IS NOT NULL) 
                         AS rep6
                         ON tab2.sample_id=rep6.sample_id AND tab2.species=rep6.species
                         WHERE
                          tab2.species IS NOT NULL
                         ORDER BY
                          tab2.sample_id,
                          tab2.species,
                          \"reads 1\" DESC
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
  
  rm("q_1", "q_tab", "db")
  gc()
}