db_connect = function(db_user = NA,
                      db_name = "EcoDynDB",
                      db_host = "192.168.2.83",
                      db_port = "5432"){
 
  # ask for DB user and password
  if(is.na(db_user)){
    db_user = readline("database user: ")
  }
  
  # db password query
  db_password = getPass::getPass("database password: ")
  
  # establish connection to EcoDynDB
  writeLines(paste0("Connecting to ", db_name, " database.\nTo disconnect from database use 'db_disconnect()'."))
  db_con <<- DBI::dbConnect(RPostgres::Postgres(), 
                        user = db_user,
                        password = db_password,
                        dbname = db_name,        
                        host = db_host,
                        port = db_port)
                        }


db_disconnect = function(db_con = db_con){
  if(exists("db_con")){
    DBI::dbDisconnect(conn = db_con)
    rm("db_con")
    cat("Disconnected from database")
  } else {
    message("No database connection exists.")
  }
}

# create an object 'db_structure' that shows the available schema and tables of the database
db_tables = function(){
  if(exists("db_con") == T){
    db_structure = DBI::dbGetQuery(db_con, "SELECT * FROM information_schema.tables WHERE table_schema NOT IN ('pg_catalog', 'information_schema') AND table_schema NOT LIKE 'pg_toast%'")[, 1:3]
    db_structure = db_structure[order(db_structure$table_schema, db_structure$table_name),]
    row.names(db_structure) = c(1 : nrow(db_structure))
    colnames(db_structure) = c("db_name", "schema_name", "table_name")
    db_structure <<- db_structure
    writeLines("created R dataframe object 'db_structure':\n")
    print(db_structure)
  } else {
    message("No database connection exists.")
  }
}
    
