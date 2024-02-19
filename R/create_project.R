# create new project


new_project = function(
    db_user = NA
    ){
 
  db_connect(db_user = db_user) 
 
  keywords = DBI::dbReadTable(db_con, DBI::Id(schema = "projects", table = "keywords"))
  projects = DBI::dbReadTable(db_con, DBI::Id(schema = "projects", table = "proj_info"))
  
  # project details
  proj_name = readline("Project name: ")
  proj_year = readline("Project start year: ")

  # check if project name already exists
  while(tolower(proj_name) %in% tolower(projects$proj_name) == T){
    message(paste0("A project with name '", proj_name, "' already exists. Please check or chose a different name.\nYou can use RHelperDB::db_projects_info(proj_name = \"", proj_name, "\") to see more information on this project.\nPlease chose if you want to change project name an continue or if you want to stop."))
    choice_1 = select.list(c("Chose new project name", "Exit"), title = "Please chose by typing '1' or '2':")
    if(choice_1 == "Chose new project name"){
      proj_name <<- readline("Project name: ")
    } else {
      stop("You decided to stop the process.")
    }
  }
  # write new projec to DB
  DBI::dbWriteTable(db_con, DBI::Id(schema = "projects", table="proj_info"), data.frame(cbind(proj_name, proj_year)), append = T)
  # renew projects to get the new project ID
  projects = DBI::dbReadTable(db_con, DBI::Id(schema = "projects", table = "proj_info"))
  
  
  proj_keywords = select.list(c(sort(keywords$keyword), "other"), graphics = T, multiple = T)
  
  if(any(proj_keywords == "other")){
    other_keys = readline("Enter new keywords seperated by ';':")
    other_keys = gsub("; ", ";", other_keys)
    other_keys <<- strsplit(other_keys, ";")[[1]]
    if(length(other_keys)!=0){
      while (any(grepl("\\,", other_keys) == T)) {
        message("Keywords should not contain ','. Please seperate them by ';'!")
        other_keys = readline("Enter new keywords seperated by ';':")
        other_keys = gsub("; ", ";", other_keys)
        other_keys <<- strsplit(other_keys, ";")[[1]]
        proj_keywords <<- proj_keywords[proj_keywords != "other"]
      }
    }
    
    if(length(other_keys) == 0){
      while(length(other_keys) == 0){
        message("When the option 'other' was chosen you must define new keywords or deselect 'other'.")
        rechoice = select.list(c("Deselect other", "Define new keywords"), title = "Please chose by typing '1' or '2':")
      
        if(rechoice == "Define new keywords"){
          other_keys = readline("Enter new keywords seperated by ';':")
          other_keys = gsub("; ", ";", other_keys)
          other_keys <<- strsplit(other_keys, ";")[[1]]
          proj_keywords <<- proj_keywords[proj_keywords != "other"]
        
          while (any(grepl("\\,", other_keys) == T)) {
            message("Keywords should not contain ','. Please seperate them by ';'!")
            other_keys = readline("Enter new keywords seperated by ';':")
            other_keys = gsub("; ", ";", other_keys)
            other_keys <<- strsplit(other_keys, ";")[[1]]
            proj_keywords <<- proj_keywords[proj_keywords != "other"]
            }
          } else {
          proj_keywords <<- unique(proj_keywords[proj_keywords != "other"])
          other_keys = NA
        }
      } 
    } 
  }
  if(is.na(other_keys) == F){
    DBI::dbWriteTable(db_con, DBI::Id(schema="projects", table="keywords"), data.frame(keyword = other_keys), append = T)
    proj_keywords <<- c(proj_keywords[proj_keywords != "other"], other_keys)
      } else{
    proj_keywords <<- proj_keywords[proj_keywords != "other"]
  }
  DBI::dbWriteTable(db_con, DBI::Id(schema="projects", table="proj_keywords"), data.frame(proj_id = projects$proj_id[projects$proj_name == proj_name], keyword = proj_keywords), append = T)
 
}
