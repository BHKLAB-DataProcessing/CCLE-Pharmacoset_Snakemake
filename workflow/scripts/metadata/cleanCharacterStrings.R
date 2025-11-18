cleanCharacterStrings <- function(name) {
  # # make sure name is a string
  name <- as.character(name)

  # if there is a colon like in "Cisplatin: 1 mg/mL (1.5 mM); 5 mM in DMSO"
  # remove everything after the colon
  name <- gsub(":.*", "", name)

  # remove ,  ;  -  +  *  $  %  #  ^  _  as well as any spaces
  name <- gsub("[\\,\\;\\-\\+\\*\\$\\%\\#\\^\\_\\s]", "", name, perl = TRUE)

  # remove hyphen
  name <- gsub("-", "", name)

  # # remove substring of round brackets and contents
  name <- gsub("\\s*\\(.*\\)", "", name)

  # # remove substring of square brackets and contents
  name <- gsub("\\s*\\[.*\\]", "", name)

  # # remove substring of curly brackets and contents
  name <- gsub("\\s*\\{.*\\}", "", name)

  # # convert entire string to uppercase
  name <- toupper(name)

  # dealing with unicode characters
  name <- gsub(
    "Unicode",
    "",
    iconv(name, "LATIN1", "ASCII", "Unicode"),
    perl = TRUE
  )

  name
}
