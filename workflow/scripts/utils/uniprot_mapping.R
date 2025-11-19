suppressPackageStartupMessages({
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required for UniProt mapping.", call. = FALSE)
  }
})

canonical_gene_symbol <- function(symbols) {
  clean <- trimws(symbols)
  clean[clean == ""] <- NA_character_
  clean <- sub("[|;,].*$", "", clean)
  clean <- gsub("[^A-Za-z0-9]", "", clean)
  clean <- toupper(clean)
  clean[clean == ""] <- NA_character_
  clean
}

coalesce_non_empty <- function(...) {
  args <- list(...)
  if (!length(args)) {
    return(NULL)
  }
  out <- args[[1]]
  if (length(args) == 1) {
    return(out)
  }
  for (i in 2:length(args)) {
    current <- args[[i]]
    mask <- is.na(out) | out == ""
    out[mask] <- current[mask]
  }
  out
}

query_uniprot_symbol <- function(symbol,
                                 organism_id = "9606",
                                 sleep_seconds = 0.1) {
  base_url <- "https://rest.uniprot.org/uniprotkb/search"
  params <- list(
    query = sprintf("gene_exact:%s AND organism_id:%s", symbol, organism_id),
    format = "tsv",
    fields = "accession,id,gene_primary,protein_name",
    size = 1
  )
  resp <- httr::GET(base_url, query = params, httr::timeout(30))
  if (httr::status_code(resp) != 200) {
    warning(
      sprintf(
        "UniProt query failed for %s (status %s)",
        symbol,
        httr::status_code(resp)
      ),
      call. = FALSE
    )
    return(
      data.table::data.table(
        gene_symbol = symbol,
        uniprot_accession = NA_character_,
        uniprot_entry = NA_character_,
        uniprot_protein_name = NA_character_
      )
    )
  }
  txt <- httr::content(resp, as = "text", encoding = "UTF-8")
  lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  if (length(lines) < 2) {
    return(
      data.table::data.table(
        gene_symbol = symbol,
        uniprot_accession = NA_character_,
        uniprot_entry = NA_character_,
        uniprot_protein_name = NA_character_
      )
    )
  }
  tbl <- utils::read.delim(
    textConnection(txt),
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!nrow(tbl)) {
    return(
      data.table::data.table(
        gene_symbol = symbol,
        uniprot_accession = NA_character_,
        uniprot_entry = NA_character_,
        uniprot_protein_name = NA_character_
      )
    )
  }
  names(tbl) <- gsub("[^A-Za-z0-9.]+", ".", names(tbl))
  result <- data.table::data.table(
    gene_symbol = symbol,
    uniprot_accession = tbl$Entry,
    uniprot_entry = tbl$Entry.Name,
    uniprot_protein_name = tbl$Protein.names
  )
  if (!is.null(sleep_seconds) && sleep_seconds > 0) {
    Sys.sleep(sleep_seconds)
  }
  result
}

map_uniprot_symbols <- function(symbols,
                                cache_file = "metadata/uniprot_gene_cache.rds",
                                organism_id = "9606") {
  symbols <- canonical_gene_symbol(symbols)
  symbols <- unique(symbols[!is.na(symbols)])
  if (!length(symbols)) {
    return(data.table::data.table(
      gene_symbol = character(),
      uniprot_accession = character(),
      uniprot_entry = character(),
      uniprot_protein_name = character()
    ))
  }
  cache <- if (file.exists(cache_file)) {
    data.table::as.data.table(readRDS(cache_file))
  } else {
    data.table::data.table(
      gene_symbol = character(),
      uniprot_accession = character(),
      uniprot_entry = character(),
      uniprot_protein_name = character()
    )
  }
  missing_symbols <- setdiff(symbols, cache$gene_symbol)
  if (length(missing_symbols)) {
    message(
      sprintf(
        "Querying UniProt for %d gene symbols...",
        length(missing_symbols)
      )
    )
    queried <- data.table::rbindlist(
      lapply(
        missing_symbols,
        function(sym) query_uniprot_symbol(sym, organism_id = organism_id)
      ),
      fill = TRUE,
      use.names = TRUE
    )
    cache <- unique(
      data.table::rbindlist(
        list(cache, queried),
        fill = TRUE,
        use.names = TRUE
      ),
      by = "gene_symbol"
    )
    dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(cache, cache_file)
  }
  cache[gene_symbol %in% symbols]
}
