
library(tidyverse)
library(httr)

# ----------------------------- Read complex data ------------------------------

# Complex tables (complextab) from Complex Portal 2026-01-14 12:01 update

# Read individual species complexes
complexes_Ce <- read_delim("Ce_6239.tsv", delim = "\t") %>% 
  mutate(Species = "C. elegans")
complexes_Hs <- read_delim("Hs_9606.tsv", delim = "\t") %>% 
  mutate(Species = "H. sapiens")
complexes_Sc <- read_delim("Sc_559292.tsv", delim = "\t") %>% 
  mutate(Species = "S. cerevisiae")

# Join data
complexes <- full_join(complexes_Ce, complexes_Hs) %>% full_join(complexes_Sc)

# Get only complexes with known stoichiometry
complexes_stoich <- complexes %>% filter(`Complex assembly` != "-")

# ------------------------- Retrieve complexes ---------------------------------

# Minimum and maximum number of subunits
N_min = 3 # Trimer
N_max = 6 # Hexamer

# Avoid homooligomers? (only one protein entity)
avoid_homo = TRUE

# How many complexes for each size?
samplex_per_N_size = 1

# AF3 size restrictions
max_single_protein_aa <- 1000  # No single protein chain can exceed this
max_complex_total_aa  <- 4000  # Total complex size cannot exceed this

# ----------------------- Count subunits and filter ---------------------------

complexes_filtered <- complexes_stoich %>%
  mutate(
    parsed = map(
      `Identifiers (and stoichiometry) of molecules in complex`,
      ~ {
        entries <- str_split(.x, "\\|")[[1]]
        tibble(
          id = str_extract(entries, "^[^(]+"),
          n  = as.integer(str_extract(entries, "(?<=\\()\\d+(?=\\))"))
        )
      }
    ),
    
    # Restriction 3: complexes must contain ONLY protein (UniProt) entities
    has_nonprot_entity   = map_lgl(parsed, ~ {
      any(!str_detect(.x$id, "^[A-Z][0-9A-Z]{5}$"))
    }),
    
    # Restriction 4: N_subunits = sum of stoichiometric coefficients
    # (implicitly enforced by summing stoich of ALL entries;
    #  if has_nonprot_entity is FALSE, all entries are proteins, so the sum
    #  equals the true oligomeric state)
    has_unknown_stoich   = map_lgl(parsed, ~ any(.x$n == 0)),
    
    protein_entries      = map(parsed, ~ filter(.x, str_detect(id, "^[A-Z][0-9A-Z]{5}$"))),
    N_subunits           = map_int(protein_entries, ~ sum(.x$n)),
    N_unique_proteins    = map_int(protein_entries, ~ n_distinct(.x$id))
    
  ) %>%
  # Restriction 3: only protein entities
  filter(!has_nonprot_entity) %>%
  # Restriction 4: known stoichiometry for all subunits
  filter(!has_unknown_stoich) %>%
  # Size range
  filter(N_subunits >= N_min, N_subunits <= N_max) %>%
  # Homooligomer filter
  { if (avoid_homo) filter(., N_unique_proteins > 1) else . }

# ----------------------- Fetch sequence from UniProt -------------------------

fetch_uniprot_sequence <- function(uniprot_id) {
  url <- paste0("https://rest.uniprot.org/uniprotkb/", uniprot_id, ".fasta")
  response <- GET(url, timeout(10))
  
  if (status_code(response) != 200) {
    warning(paste("Could not fetch sequence for:", uniprot_id))
    return(NULL)
  }
  
  fasta_text <- content(response, as = "text", encoding = "UTF-8")
  lines <- str_split(fasta_text, "\n")[[1]]
  sequence <- paste(lines[!str_starts(lines, ">") & lines != ""], collapse = "")
  return(sequence)
}

# Cache to avoid fetching the same UniProt ID multiple times
sequence_cache <- list()

fetch_cached <- function(uniprot_id) {
  if (is.null(sequence_cache[[uniprot_id]])) {
    message(paste("  Fetching:", uniprot_id))
    sequence_cache[[uniprot_id]] <<- fetch_uniprot_sequence(uniprot_id)
  }
  sequence_cache[[uniprot_id]]
}

# ----------------------- Size-aware sampling ---------------------------------

check_complex_size <- function(protein_ids_and_stoich) {
  # protein_ids_and_stoich: tibble with columns id (UniProt) and n (stoichiometry)
  
  sequences <- map(protein_ids_and_stoich$id, fetch_cached)
  
  # Restriction 5: all sequences must be fetchable
  if (any(map_lgl(sequences, is.null))) {
    failed <- protein_ids_and_stoich$id[map_lgl(sequences, is.null)]
    message("  REJECTED: could not fetch sequence(s) for: ", paste(failed, collapse = ", "))
    return(list(pass = FALSE, sequences = NULL))
  }
  
  seq_lengths <- map_int(sequences, nchar)
  
  # Restriction 2: no single protein > max_single_protein_aa
  if (any(seq_lengths > max_single_protein_aa)) {
    failed <- protein_ids_and_stoich$id[seq_lengths > max_single_protein_aa]
    message("  REJECTED: protein(s) exceed ", max_single_protein_aa, " aa: ",
            paste(sprintf("%s (%d aa)", failed, seq_lengths[seq_lengths > max_single_protein_aa]),
                  collapse = ", "))
    return(list(pass = FALSE, sequences = NULL))
  }
  
  # Restriction 1: total complex size (each sequence length * stoichiometry) 
  total_aa <- sum(seq_lengths * protein_ids_and_stoich$n)
  if (total_aa > max_complex_total_aa) {
    message("  REJECTED: total complex size ", total_aa, " aa exceeds ",
            max_complex_total_aa, " aa")
    return(list(pass = FALSE, sequences = NULL))
  }
  
  message("  ACCEPTED: ", nrow(protein_ids_and_stoich), " unique proteins | ",
          paste(sprintf("%s=%d aa x%d", protein_ids_and_stoich$id, seq_lengths,
                        protein_ids_and_stoich$n), collapse = ", "),
          " | total = ", total_aa, " aa")
  
  list(pass = TRUE, sequences = setNames(sequences, protein_ids_and_stoich$id))
}

# ----------------------- Sampling loop ---------------------------------------

set.seed(42)

groups <- complexes_filtered %>%
  group_by(Species, N_subunits) %>%
  group_keys()

complexes_sampled <- list()
valid_seqs        <- list()

for (g in seq_len(nrow(groups))) {
  
  sp <- groups$Species[g]
  ns <- groups$N_subunits[g]
  
  # Shuffle the entire pool for this species x N_subunits group
  pool <- complexes_filtered %>%
    filter(Species == sp, N_subunits == ns) %>%
    slice_sample(n = nrow(.))
  
  message(sprintf("\n[%s | N=%d] Pool size: %d complexes", sp, ns, nrow(pool)))
  
  found      <- 0
  valid_rows <- list()
  
  for (i in seq_len(nrow(pool))) {
    if (found >= samplex_per_N_size) break
    
    candidate  <- pool[i, ]
    complex_ac <- candidate$`#Complex ac`
    message(sprintf("  Trying %s (%s)...", complex_ac, candidate$`Recommended name`))
    
    prot_entries <- candidate$protein_entries[[1]]
    size_check   <- check_complex_size(prot_entries)
    
    if (size_check$pass) {
      found <- found + 1
      valid_rows[[found]]          <- candidate
      valid_seqs[[complex_ac]]     <- size_check$sequences
    }
  }
  
  if (found < samplex_per_N_size) {
    warning(sprintf(
      "Only %d/%d valid complexes found for %s N=%d after exhausting pool.",
      found, samplex_per_N_size, sp, ns
    ))
  }
  
  if (length(valid_rows) > 0) {
    complexes_sampled <- c(complexes_sampled, valid_rows)
  }
}

complexes_sampled <- bind_rows(complexes_sampled)

message("\n--- Final sampled complexes ---")
complexes_sampled %>%
  select(`#Complex ac`, `Recommended name`, Species, N_subunits,
         `Identifiers (and stoichiometry) of molecules in complex`) %>%
  print()

# NOTE: There is only one 6-mer C. elegans complex, but it doesn't satisfy the
#       limits of 4000 and 1000 residues as maximum. So, it wasn't included:
# complexes %>% filter(`Complex assembly` == "Heterohexamer") %>% filter(Species == "C. elegans")

# We added an extra block at the end to fill-in these types of gaps

# ----------------------- Create results directory ----------------------------

results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE)

write_tsv(complexes_sampled, file.path(results_dir, "complexes_sampled.tsv"))

# ----------------------- Write per-complex FASTA files -----------------------

walk(seq_len(nrow(complexes_sampled)), function(i) {
  
  row         <- complexes_sampled[i, ]
  complex_ac  <- row$`#Complex ac`
  complex_dir <- file.path(results_dir, complex_ac)
  dir.create(complex_dir, showWarnings = FALSE)
  
  prot_entries <- row$protein_entries[[1]]
  sequences    <- valid_seqs[[complex_ac]]
  
  fasta_lines <- c()
  
  for (j in seq_len(nrow(prot_entries))) {
    uniprot_id  <- prot_entries$id[j]
    stoich      <- prot_entries$n[j]
    entry_str   <- paste0(uniprot_id, "(", stoich, ")")
    header      <- paste0(">", uniprot_id, "_", stoich, "|", entry_str)
    seq         <- sequences[[uniprot_id]]
    fasta_lines <- c(fasta_lines, header, seq)
  }
  
  writeLines(fasta_lines, file.path(complex_dir, "database.fasta"))
  message("Written: ", complex_ac)
})

message("\nAll done! Results saved to: ", results_dir)


################################################################################
################### Also sample D. melanogaster & E. coli ######################
################################################################################

complexes_Dm <- read_delim("Dm_7227.tsv", delim = "\t") %>% 
  mutate(Species = "D. melanogaster")
complexes_Ec <- read_delim("Ec_83333.tsv", delim = "\t") %>% 
  mutate(Species = "E. coli")

# ----------------------- Reusable per-species pipeline ------------------------

filter_complexes <- function(complexes_df) {
  complexes_df %>%
    filter(`Complex assembly` != "-") %>%
    mutate(
      parsed = map(
        `Identifiers (and stoichiometry) of molecules in complex`,
        ~ {
          entries <- str_split(.x, "\\|")[[1]]
          tibble(
            id = str_extract(entries, "^[^(]+"),
            n  = as.integer(str_extract(entries, "(?<=\\()\\d+(?=\\))"))
          )
        }
      ),
      has_nonprot_entity   = map_lgl(parsed, ~ any(!str_detect(.x$id, "^[A-Z][0-9A-Z]{5}$"))),
      has_unknown_stoich   = map_lgl(parsed, ~ any(.x$n == 0)),
      protein_entries      = map(parsed, ~ filter(.x, str_detect(id, "^[A-Z][0-9A-Z]{5}$"))),
      N_subunits           = map_int(protein_entries, ~ sum(.x$n)),
      N_unique_proteins    = map_int(protein_entries, ~ n_distinct(.x$id))
    ) %>%
    filter(!has_nonprot_entity) %>%
    filter(!has_unknown_stoich) %>%
    filter(N_subunits >= N_min, N_subunits <= N_max) %>%
    { if (avoid_homo) filter(., N_unique_proteins > 1) else . }
}

sample_species_complexes <- function(complexes_filtered_df, valid_seqs_env) {
  # valid_seqs_env: environment holding the global valid_seqs list, so we can
  # update it by reference across calls
  
  groups <- complexes_filtered_df %>%
    group_by(Species, N_subunits) %>%
    group_keys()
  
  complexes_sampled_sp <- list()
  
  for (g in seq_len(nrow(groups))) {
    
    sp <- groups$Species[g]
    ns <- groups$N_subunits[g]
    
    pool <- complexes_filtered_df %>%
      filter(Species == sp, N_subunits == ns) %>%
      slice_sample(n = nrow(.))
    
    message(sprintf("\n[%s | N=%d] Pool size: %d complexes", sp, ns, nrow(pool)))
    
    found      <- 0
    valid_rows <- list()
    
    for (i in seq_len(nrow(pool))) {
      if (found >= samplex_per_N_size) break
      
      candidate  <- pool[i, ]
      complex_ac <- candidate$`#Complex ac`
      message(sprintf("  Trying %s (%s)...", complex_ac, candidate$`Recommended name`))
      
      prot_entries <- candidate$protein_entries[[1]]
      size_check   <- check_complex_size(prot_entries)
      
      if (size_check$pass) {
        found <- found + 1
        valid_rows[[found]] <- candidate
        assign("valid_seqs", 
               modifyList(get("valid_seqs", envir = valid_seqs_env), 
                          setNames(list(size_check$sequences), complex_ac)),
               envir = valid_seqs_env)
      }
    }
    
    if (found < samplex_per_N_size) {
      warning(sprintf(
        "Only %d/%d valid complexes found for %s N=%d after exhausting pool.",
        found, samplex_per_N_size, sp, ns
      ))
    }
    
    if (length(valid_rows) > 0) {
      complexes_sampled_sp <- c(complexes_sampled_sp, valid_rows)
    }
  }
  
  bind_rows(complexes_sampled_sp)
}

write_fasta_for_sampled <- function(sampled_df, valid_seqs_list, results_dir) {
  walk(seq_len(nrow(sampled_df)), function(i) {
    
    row         <- sampled_df[i, ]
    complex_ac  <- row$`#Complex ac`
    complex_dir <- file.path(results_dir, complex_ac)
    dir.create(complex_dir, showWarnings = FALSE)
    
    prot_entries <- row$protein_entries[[1]]
    sequences    <- valid_seqs_list[[complex_ac]]
    
    fasta_lines <- c()
    
    for (j in seq_len(nrow(prot_entries))) {
      uniprot_id  <- prot_entries$id[j]
      stoich      <- prot_entries$n[j]
      entry_str   <- paste0(uniprot_id, "(", stoich, ")")
      header      <- paste0(">", uniprot_id, "_", stoich, "|", entry_str)
      seq         <- sequences[[uniprot_id]]
      fasta_lines <- c(fasta_lines, header, seq)
    }
    
    writeLines(fasta_lines, file.path(complex_dir, "database.fasta"))
    message("Written: ", complex_ac)
  })
}

# ----------------------- Run pipeline for Dm and Ec ---------------------------

complexes_Dm_filtered <- filter_complexes(complexes_Dm)
complexes_Ec_filtered <- filter_complexes(complexes_Ec)

complexes_sampled_Dm <- sample_species_complexes(complexes_Dm_filtered, environment())
complexes_sampled_Ec <- sample_species_complexes(complexes_Ec_filtered, environment())

message("\n--- Sampled Drosophila complexes ---")
complexes_sampled_Dm %>%
  select(`#Complex ac`, `Recommended name`, Species, N_subunits,
         `Identifiers (and stoichiometry) of molecules in complex`) %>%
  print()

message("\n--- Sampled E. coli complexes ---")
complexes_sampled_Ec %>%
  select(`#Complex ac`, `Recommended name`, Species, N_subunits,
         `Identifiers (and stoichiometry) of molecules in complex`) %>%
  print()

# ----------------------- Append to main dataframe and re-save ----------------

complexes_sampled <- bind_rows(complexes_sampled, complexes_sampled_Dm, complexes_sampled_Ec)

write_tsv(complexes_sampled, file.path(results_dir, "complexes_sampled.tsv"))

# ----------------------- Write FASTA files for the new complexes only --------

write_fasta_for_sampled(complexes_sampled_Dm, valid_seqs, results_dir)
write_fasta_for_sampled(complexes_sampled_Ec, valid_seqs, results_dir)

message("\nAll done! Drosophila and E. coli results appended and saved to: ", results_dir)


# ============================================================================
# ------------------------- FILL PER-SPECIES GAPS ---------------------------
# ============================================================================
# Some species are missing complexes for specific N sizes because no
# candidate in that size class passes the AF3 length restrictions or there is no
# candidate for that specific size. This block tops up each species to the same
# total target count by randomly sampling additional complexes from ANY N within
# [N_min, N_max] that hasn't already been picked, but rejecting any
# candidate that shares a UniProt ID with a complex already selected for that
# species, to avoid near-duplicate subcomplexes.


# Combined pool across all species' filtered candidates — needed so Stage 2
# can search any N size within a species, not just the one that was missing.
all_filtered <- bind_rows(
  complexes_filtered,       # Ce, Hs, Sc
  complexes_Dm_filtered,
  complexes_Ec_filtered
)

target_n_per_species <- samplex_per_N_size * (N_max - N_min + 1)

current_counts_2 <- complexes_sampled %>% count(Species, name = "have")

species_deficits <- current_counts_2 %>%
  mutate(need = target_n_per_species - have) %>%
  filter(need > 0)

get_used_proteins <- function(species_name) {
  complexes_sampled %>%
    filter(Species == species_name) %>%
    pull(protein_entries) %>%
    map(~ .x$id) %>%
    unlist() %>%
    unique()
}

fill_cross_size <- function(species_name, need, valid_seqs_env) {
  
  message(sprintf("\n[%s] Stage 2: still need %d more (protein-disjoint, any N)...",
                  species_name, need))
  
  used_proteins <- get_used_proteins(species_name)
  already_ac    <- complexes_sampled %>% filter(Species == species_name) %>% pull(`#Complex ac`)
  
  pool <- all_filtered %>%
    filter(Species == species_name, !(`#Complex ac` %in% already_ac)) %>%
    mutate(overlaps = map_lgl(protein_entries, ~ any(.x$id %in% used_proteins))) %>%
    filter(!overlaps) %>%
    slice_sample(n = nrow(.))
  
  found    <- 0
  new_rows <- list()
  
  for (i in seq_len(nrow(pool))) {
    if (found >= need) break
    
    candidate  <- pool[i, ]
    complex_ac <- candidate$`#Complex ac`
    message(sprintf("  Trying %s (%s, N=%d)...",
                    complex_ac, candidate$`Recommended name`, candidate$N_subunits))
    
    size_check <- check_complex_size(candidate$protein_entries[[1]])
    
    if (size_check$pass) {
      found <- found + 1
      new_rows[[found]] <- candidate
      used_proteins <- c(used_proteins, candidate$protein_entries[[1]]$id)  # keep growing exclusion set
      assign("valid_seqs",
             modifyList(get("valid_seqs", envir = valid_seqs_env),
                        setNames(list(size_check$sequences), complex_ac)),
             envir = valid_seqs_env)
    }
  }
  
  if (found < need) {
    warning(sprintf("[%s] Stage 2 could only fill %d/%d — no protein-disjoint candidates left.",
                    species_name, found, need))
  }
  
  if (length(new_rows) > 0) bind_rows(new_rows) else NULL
}

stage2_rows <- list()

for (r in seq_len(nrow(species_deficits))) {
  sp <- species_deficits$Species[r]
  nd <- species_deficits$need[r]
  
  result <- fill_cross_size(sp, nd, environment())
  if (!is.null(result)) stage2_rows[[length(stage2_rows) + 1]] <- result
}

stage2_rows <- if (length(stage2_rows) > 0) bind_rows(stage2_rows) else tibble()

if (nrow(stage2_rows) > 0) {
  complexes_sampled <- bind_rows(complexes_sampled, stage2_rows)
  
  message("\n--- Stage 2 cross-size replacements added ---")
  stage2_rows %>%
    select(`#Complex ac`, `Recommended name`, Species, N_subunits,
           `Identifiers (and stoichiometry) of molecules in complex`) %>%
    print()
  
  write_tsv(complexes_sampled, file.path(results_dir, "complexes_sampled.tsv"))
  write_fasta_for_sampled(stage2_rows, valid_seqs, results_dir)
} else {
  message("\nNo Stage 2 fallback was necessary.")
}

message("\nFinal per-species totals:")
complexes_sampled %>% count(Species) %>% print(n = Inf)

message("\nFinal per-species, per-size breakdown:")
complexes_sampled %>% count(Species, N_subunits) %>% print(n = Inf)