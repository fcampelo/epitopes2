process_individual_epitope_B <- function(idx, list_data){

  ep <- list_data$Reference$Epitopes[[idx]]

  # -> ASSUMPTION: Linear B-Cell epitopes appear as
  # "FragmentOfANaturalSequenceMolecule" and contain a field named
  # "LinearSequence"
  Assays <- ep[which(names(ep) == "Assays")]
  if (length(Assays) == 0) return(NULL)

  BCell  <- sapply(Assays, function(x){"BCell" %in% names(x)})
  Assays <- Assays[which(BCell)]
  not_B  <- !any(BCell)
  not_L  <- is.null(ep$EpitopeStructure$FragmentOfANaturalSequenceMolecule$LinearSequence)

  # If it is not a linear B-Cell epitope, return NULL
  if(not_B | not_L) return(NULL)

  # ============= ONLY LINEAR B-CELL EPITOPES CROSS THIS LINE ============= #

  # Extract relevant fields.
  out <- data.table::data.table(
    pubmed_id      = nullcheck(list_data$Reference$Article$PubmedId),
    year           = nullcheck(list_data$Reference$Article$ArticleYear),
    epit_name      = nullcheck(ep$EpitopeName),
    epitope_id     = nullcheck(ep$EpitopeId),
    evid_code      = nullcheck(ep$EpitopeEvidenceCode),
    epit_struc_def = nullcheck(ep$EpitopeStructureDefines),
    sourceOrg_id   = nullcheck(ep$EpitopeStructure$FragmentOfANaturalSequenceMolecule$SourceOrganismId),
    protein_id     = nullcheck(ep$EpitopeStructure$FragmentOfANaturalSequenceMolecule$SourceMolecule$GenBankId),
    epit_seq       = nullcheck(ep$EpitopeStructure$FragmentOfANaturalSequenceMolecule$LinearSequence),
    start_pos      = nullcheck(ep$EpitopeStructure$FragmentOfANaturalSequenceMolecule$StartingPosition),
    end_pos        = nullcheck(ep$EpitopeStructure$FragmentOfANaturalSequenceMolecule$EndingPosition))

  # Double check start and end positions
  if(is.na(out$start_pos)) out$start_pos <- nullcheck(ep$ReferenceStartingPosition)
  if(is.na(out$end_pos))   out$end_pos   <- nullcheck(ep$ReferenceEndingPosition)
  out$start_pos <- as.numeric(out$start_pos)
  out$end_pos   <- as.numeric(out$end_pos)

  c1 <- !is.na(out$epit_seq)
  c2 <- !is.na(out$end_pos - out$start_pos)
  if(c1 && c2){
    # If the epitope length does not agree with its declared positions:
    if(out$end_pos - out$start_pos + 1 != nchar(out$epit_seq)){
      out$start_pos <- NA
      out$end_pos   <- NA
    }
  }

  # Get information from Assays
  host_id    <- character(length(Assays))
  class      <- host_id
  bcell_id   <- host_id
  assay_type <- host_id

  for (i in seq_along(Assays)){
    host_id[i]    <- nullcheck(Assays[[i]]$BCell$Immunization$HostOrganism$OrganismId)
    bcell_id[i]   <- nullcheck(Assays[[i]]$BCell$BCellId)
    class[i]      <- nullcheck(Assays[[i]]$BCell$AssayInformation$QualitativeMeasurement)
    assay_type[i] <- nullcheck(Assays[[i]]$BCell$AssayInformation$AssayTypeId)
  }

  out$n_assays    <- length(Assays)
  out$host_id     <- paste(host_id, collapse = ",")
  out$bcell_id    <- paste(bcell_id, collapse = ",")
  out$assay_type  <- paste(assay_type, collapse = ",")
  out$n_Positive  <- length(grep("Positive", class,ignore.case = TRUE))
  out$n_Negative  <- out$n_assays - out$n_Positive
  class           <- as.numeric(grepl("Positive", class,ignore.case = TRUE))
  out$assay_class <- paste(-1 + 2 * class, collapse = ",")

  return(out)
}
