#' SET UP ESM FEATURE CALCULATION PARAMETER
#'
#' **Parameter descriptions**:
#' @param prot.df dataframe containing a column with proteins
#' sequences
#' @param py_script_path path to the ESM Python script (NOTE: this is
#' provided as inst/utils/extract.py in this package's folder structure. The
#' script is a slightly modified version of the one available in
#' <https://github.com/facebookresearch/esm/blob/main/scripts/extract.py>).
#' @param id_column name of column in `prot.df` containing the unique protein
#' IDs.
#' @param seqs_column name of column in `prot.df` containing the protein
#' sequences
#' @param save_folder path to folder for saving the output.
#' @param model string with the full model name
#' (see <https://github.com/facebookresearch/esm>)
#' @param model.params string with model options to be used
#' @param chunk_size size of chunk to be used when proteins are
#' longer than 1024 residues (usually as large as possible and < 1024)
#' @param step_size step size to use for processing long proteins. Smaller
#' values are better, but more computationally intensive.
#' @param ncpus positive integer, number of cores to use
#' @param feat_prefix prefix to be added to the feature names.
#'
#' ESM-1b returns errors when trying to calculate features for proteins longer
#' than 1024 residues (see https://github.com/facebookresearch/esm/issues/49).
#' A possible solution (suggested in, e.g., https://github.com/brianhie/evolocity/issues/2)
#' is to break longer sequences in smaller overlapping windows, predict
#' independently and then aggregate in post-processing, averaging where needed.
#' We follow this strategy here.
#'
#' For the feature calculation, any non-standard AA character is replaced by
#' the `<unk>` placeholder.

prot.df = prot.df # <---- data frame with sequences
ncpus = 1
id_column  = "Info_protein_id"
seqs_column = "Info_protein_sequence"
save_folder = "output/esm1b_features"
feat_prefix = "feat_esm1b_"
chunk_size = 1000
step_size = 250
EXTRACT_SCRIPT_PATH = "extract.py"
MODEL_SPEC  = "esm1b_t33_650M_UR50S"
MODEL_OPTS  = "--include per_tok --nogpu --repr_layers 33"
