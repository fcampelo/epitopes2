FASTA_PATH  = paste0(save_folder, "/proteins_masked.fa")
SAVE_FOLDER = paste0(save_folder, "/out")

cmdline <- paste0("python3", " ",
                  EXTRACT_SCRIPT_PATH, " ",
                  MODEL_SPEC, " ",
                  FASTA_PATH, " ",
                  SAVE_FOLDER, " ",
                  MODEL_OPTS)

system(cmdline)

# Alternatively, run the line below (adapt as needed)
# on the terminal
# python3 extract.py esm1b_t33_650M_UR50S output/esm1b_features/proteins_masked.fa output/esm1b_features/out --include per_tok --nogpu --repr_layers 33
# python3 extract.py esm1b_t33_650M_UR50S output/esm1b_features/proteins_masked.fa output/esm1b_features/out --include per_tok --nogpu --repr_layers 33
