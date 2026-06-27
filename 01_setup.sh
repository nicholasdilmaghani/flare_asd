repo_url=https://github.com/drewmard/FLARE.git
repo_dir=FLARE
if [[ ! -d "$repo_dir/.git" ]]; then
    git clone --depth 1 "$repo_url" "$repo_dir"
fi

synapse get -r syn64717038
synapse get -r syn64713923

# 1. Replace the two hardcoded snp_identifier lines with a comment
sed -i '/^snp_identifier = "snp_id"$/{N; s/snp_identifier = "snp_id"\nsnp_identifier = "variant_id"/# snp_identifier auto-detected from data columns/}' FLARE/scripts/FLARE_Preprocess.R
# 2. Insert auto-detection right after df is loaded
sed -i '/df = initial_data_load(opt$input_file)/a snp_identifier = if ("snp_id" %in% colnames(df)) "snp_id" else "variant_id"' FLARE/scripts/FLARE_Preprocess.R



# 1. Replace the hardcoded snp_identifier line with a comment
sed -i 's/^snp_identifier = "variant_id"$/# snp_identifier auto-detected from data columns/' FLARE/scripts/FLARE_Predict.R

# 2. Insert detection right after df is loaded inside FLARE_Predict
sed -i '/df = fread(f\.input,data\.table = F,stringsAsFactors = F)/a\  snp_identifier = if ("snp_id" %in% colnames(df)) "snp_id" else "variant_id"' FLARE/scripts/FLARE_Predict.R

sed -i 's/c("chr", snp_identifier)/c("chr", snp_identifier, "phylop")/' ./FLARE/scripts/FLARE_Predict.R

sed -i '/x <- as.matrix(df\[ind_chr_include, cols_exclude\])/a\    if (nrow(x) == 0) next' ./FLARE/scripts/FLARE_Predict.R
