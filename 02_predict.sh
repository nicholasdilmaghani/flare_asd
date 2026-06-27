module load r/4.5.0
./FLARE/scripts/FLARE_Preprocess.R \
  -i  ./asd.all_dataset.K562_bias.annot2.txt.gz \
  -o  ./ASD.FLARE-fb.txt \
  -m  fetal_brain

# Train: lasso, 4-fold CV for lambda, one model per held-out chromosome
mkdir models
./FLARE/scripts/FLARE_Training.R -i ./ASD.FLARE-fb.txt -o ./models/ASD.FLARE-fb

# Predict FLARE-fb scores for the ASD variants
./FLARE/scripts/FLARE_Predict.R \
  -i ./ASD.FLARE-fb.txt \
  -o ./ASD.FLARE-fb.ASD.txt \
  -m ./models/ASD.FLARE-fb

