normalization:
	mkdir -p bin/NORMALIZED_DATASETS
	RScript src/normalization/NormalizationScript.R TMM ALL

clean dataset:
	rm -R bin/NORMALIZED_DATASETS
