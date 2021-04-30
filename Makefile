normalization:
	mkdir -p bin/NORMALIZED_DATASETS
	RScript src/normalization/NormalizationScript.R TMM ALL
	RScript src/normalization/NormalizationScript.R SUBSAMPLING ALL
	#RScript src/normalization/NormalizationScript.R SCTRANSFORM ALL
	#RScript src/normalization/NormalizationScript.R CLR_COMPOSITIONS ALL

clean dataset:
	rm -R bin/NORMALIZED_DATASETS
