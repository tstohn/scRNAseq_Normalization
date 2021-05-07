#creates a directory with all normalizations for evry dataset at bin/NORMALIZED_DATASETS
normalization:
	mkdir -p bin/NORMALIZED_DATASETS
	RScript --quiet src/normalization/NormalizationScript.R TMM ALL
	RScript --quiet src/normalization/NormalizationScript.R SUBSAMPLING ALL
	#RScript --quiet src/normalization/NormalizationScript.R SCTRANSFORM ALL
	RScript --quiet src/normalization/NormalizationScript.R CLR_COMPOSITIONS ALL
#creates a directory with all benchmarks for every dataset found in bin/NORMALIZED_DATASETS
benchmark:
	mkdir -p bin/BENCHMARKED_DATASETS
	. ./scRNAseq/bin/activate && \
	python3 src/benchmark/benchmark.py bin/NORMALIZED_DATASETS

make all:
	make normalization
	make benchmark

clean_dataset:
	rm -R bin/NORMALIZED_DATASETS
clean_benchmark:
	rm -R bin/BENCHMARKED_DATASETS

