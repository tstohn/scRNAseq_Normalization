#creates a directory with all normalizations for evry dataset at bin/NORMALIZED_DATASETS
normalization:
	mkdir -p bin/NORMALIZED_DATASETS

	RScript --quiet src/normalization/NormalizationScript.R TMM ALL
	RScript --quiet src/normalization/NormalizationScript.R SUBSAMPLING ALL
	RScript --quiet src/normalization/NormalizationScript.R CLR_COMPOSITIONS ALL

make graph_norm:
	python3 src/methods/GraphNormalization/main.py /Users/t.stohn/Desktop/Normalization/PIPELINE/scRNAseq_Normalization/bin/NORMALIZED_DATASETS/scIDseq-vanEijl-raw-counts/TMM.tsv
	python3 src/methods/GraphNormalization/main.py /Users/t.stohn/Desktop/Normalization/PIPELINE/scRNAseq_Normalization/bin/NORMALIZED_DATASETS/scIDseq-vanEijl-raw-counts_NoEGFRemoved/TMM.tsv
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
update_requirements:
	. ./scRNAseq/bin/activate && \
	pip3 freeze > .requirements_Py.txt
	./.freeze_R.sh
install_requirements:
	. ./scRNAseq/bin/activate && \
	pip3 install -r .requirements_Py.txt
	./.install_R.sh