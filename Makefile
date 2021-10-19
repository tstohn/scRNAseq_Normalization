#creates a directory with all normalizations for evry dataset at bin/NORMALIZED_DATASETS
normalization:
	$(eval DATA := $(if $(DATA),$(DATA),"ALL"))

	mkdir -p bin/NORMALIZED_DATASETS
	mkdir -p bin/FILTERED_DATASETS

	Rscript --quiet src/normalization/NormalizationScript.R TMM $(DATA)
	Rscript --quiet src/normalization/NormalizationScript.R SUBSAMPLING $(DATA)
	Rscript --quiet src/normalization/NormalizationScript.R CLR_COMPOSITIONS $(DATA)
	Rscript --quiet src/normalization/NormalizationScript.R EXPORT_FILTERED_DATA $(DATA)
	python3 src/methods/GraphNormalization/main.py /Users/t.stohn/Desktop/Normalization/PIPELINE/datasets/2_preprocessed/$(DATA)


make graph_norm:
	python3 src/methods/GraphNormalization/main.py /Users/t.stohn/Desktop/Normalization/PIPELINE/scRNAseq_Normalization/bin/NORMALIZED_DATASETS/scIDseq-vanEijl-raw-counts/TMM.tsv
	python3 src/methods/GraphNormalization/main.py /Users/t.stohn/Desktop/Normalization/PIPELINE/scRNAseq_Normalization/bin/NORMALIZED_DATASETS/scIDseq-vanEijl-raw-counts_NoEGFRemoved/TMM.tsv
#creates a directory with all benchmarks for every dataset found in bin/NORMALIZED_DATASETS
benchmark:
	mkdir -p bin/BENCHMARKED_DATASETS
	. ./scRNAseq/bin/activate && \
	python3 src/normBenchmark/benchmark.py bin/NORMALIZED_DATASETS

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


test_simulation:
	python3 ./src/methods/ABseqSimulation/main.py ./examples/simulation.ini

test_simulation_benchmark:
	python3 ./src/methods/ABseqSimulationBenchmark/main.py ../datasets/Simulations/
