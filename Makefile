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
	python3 src/methods/GraphNormalization/main.py ./bin/NORMALIZED_DATASETS/scIDseq-vanEijl-raw-counts/TMM.tsv ./bin
	python3 src/methods/GraphNormalization/main.py ./bin/NORMALIZED_DATASETS/scIDseq-vanEijl-raw-counts_REMOVED_No_EGF/TMM.tsv ./bin
#creates a directory with all benchmarks for every dataset found in bin/NORMALIZED_DATASETS
benchmark:
	mkdir -p bin/BENCHMARKED_DATASETS
	. ./scRNAseq/bin/activate && \
	python3 src/benchmark/main.py bin/NORMALIZED_DATASETS

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
	python3 ./src/simulation/ABseqSimulationBenchmark/main.py ./test/data/SimulationFolder/test_simulation.ini --d 2 -k

test_simulation_benchmark:
	python3 ./src/simulation/ABseqSimulationBenchmark/main.py --t 2 ./test/data/SimulationFolder/

test_benchmarking_normalized_files:
	python3 ./src/benchmark/main.py --groundtruth --iniFile ./test/data/BenchmarkFolder/test_simulation.ini --t -1 ./test/data/BenchmarkFolder/benchmarkFiles

# RUNS for analysis of paper
#correlation runs for increasing number of high abudnance proteins (with hig variance)
run_1:
	python3 ./src/simulation/ABseqSimulationBenchmark/main.py ../Simulations/corr_simulation_1.ini --d 2 -k --t 10
run_2:
	python3 ./src/simulation/ABseqSimulationBenchmark/main.py ../Simulations/trueCorr_simulation_2.ini --d 2 -k --t 10
run_3:
	python3 ./src/simulation/ABseqSimulationBenchmark/main.py ../Simulations/strengthCorr_simulation_3.ini --d 2 -k --t 10
run_4:
	python3 ./src/simulation/ABseqSimulationBenchmark/main.py ../Simulations/combination_4.ini --d 2 -k --t 10

run_cluster_1:
	python3 ./src/simulation/ABseqSimulationBenchmark/main.py ../Simulations/increaseNumProteins_1.ini --d 2 -k --t 10
