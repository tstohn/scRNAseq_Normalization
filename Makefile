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


#Runs to re-create paper results










# PERSONAL TEST RUNS
#correlation runs for increasing number of high abudnance proteins (with hig variance)
run_highProteinCorr:
	python3 ./src/simulation/ABseqSimulationBenchmark/main.py ./simulation_inis/1_highProteinCorr.ini --d 5 -k --t 10
#50 proteins with 10 highly inter-correlated -> leads to reduced true corr, and anti-corr between random prot and inter-corr prot
run_highCorr:
	python3 ./src/simulation/ABseqSimulationBenchmark/main.py ./simulation_inis/2_highCorr.ini --d 5 -k --t 10

run_increasedCorr:
	python3 ./src/simulation/ABseqSimulationBenchmark/main.py ./simulation_inis/increasingCorrelations.ini --d 10 -k --t 30

#here we have two distributions (around zero and a high one)
#and we test how well we recover corr as we increase the width of the zero dist

#30 crushes TMM/ CLR
#20 still TMM is best, (also clr shot but maybe an outlier we don t know, but LIB definitely shit)
run_definedCorr:
	python3 ./src/simulation/ABseqSimulationBenchmark/main.py ./simulation_inis/definedCorrelations.ini --d 2 -k --t 5

run_OneDefinedCorr:
	python3 ./src/simulation/ABseqSimulationBenchmark/main.py ./simulation_inis/OneDefinedCorrelations.ini --d 1 -k --t 5
run_RandomCorr:
	python3 ./src/simulation/ABseqSimulationBenchmark/main.py ./simulation_inis/RandomHighCorrelations.ini --d 1 -k --t 5

run_ClusterSpecificRandomCorr:
	python3 ./src/simulation/ABseqSimulationBenchmark/main.py ./simulation_inis/ClusterSpecificRandomCorr.ini --d 1 -k --t 5

run_foldChangeRange:
	python3 ./src/simulation/ABseqSimulationBenchmark/main.py ./simulation_inis/logFoldChangeRange.ini --d 1 -k --t 10 --knn 20 --metric euclidean

run_highPosCorr:
	python3 ./src/simulation/ABseqSimulationBenchmark/main.py ./simulation_inis/highProCorr.ini --d 5 -k --t 10

run_rarePop:
	python3 ./src/simulation/ABseqSimulationBenchmark/main.py ./simulation_inis/rareCellType.ini --d 1 -k --t 10

run_sameProteinChange:
	python3 ./src/simulation/ABseqSimulationBenchmark/main.py ./simulation_inis/sameProteinFoldChange.ini --d 1 -k --t 10

run_linearTraj:
	python3 ./src/simulation/ABseqSimulationBenchmark/main.py ./simulation_inis/linearTrajectory.ini --d 1 -k --t 10
run_twoClust:
	python3 ./src/simulation/ABseqSimulationBenchmark/main.py ./simulation_inis/twoClust.ini --d 1 -k --t 10
	
	
run_cycleTraj:
	python3 ./src/simulation/ABseqSimulationBenchmark/main.py ./simulation_inis/cyclicTrajectory.ini --d 1 -k --t 10


#BioSB Runs:
run_foldChangeRange_biosb:
	python3 ./src/simulation/ABseqSimulationBenchmark/main.py ./simulation_inis/logFoldChangeRange.ini --d 1 -k --t 10 --knn 20 --metric euclidean

run_highPosCorr_biosb:
	python3 ./src/simulation/ABseqSimulationBenchmark/main.py ./simulation_inis/highCorrBioSB.ini --d 1 -k --t 10

run_linearTraj_biosb:
	python3 ./src/simulation/ABseqSimulationBenchmark/main.py ./simulation_inis/linearTrajectoryBioSB.ini --d 1 -k --t 10

run_allProteinsFoldChange:
	python3 ./src/simulation/ABseqSimulationBenchmark/main.py ./simulation_inis/allProtFoldChange.ini --d 1 -k --t 10

run_linearTraj_small:
	python3 ./src/simulation/ABseqSimulationBenchmark/main.py ./simulation_inis/linearTrajectorySmall.ini --d 1 -k --t 1
