using CUDA, FFTW, LinearAlgebra, Random, 
      DoubleFloats, JLD2, DataStructures, 
      HypothesisTests, Distributions, CSV, 
      DataFrames, ArgParse, Mustache

      
includet("src/SETUP/TYPE_setup.jl")
includet("src/SETUP/GPU_setup.jl")
includet("src/SETUP/SIM_setup.jl")
includet("src/SETUP/LOAD_setup.jl")
includet("src/SETUP/SEARCH_setup.jl")
includet("src/SETUP/WRAPPER_setup.jl")
includet("src/SETUP/CDL_setup.jl")
includet("src/LOAD/read.jl")
includet("src/LOAD/save.jl")
includet("src/LOAD/to_dummy.jl")
includet("src/LOAD/sample.jl")
includet("src/CDL/csc_kernels.jl")
includet("src/CDL/cdl_kernels.jl")
includet("src/CDL/train.jl")
includet("src/SEARCH/store.jl")
includet("src/SEARCH/extractions.jl")
includet("src/SEARCH/helpers.jl")
includet("src/SEARCH/PWM_Touzet.jl")
includet("src/SEARCH/scan.jl")
includet("src/SEARCH/extend.jl")
includet("src/SEARCH/trim.jl")
includet("src/SEARCH/e_value_filter.jl")
includet("src/SEARCH/allr.jl")
includet("src/WRAPPER/wrap.jl")
includet("src/SAVE/cover_calculations.jl")
includet("src/SAVE/template.jl")
includet("src/simulate_motif_tests/simulation_wrap.jl")

JASPAR_c2h2 = "/home/shane/Desktop/Motif_package_tests/expr_data/JASPAR/c2h2/MA0463.1.sites"
JASPAR_c2h2 = "/home/shane/Desktop/Motif_package_tests/expr_data/JASPAR/c2h2/MA0095.2.sites"
JASPAR_c2h2_out = "/home/shane/Desktop/Motif_package_tests/expr_data/JASPAR_out/c2h2"

data = FASTA_DNA{int_t, dat_t}(JASPAR_c2h2, target_folder);
search_setup = SEARCH_setup{int_t, dat_t}();
cdl_setup = CDL_setup{int_t, dat_t, cdl_c}();
g = good_stuff{int_t,dat_t,cdl_c}(data, cdl_setup, search_setup);
@time find_motif(g)


