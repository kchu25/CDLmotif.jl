module CDLmotif

using CUDA, FFTW, LinearAlgebra, Random, 
      DoubleFloats, JLD2, DataStructures, 
      HypothesisTests, Distributions, CSV, 
      DataFrames, ArgParse, Mustache, StatsBase,
      Gadfly, Cairo



import Base.view

export test_mixture,
       test_single,
       test_gapped,
       test_kg,
       find_motifs_fasta,
       test_one_dataset,
       test_mixture_triple,
       find_motifs_fasta_folder
    
include("SETUP/TYPE_setup.jl")
include("SETUP/GPU_setup.jl")
include("SETUP/SIM_setup.jl")
include("SETUP/LOAD_setup.jl")
include("SETUP/SEARCH_setup.jl")
include("SETUP/WRAPPER_setup.jl")
include("SETUP/CDL_setup.jl")
include("LOAD/read.jl")
include("LOAD/save.jl")
include("LOAD/to_dummy.jl")
include("LOAD/sample.jl")
include("CDL/csc_kernels.jl")
include("CDL/cdl_kernels.jl")
include("CDL/train.jl")
include("SEARCH/store.jl")
include("SEARCH/extractions.jl")
include("SEARCH/helpers.jl")
include("SEARCH/PWM_Touzet.jl")
include("SEARCH/scan.jl")
include("SEARCH/extend.jl")
include("SEARCH/trim.jl")
include("SEARCH/e_value_filter.jl")
include("SEARCH/allr.jl")
include("WRAPPER/wrap.jl")
include("SAVE/cover_calculations.jl")
include("SAVE/template.jl")
include("SAVE/save_learned_results.jl")
include("simulate_motif_tests/simulation_wrap.jl")

"""
    find_motifs_fasta(fasta_path, output_folder)
Do a motif discovery to the DNA sequences in the fasta file specified by the fastapath.

Input:
* `fasta_path`:  A string that's the input fasta file's absolute filepath.
* `output_folder`: A string that's the output folder's absolute filepath; 
                   all the motif discovery results will be stored here.
Output:
    The motif discovery results stored in the folder specified in `output_folder`
"""
function find_motifs_fasta(fasta_path::String, output_folder::String)
    data = FASTA_DNA{int_t, dat_t}(fasta_path, output_folder);
    g = try_to_find_motif(data);
    save_result_fasta(g, output_folder);
    # !isnothing(g) && save_jaspar(output_folder*"/jaspar", g)
end

"""
    find_motifs_fasta_folder(fasta_path, output_folder)
Do a motif discovery to the DNA sequences in each fasta files in the input_folder.

Input:
* `input_folder`:  A string that's the input folder's absolute file path; 
                   the input folder contains multiple fasta files. This folder
                   *must* contain only fasta files.
* `output_folder`: A string that's the output folder's absolute file path;
                   this output folder will contain multiple folders. Each 
                   folder will store the motif discovery results that 
                   correspond to an input fasta file.
Output:
   The motif discovery results stored in the folder specified in `output_folder`
"""
function find_motifs_fasta_folder(input_folder::String, output_folder::String)
    for fasta in readdir(input_folder)
        fasta_input = input_folder*"/"*fasta;
        output_loc = output_folder*"/"*fasta;
        find_motifs_fasta(fasta_input, output_loc);
        GC.gc(true);
        GC.gc(false);
        CUDA.reclaim()
    end
end

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--input", "-i"
            help = "input folder"
            required = true
        "-output", "-o"
            help = "output folder"
            required = true
    end
    return parse_args(s)
end

function julia_main()
    parsed_args = parse_commandline()
    find_motifs_fasta_folder(parsed_args["input"], parsed_args["output"]);
end



end
