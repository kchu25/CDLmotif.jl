# setup file that preceeds all the setup files
const dat_t = Float32;          # data_matrix_type
const cdl_c = ComplexF32;       # cdl_complex_type
const int_t = Int32;            # integer type
promote_i(x...) = Int.(x);


######### constants for simulated data #########
# concentration parameter for Dirichlet that generate Categorical Distributions
const motif_alpha_pc = 0.1;             
const DNA_char_cap = ['A','C','G','T'];
const DNA_char_nocap = ['a','c','g','t'];
const mixture_weights_alpha = 50;               # α for generating mixture weight using Dirichlet
                                                # assume for now mixture weights are roughly the same            
                                                # Note that modes are vector of unit ranges
const percent_complement = .5;                  # probability that a motif is represented in its reverse complement direction
const DNA_complement = Dict('A'=>'T','C'=>'G','G'=>'C','T'=>'A');                     

const background = Categorical(0.25*ones(4));   # simple i.i.d flat background for now
const bg_array = [.25, .25, .25, .25];
const sample_near_center = 0;                   #= increment and decrement the avail. position 
                                                    from the start and the end, resp. =#

const v_f64 = Vector{Float64};
const cat_d = Categorical{Float64, v_f64};
const u_conc = Union{Float64, v_f64, Nothing};

# tests
const num_seq = 800;
const num_expr_each = 1;

# dummy reference # might change the type later for more flexibility
const dummy = Dict('A'=>Array{Float32}([1, 0, 0, 0]), 
                   'C'=>Array{Float32}([0, 1, 0, 0]),
                   'G'=>Array{Float32}([0, 0, 1, 0]), 
                   'T'=>Array{Float32}([0, 0, 0, 1]));

################################################

######## fasta load settings: ##################
const max_num_read_fasta = 10000;
const nmbr_regex = r"\-?[0-9]+\.?[0-9]*"; 
################################################

################ code tests ####################
const target_folder = "tmp";
const fasta_loc = "src/public_test_data/MA1374.1.sites";
const jaspar_mat_loc = "src/public_test_data/MA1374.1.jaspar";
################################################

############# for PWM_Touzet.jl ################
const _granularity_ = 1e-1; # initial granularity for score2pvalue and pval2score
const _k_ = 100; # decreasing factor for finer granularity in each iteration
const _bg_ = [.25,.25,.25,.25]; # default background
################################################