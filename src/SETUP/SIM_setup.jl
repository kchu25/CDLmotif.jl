#=
e.g. max([1:5, 2:9]) returns 9
=#
function max_v(v::Vector{UnitRange{Int}})
    _max = -Inf;
    for vec in v
        η = maximum(vec);
        _max = η > _max ? η : _max;
    end
    return _max
end

# product categorical 
struct product_categorical
    pc::Vector{cat_d}      # product categorical 
    len::Int              # length of the product
    dim::Int              # dimension for each categorical in the product
    alpha::u_conc           # concentration parameter; 
             
    function product_categorical(L::Int, 
                                 dim::Int, 
                                 alpha::Union{Real, v_f64})             
        if typeof(alpha) == v_f64 
            length(alpha) != dim && error("alpha must be of length same as the dimension");                            
        end
        D = typeof(alpha) == v_f64 ? Dirichlet(alpha) : Dirichlet(dim, alpha);
        new([Categorical(rand(D)) for _ in 1:L], L, dim, alpha)
    end

    function product_categorical(mat::Matrix{T}) where T <: Real
        L, width = size(mat);
        @assert width == 4 "the input matrix must have width 4"
        new([Categorical(mat[i,:]) for i = 1:L], L, 4, nothing)
    end
end

# view(v::product_categorical) = [round(v.pc[j].p[i], digits=2) for j = 1:v.len, i = 1:v.dim];
Base.rand(v::product_categorical) = [rand(v.pc[i]) for i = 1:v.len];


# motifs definitions
struct single_part_motif
    P::product_categorical # P stands for probability measure
    len::Int64
    function single_part_motif(L::T) where T <: Integer
        new(product_categorical(L, 4, motif_alpha_pc), L)
    end
    function single_part_motif(mat::Matrix{T}) where T <: Real
        new(product_categorical(mat), size(mat,1));
    end
end

# example ----------------------------------
# a = single_part_motif(5);
# view(a)
# ------------------------------------------

# an intermediate type for gap motifs
struct k_parts_motif
    P::Vector{single_part_motif}
    K::Int64
    """
    K: number of parts in the motif
    L_vec: vector of lengths L's for each product categorical
    """
    function k_parts_motif(length_vec::Vector{Int})
        K = length(length_vec);
        new([single_part_motif(length_vec[i]) for i = 1:K], K)
    end
end

#= 
motif with k parts and gaps in between every adjacent parts
Input:
    L_vec:      Length of each of the parts of the motif
    gap_vec:    Length of each of the gap in between parts --
                gap_vecᵢ specifies the gap in between ith and 
                (i+1)th part        
    gap_dist:   Multinomial that generates the alphabets in an i.i.d fashion
=# 

struct gapped_k_parts_motif
    P_motifs::k_parts_motif
    gap_len::Vector{Int}    # the maximal gap length that can occur in between each parts
    P_gap::Categorical{Float64, Vector{Float64}}
    K::Int
    function gapped_k_parts_motif(L_vec::Vector{T}, gap_vec::Vector{T}, gap_dist=bg_array) where T <: Integer       
        K = length(L_vec);
        @assert length(gap_vec) == K-1 "Length of the gap_vec (2nd argument) must be length(L_vec)-1"
        new(k_parts_motif(L_vec), gap_vec, Categorical(gap_dist), K)
    end
end


struct mixture_gapped_k_parts_motifs
    modes::Vector{UnitRange{Int}} # inclusive
    mixture_weights::cat_d
    motif::gapped_k_parts_motif
    num_modes::Int
    function mixture_gapped_k_parts_motifs(
        L_vec::Vector{Int}, 
        gap_vec::Vector{Int}, 
        modes::Vector{UnitRange{Int64}},
        gap_dist=bg_array
    )   
        num_modes = length(modes);
        new(modes, 
            Categorical(rand(Dirichlet(num_modes, mixture_weights_alpha))),
            gapped_k_parts_motif(L_vec, gap_vec, gap_dist),         
            num_modes
        )
    end
end

# a simpler mixture motif than the one above

struct mixture_k_parts_motifs
    modes::Vector{UnitRange{Int}} # inclusive
    mixture_weights::cat_d
    motif::single_part_motif
    num_modes::Int
    function mixture_k_parts_motifs(modes::Vector{UnitRange{Int64}})   
        K = max_v(modes);
        num_modes = length(modes);
        new(modes, 
            Categorical(rand(Dirichlet(num_modes, mixture_weights_alpha))),
            single_part_motif(K),            
            num_modes
        )
    end
end

# motif type
const motif_type = Union{single_part_motif, 
                         gapped_k_parts_motif, 
                         mixture_k_parts_motifs,
                         mixture_gapped_k_parts_motifs
};

# type for a (simulated) single dna background string with a motif in it 
const sim_dna_str_w_motif = NamedTuple{(:str, :motif_where, :mode, :complement), 
                    Tuple{String, UnitRange{Int64}, Int64, Bool}
                    }; # mode = 0 => no motif
const public_dna_str_w_motif = NamedTuple{(:str, :motif_where, :mode), 
                    Tuple{String, UnitRange{Int64}, Int64}
                    };
