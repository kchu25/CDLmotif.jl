#=
A record data structure for each dna sequence in the dataset
    T1: whatever Int type
    T2: smaller Int type to save space; e.g. Int16, Int32    
=#
mutable struct activation_record{T <: Integer,                             
                                 S <: Real}    
    K::T # number of filters to be recorded
    N::T # number of sequences 
    T::T # number of trainings
    C::T # width of the cover area
    uinds::Vector{T}

    #= start_positions: 
    A vector of N matrices; N is the number of data points
        row indicates the filter indices, 
        columns are the starting positions at each dna sequence
    Each matrix is of dimension C x K, here C is the "covered area"
        set by the reduce-search-space procedure by the ssk kernel
    =#
    start_positions::Vector{BitMatrix}
    #=
    Filters' lengths: if there are K filters, then this is a 
        K dimension vector that record each filter's length
    =#    
    filter_lengths::Vector{T}
    #=
    Training number: A K-dimensional vector that records which 
        training (labelled by positive integers) each filter 
        belongs to
    =#
    training_num::Vector{T} # the only redundant part (at least for now for easier debugging)
    #=
    cover_starts record: needed b/c of the reduce-search-space procedure
    =#
    cover_starts::Array{T, 2}
    #= 
    codes: store the strength of the code coefficients
    =#
    codes::CuArray{S,4}

    function activation_record{T,S}(N::Integer, 
                                        C::Integer, 
                                        M::Integer
            ) where {T <: Integer, S <: Real}
        # note that here there's zero columns, so no uncertainty about the entries being initialized
        bitMatrices = [BitMatrix(undef, (C,0)) for _=1:N]; 
        new(0,T(N),0,T(C),
            Vector{T}(), 
            bitMatrices,
            Vector{T}(),
            Vector{T}(),
            Matrix{T}(undef, promote_i(N,0)),
            CuArray{S,4}(undef, promote_i(2,C,M,N))
            );
    end
end

#=
struct for inferring the topology of the activations
    construction:
        K: total number of filters (from all of the trainings)
        max_fil_len: maximum length of the filters
        CminusFsminus1: cover_length minus the smallest filter length minus 1
                        which is C-|fₛ|-1; this is the maximum length that the 
                        filters can be apart from each other in the cover area

Note: 1. I am using "topology" in the loose sense here
      2. T should be Int16 or Int32 since 
                - K: total number of filters is typically less than 100
                - Use Int16 if overlap/non_overlap/next_to count is less than 32,768 
=#
mutable struct activations_topology{T<: Integer}
    K::T                           # number of filters recorded
    overlap::CuArray{T, 3}
    non_overlap::CuArray{T, 3}
    next_to::CuArray{T, 2}
    function activations_topology{T}(
                                  K::Integer, 
                                  max_fil_len::Integer, 
                                  CminusFsminus1::Integer
                                  ) where {T <: Integer}
        new(K, CUDA.fill(T(0), promote_i(K,K,max_fil_len)), 
               CUDA.fill(T(0), promote_i(K,K,CminusFsminus1)),
               CUDA.fill(T(0), promote_i(K,K))
            );
    end
end

mutable struct motifs{T <: Integer, S <: Real}
    pfms::Vector{Matrix{S}}
    pwms::Union{Nothing,Vector{Matrix{S}}}
    thresh::Union{Nothing, Vector{S}}
    lens::Vector{T}
    num_motifs::T
    positions::Union{Nothing, Vector{Dict{T, T}}}
    scores::Union{Nothing, Vector{Dict{T, S}}}
    use_comp::Union{Nothing, Vector{Dict{T, Bool}}}
end

findparam_int(ms::motifs{T,S}) where {T,S} = T;
findparam_real(ms::motifs{T,S}) where {T,S} = S;

