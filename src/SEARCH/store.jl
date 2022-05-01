function record_training_result!(c::cdl, 
                                 thresh::Real, 
                                 filter_spec, 
                                 cover_starts, 
                                 ar::activation_record, 
                                 train_num
                                 )
    @assert thresh > 0 "threshold must be greater than 0"
    # record the cover_starts for this training result
    # TODO eliminate cover_starts or make type compatible
    ar.cover_starts = cat(ar.cover_starts, cover_starts, dims=2);
    ar.T += 1;
    
    # record the activations    
    ar.codes[train_num,:,:,:] = copy(c.X_best[c.gamma+1:4:c.gamma+c.L,:,:]);
    # active_indicates =  Array(c.X_best) .> thresh;
    active_indicates = Array(ar.codes[train_num,:,:,:] .> thresh);

    indices = findall(active_indicates);

    # get all the indices
    uinds = sort(unique([i[2] for i in indices]));  # what is uinds?
    ar.uinds = append!(ar.uinds, uinds);
    uinds_len = length(uinds);
    ar.training_num = append!(ar.training_num, [train_num for _ = 1:uinds_len]);
    ar.filter_lengths = append!(ar.filter_lengths, filter_spec[uinds]);    
    println("$(length(uinds)) filters added!");
    ar.K += uinds_len;

    # for loop here...
    for n = 1:c.N
       #=   construct the indicator indices for the spatial part of the sequence
            1. first part of the indexing leaves the masked part out
            2. second part of the indexing tranform the code back to spatial 
                indexing of ACGT (instead of dummy) 
            =#
        # ar.start_positions[n] = hcat(ar.start_positions[n], active_indicates[(c.γ+1):4:c.γ+c.L, uinds, n]);            
        ar.start_positions[n] = hcat(ar.start_positions[n], active_indicates[:, uinds, n]);            
    end
end

"""
Kernel to store the "topological" results

Note that we adopt the convention that "row always occur before the column"
i.e. as in over, non_over, next_to arrays, the row (first dimension) and the 
column (second dimension) 

Inputs:
    start_pos_n: C x K matrix, C is cover length, K is number of filters (Bool)
                starting positions wrt nth pt.
    train_num: K-dimensional vector that records each filter is from which training  (UInt8)
    over: K x K x max_filter_len; store the overlap results (UInt16)
    non_over: K x K x (C-|fₛ|-1); |fₛ| is the smallest filter length. store the non-overlap results (UInt16)
    next_to: K x K; store the pairs that are exact next to each other (no gap) (UInt16)
    fil_lens: K-dimensional vector that records each filter's length (UInt8)
    cover_starts: T-dimensional vector that records where the cover actually starts for this data-point,
                   for each training; T is the total number of training (positive integer) (UInt8)

    note that the variable r,c corresponds to row, column respectively.
        The variable C is the cover size; it's doesn't have anything to do with variable c.
"""

function topology_at_nth_pt!(start_pos, train_num, over, non_over, next_to, fil_lens, cover_starts, cur_n)
    r = (blockIdx().x - 1) * blockDim().x + threadIdx().x;
    c = (blockIdx().y - 1) * blockDim().y + threadIdx().y;
    n = (blockIdx().z - 1) * blockDim().z + threadIdx().z;
    C, K, batch_size = size(start_pos);
    _,_,CmFm1 = size(non_over);
    CK = C*K;

    # TODO check: staring and end positions are inclusive, right? think so, since = 0 implies next_to
    @inbounds if r ≤ CK && c ≤ CK && n ≤ batch_size
        n_pt = cur_n + n -1;                              # the current nth-data-point        
        pr_s_cover = (r-1)%C+1;                           # starting position (in the cover)
        fr = CUDA.Int(((r-1)-(pr_s_cover-1))/C)+1;        # filter number
        tr = train_num[fr];                               # training number of this filter
        pr_s = pr_s_cover+cover_starts[n_pt,tr]-1;        # start position (in the actual data-set) "position-right-start"
        pr_e = pr_s+fil_lens[fr]-1;                       # end position (in the actual data-set)   "position-right-end"

        pc_s_cover = (c-1)%C+1;                           # similar to the above...
        fc = CUDA.Int(((c-1)-(pc_s_cover-1))/C)+1;
        tc = train_num[fc];
        pc_s = pc_s_cover+cover_starts[n_pt,tc]-1;
        
        #= process the pair only if:
         1. the pair of filters must be distinct         
         2. the order is correct: "row always occur before the column"
         3. both are actually activated
         4. both are not exact matches
        =#
        if (fr != fc) && (pr_s ≤ pc_s) && start_pos[pr_s_cover, fr, n] && start_pos[pc_s_cover, fc, n] && tr != tc 
            distance_bt = pc_s-pr_e-1;
            if !(fil_lens[fr] == fil_lens[fc] && pr_s == pc_s) # exculde the exact matches
                if distance_bt < 0 # distance_bt ≤ maxFilLen shouldn't be needed; check later
                    distance_bt = CUDA.Int(-1*distance_bt);
                    CUDA.@atomic over[fr,fc, distance_bt] += 1;
                elseif distance_bt > 0 && distance_bt ≤ CmFm1 # ignore when distance_bt > CmFm1 
                    CUDA.@atomic non_over[fr,fc,distance_bt] += 1;
                else
                    CUDA.@atomic next_to[fr,fc] += 1;
                end
            end
        end
    end
    return nothing
end

function get_batch_size_for_gpu_spatial(N::Integer, 
                                        K::Integer, 
                                        L::Integer; 
                                        use_free=0.4
                                        )
    free_bytes, _ = Mem.info();    
    bytes_needed = K*N*L;
    space_avail = use_free*free_bytes;
    (space_avail > bytes_needed) && (return N);
    return Int(floor(space_avail / (N*L)))
end

function get_spatial_info!(ar::activation_record)
    batch_size = get_batch_size_for_gpu_spatial(ar.N, ar.K, ar.C);
    CK = ar.C * ar.K; b_size_ = ceil.(Int, (CK, CK, batch_size) ./ threads_3d);
    at = activations_topology{findparam_int(ar)}(
                              ar.K, 
                              maximum_filter_length(ar), 
                              maximum_apart_length(ar));
    train_nums = cu(ar.training_num);
    fil_lens = cu(ar.filter_lengths);
    cover_starts = cu(ar.cover_starts);

    for n = 1:batch_size:ar.N
        # if n % batch_size == 1
        #     println("at $n to $(n+batch_size-1)")
        # end
        # TODO fix batch_size index out of bound error (try N=3000 to see this error) -- test this

        start_pos = cu(reduce((x,y)->cat(x,y;dims=3), ar.start_positions[n:min(n+batch_size-1, ar.N)]));
        CUDA.@sync @cuda threads=ker_3d blocks=b_size_ topology_at_nth_pt!(start_pos, train_nums, at.overlap, at.non_overlap, at.next_to, fil_lens, cover_starts, n);
        start_pos = nothing; 
    end    
    train_nums = nothing; fil_lens = nothing; cover_starts = nothing;
    return at;
end