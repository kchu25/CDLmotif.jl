function get_batch_size_for_gpu_pair_merge(num_pairs::Integer, 
                                           N::Integer, 
                                           L::Integer,
                                           type_real; 
                                           use_free=0.25                                           
                                           )::Integer
    # TODO make a description about this
    @assert type_real == Float32 || type_real == FLoat16 "type_real must be Float32 or FLoat16"
    multiplier = type_real == Float32 ? 4 : 2;
    free_bytes, _ = Mem.info();    
    bytes_needed = multiplier*num_pairs*N*L;    # 
    space_avail = use_free*free_bytes; # use 25% of the free space
    (space_avail > bytes_needed) && (return num_pairs);
    return Int(floor(space_avail / (multiplier*N*L))) 
end

####################### use gpu to get pairwise information #######################
function get_pair_pfms!(at::activations_topology, 
                       ar::activation_record, 
                       g::good_stuff)
    type_real = findparam_real(ar);
    type_int = findparam_int(ar);
    max_span_len = 2*g.cdl.filter_size-1;
    all_pfms = CuArray{type_real, 3}(undef, (0, 4, max_span_len)); 
    span_lens = CuArray{type_int, 1}(undef, 0);
    # setup kernels     
    z = at.overlap .> g.cdl.count_threshold;

    indices = Array(findall(z)); num_pairs = length(indices);    
    batch_size = get_batch_size_for_gpu_pair_merge(num_pairs, ar.N, ar.C, type_real);
    ar_start_positions = cu(reduce((x,y)->cat(x,y;dims=3), ar.start_positions));
    ar_fil_lens = cu(ar.filter_lengths);
    ar_cover_starts = cu(ar.cover_starts);
    ar_codes = ar.codes;
    ar_train_num = cu(ar.training_num);
    ar_uinds = cu(ar.uinds)

    for i = 1:batch_size:num_pairs
        # prepare batch size information
        batch_range = i:min(i+batch_size-1,num_pairs);
        batch_range_size = batch_range[end] == num_pairs ? batch_range[end]-batch_range[1]+1 : batch_size;
        # init data structure to contain pair merge information
        indices_gpu = cu(reduce(vcat, [type_int.([i[1] i[2] i[3]]) for i in indices[batch_range]]));
        filter_pair_starts = CUDA.fill(false, promote_i(batch_range_size, ar.N, ar.C));
        span_len = CUDA.zeros(type_int, batch_range_size);
        b_size_ = ceil.(Int, promote_i(batch_range_size, ar.N) ./ threads_2d);

        CUDA.@sync @cuda threads=ker_2d blocks=b_size_ get_f_pair_starts(indices_gpu, 
                                                                         ar_start_positions, 
                                                                         ar_fil_lens, 
                                                                         filter_pair_starts);

        CUDA.@sync @cuda threads=ker_1d blocks=ceil(Int,batch_range_size/threads_1d) get_span_len(span_len, 
                                                                                       ar_fil_lens, 
                                                                                       indices_gpu);        
        # MSA
        msas = CUDA.zeros(type_real, promote_i(batch_range_size, 4*max_span_len, ar.N)); # TODO: change this to Int16? (faster)
        CUDA.@sync @cuda threads=ker_3d blocks=b_size_3d(msas) msa_fill(msas, 
                                                                     filter_pair_starts, 
                                                                     span_len, 
                                                                     g.data.data_matrix_gpu, 
                                                                     ar_codes, 
                                                                     ar_train_num, 
                                                                     ar_fil_lens, 
                                                                     indices_gpu, 
                                                                     ar_uinds);
        # sum up the counts
        msas = dropdims(sum(msas, dims=3),dims=3); 
        CUDA.@sync @cuda threads=ker_2d blocks=b_size_2d(msas) msa_pseudocount(msas, span_len);

        # divide by the total number that adds up the counts
        msas = reshape(msas, (size(msas,1), 4, Int(size(msas,2)/4)));          
        @assert length(findall(msas .< 0)) == 0 "found negative entries in msa"     

        pfms = msas ./ dropdims(sum(msas[:,:,1],dims=2),dims=2); 

        span_lens = cat(span_lens, span_len, dims=1);
        all_pfms = cat(all_pfms, pfms, dims=1);

        msas = nothing; 
        span_len = nothing; 
        filter_pair_starts = nothing; 
        indices_gpu = nothing; 
        pfms = nothing;
    end

    # free up memory 
    ar_start_positions = nothing; 
    ar_fil_lens = nothing; 
    ar_cover_starts = nothing; 
    # GC.gc(true); GC.gc(false); CUDA.reclaim();    

    # # return a motifs' struct
    pfms_cpu = Array(all_pfms);
    num_pfms = size(pfms_cpu,1);
    span_len_cpu = Array(span_lens);
    pfms = [pfms_cpu[i,1:4,1:span_len_cpu[i]] for i = 1:num_pfms];
    pwms =  get_pwms(pfms, g.search.bg);
    thresh = get_thresh(pwms, g.search.pval_thresh, findparam_real(g));
    lens = type_int.([size(pfm,2) for pfm in pfms]);

    g.ms = motifs(
            pfms,
            pwms,
            thresh,
            lens,
            type_int(num_pfms),
            nothing, 
            nothing,
            nothing 
            # Dict{type_int, Bool}(i=>false for i = 1:num_pfms),
            # lens
        );
end


# function first_scan_pfms_return_motifs(all_pfms, span_lens, pval_thresh, N, L, data_dat_gpu)
#     # back to cpu for the prep of p-value/score computations
#     pfms_cpu = Array(all_pfms);
#     num_pfms = size(pfms_cpu,1);
#     span_len_cpu = Array(span_lens);
#     pfms_batch = [pfms_cpu[i,1:4,1:span_len_cpu[i]] for i = 1:num_pfms];
#     pwms_cpu = [pfm2pwm(pfm) for pfm in pfms_batch];
#     thresh = cu([pvalue2score(pwm, pval_thresh) for pwm in pwms_cpu]);
#     # prep for the scan
#     pwms = CUDA.zeros(Float16, size(all_pfms)); # TODO change type
#     CUDA.@sync @cuda threads=t1 blocks=ceil(Int, num_pfms/t1) pfm_arr2pwm_arr(all_pfms,pwms);
#     # scan
#     pos_scores = CUDA.zeros(Float16, num_pfms, N, L); # TODO change type
#     CUDA.@sync @cuda threads=ker_3d blocks=b_size_3d(pos_scores) greedy_search!(pwms, 
#                                                                  data_dat_gpu, 
#                                                                  span_lens, 
#                                                                  pos_scores,
#                                                                  thresh)
#     pos_scores_pos = dropdims(map(x-> x[3]==1 ? Int32(0) : Int32(x[3]), argmax(pos_scores, dims=3)),dims=3);
#     # pos_scores_cpu = Array(pos_scores); # don't need scores (for now)
#     pos_scores_pos_cpu = Array(pos_scores_pos);
#     found = Array(findall(pos_scores_pos_cpu .> 0));
#     # record the positions
#     positions = [Dict{Int16, Int16}() for _ = 1:num_pfms]; # TODO change type
#     # scores = [Dict{Int16, Float64}() for _ = 1:num_pfms]; # don't need scores (for now)
#     for f in found
#         positions[f[1]][f[2]] = pos_scores_pos_cpu[f[1],f[2]];
#     end

#     # free up memory
#     pwms = nothing; pos_scores = nothing; thresh = nothing; pos_scores_pos = nothing;
#     GC.gc(true); CUDA.reclaim(); GC.gc(true); GC.gc(false); CUDA.reclaim();   

#     return motifs(
#         pfms_batch,
#         Int16.([size(pfm,2) for pfm in pfms_batch]),
#         num_pfms,
#         positions, 
#         nothing, 
#         Dict{Int16, Bool}(i=>false for i = 1:num_pfms), # TODO change type
#         Int16.([size(pfm,2) for pfm in pfms_batch])
#     )
# end

function get_f_pair_starts(indices_gpu, 
                           ar_start_positions, 
                           ar_fil_lens, 
                           filter_pair_starts
                           )
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x;
    n = (blockIdx().y - 1) * blockDim().y + threadIdx().y;
    i_max = size(indices_gpu, 1);
    C, _, N = size(ar_start_positions);
    if i ≤ i_max && n ≤ N
        fil_ind_1 = indices_gpu[i,1]; fil_ind_2 = indices_gpu[i,2]; o_len = indices_gpu[i,3];
        @inbounds for c = 1:C
            #= if fil2 is activated at c, it's now at an "allowed" place to examine how fil1 is activated, 
            and fil1 is activated at o_len+c-ar_fil_lens[fil_ind_1] =#
            fil_1_at = o_len+c-ar_fil_lens[fil_ind_1];
            if ar_start_positions[c,fil_ind_2,n]==true && fil_1_at > 0 && ar_start_positions[fil_1_at,fil_ind_1,n]==true
                filter_pair_starts[i,n,fil_1_at] = true;
            end
        end 
    end
    return nothing
end

function get_span_len(span_len, ar_fil_lens, indices)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x;
    i_max = size(indices, 1);
    if i ≤ i_max
        # span_len = overlap_len ≥ f2_len ? overlap_len : f1_len - overlap_len + f2_len;
        span_len[i] = indices[i,3] ≥ ar_fil_lens[indices[i,2]] ? indices[i,3] : ar_fil_lens[indices[i,1]]-indices[i,3]+ar_fil_lens[indices[i,2]];
    end
    return nothing
end

function msa_fill(msas, filter_pair_starts, 
                  span_len, data_dat_gpu, 
                  ar_codes, ar_train_num, 
                  ar_fil_lens, indices_gpu, 
                  ar_uinds
                  )
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x; # ith pair 
    l = (blockIdx().y - 1) * blockDim().y + threadIdx().y; # lth position in msa (4*l)
    n = (blockIdx().z - 1) * blockDim().z + threadIdx().z; # nth string 
    i_max, l_max, N = size(msas);
    C = size(filter_pair_starts, 3);
    if i ≤ i_max && l ≤ l_max && n ≤ N
        # training number for both filters
        fil_1_t = ar_train_num[indices_gpu[i,1]]; 
        fil_2_t = ar_train_num[indices_gpu[i,2]];
        # the filter number
        fil_1 = ar_uinds[indices_gpu[i,1]];
        fil_2 = ar_uinds[indices_gpu[i,2]];
        # if it actually is activate here and does not cross the boundary
        # indices_gpu[i,3] is overlap between the two filters
        @inbounds for c = 1:(C-span_len[i]+1)
            # if we have a start right here, and l is within the span len
            if filter_pair_starts[i,n,c] != 0 && l ≤ 4*span_len[i] 
                offset = (c-1)*4;            
                fil_2_start = (c + ar_fil_lens[indices_gpu[i,1]]) - indices_gpu[i,3];
                code_coeff_avg = (ar_codes[fil_1_t,c,fil_1,n] + ar_codes[fil_2_t,fil_2_start,fil_2,n])/2;
                # code_coeff_avg < 0 && (println(c); break)
                msas[i,l,n] += code_coeff_avg*(data_dat_gpu[offset+l,n]);
            end
        end
    end
    return nothing
end

function msa_pseudocount(msas_2d, span_len)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x; # ith pair 
    l = (blockIdx().y - 1) * blockDim().y + threadIdx().y; # lth position
    i_max, l_max = size(msas_2d);
    if i ≤ i_max && l ≤ l_max && l ≤ 4*span_len[i]
        @inbounds msas_2d[i,l] += 0.1;
    end
    return nothing
end

# function pfm_arr2pwm_arr(pfms, pwms)
#     i = (blockIdx().x - 1) * blockDim().x + threadIdx().x; # ith pair 
#     I, A, maxlen = size(pfms);
#     if i ≤ I
#         for a = 1:A
#             for l = 1:maxlen
#                 pwms[i,a,l] = pfms[i,a,l] != 0 ? log2(pfms[i,a,l]/.25) : pwms[i,a,l];
#             end
#         end
#     end
#     return nothing
# end

