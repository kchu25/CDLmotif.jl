struct ground_truth{T <: Integer, T2 <: Integer}
    mode_parts::Vector{Tuple{T,T}}              # (mode, k)
    num_mp::T
    covering::Dict{Tuple{T,T}, Dict{T2,T2}}     # key=(mode,k), val=Dict(n=>pos); val[n] may be empty, and that's fine
    lens::Vector{T}
end

function get_gt(motif::gapped_k_parts_motif, data::Simulated_DNA_Dataset)
    mode_parts = [(1,k) for k = 1:motif.K]
    lens = [i.len for i in motif.P_motifs.P]
    covering = Dict(mp=>Dict{Int16,Int16}() for mp in mode_parts);
    for n = 1:data.N
        if data.raw_data[n].mode > 0
            start_ = data.raw_data[n].motif_where[1];
            end_ = data.raw_data[n].motif_where[end];
            motif_where = start_:end_;
            str = data.raw_data[n].str[motif_where];
            upper_on = false; 
            pos = [];
            cur_motif = 1;
            counter = 0;
            for (ind,s) in enumerate(str)        
                if isuppercase(s) 
                    if !upper_on && counter == 0   
                        push!(pos, start_+ind-1); upper_on = true;                    
                    end      
                    counter+=1;          
                    if counter == lens[cur_motif]
                        cur_motif += 1;
                        upper_on = false;
                        counter = 0;
                    end                
                else # if it's smaller case
                    upper_on = false;
                end
            end
            if data.raw_data[n].complement
                for (ind,p) in enumerate(reverse(pos))
                    covering[mode_parts[ind]][n] = p;
                end
            else
                for (ind,p) in enumerate(pos)
                    covering[mode_parts[ind]][n] = p;
                end
            end                
        end
    end
    return ground_truth(mode_parts, length(mode_parts), covering, lens)
end

function get_gt(motif::mixture_gapped_k_parts_motifs, data::Simulated_DNA_Dataset)
    # use 0 to disregard the mode index.. 
    # TODO: maybe make it more readable later
    mode_parts = [(0,k) for k = 1:motif.motif.K];
    lens = [i.len for i in motif.motif.P_motifs.P];
    covering = Dict(mp=>Dict{Int16,Int16}() for mp in mode_parts);
    for n = 1:data.N
        if data.raw_data[n].mode > 0
            # println(n)
            start_ = data.raw_data[n].motif_where[1];
            end_ = data.raw_data[n].motif_where[end];
            motif_where = start_:end_;
            str = data.raw_data[n].str[motif_where];

            mode_ = data.raw_data[n].mode;
            mode_UnitRange = motif.modes[mode_];

            upper_on = false; 
            pos = [];
            cur_motif = minimum(mode_UnitRange);
            counter = 0;

            for (ind,s) in enumerate(str)                
                if isuppercase(s) 
                    if !upper_on && counter == 0 
                        push!(pos, start_+ind-1); upper_on = true;                    
                    end      
                    counter+=1;       
                    if counter == lens[cur_motif]                    
                        cur_motif += 1;
                        upper_on = false;
                        counter = 0;
                    end                
                else # if it's smaller case
                    upper_on = false;
                end
            end        
            
            # note that length(pos) == length(mode_UnitRange) is true
            @assert length(pos) == length(mode_UnitRange)
            if data.raw_data[n].complement
                for (ind,p) in enumerate(reverse(pos))
                    covering[(0,mode_UnitRange[ind])][n] = p;
                end
            else
                for (ind,p) in enumerate(pos)
                    covering[(0,mode_UnitRange[ind])][n] = p;
                end
            end        
        end
    end
    ground_truth(mode_parts, length(mode_parts), covering, lens)
end

function get_gt(motif::mixture_k_parts_motifs, data::Simulated_DNA_Dataset)
    mode_parts = [(i,1) for i = 1:motif.num_modes];
    lens = [i[end]-i[1]+1 for i in motif.modes];
    covering = Dict(mp=>Dict{Int16,Int16}() for mp in mode_parts);
    for n = 1:data.N
        if data.raw_data[n].mode > 0
            start_ = data.raw_data[n].motif_where[1];
            mode_ = data.raw_data[n].mode;
            covering[mode_parts[mode_]][n] = start_;
        end
    end
    return ground_truth(mode_parts, length(mode_parts), covering, lens)
end

function get_gt(motif::single_part_motif, data::Simulated_DNA_Dataset)
    mode_parts = [(1,1)];
    lens = [motif.len];
    covering = Dict(mp=>Dict{Int16,Int16}() for mp in mode_parts);
    for n = 1:data.N
        if data.raw_data[n].mode > 0
            start_ = data.raw_data[n].motif_where[1];
            mode_ = data.raw_data[n].mode;
            covering[mode_parts[mode_]][n] = start_;
        end
    end
    return ground_truth(mode_parts, length(mode_parts), covering, lens)
end

#= 
return an ordered (according to mode_parts) array of 
covering starts for each mode_parts in seqeunce n 
=#
function covering_at(g::ground_truth, n::Integer) 
    starts_at = Vector{Integer}();
    # for all the (mode, k) in the ground truth motif
    for (mode,k) in g.mode_parts
        if haskey(g.covering[(mode,k)], n)
            push!(starts_at, g.covering[(mode,k)][n]);
        else
            push!(starts_at, -1);
        end
    end
    return starts_at
end

overlap(r1::UnitRange, r2::UnitRange) = r1[1] ≤ r2[1] ≤ r1[end] || r1[1] ≤ r2[end] ≤ r1[end] || r2[1] ≤ r1[1] ≤ r2[end] || r2[1] ≤ r1[end] ≤ r2[end];
not_found(g_start::Integer, f_start::Integer) = g_start == -1 || f_start == -1;
appeared(r::UnitRange) = r[1] != -1 && r[end] != -1;
sum_range_vec(rs::Vector{UnitRange}) = [r[1] != -1 ? r[end]-r[1]+1 : 0 for r in rs];
sum_range(rs::Vector{UnitRange}) = sum(sum_range_vec(rs));

function num_overlap(r1::UnitRange, r2::UnitRange)
    @assert r1[1] ≤ r1[end] && r2[1] ≤ r2[end] "range is not valid"
    if r1[1] ≤ r2[1] ≤ r1[end]
        return min(r1[end],r2[end])-r2[1]+1;
    elseif  r2[1] ≤ r1[1] ≤ r2[end]
        return min(r1[end],r2[end])-r1[1]+1;
    else
        return 0;
    end
end

#=
tests:
    num_overlap(1:1,1:5) == 1
    num_overlap(13:16,1:5) == 0
    num_overlap(1:1,1:5) == 1
    num_overlap(1:1,1:5) == 1

=#

#= 
Return the unit ranges in bg that does not intersect with ground truths
e.g. r_bg = 1:100
     r_gts = [2:6, 8:15, 55:60]
     return [1:1, 7:7, 16:54, 61:100]     

    test:
    a = 1:100;
    r1 = 2:6;
    r2 = 8:15;
    r3 = 55:60;
    rs = [r1,r2,r3];
    get_bg_UnitRanges(a,rs)
=#
function get_bg_UnitRanges(r_bg::UnitRange, r_gts::Vector{UnitRange}) 
    set_r_gts = Set{Int}();
    for r in r_gts union!(set_r_gts, Set(r)) end
    bg_wo_r_gts = sort(collect(setdiff(Set(r_bg), set_r_gts)));

    last=0; bg = Vector{Vector{Int}}(); push!(bg,Vector{Int}()); counter = 1;
    for i in bg_wo_r_gts        
        if last == i-1 || last == 0
            push!(bg[counter], i)
        else
            push!(bg, Vector{Int}())
            counter += 1;
            push!(bg[counter], i)
        end
        last = i;
    end
    return [i[1]:i[end] for i in bg]
end

#=
increment the overlap count between each mode of the ground truth
and the pfms in a sequence 
=#
function overlap_gt_vs_pfm!(gt_overlap_matrix::Matrix{Q},
                                 r_gts,
                                 r_fs
                                 ) where Q <: Real
    for (g_ind, rg) in enumerate(r_gts)
        for (f_ind, rf) in enumerate(r_fs)        
            if rf[1] != -1
                gt_overlap_matrix[g_ind, f_ind] += num_overlap(rg,rf);            
            end
        end
    end
end

#=
increment the overlap count between the background 
and each of the pfms in a sequence 
=#
function overlap_bg_vs_pfm!(overlap_vec::Vector{Q},
                                 r_bgs,       
                                 r_fs                         
                                ) where Q <: Real
    for (_, rb) in enumerate(r_bgs)
        for (f_ind, rf) in enumerate(r_fs)
            overlap_vec[f_ind] += num_overlap(rb,rf);
        end
    end
end

function get_overlaps(gt::ground_truth, ms, data)
    L = Int(size(data.data_matrix,1)/4);
    num_fs = length(ms.pfms);
    gt_overlaps = zeros((gt.num_mp, num_fs)); # true positive 
    bg_overlaps = zeros(num_fs);
    activated_f_area = 0;
    activated_f_individial = zeros(num_fs);
    gt_area = 0;
    gt_area_each = zeros(gt.num_mp);
    bg_area = 0;

    # convention: if no key (at nth string) for pfm f, then return -1
    @inbounds for n = 1:data.N
        if data.raw_data[n].mode > 0
            gt_covering_at_n = covering_at(gt, n);
            sum_g_len = 0.0;
            region_gs = Vector{UnitRange}();
            activated_fs = Vector{UnitRange}();

            for g = 1:gt.num_mp 
                g_s = gt_covering_at_n[g]; # if absent, this is -1
                g_e = g_s != -1 ? g_s+gt.lens[g]-1 : -1;
                push!(region_gs, g_s:g_e);
                g_area = g_e-g_s+1;
                gt_area_each[g] += g_s != -1 ? g_area : 0;
                sum_g_len += g_s != -1 ? g_area : 0;
            end
            for f = 1:num_fs                            
                f_s = haskey(ms.positions[f], n) ? ms.positions[f][n] : -1;
                f_e = f_s != -1 ? f_s + ms.lens[f]-1 : -1;
                push!(activated_fs,  f_s:f_e);
            end
            region_bg = get_bg_UnitRanges(1:L, region_gs);
            overlap_gt_vs_pfm!(gt_overlaps, region_gs, activated_fs);
            overlap_bg_vs_pfm!(bg_overlaps, region_bg, activated_fs);       
            activated_f_individial = activated_f_individial .+ sum_range_vec(activated_fs);
            activated_f_area += sum_range(activated_fs);
            gt_area += sum_g_len; 
            bg_area += L-sum_g_len;
        end
    end
    return gt_area, 
           gt_area_each, 
           bg_area, 
           activated_f_area, 
           activated_f_individial, 
           gt_overlaps, 
           bg_overlaps
end

function get_overlaps_c(gt::ground_truth, ms, data)
    L = Int(size(data.data_matrix,1)/4);
    num_fs = length(ms.pfms);
    gt_overlaps = zeros((gt.num_mp, num_fs)); # true positive 
    gt_overlaps_c = zeros((gt.num_mp, num_fs)); # true positive 
    bg_overlaps = zeros(num_fs);
    activated_f_area = 0;
    activated_f_individial = zeros(num_fs);
    gt_area = 0;
    gt_area_each = zeros(gt.num_mp);
    gt_area_each_c = zeros(gt.num_mp);
    bg_area = 0;

    # convention: if no key (at nth string) for pfm f, then return -1
    @inbounds for n = 1:data.N
        if data.raw_data[n].mode > 0
            gt_covering_at_n = covering_at(gt, n);
            sum_g_len = 0.0;
            region_gs = Vector{UnitRange}();
            activated_fs = Vector{UnitRange}();

            for g = 1:gt.num_mp 
                g_s = gt_covering_at_n[g]; # if absent, this is -1
                g_e = g_s != -1 ? g_s+gt.lens[g]-1 : -1;
                push!(region_gs, g_s:g_e);
                g_area = g_e-g_s+1;
                if data.raw_data[n].complement
                    gt_area_each_c[g] += g_s != -1 ? g_area : 0;
                else
                    gt_area_each[g] += g_s != -1 ? g_area : 0;
                end                                                 
                sum_g_len += g_s != -1 ? g_area : 0;
            end
            for f = 1:num_fs                            
                f_s = haskey(ms.positions[f], n) ? ms.positions[f][n] : -1;
                f_e = f_s != -1 ? f_s + ms.lens[f]-1 : -1;
                push!(activated_fs,  f_s:f_e);
            end
            region_bg = get_bg_UnitRanges(1:L, region_gs);
            if data.raw_data[n].complement
                overlap_gt_vs_pfm!(gt_overlaps_c, region_gs, activated_fs);
            else
                overlap_gt_vs_pfm!(gt_overlaps, region_gs, activated_fs);
            end
            overlap_bg_vs_pfm!(bg_overlaps, region_bg, activated_fs);   
            activated_f_individial = activated_f_individial .+ sum_range_vec(activated_fs);
            activated_f_area += sum_range(activated_fs);
            gt_area += sum_g_len; 
            bg_area += L-sum_g_len;
        end
    end
    return gt_area, 
           gt_area_each, 
           gt_area_each_c,
           bg_area, 
           activated_f_area, 
           activated_f_individial, 
           gt_overlaps, 
           gt_overlaps_c,
           bg_overlaps
end

# return the pfm selected for each atomic motif
majority_vote(gt_overlaps::Matrix{T}) where T <: Real = [argmax(gt_overlaps[i,:]) for i = 1:size(gt_overlaps,1)];

gt_cover_count(gt, gt_overlaps, maj_vote_vec) = sum([gt_overlaps[i,m] for (i,m) in zip(1:gt.num_mp ,maj_vote_vec)]);

# the ratio of ground truth motif being covered by the majority vote chosen pfms
gt_cover_ratio(gt, gt_overlaps, gt_area, maj_vote_vec) = gt_cover_count(gt, gt_overlaps, maj_vote_vec) / gt_area;

# the ratio of incorrect activation (majority vote filters activations in the ground truth) vs all the activations
function false_cover_ratio(gt, bg_overlaps, gt_overlaps, maj_vote_vec, activated_f_area)
    (sum(bg_overlaps)+sum(gt_overlaps[:])-gt_cover_count(gt, gt_overlaps, maj_vote_vec))/activated_f_area
end

# ratio of ground truth being covered by its respective pfm
function maj_cover_ratio(mag_vec, gt_overlaps, activated_f_individial)
    [gt_overlaps[ind,m] / activated_f_individial[m] for (ind,m) in enumerate(mag_vec)]
end

false_backgroun_ratio(bg_overlaps, activated_f_individial) = bg_overlaps ./ activated_f_individial

function get_summary(ms, motif, data)
    gt = get_gt(motif, data);
    gt_area, _, _, activated_f_area, activated_f_individial, gt_overlaps, bg_overlaps = get_overlaps(gt, ms, data)
    summarize(gt, gt_overlaps, gt_area, activated_f_area, activated_f_individial, bg_overlaps);
end

function activate_cover_ratio(gt_overlaps, bg_overlaps, activated_f_individial)
    reduce(hcat, [[gt_overlaps[:,i] ./ activated_f_individial[i]; 
                bg_overlaps[i]/activated_f_individial[i]] 
                for i = 1:length(activated_f_individial)])
end

# function gt_cover_public_data(ms, motif, data::DNAdataset)
#     gt = get_gt(motif, data);
#     gt_area, gt_area_each, bg_area, activated_f_area, activated_f_individial, gt_overlaps, bg_overlaps = get_overlaps(gt, ms, data);    
#     sum_overlap = sum(gt_overlaps);
#     gt_cover_ratio = sum(gt_overlaps) ./ gt_area;
#     gt_n_cover = 1 - gt_cover_ratio;
#     a_cover = sum_overlap / activated_f_area;
#     fp_ratio = sum(bg_overlaps)/activated_f_area;
#     return gt_cover_ratio, fp_ratio, gt_n_cover, a_cover
# end

function gt2transfac(g, target_folder::String)

end

function gt_cover_by(g, target_folder::String)

    gt = get_gt(g.data.motif, g.data);
    gt_area, gt_area_each, _, 
        activated_f_area, activated_f_individial, 
            gt_overlaps, bg_overlaps = get_overlaps(gt, g.ms, g.data);

    #=
    0. Get all the performance measures    
    =#
    sum_overlap = sum(gt_overlaps);
    gt_cover_ratio = sum(gt_overlaps) ./ gt_area;
    gt_n_cover = 1 - gt_cover_ratio;
    a_cover = sum_overlap / activated_f_area;
    fp_ratio = sum(bg_overlaps)/activated_f_area;
    
    g.performance.gt_cover_perc = gt_cover_ratio;
    g.performance.false_cover_perc = gt_n_cover;
    g.performance.gt_n_cover_perc = a_cover;
    g.performance.a_cover_perc = fp_ratio;

    #=
    1. Make the CSV file that contains information on how the motif binds
        -- how each part of the ground truth is covered
    =#
    header = reshape(["by", "cover_ratio"], (1,2));
    y_axis = vcat(["D$j" for j = 1:size(gt_overlaps,2)],["NC"]);    
    for i = 1:size(gt_overlaps,1) # for each of the ground truth part
        covered_by_leaned_pfms = gt_overlaps[i,:] ./ gt_area_each[i];
        uncovered = (gt_area_each[i]-sum(gt_overlaps[i,:])) / gt_area_each[i];
        y_axis_2 = vcat(covered_by_leaned_pfms, uncovered);
        gt_cover_by_i = vcat(header, hcat(y_axis, y_axis_2))
        CSV.write(target_folder*"/gt_covered_by_$i.csv", 
            Tables.table(gt_cover_by_i), writeheader=false);
    end

    #=
    2. Plot the score contribution of of each discovered motif
    =#
    header = reshape(["discovered", "contribution"], (1,2));
    scores = [round(sum(values(g.ms.scores[i])),digits=3)
                     for i = 1:g.ms.num_motifs]; # scores of each discovered motifs
    sort_perm = sortperm(scores, rev=true);
    y_axis = ["D$j" for j = 1:g.ms.num_motifs][sort_perm];
    score_contributions = vcat(header, hcat(y_axis, scores[sort_perm]));
    CSV.write(target_folder*"/score_contributions.csv", 
                Tables.table(score_contributions), writeheader=false);

    #=
    3. Get the performance-coefficient-d for discovered motifs
    =#                
    acr = activate_cover_ratio(gt_overlaps, bg_overlaps, activated_f_individial);
    header = reshape(["target","cover_ratio"],(1,2));
    y_axis = nothing;
    if typeof(g.data.motif) == single_part_motif
        y_axis = vcat(["G1", "BG"]);
    elseif typeof(g.data.motif) == mixture_k_parts_motifs
        y_axis = vcat(["G$i" for i = 1:g.data.motif.num_modes], ["BG"])
    elseif typeof(g.data.motif) == gapped_k_parts_motif
        y_axis = vcat(["G$i" for i = 1:g.data.motif.K], ["BG"])
    elseif typeof(g.data.motif) ==  mixture_gapped_k_parts_motifs
        y_axis = vcat(["G$i" for i = 1:g.data.motif.motif.K], ["BG"])
    end
    acr_xy = [vcat(header, hcat(y_axis, acr[:,i])) for i = 1:size(acr,2)];
    for (ind,a) in enumerate(acr_xy)
        CSV.write(target_folder*"/ratio$ind.csv",  Tables.table(a), writeheader=false);
    end

    return maximum(scores);
end

function get_cover_by_fasta(g, target_folder::String)
    # -----------------------
    # header = reshape(["discovered", "contribution"], (1,2));
    scores = [round(sum(values(g.ms.scores[i])),digits=1)
                     for i = 1:g.ms.num_motifs]; # scores of each discovered motifs
    sort_perm = sortperm(scores, rev=true);

    for (ind,pfm) in enumerate(g.ms.pfms[sort_perm]) 
        CSV.write(target_folder*"/d_pfm$ind.csv",  Tables.table(pfm), writeheader=false) 
    end

    # ----------------------- evalues calculations
    evalues = get_evalues_each(g);
    evalues = [round(i, sigdigits=2) for i in evalues]
    # evalues = evalues[sort_perm];

    # -----------------------
    header = reshape(["discovered", "msa_count", "ssc", "evalues"], (1,4));
    y_axis = ["D$j" for j = 1:g.ms.num_motifs][sort_perm];
    msa_count = [length(g.ms.positions[j]) for j = 1:g.ms.num_motifs][sort_perm];
    score_contributions = vcat(header, hcat(y_axis, msa_count, scores[sort_perm], evalues[sort_perm]));
    CSV.write(target_folder*"/info.csv", 
                Tables.table(score_contributions), writeheader=false);
    

    # sort the motifs

    return maximum(scores);
end


function make_motif_type_file(motif, target_folder)
    mode_tuple_strings = [];  gaps_count = []; mixture_weights = [];
    if typeof(motif) == single_part_motif
        push!(mode_tuple_strings, "G1")
        push!(gaps_count, "-")
    elseif typeof(motif) == mixture_k_parts_motifs
        for i = 1:motif.num_modes 
            push!(mode_tuple_strings, "G$i") 
            push!(gaps_count, "-")
            push!(mixture_weights, round(motif.mixture_weights.p[i], digits=3))
        end    
    elseif typeof(motif) == gapped_k_parts_motif
        push!(mode_tuple_strings, join(["G$i" for i = 1:motif.K],","))
        gaps_count = [join(["$g" for g in motif.gap_len], ",")];
    elseif typeof(motif) == mixture_gapped_k_parts_motifs
        for i = 1:motif.num_modes 
            s = motif.modes[i][1]; e = motif.modes[i][end];
            push!(mode_tuple_strings, join(["G$j" for j = s:e],","))
            push!(gaps_count, join(["$(motif.motif.gap_len[g])" for g = s:e-1], ","))
            push!(mixture_weights, round(motif.mixture_weights.p[i], digits=3))
        end  
    end
    motif_type = isempty(mixture_weights) ? 
                    hcat(mode_tuple_strings, gaps_count,) : 
                    hcat(mode_tuple_strings, gaps_count, mixture_weights);
    CSV.write(target_folder*"/motif_type.csv",  Tables.table(motif_type), writeheader=false, delim="::");
end

function get_perform_coeff(ms,
                           motif::Union{single_part_motif, 
                                         mixture_k_parts_motifs, 
                                        gapped_k_parts_motif, mixture_gapped_k_parts_motifs},
                           data::Simulated_DNA_Dataset;
                           digits=3
                            )
    gt = get_gt(motif, data);
    gt_area, _, _, _, activated_f_area, _, gt_overlaps, gt_overlaps_c, _ = get_overlaps_c(gt, ms, data);    

    overlap_ = [maximum(gt_overlaps[i,:]) for i = 1:size(gt_overlaps,1)]
    overlap_c = [maximum(gt_overlaps_c[i,:]) for i = 1:size(gt_overlaps_c,1)]
    
    performance_coeff = (sum(overlap_)+sum(overlap_c)) / 
                    (gt_area+activated_f_area-sum(gt_overlaps)-sum(gt_overlaps_c));
    return round(performance_coeff, digits=digits)
end

function get_perform_coeff_test(ms,
                           motif::Union{single_part_motif, 
                                         mixture_k_parts_motifs, 
                                        gapped_k_parts_motif, mixture_gapped_k_parts_motifs},
                           data::Simulated_DNA_Dataset;
                           digits=3
                            )
    gt = get_gt(motif, data);
    gt_area, _, _, _, activated_f_area, _, gt_overlaps, gt_overlaps_c, _ = get_overlaps_c(gt, ms, data);    

    overlap_ = [maximum(gt_overlaps[i,:]) for i = 1:size(gt_overlaps,1)]
    overlap_c = [maximum(gt_overlaps_c[i,:]) for i = 1:size(gt_overlaps_c,1)]
    
    performance_coeff = (sum(overlap_)+sum(overlap_c)) / 
                    (gt_area+activated_f_area-sum(gt_overlaps)-sum(gt_overlaps_c));
    
    ###############
    gt_area, _, _, activated_f_area, _, gt_overlaps, bg_overlaps = get_overlaps(gt, ms, data);    
    sum_overlap = sum(gt_overlaps); gt_cover_ratio = sum(gt_overlaps) ./ gt_area;
    gt_n_cover = 1 - gt_cover_ratio; a_cover = sum_overlap / activated_f_area;
    fp_ratio = sum(bg_overlaps)/activated_f_area;

    return round(performance_coeff, digits=digits), gt_cover_ratio, fp_ratio, gt_n_cover, a_cover
end




# function get_perform_coeff_u(ms,
#     motif::Union{single_part_motif, 
#                 mixture_k_parts_motifs, 
#                 gapped_k_parts_motif, mixture_gapped_k_parts_motifs},
#     data::Simulated_DNA_Dataset;
#     digits=3
#     )
#     gt = get_gt(motif, data);
#     gt_area, gt_area_each, bg_area, activated_f_area, activated_f_individial, gt_overlaps, bg_overlaps = get_overlaps(gt, ms, data);    
    
#     performance_coeff = sum(gt_overlaps) / 
#                     (gt_area+activated_f_area-sum(gt_overlaps));
#     return round(performance_coeff, digits=digits)
# end

