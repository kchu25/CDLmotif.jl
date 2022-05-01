function save_gt_motifs(target_folder::String, motif::mixture_gapped_k_parts_motifs)    
    for j = 1:length(motif.motif.P_motifs.P)
        gt = reduce(hcat, [motif.motif.P_motifs.P[j].P.pc[i].p for i = 1:length(motif.motif.P_motifs.P[j].P.pc)]);
        CSV.write(target_folder*"/pfm_gt$j.csv",  Tables.table(gt), writeheader=false);
    end     
end

function save_gt_motifs(target_folder::String, motif::gapped_k_parts_motif)    
    for j = 1:motif.K
        gt = reduce(hcat, [motif.P_motifs.P[j].P.pc[i].p for i = 1:length(motif.P_motifs.P[j].P.pc)]);
        CSV.write(target_folder*"/pfm_gt$j.csv",  Tables.table(gt), writeheader=false);
    end
end

function save_gt_motifs(target_folder::String, motif::mixture_k_parts_motifs)    
    for j = 1:motif.num_modes
        gt = reduce(hcat, [motif.motif.P.pc[i].p for i = motif.modes[j]]);
        CSV.write(target_folder*"/pfm_gt$j.csv",  Tables.table(gt), writeheader=false);
    end
end

function save_gt_motifs(target_folder::String, motif::single_part_motif)    
    gt = reduce(hcat, [motif.P.pc[i].p for i = 1:length(motif.P.pc)]);
    CSV.write(target_folder*"/pfm_gt1.csv",  Tables.table(gt), writeheader=false);
end

# _m_ = single_part_motif(12);
# motif=_m_;
# gt = reduce(hcat, [motif.P.pc[i].p for i = 1:length(motif.P.pc)]);
# gt = reduce(hcat, [motif.P.pc[i].p for i = 5:11]);
# q = (floor.(gt .* 10000));
# q = (floor.(gt .* 10000)).+900;
# io = open(target_folder*"/pfm_gt123.csv", "w")
# # println(io, "ID\t")
# # println(io, "XX\t")
# # println(io, "BF\t")
# # println(io, "XX\t")
# println(io, "P0\tA\tC\tG\tT")
# for i = 1:size(q,2)
#     if i < 10
#         println(io, "0"*"$i\t$(q[1,i])\t$(q[2,i])\t$(q[3,i])\t$(q[4,i])")
#     else
#         println(io, "1"*"$i\t$(q[1,i])\t$(q[2,i])\t$(q[3,i])\t$(q[4,i])")
#     end
# end
# println(io, "XX\t")
# close(io)



function save_sim_data_as_fasta(target_folder::String, 
                                raw_data, 
                                N::Integer,
                                motif)
    !isdir(target_folder) && mkpath(target_folder);
    save_gt_motifs(target_folder, motif);

    strings = [raw_data[i].str for i = 1:N];
    strings_all_upcase = [uppercase.(raw_data[i].str) for i = 1:N];
    
    open(target_folder*"/data_w_answer.fa","w") do file
        for (ind,s) in enumerate(strings)
            write(file, string(">sequence_", string(ind),"\n"));
            write(file, string(s,"\n"))
        end
    end
    open(target_folder*"/data.fa","w") do file
        for (ind,s) in enumerate(strings_all_upcase)
            write(file, string(">sequence_", string(ind),"\n"));
            write(file, string(s,"\n"))
        end
    end
end





