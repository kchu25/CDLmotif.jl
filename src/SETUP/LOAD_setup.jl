abstract type Abstract_DNA_dataset end
abstract type Simulated_DNA_Dataset <: Abstract_DNA_dataset end
abstract type Real_DNA_Dataset <: Abstract_DNA_dataset end

function get_data_matrices(dna_read, fasta_location::String, target_folder::String, S)
    # run(`./scripts_bash/create_bg.sh  "$(fasta_loc)" "$target_folder"`);
    out_str = "$target_folder/shuffled_bg.fa"
    run(`fasta-shuffle-letters -kmer 3 "$fasta_location" "$out_str"`);
    # TODO find a way to replace the above later?
    # println("$target_folder/shuffled_bg.fa")
    data_matrix = data_2_dummy(dna_read; F=S);
    data_matrix_bg = data_2_dummy(read_fasta("$target_folder/shuffled_bg.fa"));
    return data_matrix, data_matrix_bg
end

# function create_data_matrix(motif, num_data_pts, data_pt_len, bern_prob, S)
#     raw_data = sample_backgound_with_motif_multiple(motif,num_data_pts,data_pt_len,bern_prob);
#     data_2_dummy(raw_data; F=S);
# end


struct Sim_DNA{T <: Integer, S <: Real} <: Simulated_DNA_Dataset
    raw_data::Vector{sim_dna_str_w_motif}
    motif::motif_type
    N::T
    L::T
    data_matrix::Matrix{S}
    data_matrix_gpu::Union{Nothing, CuArray{S,2}}
    data_matrix_bg::Union{Nothing, Matrix{S}}
    target_folder::Union{Nothing, String}
    prob_per_seq::S

    function Sim_DNA{T, S}(motif::motif_type, 
                        num_data_pts::Integer,
                        output_folder::String,
                        data_pt_len::Integer=100,
                        bern_prob::Real=1.0
        ) where {T <: Integer, S <: Real}        
        # simulate the data
        raw_data = sample_backgound_with_motif_multiple(motif,num_data_pts,data_pt_len,bern_prob);
        # save the simulated data as a fasta file in the target folder
        save_sim_data_as_fasta(output_folder, raw_data, num_data_pts, motif);
        data_matrix, data_matrix_bg = get_data_matrices(raw_data, 
                                                        "$output_folder/data.fa", 
                                                        output_folder, S);
        new(
            raw_data,
            motif,
            T(num_data_pts),
            T(size(data_matrix,1)/4),
            data_matrix,
            cu(data_matrix),
            data_matrix_bg,
            output_folder,
            S(bern_prob)            
        )
    end
    
    function Sim_DNA{T,S}(motif::motif_type, 
                          num_data_pts::Integer, 
                          data_pt_len::Integer=100,
                          bern_prob::Real=1.0) where {T <: Integer, S <: Real}     
        raw_data = sample_backgound_with_motif_multiple(motif,num_data_pts,data_pt_len,bern_prob);
        data_matrix = data_2_dummy(raw_data; F=S);
        new(
            raw_data, motif, T(num_data_pts), T(size(data_matrix,1)/4),
            data_matrix, nothing, nothing, nothing, S(bern_prob)
        )                          
    end
end

struct FASTA_DNA{T <: Integer, S <: Real} <: Real_DNA_Dataset
    N::T
    L::T
    raw_data::Vector{String}
    data_matrix::Matrix{S}
    data_matrix_gpu::CuArray{S,2}
    data_matrix_bg::Matrix{S}
    target_folder::String
    function FASTA_DNA{T, S}(fasta_location::String, 
                        output_folder::String;
                        max_entries=max_num_read_fasta
                        ) where {T <: Integer, S <: Real}       
    !isdir(output_folder) && mkpath(output_folder); 
    dna_read = read_fasta(fasta_location; max_entries);
    data_matrix, data_matrix_bg = get_data_matrices(dna_read, 
                                                    fasta_location, 
                                                    output_folder, S);
    new(        
        T(length(dna_read)),
        T(size(data_matrix,1)/4),
        dna_read,
        data_matrix,
        cu(data_matrix),
        data_matrix_bg,
        output_folder)          
    end
end

get_N(d::Abstract_DNA_dataset) = d.N;
get_L(d::Abstract_DNA_dataset) = d.L;
get_data_matrix(d::Abstract_DNA_dataset) = d.data_matrix;
get_data_matrix_bg(d::Abstract_DNA_dataset) = d.data_matrix_bg;
get_target_folder(d::Abstract_DNA_dataset) = d.target_folder;
is_simulated_data(d::Abstract_DNA_dataset) = typeof(d) <: Simulated_DNA_Dataset ? true : false;
findparam_real(d::FASTA_DNA{T,S}) where {T,S} = S;
findparam_real(d::Sim_DNA{T,S}) where {T,S} = S;

