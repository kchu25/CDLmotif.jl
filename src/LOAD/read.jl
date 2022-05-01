###### function that reads a fasta file ################
function reading(filepath::String, max_entries=max_num_read_fasta)
    # read the file
    reads = read(filepath, String);
    # process the fasta file to get the DNA part
    # rid of sequences that contains "n"
    dna_reads = Vector{String}();
    for i in split(reads, '>')
        if !isempty(i)
            this_read = join(split(i, "\n")[2:end]);        
            if !occursin("N", this_read) && !occursin("n", this_read)
                push!(dna_reads, this_read);
            end
        end
    end    
    # dna_reads = [join(split(i, "\n")[2:end]) for i in split(reads, '>') if !isempty(i) && !occursin("N", i)]; 
    dna_reads = length(dna_reads) > max_entries ? dna_reads[1:max_entries] : dna_reads;
    return dna_reads
end

function read_fasta(filepath::String; 
                    max_entries=max_num_read_fasta
                    )::Vector{String}
    #= read a fasta file =#
    dna_reads = reading(filepath, max_entries);   
    # convert all DNA seqeunce to uppercase
    return [uppercase(i) for i in dna_reads]
end

function read_fasta_that_has_ground_truth(filepath::String, 
                                        jaspar_mat_path::String; 
                                        max_entries=max_num_read_fasta
    )
    #= read the fasta file that contains the ground truth 
    (ground truths are in capital letters, as fasta files 
    from public databases typically contains it), return 
    1. dna reads and 2. ground truth 

    note: ok maybe it's not really ground truth since 
    by defn ground truth is unknown but we still keep that 
    labeled info
    =#
    
    #= read a fasta file =#
    dna_reads = reading(filepath, max_entries);   
    read_len = length(dna_reads[1]);
    @assert all(length.(dna_reads) .== read_len) "seqs must have same length"
    raw_data = Vector{NamedTuple{(:str, :motif_where, :mode), 
                    Tuple{String, UnitRange{Int64}, Int64}}}();
    num_ground_truth = 0;
    for str in dna_reads
        characters = [i for i in str];
        motif_found = findall(isuppercase.(characters) .== 1);
        num_ground_truth += isempty(motif_found) ? 0 : 1;
        motif_start = minimum(motif_found); motif_end = maximum(motif_found);
        motif_range =  motif_start:motif_end;        
        push!(raw_data, (str=str, motif_where=motif_range, mode=1));
    end

    OOPS = num_ground_truth == length(dna_reads);
    motif = single_part_motif(get_JASPAR_PFM(jaspar_mat_path));

    dna_reads = uppercase.(dna_reads);
    
    return uppercase.(dna_reads), motif, raw_data, OOPS
end

function get_JASPAR_PFM(jaspar_path::String)
    read_jaspar = read(jaspar_path, String);    
    s = split(read_jaspar, "\n");
    parsed = [[parse(Float64, i.match) for i in eachmatch(nmbr_regex, s[j])] for j = 2:length(s)-1];
    mat = reduce(hcat, parsed);
    return mat ./ sum(mat, dims=2)
end
########################################################

###### reverse complement of a string ##################
function reverse_complement(s::String)    
    join(islowercase(s[si]) ? s[si] : DNA_complement[s[si]] for si = length(s):-1:1)
end
########################################################

