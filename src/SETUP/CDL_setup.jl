#=
Input:
    Lᵧ:padded signal length  
    M: number of filters
    filter_len: length of each filter (must be multiple of 3s);
Return: 
    a "projector" matrix such that when elementwise
    multiply with the filter matrix D, it zeros out 
    the components that should be zero, as specified 
    by the constraint on D. 

Note: This should be only executed once.
=#

function make_D_projector(Lᵧ::Integer, 
                          M::Integer, 
                          filter_len::Vector{T}, 
                          S::DataType
                          ) where {T <: Integer}
    @assert S <: Real "type S must be subtype of Real"
    P = Array{S, 2}(undef, (Lᵧ, M));
    for i = 1:M
        P[:,i] = vcat([S(1.0) for i in 1:filter_len[i]], 
                [S(0f0) for j in (filter_len[i]+1):Lᵧ]);
    end    
    return P;
end

function enforce_dummy_code_constraint!(X)
    l = (blockIdx().x - 1) * blockDim().x + threadIdx().x;
    m = (blockIdx().y - 1) * blockDim().y + threadIdx().y;
    n = (blockIdx().z - 1) * blockDim().z + threadIdx().z;
    Lᵧ, M , N = size(X);
    if l ≤ Lᵧ && m ≤ M && n ≤ N
        if (l-1) % 4 != 0.0
            @inbounds X[l,m,n] = 0.0;
        end
    end
    return nothing
end

function make_d̂!(d̂, D̂)
    #= 
    It's weird that CUDA.jl does not support repeat for multidimensional arrays
    https://github.com/JuliaGPU/CUDA.jl/issues/1051
    Until this is updated, use this kernel to initialize d̂
    =#
    l = (blockIdx().x - 1) * blockDim().x + threadIdx().x;
    m = (blockIdx().y - 1) * blockDim().y + threadIdx().y;
    n = (blockIdx().z - 1) * blockDim().z + threadIdx().z;
    Lᵧ, M , N = size(d̂);
    if l ≤ Lᵧ && m ≤ M && n ≤ N
        @inbounds d̂[l,m,n] = D̂[l,m];
    end
    return nothing
end

# http://www.cs.cmu.edu/~nasmith/papers/smith+tromble.tr04.pdf
function rand_probvec(rng, d)
    unif = rand(rng, d - 1)
    T = eltype(unif)
    w = zeros(T, d + 1)
    w[2:d] .= sort(unif)
    w[d+1] = one(T)
    return diff(w)
end

rand_probvec(d) = rand_probvec(Random.GLOBAL_RNG, d);

mutable struct cdl{T <: Integer, S <: Real, Q <: Complex}
    # hyperparameters
    N::T                                              # number of data points
    M::T                                              # number of "positive and negative" filters
    filter_len::Vector{T}                             # length of each filter (must be multiple of 4s);
    filter_len_gpu::CuArray{T, 1}                       # length of each filter; gpu version;
    rho::S                                              # penalty parameter for CSC
    sigma::S                                              # penalty parameter for CDL
    λ::S                                              # sparsity regularization parameter

    # essential helper variables 
    gamma::T                   # max. filter length
    L::T
    Lᵧ::T                   # padded signal length  
    P_D::CuArray{S, 2}     # Projector; projects features to its constraint
    W::CuArray{S, 2}       # spatial mask
    WᵀW::CuArray{S, 2}     # (spatial mask)^T(spatial mask)
    WᵀS::CuArray{S, 2}

    # model parameters
    X::CuArray{S, 3}           # sparse code
    d̂::CuArray{Q, 3}        # dictionary (filter) for consensus update in the Fourier domain
    d::CuArray{S, 3}           # dictionary (filter) for consensus update in the spatial domain
    Td::CuArray{S, 2}

    # transformed model parameters
    X̂::CuArray{Q, 3}        # sparse code in the Fourier domain
    D̂ₒX̂::CuArray{Q, 2}      # elem-wise product of filters in code in the Fourier domain
    DX::CuArray{S, 2}          # product of filters and code in the spatial domain
    WDX::CuArray{S, 2}         # DX after mask W has been applied
    
    # store the parameters that achieves the best objective value
    obj_best::S
    D_best::CuArray{S, 2}
    X_best::CuArray{S, 3}

    # auxiliary variables for csc
    Y_1::CuArray{S, 2}
    Y_2::CuArray{S, 3}
    Y_3::CuArray{S, 3}
    Ŷ_1::CuArray{Q, 2}
    Ŷ_2::CuArray{Q, 3}
    Ŷ_3::CuArray{Q, 3}

    # lagragian multiplier for csc
    y_1::CuArray{S, 2}
    y_2::CuArray{S, 3}
    y_3::CuArray{S, 3}
    ŷ_1::CuArray{Q, 2}
    ŷ_2::CuArray{Q, 3}
    ŷ_3::CuArray{Q, 3}

    #  auxiliary variables for cdl
    G_1::CuArray{S, 2}
    G_2::CuArray{S, 2}
    G_3::Union{CuArray{S, 2}, Nothing}
    Ĝ_1::CuArray{Q, 2}
    Ĝ_2::CuArray{Q, 2}
    Ĝ_3::Union{CuArray{Q, 2}, Nothing}
    # lagragian multiplier for cdl
    g_1::CuArray{S, 3}
    g_2::CuArray{S, 2}
    g_3::Union{CuArray{S, 3}, Nothing}
    ĝ_1::CuArray{Q, 3}
    ĝ_2::CuArray{Q, 2}
    ĝ_3::Union{CuArray{Q, 3}, Nothing}

    # helper variables
    dots_csc::CuArray{Q, 2}
    sum_dots_csc::CuArray{Q, 2}
    Ã_csc::CuArray{Q, 2}
    RHS_csc::CuArray{Q, 3}
    alpha_csc::CuArray{Q, 2}
    dots_cdl::CuArray{Q, 3}
    sum_dots_cdl::CuArray{Q, 3}
    X_cdl::CuArray{Q, 3}
    RHS_cdl::CuArray{Q, 3}
    ifft_elem_mult::CuArray{Q, 2}
    beta_cdl::CuArray{Q, 2}
    v_avg_2::CuArray{S, 2}
    v_avg_1::CuArray{S, 2}

    #= 
    fft:
    fft_l_LᵧMN  apply fft along dimension indexed by l from an input with shape (Lᵧ, M, N)
    fft_l_LᵧM   apply fft along dimension indexed by l from an input with shape (Lᵧ, M)
    fft_l_LᵧN   apply fft along dimension indexed by l from an input with shape (Lᵧ, N) 
    =#    
    fft_l_LᵧM::CUDA.CUFFT.cCuFFTPlan{Q, -1, false, 2}  # FFT to be applied on columns of (Lᵧ,M) matrices
    fft_l_LᵧN::CUDA.CUFFT.cCuFFTPlan{Q, -1, false, 2} 
    fft_l_LᵧMN::CUDA.CUFFT.cCuFFTPlan{Q, -1, false, 3} # FFT to be applied on first dimension of (Lᵧ,M,N) tensors
    ifft_l_LᵧM::AbstractFFTs.ScaledPlan{Q, CUDA.CUFFT.cCuFFTPlan{Q, 1, false, 2}, S}
    ifft_l_LᵧN::AbstractFFTs.ScaledPlan{Q, CUDA.CUFFT.cCuFFTPlan{Q, 1, false, 2}, S}
    ifft_l_LᵧMN::AbstractFFTs.ScaledPlan{Q, CUDA.CUFFT.cCuFFTPlan{Q, 1, false, 3}, S}

    function cdl{T,S,Q}(
            M::T,                                     # number of filters
            filter_len::Vector{T},                    # length of each filter (must be multiple of 4s);
            rho::S,                               # penalty parameter for CSC
            sigma::S,                               # penalty parameter for CDL
            λ::S,                               # sparsity regularization parameter
            dat::CuArray{S, 2}                 # data
            ) where {T <: Integer, S <: Real, Q <: Complex}

        @assert size(dat,1) % 4 == 0 "first dimension of dat must be multiple of 4s"

        filter_len_gpu = cu([T(i) for i in filter_len]);
        L = convert(T, size(dat,1));
        N = convert(T, size(dat,2));
        gamma = convert(T, maximum(filter_len));                                       # max filter length
        Lᵧ = T(L + 2*gamma);                                                            # padded signal length  

        P_D = cu(make_D_projector(Lᵧ, M, filter_len, S));                              # Projector; projects features to its constraint
        D =  P_D .* cu(reduce(hcat,[reduce(vcat, [S.(rand_probvec(4)) for _ = 1:Int(Lᵧ/4)]) for _ = 1:M]));         # Initialize filters as columns 

        # since the shape input seems to have to be in Int64.. 
        LMN_shape = promote_i(Lᵧ, M, N);
        LN_shape = promote_i(Lᵧ, N);
        LL_shape = promote_i(L, Lᵧ);
        LM_shape = size(D);

        X = CUDA.randn(S, LMN_shape)             # initialize code X    
        CUDA.@sync @cuda threads=ker_3d blocks=b_size_3d(X) enforce_dummy_code_constraint!(X); # constraint X to its spatial support

        # fft
        fft_l_LᵧM = plan_fft(D, 1);                     # FFT to be applied on columns of (Lᵧ,M) matrices
        fft_l_LᵧMN = plan_fft(X, 1);                    # FFT to be applied on first dimension of (Lᵧ,M,N) tensors
        D̂ = fft_l_LᵧM * D;                              # FFT on dictionary
        X̂ = fft_l_LᵧMN * X;                             # FFT on the code
        ifft_l_LᵧM = plan_ifft(D̂, 1);
        ifft_l_LᵧMN = plan_ifft(X̂, 1);
        D̂ₒX̂ = dropdims(sum(D̂ .* X̂, dims=2),dims=2); 
        ifft_l_LᵧN = plan_ifft(D̂ₒX̂, 1);
        DX = real.(ifft_l_LᵧN * D̂ₒX̂);
        fft_l_LᵧN = plan_fft(DX, 1);

        # consensus formulation for dictionary
        d̂ = CuArray{Q}(undef, LMN_shape);
        CUDA.@sync @cuda threads=ker_3d blocks=b_size_3d(d̂) make_d̂!(d̂, D̂);
        Td = CUDA.zeros(S, LN_shape);
        d = real.(ifft_l_LᵧMN * d̂);

        # mask
        W = CUDA.zeros(S, LL_shape); 
        W[:,gamma+1:gamma+L] = Matrix(I,L,L);
        WᵀW = transpose(W)W;
        WDX = W*DX;                             # trim out the parts that aren't covered
        WᵀS = transpose(W)dat;

        # auxiliary variables for csc
        Y_1 = CUDA.zeros(S, LN_shape); Y_2 = deepcopy(X); Y_3 = deepcopy(X); 
        Ŷ_1 = fft_l_LᵧN * Y_1; Ŷ_2 = fft_l_LᵧMN * Y_2;  Ŷ_3 = deepcopy(Ŷ_2); 
        # Lagrangian multipliers for csc
        y_1 = deepcopy(Y_1); y_2 = CUDA.zeros(S, LMN_shape); y_3 = deepcopy(y_2);
        ŷ_1 = deepcopy(Ŷ_1); ŷ_2 = fft_l_LᵧMN * y_2; ŷ_3 = deepcopy(ŷ_2); 
        # axuiliary variables for CDL
        G_1 = deepcopy(D); G_2 = deepcopy(DX); G_3 = deepcopy(D);
        Ĝ_1 = fft_l_LᵧM * G_1; Ĝ_2 = fft_l_LᵧN * G_2; Ĝ_3 = deepcopy(Ĝ_1);
        # Lagrangian multipliers for cdl
        g_1 = deepcopy(y_2); g_2 = deepcopy(Y_1); g_3 = deepcopy(y_2);
        ĝ_1 = deepcopy(ŷ_2); ĝ_2 = deepcopy(Ŷ_1); ĝ_3 = deepcopy(ŷ_2);

        # store the parameters that achieves the best objective value
        D_best = deepcopy(G_1);
        X_best = deepcopy(Y_3);

        obj_best = 0.5*norm(W * DX .- dat)^2 + λ*reduce(+, abs.(X));

        # helper variables
        dots_csc = CuArray{Q}(undef, LM_shape);
        sum_dots_csc = CuArray{Q}(undef, (size(D,1),1));
        Ã_csc = CuArray{Q}(undef, LM_shape);
        RHS_csc = CuArray{Q}(undef, LMN_shape);
        alpha_csc = CuArray{Q}(undef, LN_shape);
        dots_cdl = CuArray{Q}(undef, LMN_shape);
        sum_dots_cdl = CuArray{Q}(undef, (LMN_shape[1],1,LMN_shape[3]));
        X_cdl = CuArray{Q}(undef, LMN_shape);
        RHS_cdl = CuArray{Q}(undef, LMN_shape); 
        ifft_elem_mult = CuArray{Q}(undef, LN_shape);
        beta_cdl = CuArray{Q}(undef, LN_shape);
        v_avg_2 = CuArray{S}(undef, LM_shape);
        v_avg_1 = CuArray{S}(undef, LM_shape);

        # free up memory
        CUDA.unsafe_free!(D); CUDA.unsafe_free!(D̂);

        new(N, M, filter_len, filter_len_gpu, 
            rho, sigma, λ, gamma, L, Lᵧ, P_D, W, WᵀW, WᵀS, 
            X, d̂, d, Td, X̂, D̂ₒX̂, DX, WDX, obj_best, 
            D_best, X_best, Y_1, Y_2, Y_3, Ŷ_1, 
            Ŷ_2, Ŷ_3, y_1, y_2, y_3, ŷ_1, ŷ_2, 
            ŷ_3, G_1, G_2, G_3, Ĝ_1, Ĝ_2, Ĝ_3, 
            g_1, g_2, g_3, ĝ_1, ĝ_2, ĝ_3, 
            dots_csc, sum_dots_csc, Ã_csc,
            RHS_csc, alpha_csc, dots_cdl,
            sum_dots_cdl, X_cdl, RHS_cdl,
            ifft_elem_mult, beta_cdl, 
            v_avg_2, v_avg_1, fft_l_LᵧM, 
            fft_l_LᵧN, fft_l_LᵧMN, ifft_l_LᵧM, 
            ifft_l_LᵧN, ifft_l_LᵧMN);
    end
end

function init_cdl(cdl_setup::CDL_setup{T,S,Q}, 
                  dna_dataset::Union{Sim_DNA{T,S}, FASTA_DNA{T,S}}
                  ) where {T <: Integer, S <: Real, Q <: Complex}
    cdl{T,S,Q}(cdl_setup.M, 
               [T(4*cdl_setup.filter_size) for _ = 1:cdl_setup.M],
               cdl_setup.rho,
               cdl_setup.sigma,
               cdl_setup.lambda,
               dna_dataset.data_matrix_gpu
               );
end

cdl_obj(c::cdl, signals) = 0.5*norm(c.W * c.Td .- signals)^2;
objective(c::cdl, signals) = 0.5*norm(c.W * c.DX .- signals)^2 + c.λ*reduce(+, abs.(c.X));

findparam_int(c::cdl{T,S,Q}) where {T,S,Q} = T;
findparam_real(c::cdl{T,S,Q}) where {T,S,Q} = S;
findparam_complex(c::cdl{T,S,Q}) where {T,S,Q} = Q;