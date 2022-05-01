######## dict update using Sherman Morrison trick ##########
function sherman_CDL_RHS!(RHS, X̂, Ĝ_1, ĝ_1, Ĝ_2, ĝ_2, Ĝ_3, ĝ_3)
    l = (blockIdx().x - 1) * blockDim().x + threadIdx().x;
    m = (blockIdx().y - 1) * blockDim().y + threadIdx().y;
    n = (blockIdx().z - 1) * blockDim().z + threadIdx().z;
    Lᵧ, M , N = size(X̂);
    if l ≤ Lᵧ && m ≤ M && n ≤ N        
        @inbounds RHS[l,m,n] = conj(X̂[l,m,n])*(Ĝ_2[l,n]-ĝ_2[l,n])+(Ĝ_1[l,m]-ĝ_1[l,m,n])+(Ĝ_3[l,m]-ĝ_3[l,m,n]);
    end
    return nothing
end

function sherman_CDL_dot!(β, X̃, RHS)
    l = (blockIdx().x - 1) * blockDim().x + threadIdx().x;
    m = (blockIdx().y - 1) * blockDim().y + threadIdx().y;
    n = (blockIdx().z - 1) * blockDim().z + threadIdx().z;
    Lᵧ, M , N = size(X̃);
    if l ≤ Lᵧ && m ≤ M && n ≤ N
        @inbounds β[l,n] = sum(X̃[l,:,n] .* RHS[l,:,n]);
    end
    return nothing
end

function sherman_CDL_combine(d̃, RHS, β, X)
    l = (blockIdx().x - 1) * blockDim().x + threadIdx().x;
    m = (blockIdx().y - 1) * blockDim().y + threadIdx().y;
    n = (blockIdx().z - 1) * blockDim().z + threadIdx().z;
    Lᵧ, M , N = size(d̃);
    if l ≤ Lᵧ && m ≤ M && n ≤ N
        @inbounds d̃[l,m,n] = RHS[l,m,n]-β[l,n]*conj(X[l,m,n]);
    end
    return 
end
############################################################

############################################################
function projection2box!(G_3, v̄, P_D)
    l = (blockIdx().x - 1) * blockDim().x + threadIdx().x;
    m = (blockIdx().y - 1) * blockDim().y + threadIdx().y;
    Lᵧ, M = size(G_3);
    @inbounds if l ≤ Lᵧ && m ≤ M
        if v̄[l,m] ≥ 1
            G_3[l,m] = P_D[l,m] != 0 ? 1 : 0;
        elseif v̄[l,m] ≤ 0
            G_3[l,m] = 0;
        else
            G_3[l,m] = P_D[l,m] != 0 ? v̄[l,m] : 0;
        end
    end
    return nothing;
end

function prox_simplex_zeros!(g_1, v̄, P_D)
    l = (blockIdx().x - 1) * blockDim().x + threadIdx().x;
    m = (blockIdx().y - 1) * blockDim().y + threadIdx().y;
    Lᵧ, M = size(g_1);
    if l ≤ Lᵧ && m ≤ M
        if P_D[l,m] != 0
            i = 4*CUDA.Int32(floor(l/4));
            @inbounds g_1[l,m] = v̄[l,m] - 0.25*(v̄[i+1,m]+v̄[i+2,m]+v̄[i+3,m]+v̄[i+4,m]-1);
        else
            @inbounds g_1[l,m] = 0;
        end
    end
    return nothing
end

function update_g₂!(g₂, WᵀW, WᵀS, Td̃, η₂, sigma)
    l = (blockIdx().x - 1) * blockDim().x + threadIdx().x;
    n = (blockIdx().y - 1) * blockDim().y + threadIdx().y;
    Lᵧ, N = size(Td̃);
    if l ≤ Lᵧ && n ≤ N
        @inbounds g₂[l,n] = (WᵀS[l,n]+sigma*(Td̃[l,n] + η₂[l,n])) / (WᵀW[l,l] + sigma);
    end
    # as we can see that components outside of the support better be small
    return nothing
end
############################################################