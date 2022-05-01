
######## Code update using Sherman Morrison trick ##########
#=  This part implments Brendt Wohlberg's 
    Efficient convolutional sparse coding
https://ieeexplore.ieee.org/abstract/document/6854992 
=#

function sherman_CSC_RHS!(result, D̂, Y1, y1, Y2, y2, Y3, y3)
    l = (blockIdx().x - 1) * blockDim().x + threadIdx().x;
    m = (blockIdx().y - 1) * blockDim().y + threadIdx().y;
    n = (blockIdx().z - 1) * blockDim().z + threadIdx().z;
    Lᵧ, M, N = size(Y2); # V₂ is of size (Lᵧ,M,N)
    if l ≤ Lᵧ && m ≤ M && n ≤ N
        @inbounds result[l,m,n] = conj(D̂[l,m])*(Y1[l,n]-y1[l,n])+(Y2[l,m,n]-y2[l,m,n])+(Y3[l,m,n]-y3[l,m,n]);
    end
    return nothing
end

function sherman_CSC_dot!(re, left, right)
    #= 
    re is the result matrix that stores the result of dot product. (Lᵧ,N)
    left is the precomputed sherman_left horizontal vectors (Lᵧ,M)
    right is the precomputed sherman_right vectors (Lᵧ,M,N)
    =# 
    l = (blockIdx().x - 1) * blockDim().x + threadIdx().x;
    n = (blockIdx().y - 1) * blockDim().y + threadIdx().y;
    Lᵧ, M, N = size(right);
    if l ≤ Lᵧ && n ≤ N
        re[l,n] = 0;
        for m in 1:M
            @inbounds re[l,n] += left[l,m]*right[l,m,n];
        end
    end
    return nothing
end

function sherman_CSC_combine!(solution, dots, D̂, right, rho)
    #= 
    dots is (Lᵧ,N)
    =#
    l = (blockIdx().x - 1) * blockDim().x + threadIdx().x;
    m = (blockIdx().y - 1) * blockDim().y + threadIdx().y;
    n = (blockIdx().z - 1) * blockDim().z + threadIdx().z;
    Lᵧ, M , N = size(solution);
    if l ≤ Lᵧ && m ≤ M && n ≤ N
        @inbounds solution[l,m,n] = (right[l,m,n]-dots[l,n]*conj(D̂[l,m]))/rho;
    end
    return nothing
end
############################################################

############# l1-norm (sparsity) update ####################
function shrink!(Y₂, X, y₂, λ)
    # source is a "big" matrix (ML x N) 
    l = (blockIdx().x - 1) * blockDim().x + threadIdx().x;
    m = (blockIdx().y - 1) * blockDim().y + threadIdx().y;   
    n = (blockIdx().z - 1) * blockDim().z + threadIdx().z; 
    Lᵧ, M , N = size(Y₂)
    if l ≤ Lᵧ && m ≤ M && n ≤ N 
            tmp = X[l,m,n]+y₂[l,m,n];
            @inbounds Y₂[l,m,n] = sign(tmp) * max(abs(tmp) - λ, 0);
    end
    return nothing
end
############################################################

############# projection update ############################
function projectionX!(Y₃, X, y₃)
    l = (blockIdx().x - 1) * blockDim().x + threadIdx().x;
    m = (blockIdx().y - 1) * blockDim().y + threadIdx().y;   
    n = (blockIdx().z - 1) * blockDim().z + threadIdx().z; 
    Lᵧ, M, N = size(Y₃);

    if l ≤ Lᵧ && m ≤ M && n ≤ N
        @inbounds if (l-1) % 4 != 0f0
            Y₃[l,m,n] = 0f0;
        else
            Y₃[l,m,n] = X[l,m,n] + y₃[l,m,n];
        end
    end
    return nothing
end
############################################################

############################################################
function update_Y₁!(Y₁, WᵀW, WᵀS, DX, y₁, rho)
    l = (blockIdx().x - 1) * blockDim().x + threadIdx().x;
    n = (blockIdx().y - 1) * blockDim().y + threadIdx().y;
    Lᵧ, N = size(Y₁);
    if l ≤ Lᵧ && n ≤ N
        @inbounds Y₁[l,n] = (WᵀS[l,n]+rho*(DX[l,n] + y₁[l,n])) / (WᵀW[l,l] + rho);
    end
    # as we can see that components outside of the support better be small
    return nothing
end
############################################################






