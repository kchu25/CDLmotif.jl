function code_update!(c::cdl, 
                      num_iter::Integer,
                      type_complex::DataType; 
                      sherman_const_csc=2f0)
    for _ = 1:num_iter
        ############################ proximal ℓ-1 update ########################################
        CUDA.@sync @cuda threads=ker_3d blocks = b_size_3d(c.Y_2) shrink!(c.Y_2, c.X, c.y_2, c.λ);       
        c.Ŷ_2 = c.fft_l_LᵧMN * c.Y_2;
        #########################################################################################

        ############################ lagrange multiplier update #################################
        c.y_1 .= c.y_1 .+ c.DX .- c.Y_1; c.ŷ_1 = c.fft_l_LᵧN * c.y_1; 
        c.y_2 .= c.y_2 .+ c.X .- c.Y_2; c.ŷ_2 = c.fft_l_LᵧMN * c.y_2;
        c.y_3 .= c.y_3 .+ c.X .- c.Y_3; c.ŷ_3 = c.fft_l_LᵧMN * c.y_3;
        #########################################################################################

        ############## Code update with the Sherman Morrison formula trick ######################
        c.dots_csc =  c.Ĝ_1 .* conj.(c.Ĝ_1);
        c.sum_dots_csc = reduce(+, c.dots_csc; dims=2) .+ sherman_const_csc;
        c.Ã_csc = c.Ĝ_1 ./ c.sum_dots_csc;  
        c.RHS_csc = CuArray{type_complex, 3}(undef, promote_i(c.Lᵧ, c.M, c.N));      
        c.alpha_csc = CuArray{type_complex, 2}(undef, promote_i(c.Lᵧ, c.N));       
        CUDA.@sync @cuda threads=ker_3d blocks=b_size_3d(c.RHS_csc) sherman_CSC_RHS!(c.RHS_csc, c.Ĝ_1, c.Ŷ_1, c.ŷ_1, c.Ŷ_2, c.ŷ_2, c.Ŷ_3, c.ŷ_3);                    
        CUDA.@sync @cuda threads=ker_2d blocks=b_size_2d(c.alpha_csc) sherman_CSC_dot!(c.alpha_csc, c.Ã_csc, c.RHS_csc);
        CUDA.@sync @cuda threads=ker_3d blocks=b_size_3d(c.X̂) sherman_CSC_combine!(c.X̂, c.alpha_csc, c.Ĝ_1, c.RHS_csc, sherman_const_csc);
        c.X = real.(c.ifft_l_LᵧMN * c.X̂); # iFFT
        c.DX = real.(c.ifft_l_LᵧN * dropdims(sum(c.Ĝ_1.* c.X̂, dims=2), dims=2)); 
        c.WDX = c.W * c.DX;

        # CUDA.unsafe_free!(dots_csc); CUDA.unsafe_free!(sum_dots_csc); CUDA.unsafe_free!(Ã_csc); 
        # CUDA.unsafe_free!(RHS_csc); CUDA.unsafe_free!(alpha_csc); 
        #########################################################################################

        ############################ constraint update ##########################################
        CUDA.@sync @cuda threads=ker_2d blocks=b_size_2d(c.Y_1) update_Y₁!(c.Y_1, c.WᵀW, c.WᵀS, c.DX, c.y_1, c.rho);
        c.Ŷ_1 = c.fft_l_LᵧN * c.Y_1;
        #########################################################################################

        ############################ project X to its spatial constraint ########################
        CUDA.@sync @cuda threads=ker_3d blocks=b_size_3d(c.Y_3) projectionX!(c.Y_3, c.X, c.y_3);
        c.Ŷ_3 = c.fft_l_LᵧMN * c.Y_3;
        #########################################################################################        
    end
end

function dict_update!(c::cdl, 
                      num_iter::Integer,
                      type_complex::DataType; 
                      sherman_const_cdl=2f0)
    for _ = 1:num_iter
        ##################### dict update using Sherman Morrison ################################
        c.dots_cdl = c.Ŷ_3 .* conj.(c.Ŷ_3);
        c.sum_dots_cdl = reduce(+, c.dots_cdl; dims=2) .+ sherman_const_cdl;
        c.X_cdl = c.Ŷ_3 ./ c.sum_dots_cdl;
        c.RHS_cdl = CUDA.zeros(type_complex, promote_i(c.Lᵧ, c.M, c.N));        
        c.ifft_elem_mult = c.ifft_l_LᵧN * dropdims(sum(c.Ŷ_3 .* c.d̂, dims=2), dims=2);
        c.Td = real.(c.ifft_elem_mult);
        CUDA.@sync @cuda threads=ker_3d blocks=b_size_3d(c.X̂) sherman_CDL_RHS!(c.RHS_cdl, c.Ŷ_3, c.Ĝ_1, c.ĝ_1, c.Ĝ_2, c.ĝ_2, c.Ĝ_3, c.ĝ_3);
        c.beta_cdl = dropdims(sum(c.X_cdl .* c.RHS_cdl; dims=2), dims=2);
        CUDA.@sync @cuda threads=ker_3d blocks=b_size_3d(c.d̂) sherman_CDL_combine(c.d̂, c.RHS_cdl, c.beta_cdl, c.Ŷ_3);
        c.d = real.(c.ifft_l_LᵧMN * c.d̂); # iFFT

        # CUDA.unsafe_free!(dots_cdl); CUDA.unsafe_free!(sum_dots_cdl); CUDA.unsafe_free!(X_cdl); 
        # CUDA.unsafe_free!(RHS_cdl); CUDA.unsafe_free!(ifft_elem_mult); CUDA.unsafe_free!(beta_cdl); 
        #########################################################################################

        ######################################################################################### 
        CUDA.@sync @cuda threads = ker_2d blocks = b_size_2d(c.G_2) update_g₂!(c.G_2, c.WᵀW, c.WᵀS, c.Td, c.g_2, c.sigma);
        c.Ĝ_2 = c.fft_l_LᵧN * c.G_2;
        c.g_1 = c.g_1 .+ c.d .- c.G_1; c.ĝ_1 = c.fft_l_LᵧMN * c.g_1;
        c.g_2 = c.g_2 .+ c.Td .- c.G_2; c.ĝ_2 = c.fft_l_LᵧN * c.g_2;
        c.g_3 = c.g_3 .+ c.d .- c.G_3; c.ĝ_3 = c.fft_l_LᵧMN * c.g_3;
        #########################################################################################

        ######################### projection to non-neg orthant#################################
        c.v_avg_2 = dropdims(sum(c.d .+ c.g_3, dims=3) ./ c.N, dims=3);
        CUDA.@sync @cuda threads = ker_2d blocks = b_size_2d(c.G_3) projection2box!(c.G_3, c.v_avg_2, c.P_D);        
        c.Ĝ_3 = c.fft_l_LᵧM * c.G_3;
        # CUDA.unsafe_free!(v_avg_2);  
        ########################################################################################

        ########################## proximal update ##############################################
        c.v_avg_1 = dropdims(sum(c.d .+ c.g_1, dims=3) ./ c.N, dims=3);
        CUDA.@sync @cuda threads = ker_2d blocks = b_size_2d(c.G_1) prox_simplex_zeros!(c.G_1, c.v_avg_1, c.P_D);
        c.Ĝ_1 = c.fft_l_LᵧM * c.G_1;
        # CUDA.unsafe_free!(v_avg_1); 
        #########################################################################################      
    end
end

function update_best!(c::cdl, dat)
    obj_val = objective(c, dat);
    if obj_val < c.obj_best
        c.obj_best = obj_val;        
        c.D_best = dropdims(sum(c.d, dims=3) ./ c.N, dims=3);
        c.X_best = copy(c.X);
    end
end
