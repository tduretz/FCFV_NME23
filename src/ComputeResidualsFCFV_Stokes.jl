@doc """
Computes residual of local and global equations discretised with FCFV.
""" ComputeResidualsFCFV_Stokes_o1

function ComputeResidualsFCFV_Stokes_o1(Vxh, Vyh, Pe, mesh, ae, be, ze, sex, sey, VxDir, VyDir, SxxNeu, SyyNeu, SxyNeu, SyxNeu, Formulation)

    # reconstruct element value and flux
    nel = mesh.nel
    ue  = zeros(nel,2)
    Le  = zeros(nel,2,2)

    for e=1:nel    
        nfac      = mesh.nf_el
        vole      = mesh.Ω[e]
        ue[e,:]   = be[e,:]/ae[e]
        Le[e,:,:] = -1.0/vole * ze[e,:,:]
        τi        = mesh.τ[e]
        for i=1:nfac
            nodei = mesh.e2f[e,i]
            if mesh.bc[nodei] != 1
                Γi        = mesh.Γ[e,i]
                n         = [mesh.n_x[e,i]; mesh.n_y[e,i]]
                u         = [Vxh[nodei]; Vyh[nodei]]
                ue[e,:]   = ue[e,:]   .+ Γi*τi*u'[:]/ae[e]
                if Formulation == :Gradient
                    Le[e,:,:] = Le[e,:,:] .- 1.0/vole * Γi * (n*u')
                elseif Formulation == :SymmetricGradient
                    Le[e,:,:] = Le[e,:,:] .- 1.0/vole * Γi * (n*u' .+ u*n')
                end
            end
        end
    end

    # Check out the 3 residuals: 3 for local equations and 2 for the global
    F_glob1  = zeros(nel,mesh.nf_el,2);
    F_glob2  = zeros(nel);
    F_eq1    = zeros(nel,2,2);
    F_eq2    = zeros(nel,2);
    F_eq3    = zeros(nel,1);
    D        = zeros(2,2)

    for e=1:nel
        nfac         = mesh.nf_el
        vole         = mesh.Ω[e]
        η            = mesh.ke[e]
        D           .= [η  η/mesh.δ[e]; η/mesh.δ[e] η]
        F_eq1[e,:,:] = vole*Le[e,:,:]
        F_eq2[e,1]   = -sex[e]*vole
        F_eq2[e,2]   = -sey[e]*vole
        F_eq3[e]     = 0.0
        for i=1:nfac
            nodei  = mesh.e2f[e,i]
            dAi    = mesh.Γ[e,i]
            taui   = mesh.τ[nodei]  
            n      = [mesh.n_x[e,i]; mesh.n_y[e,i]]
            bci    = mesh.bc[nodei]
            Xi     = 0.0 + (bci== 2)*1.0
            Ji     = 0.0 + (bci==-1)*1.0
            tix    = n[1]*SxxNeu[nodei] + n[2]*SxyNeu[nodei]
            tiy    = n[1]*SyxNeu[nodei] + n[2]*SyyNeu[nodei] 
            ti     = [tix; tiy]
            if mesh.bc[nodei] == 1
                u = [VxDir[nodei]; VyDir[nodei]]
            else
                u = [Vxh[nodei]; Vyh[nodei]]
            end
            # Global residual 1 (momentum) 
            if Formulation == :Gradient
                F_glob1[e,i,:] .+=  dAi .* ( (n'*(D.*Le[e,:,:])) .+ Pe[e]*n' .+ taui*ue[e,:]' .- taui*u' .+ ti'*Xi .+ Ji*(n'*(D.*Le[e,:,:])') )'
            elseif Formulation == :SymmetricGradient
                F_glob1[e,i,:] .+=  dAi .* ( (n'*(D.*Le[e,:,:])) .+ Pe[e]*n' .+ taui*ue[e,:]' .- taui*u' .+ ti'*Xi  )'
            end
            F_eq2[e,:]       = F_eq2[e,:] + dAi*taui*ue[e,:]
            # Global residual 2 (continuity)
            F_glob2[e]   = F_glob2[e] .+ dAi * u'*n
            if Formulation == :Gradient
                F_eq1[e,:,:] = F_eq1[e,:,:] .+ dAi * n * u'
            elseif Formulation == :SymmetricGradient
                F_eq1[e,:,:] = F_eq1[e,:,:] .+ dAi * (n*u' .+ u*n')
            end
            F_eq2[e,:] = F_eq2[e,:] .- dAi*taui*u[:]
            F_eq3[e,:] = F_eq3[e,:] .+ dAi*u'*n
        end
    end

    # Loop through elements and add global residual contribution of eac
    # individual faces to global nodal residual
    nnod = mesh.nf
    F_nodes_x = zeros(nnod,1);
    F_nodes_y = zeros(nnod,1);
    for e=1:nel
        nfac = mesh.nf_el
        for i=1:nfac
            nodei = mesh.e2f[e,i]
                if mesh.bc[nodei] != 1
                F_nodes_x[nodei] += F_glob1[e,i,1] 
                F_nodes_y[nodei] += F_glob1[e,i,2] 
            end
        end
    end
    @printf("Residual of local  equation 20a: %2.2e\n", norm(F_eq1)/sqrt(length(F_eq1)))
    @printf("Residual of local  equation 20b: %2.2e\n", norm(F_eq2)/sqrt(length(F_eq2)))
    @printf("Residual of local  equation 20c: %2.2e\n", norm(F_eq3)/sqrt(length(F_eq3)))
    @printf("Residual of global equation 19a: %2.2e\n", norm(F_nodes_x)/sqrt(length(F_nodes_x)))
    @printf("Residual of global equation 19a: %2.2e\n", norm(F_nodes_y)/sqrt(length(F_nodes_y)))
    @printf("Residual of global equation 19b: %2.2e\n", norm(F_glob2)/sqrt(length(F_glob2)))
    return
end