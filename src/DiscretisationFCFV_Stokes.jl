# A = 0.9
A = 1.0

@doc """
Computes FCFV dicretisation elements `α`, `β` and `Ζ`. The discretisation is based on first order scheme (see Sevilla et al., 2018).
""" ComputeFCFV

function ComputeFCFV(mesh, sex, sey, VxDir, VyDir, SxxNeu, SyyNeu, SxyNeu, SyxNeu, τr, Formulation)

    α = zeros(mesh.nel)
    β = zeros(mesh.nel,2)
    Ζ = zeros(mesh.nel,2,2)

    # Loop through elements
     for e=1:mesh.nel 
        β[e,1] += mesh.Ω[e]*sex[e] # Source term
        β[e,2] += mesh.Ω[e]*sey[e] # Source term
        # Loop through faces
        for i=1:mesh.nf_el 
            nodei = mesh.e2f[e,i]
            bc    = mesh.bc[nodei]
            Γi    = mesh.Γ[e,i]
            ni_x  = mesh.n_x[e,i]
            ni_y  = mesh.n_y[e,i]
            τe    = mesh.τe[e]  # Stabilisation parameter for the face
            if Formulation==:Gradient
                Ζ[e,1,1] += (bc==1) * Γi*ni_x*VxDir[nodei] # Dirichlet
                Ζ[e,1,2] += (bc==1) * Γi*ni_x*VyDir[nodei] # Dirichlet
                Ζ[e,2,1] += (bc==1) * Γi*ni_y*VxDir[nodei] # Dirichlet
                Ζ[e,2,2] += (bc==1) * Γi*ni_y*VyDir[nodei] # Dirichlet
            elseif Formulation==:SymmetricGradient
                Ζ[e,1,1] += (bc==1) * Γi*(ni_x*VxDir[nodei] + A*ni_x*VxDir[nodei]) # Dirichlet
                Ζ[e,1,2] += (bc==1) * Γi*(ni_x*VyDir[nodei] + A*ni_y*VxDir[nodei]) # Dirichlet
                Ζ[e,2,1] += (bc==1) * Γi*(ni_y*VxDir[nodei] + A*ni_x*VyDir[nodei]) # Dirichlet
                Ζ[e,2,2] += (bc==1) * Γi*(ni_y*VyDir[nodei] + A*ni_y*VyDir[nodei]) # Dirichlet
            end
            β[e,1]       += (bc==1) * Γi*τe*VxDir[nodei]   # Dirichlet
            β[e,2]       += (bc==1) * Γi*τe*VyDir[nodei]   # Dirichlet
            α[e]         +=           Γi*τe
            mesh.τ[nodei] =  τe 
        end
    end
    return α, β, Ζ
end

#--------------------------------------------------------------------#

@doc """
Computes element values of velocity vector and deviatoric stress tensor components. 
""" ComputeElementValues

function ComputeElementValues(mesh, Vxh, Vyh, Pe, α, β, Ζ, VxDir, VyDir, Formulation)

    Vxe         = zeros(mesh.nel)
    Vye         = zeros(mesh.nel)
    τxxe        = zeros(mesh.nel)
    τyye        = zeros(mesh.nel)
    τxye        = zeros(mesh.nel) 
    D           = zeros(2,2)

     for e=1:mesh.nel
    
        η       =  mesh.ke[e]
        D       = [η  η/mesh.δ[e]; η/mesh.δ[e] η]
        Ω       =  mesh.Ω[e]
        Vxe[e]  =  β[e,1]/α[e]
        Vye[e]  =  β[e,2]/α[e]
        τxxe[e] =  η/Ω*Ζ[e,1,1]
        τyye[e] =  η/Ω*Ζ[e,2,2] 
        τxye[e] =  η/Ω*0.5*(Ζ[e,1,2] + Ζ[e,2,1])
        
        for i=1:mesh.nf_el
            
            # Face
            nodei = mesh.e2f[e,i]
            bc    = mesh.bc[nodei]
            Γi    = mesh.Γ[e,i]
            ni_x  = mesh.n_x[e,i]
            ni_y  = mesh.n_y[e,i]
            τi    = mesh.τ[nodei]  # Det face stabilisation equal to element stabilisation

            # Assemble
            Vxe[e]  += (bc!=1) *  Γi*τi*Vxh[nodei]/α[e]
            Vye[e]  += (bc!=1) *  Γi*τi*Vyh[nodei]/α[e]
            τxxe[e] += (bc!=1) *  D[1,1]/Ω*Γi*ni_x*Vxh[nodei]
            τyye[e] += (bc!=1) *  D[2,2]/Ω*Γi*ni_y*Vyh[nodei]
            τxye[e] += (bc!=1) *  D[1,2]*0.5*( 1.0/Ω*Γi*( ni_x*Vyh[nodei] + ni_y*Vxh[nodei] ) )
         end
        τxxe[e] *= 2.0
        τyye[e] *= 2.0
        τxye[e] *= 2.0
    end
    return Vxe, Vye, τxxe, τyye, τxye
end

#--------------------------------------------------------------------#

@doc """
Assembles linear system of equations 
""" ElementAssemblyLoop

function ElementAssemblyLoop(mesh, α, β, Ζ, VxDir, VyDir, σxxNeu, σyyNeu, σxyNeu, σyxNeu, Formulation) 

    # Assemble element matrices and rhs
    Kuui = zeros(2*mesh.nf_el, 2*mesh.nf_el, mesh.nel)
    Kuuj = zeros(2*mesh.nf_el, 2*mesh.nf_el, mesh.nel)
    Kuuv = zeros(2*mesh.nf_el, 2*mesh.nf_el, mesh.nel)
    Muuv = zeros(2*mesh.nf_el, 2*mesh.nf_el, mesh.nel)
    fu   = zeros(2*mesh.nf_el, mesh.nel)
    Kupi = zeros(2*mesh.nf_el, mesh.nel)
    Kupj = zeros(2*mesh.nf_el, mesh.nel)
    Kupv = zeros(2*mesh.nf_el, mesh.nel)
    fp   = zeros(mesh.nel)
    nf   = mesh.nf_el

     for e=1:mesh.nel 

        # Element properties
        Ωe = mesh.Ω[e]
        ηe = mesh.ke[e]

        for i=1:mesh.nf_el 

            ni_x, ni_y = mesh.n_x[e,i], mesh.n_y[e,i]
            nodei = mesh.e2f[e,i]
            bci   = mesh.bc[nodei]
            ȷ     = 0.0 + (bci==-1)*1.0 # indicates interface
            Γi    = mesh.Γ[e,i]
            τi    = mesh.τ[nodei]  
                
            for j=1:mesh.nf_el

                nj_x, nj_y  = mesh.n_x[e,j], mesh.n_y[e,j]
                nodej = mesh.e2f[e,j]
                bcj   = mesh.bc[nodej]   
                Γj    = mesh.Γ[e,j]
                τj    = mesh.τ[nodej]  
                δ     = 0.0 + (i==j)*1.0    # Delta operator
                on    = (bci!=1) & (bcj!=1) # Activate nodal connection if not Dirichlet!
                        
                # Element matrix components
                ninj = ni_x*nj_x + ni_y*nj_y

                # Element matrix 
                if Formulation==:Gradient
                    Kuuv[j   , i   , e] = on * -Γi * (α[e]^-1*τi*τj*Γj - ηe*Ωe^-1*Γj*(ninj + ȷ*ni_x*nj_x) - τi*δ) # u1u1
                    Kuuv[j+nf, i   , e] = on * -Γi * (                 - ηe*Ωe^-1*Γj*(       ȷ*ni_y*nj_x)       ) # u1u2
                    Kuuv[j   , i+nf, e] = on * -Γi * (                 - ηe*Ωe^-1*Γj*(       ȷ*ni_x*nj_y)       ) # u2u1
                    Kuuv[j+nf, i+nf, e] = on * -Γi * (α[e]^-1*τi*τj*Γj - ηe*Ωe^-1*Γj*(ninj + ȷ*ni_y*nj_y) - τi*δ) # u2u2
                elseif Formulation==:SymmetricGradient
                    Kuuv[j   , i   , e] = on * -Γi * (α[e]^-1*τi*τj*Γj - ηe*Ωe^-1*Γj*(ninj + A*ni_x*nj_x) - τi*δ) # u1u1
                    Kuuv[j+nf, i   , e] = on * -Γi * (                 - ηe*Ωe^-1*Γj*(       A*ni_y*nj_x)       ) # u1u2
                    Kuuv[j   , i+nf, e] = on * -Γi * (                 - ηe*Ωe^-1*Γj*(       A*ni_x*nj_y)       ) # u2u1
                    Kuuv[j+nf, i+nf, e] = on * -Γi * (α[e]^-1*τi*τj*Γj - ηe*Ωe^-1*Γj*(ninj + A*ni_y*nj_y) - τi*δ) # u2u2
                end
                # PC - deactivate terms from new interface implementation
                Muuv[j   , i   , e] = on * -Γi * (α[e]^-1*τi*τj*Γj - ηe*Ωe^-1*Γj*ninj - τi*δ) # u1u1
                Muuv[j+nf, i+nf, e] = on * -Γi * (α[e]^-1*τi*τj*Γj - ηe*Ωe^-1*Γj*ninj - τi*δ) # u2u2

                # Connectivity
                Kuui[j   , i   , e]  = nodei;         Kuui[j+nf, i   , e]  = nodei
                Kuuj[j   , i   , e]  = nodej;         Kuuj[j+nf, i   , e]  = nodej+mesh.nf
                Kuui[j   , i+nf, e]  = nodei+mesh.nf; Kuui[j+nf, i+nf, e]  = nodei+mesh.nf
                Kuuj[j   , i+nf, e]  = nodej;         Kuuj[j+nf, i+nf, e]  = nodej+mesh.nf
            end 
            # RHS
            Xi    = 0.0 + (bci==2)*1.0
            tix   = ni_x*σxxNeu[nodei] + ni_y*σxyNeu[nodei]
            tiy   = ni_x*σyxNeu[nodei] + ni_y*σyyNeu[nodei]   
            if Formulation==:Gradient
                niΖ_x = ni_x*(Ζ[e,1,1] +  ȷ*Ζ[e,1,1]) + ni_y*(Ζ[e,2,1] + ȷ*Ζ[e,1,2]) 
                niΖ_y = ni_x*(Ζ[e,1,2] +  ȷ*Ζ[e,2,1]) + ni_y*(Ζ[e,2,2] + ȷ*Ζ[e,2,2])
            elseif Formulation==:SymmetricGradient
                niΖ_x = ni_x*Ζ[e,1,1] + ni_y*Ζ[e,2,1] 
                niΖ_y = ni_x*Ζ[e,1,2] + ni_y*Ζ[e,2,2]     
            end
            feix  = (bci!=1) * -Γi * (-α[e]^-1*τi*β[e,1] + ηe*Ωe^-1*niΖ_x - tix*Xi)
            feiy  = (bci!=1) * -Γi * (-α[e]^-1*τi*β[e,2] + ηe*Ωe^-1*niΖ_y - tiy*Xi)
            # up block
            Kupv[i   , e] -= (bci!=1) * Γi*ni_x
            Kupv[i+nf, e] -= (bci!=1) * Γi*ni_y
            Kupi[i   , e]  = nodei
            Kupj[i   , e]  = e
            Kupi[i+nf, e]  = nodei + mesh.nf
            Kupj[i+nf, e]  = e
            # Dirichlet nodes - uu block
            Kuuv[i   , i   , e] += (bci==1) * 1e0
            Kuuv[i+nf, i+nf, e] += (bci==1) * 1e0
            Muuv[i   , i   , e] += (bci==1) * 1e0
            Muuv[i+nf, i+nf, e] += (bci==1) * 1e0
            fu[i   ,e]          += (bci!=1) * feix + (bci==1) * VxDir[nodei] * 1e0
            fu[i+nf,e]          += (bci!=1) * feiy + (bci==1) * VyDir[nodei] * 1e0
            # Dirichlet nodes - pressure RHS
            fp[e]               -= (bci==1) * Γi*(VxDir[nodei]*ni_x + VyDir[nodei]*ni_y)
        end
    end

    # Call sparse assembly
    tsparse = @elapsed Kuu, Muu, Kup, fu = Sparsify( Kuui, Kuuj, Kuuv, Muuv, Kupi, Kupj, Kupv, fu, mesh.nf, mesh.nel)

    return Kuu, Muu, Kup, fu, fp, tsparse
end

#--------------------------------------------------------------------#

@doc """
Assembles linear system of equations 
""" ElementAssemblyLoopNEW

function ElementAssemblyLoopNEW(mesh, α, β, Ζ, VxDir, VyDir, σxxNeu, σyyNeu, σxyNeu, σyxNeu, Formulation) 

    # Assemble element matrices and rhs
    Kuui = zeros(2*mesh.nf_el, 2*mesh.nf_el, mesh.nel)
    Kuuj = zeros(2*mesh.nf_el, 2*mesh.nf_el, mesh.nel)
    Kuuv = zeros(2*mesh.nf_el, 2*mesh.nf_el, mesh.nel)
    Muuv = zeros(2*mesh.nf_el, 2*mesh.nf_el, mesh.nel)
    fu   = zeros(2*mesh.nf_el, mesh.nel)
    Kupi = zeros(2*mesh.nf_el, mesh.nel)
    Kupj = zeros(2*mesh.nf_el, mesh.nel)
    Kupv = zeros(2*mesh.nf_el, mesh.nel)
    fp   = zeros(mesh.nel)
    nf   = mesh.nf_el
    D    = zeros(2,2)

     for e=1:mesh.nel 

        # Rheology
        # D1 = [1.0/mesh.λ[e] 1.0/mesh.δ[e]/mesh.λ[e]; 1.0/mesh.δ[e]/mesh.λ[e] 1.0/mesh.λ[e]]
        # D1 = [mesh.ke[e] mesh.ke[e]/mesh.δ[e]; mesh.ke[e]/mesh.δ[e] mesh.ke[e]] .* 1.0/mesh.Ω[e]

        if mesh.δ[e] ==1
            D .= [mesh.ke[e] mesh.ke[e]/mesh.δ[e]; mesh.ke[e]/mesh.δ[e] mesh.ke[e]] .* 1.0/mesh.Ω[e]
        else
            Dv = inv([ mesh.ke[e]  mesh.ke[e]/mesh.δ[e]; mesh.ke[e]/mesh.δ[e] mesh.ke[e]]) .* mesh.Ω[e]
            D .= inv(Dv)
        end

        # @show D1
        # @show D

        # Dv  = [mesh.ke[e] 0. 0.; 0. mesh.ke[e] 0.0; 0.0 0.0 mesh.ke[e]/ mesh.δ[e]]
        # Dv1 = inv( inv(Dv).*mesh.Ω[e]    )     
        # D   = [Dv1[1,1] Dv1[3,3]; Dv[3,3] Dv[2,2]] 
        # @show mesh.ke[e]
        # @show [ mesh.ke[e]  mesh.ke[e]/mesh.δ[e]; mesh.ke[e]/mesh.δ[e] mesh.ke[e]]



        for i=1:mesh.nf_el 

            ni_x, ni_y = mesh.n_x[e,i], mesh.n_y[e,i]
            nodei = mesh.e2f[e,i]
            bci   = mesh.bc[nodei]
            ȷ     = 0.0 + (bci==-1)*1.0 # indicates interface
            Γi    = mesh.Γ[e,i]
            τi    = mesh.τ[nodei]  
                
            for j=1:mesh.nf_el

                nj_x, nj_y  = mesh.n_x[e,j], mesh.n_y[e,j]
                nodej = mesh.e2f[e,j]
                bcj   = mesh.bc[nodej]   
                Γj    = mesh.Γ[e,j]
                τj    = mesh.τ[nodej]  
                δ     = i==j                # Delta operator
                on    = (bci!=1) & (bcj!=1) # Activate nodal connection if not Dirichlet!
                        
                # Element matrix components
                ninj = ni_x*nj_x + ni_y*nj_y

                # Element matrix 
                if Formulation==:Gradient
                    Kuuv[j   , i   , e] = on * -Γi * (α[e]^-1*τi*τj*Γj - D[1,1]*Γj*(ninj + ȷ*ni_x*nj_x) - τi*δ) # u1u1
                    Kuuv[j   , i+nf, e] = on * -Γi * (                 - D[1,2]*Γj*(       ȷ*ni_x*nj_y)       ) # u2u1
                    Kuuv[j+nf, i   , e] = on * -Γi * (                 - D[2,1]*Γj*(       ȷ*ni_y*nj_x)       ) # u1u2
                    Kuuv[j+nf, i+nf, e] = on * -Γi * (α[e]^-1*τi*τj*Γj - D[2,2]*Γj*(ninj + ȷ*ni_y*nj_y) - τi*δ) # u2u2
                    # Kuuv[j   , i   , e] = on * -Γi * (α[e]^-1*τi*τj*Γj - Γj*(D[1,1]*ninj  +   D[1,1]*A*ni_x*nj_x) - τi*δ) # u1u1
                    # Kuuv[j   , i+nf, e] = on * -Γi * (                 - Γj*(                 D[1,2]*A*ni_x*nj_y)       ) # u2u1
                    # Kuuv[j+nf, i   , e] = on * -Γi * (                 - Γj*(                 D[2,1]*A*ni_y*nj_x)       ) # u1u2
                    # Kuuv[j+nf, i+nf, e] = on * -Γi * (α[e]^-1*τi*τj*Γj - Γj*(D[2,2]*ninj  +   D[2,2]*A*ni_y*nj_y) - τi*δ) # u2u2
                elseif Formulation==:SymmetricGradient
                    Ωe = mesh.Ω[e]
                    # Kuuv[i   , j   , e] = on * -Γi * (α[e]^-1*τi*τj*Γj - D[1,1]*Γj*(ninj + A*ni_x*nj_x) - τi*δ) # u1u1
                    # Kuuv[i   , j+nf, e] = on * -Γi * (                 - D[1,2]*Γj*(       A*ni_y*nj_x)       ) # u1u2
                    # Kuuv[i+nf, j   , e] = on * -Γi * (                 - D[2,1]*Γj*(       A*ni_x*nj_y)       ) # u2u1
                    # Kuuv[i+nf, j+nf, e] = on * -Γi * (α[e]^-1*τi*τj*Γj - D[2,2]*Γj*(ninj + A*ni_y*nj_y) - τi*δ) # u2u2

                    Kuuv[i   , j   , e] = on * -Γi * (α[e]^-1*τi*τj*Γj - D[1,1]*Γj*(ninj + A*ni_x*nj_x) - τi*δ) # u1u1
                    Kuuv[i   , j+nf, e] = on * -Γi * (                 - D[1,2]*Γj*(       A*ni_y*nj_x)       ) # u1u2
                    Kuuv[i+nf, j   , e] = on * -Γi * (                 - D[2,1]*Γj*(       A*ni_x*nj_y)       ) # u2u1
                    Kuuv[i+nf, j+nf, e] = on * -Γi * (α[e]^-1*τi*τj*Γj - D[2,2]*Γj*(ninj + A*ni_y*nj_y) - τi*δ) # u2u2

                end
                # PC - deactivate terms from new interface implementation
                Muuv[j   , i   , e] = on * -Γi * (α[e]^-1*τi*τj*Γj - D[1,1]*Γj*ninj - τi*δ) # u1u1
                Muuv[j+nf, i+nf, e] = on * -Γi * (α[e]^-1*τi*τj*Γj - D[2,2]*Γj*ninj - τi*δ) # u2u2

                # Connectivity
                Kuui[i   , j   , e]  = nodei;         Kuui[i+nf, j   , e]  = nodei+mesh.nf
                Kuuj[i   , j   , e]  = nodej;         Kuuj[i+nf, j   , e]  = nodej
                Kuui[i   , j+nf, e]  = nodei;         Kuui[i+nf, j+nf, e]  = nodei+mesh.nf
                Kuuj[i   , j+nf, e]  = nodej+mesh.nf; Kuuj[i+nf, j+nf, e]  = nodej+mesh.nf

                # Kuui[j   , i   , e]  = nodei;         Kuui[j+nf, i   , e]  = nodei
                # Kuuj[j   , i   , e]  = nodej;         Kuuj[j+nf, i   , e]  = nodej+mesh.nf
                # Kuui[j   , i+nf, e]  = nodei+mesh.nf; Kuui[j+nf, i+nf, e]  = nodei+mesh.nf
                # Kuuj[j   , i+nf, e]  = nodej;         Kuuj[j+nf, i+nf, e]  = nodej+mesh.nf
            end 
            # RHS
            Xi    = 0.0 + (bci==2)*1.0
            tix   = ni_x*σxxNeu[nodei] + ni_y*σxyNeu[nodei]
            tiy   = ni_x*σyxNeu[nodei] + ni_y*σyyNeu[nodei]   
            if Formulation==:Gradient
                niΖ_x = ni_x*(Ζ[e,1,1] +  ȷ*Ζ[e,1,1]) + ni_y*(Ζ[e,2,1] + ȷ*Ζ[e,1,2]) 
                niΖ_y = ni_x*(Ζ[e,1,2] +  ȷ*Ζ[e,2,1]) + ni_y*(Ζ[e,2,2] + ȷ*Ζ[e,2,2])
            elseif Formulation==:SymmetricGradient
                niΖ_x = ni_x*Ζ[e,1,1] + ni_y*Ζ[e,2,1] 
                niΖ_y = ni_x*Ζ[e,1,2] + ni_y*Ζ[e,2,2]     
            end
            # feix  = (bci!=1) * -Γi * (-α[e]^-1*τi*β[e,1] + D[1,1] *niΖ_x - tix*Xi)
            # feiy  = (bci!=1) * -Γi * (-α[e]^-1*τi*β[e,2] + D[2,2] *niΖ_y - tiy*Xi)
            feix  = (bci!=1) * -Γi * (-α[e]^-1*τi*β[e,1] + (D[1,1]*ni_x*Ζ[e,1,1] + D[1,1]*ni_y*Ζ[e,2,1])  - tix*Xi)
            feiy  = (bci!=1) * -Γi * (-α[e]^-1*τi*β[e,2] + (D[2,2]*ni_x*Ζ[e,1,2] + D[2,2]*ni_y*Ζ[e,2,2])  - tiy*Xi)
            # up block
            Kupv[i   , e] -= (bci!=1) * Γi*ni_x
            Kupv[i+nf, e] -= (bci!=1) * Γi*ni_y
            Kupi[i   , e]  = nodei
            Kupj[i   , e]  = e
            Kupi[i+nf, e]  = nodei + mesh.nf
            Kupj[i+nf, e]  = e
            # Dirichlet nodes - uu block
            Kuuv[i   , i   , e] += (bci==1) * 1e0
            Kuuv[i+nf, i+nf, e] += (bci==1) * 1e0
            Muuv[i   , i   , e] += (bci==1) * 1e0
            Muuv[i+nf, i+nf, e] += (bci==1) * 1e0
            fu[i   ,e]          += (bci!=1) * feix + (bci==1) * VxDir[nodei] * 1e0
            fu[i+nf,e]          += (bci!=1) * feiy + (bci==1) * VyDir[nodei] * 1e0
            # Dirichlet nodes - pressure RHS
            fp[e]               -= (bci==1) * Γi*(VxDir[nodei]*ni_x + VyDir[nodei]*ni_y)
        end
    end

    # Call sparse assembly
    tsparse = @elapsed Kuu, Muu, Kup, fu = Sparsify( Kuui, Kuuj, Kuuv, Muuv, Kupi, Kupj, Kupv, fu, mesh.nf, mesh.nel)

    return Kuu, Muu, Kup, fu, fp, tsparse
end

#--------------------------------------------------------------------#

@doc """
Generate sparse matrix blocks from triplets given in COO format.
""" Sparsify

function Sparsify( Kuui, Kuuj, Kuuv, Muuv, Kupi, Kupj, Kupv, fuv, nf, nel)

    _one = ones(size(Kupi[:]))
    Kuu  =       dropzeros(sparse(Kuui[:], Kuuj[:], Kuuv[:], nf*2, nf*2))
    Muu  =       dropzeros(sparse(Kuui[:], Kuuj[:], Muuv[:], nf*2, nf*2))
    Kup  =       dropzeros(sparse(Kupi[:], Kupj[:], Kupv[:], nf*2, nel ))
    fu   = Array(dropzeros(sparse(Kupi[:],    _one,  fuv[:], nf*2,   1 )))

    return Kuu, Muu, Kup, fu
end










function ElementAssemblyLoopToy(mesh, α, β, Ζ, VxDir, VyDir, σxxNeu, σyyNeu, σxyNeu, σyxNeu, Formulation) 

    # Assemble element matrices and rhs
    Kuui = zeros(2*mesh.nf_el, 2*mesh.nf_el, mesh.nel)
    Kuuj = zeros(2*mesh.nf_el, 2*mesh.nf_el, mesh.nel)
    Kuuv = zeros(2*mesh.nf_el, 2*mesh.nf_el, mesh.nel)
    Muuv = zeros(2*mesh.nf_el, 2*mesh.nf_el, mesh.nel)
    fu   = zeros(2*mesh.nf_el, mesh.nel)
    Kupi = zeros(2*mesh.nf_el, mesh.nel)
    Kupj = zeros(2*mesh.nf_el, mesh.nel)
    Kupv = zeros(2*mesh.nf_el, mesh.nel)
    fp   = zeros(mesh.nel)
    nf   = mesh.nf_el

    nfac = mesh.nf
    nel = mesh.nel

    Kuu = spzeros(2nfac, 2nfac)
    Kup = spzeros(2nfac, nel)
    Muu = spzeros(2nfac, 2nfac)
    fu  = zeros(2nfac)
    fp  = zeros(nel)

    ni  = zeros(2)
    nj  = zeros(2)
    I2d =  sparse(I,2,2)

    for e=1:mesh.nel 

        Ωe = mesh.Ω[e]
        fp_loc = 0.

        for i=1:nf

            ni[1] = mesh.n_x[e,i]
            ni[2] = mesh.n_y[e,i]
            Γi    = mesh.Γ[e,i]
            nodei = mesh.e2f[e,i]
            τi    = mesh.τ[nodei]
            bci   = mesh.bc[nodei] 
            ȷ     = 0.0 + (bci==-1)*1.0 # indicates interface
            fp_loc = 0.

            for j=1:nf

                nj[1] = mesh.n_x[e,j]
                nj[2] = mesh.n_y[e,j]
                nodej = mesh.e2f[e,j]
                bcj   = mesh.bc[nodej]   
                Γj    = mesh.Γ[e,j]
                τj    = mesh.τ[nodej]  
                δ     = 0.0 + (i==j)*1.0    # Delta operator
                on    = (bci!=1) & (bcj!=1) # Activate nodal connection if not Dirichlet!

                if Formulation==:Gradient
                    Kuu_loc = - on * Γi * (α[e]^-1*τi*τj*Γj*I2d .- 1.0/mesh.λ[e]*Γj*(ni'*nj*I2d .+ ȷ*nj*ni') .- τi*δ*I2d )
                elseif Formulation==:SymmetricGradient
                    Kuu_loc = - on * Γi * (α[e]^-1*τi*τj*Γj*I2d .- 1.0/mesh.λ[e]*Γj*(ni'*nj*I2d .+ A*nj*ni') .- τi*δ*I2d )
                end
                
                vBC     = [VxDir[nodej]; VyDir[nodej]]
                fp_loc -= (bcj==1) * Γj*nj'*vBC 

                Kuu[nodei,      nodej     ] += Kuu_loc[1,1]
                Kuu[nodei,      nodej+nfac] += Kuu_loc[1,2]
                Kuu[nodei+nfac, nodej     ] += Kuu_loc[2,1]
                Kuu[nodei+nfac, nodej+nfac] += Kuu_loc[2,2]
            end

            if  (bci==1) 
                Kuu[nodei,      nodei     ] += 1.
                Kuu[nodei+nfac, nodei+nfac] += 1.
            end

            Xi                 = 0.0 + (bci==2)*1.0
            σBC                = [ σxxNeu[nodei] σxyNeu[nodei]; σyxNeu[nodei] σyyNeu[nodei] ]
 
            Kup_loc            = (bci!=1) * -Γi.*ni
            if Formulation==:Gradient
                fu_loc             = -Γi*(-α[e]^-1*τi*β[e,:]' .+ 1.0/mesh.λ[e] * ni'*(Ζ[e,:,:] + ȷ*Ζ[e,:,:]')  - Xi*ni'*σBC)
            elseif Formulation==:SymmetricGradient
                fu_loc             = -Γi*(-α[e]^-1*τi*β[e,:]' .+ 1.0/mesh.λ[e] * ni'*Ζ[e,:,:]  - Xi*ni'*σBC)
            end
            fu[nodei]         += (bci!=1) * fu_loc[1] + (bci==1)*VxDir[nodei]
            fu[nodei+nfac]    += (bci!=1) * fu_loc[2] + (bci==1)*VyDir[nodei]
            Kup[nodei     ,e] += Kup_loc[1]
            Kup[nodei+nfac,e] += Kup_loc[2]

        end
        fp[e] = fp_loc
    end
    t = 0.
    return Kuu, Muu, Kup, fu, fp, t
end