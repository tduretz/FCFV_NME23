A = 1.0

@doc """
Computes FCFV dicretisation elements `α`, `β` and `Ζ`. The discretisation is based on first order scheme (see Sevilla et al., 2018).
""" ComputeFCFV

function ComputeFCFV(mesh, sex, sey, VxDir, VyDir, SxxNeu, SyyNeu, SxyNeu, SyxNeu, τr, Formulation)

    α = zeros(mesh.nel)
    β = zeros(mesh.nel,2)
    Ζ = zeros(mesh.nel,2,2)

    # Loop through elements
    @inbounds for e=1:mesh.nel 
        β[e,1] += mesh.Ω[e]*sex[e] # Source term
        β[e,2] += mesh.Ω[e]*sey[e] # Source term
        # Loop through faces
        for i=1:mesh.nf_el 
            nodei = mesh.e2f[e,i]
            bc    = mesh.bc[nodei]
            Γi    = mesh.Γ[e,i]
            ni_x  = mesh.n_x[e,i]
            ni_y  = mesh.n_y[e,i]
            τi    = τr#*mesh.ke[e]  # Stabilisation parameter for the face
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
            β[e,1]       += (bc==1) * Γi*τi*VxDir[nodei]   # Dirichlet
            β[e,2]       += (bc==1) * Γi*τi*VyDir[nodei]   # Dirichlet
            α[e]         +=           Γi*τi
            mesh.τ[nodei] = τi
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

    @inbounds for e=1:mesh.nel
    
        η       =  mesh.ke[e]
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
            τi    = mesh.τ[nodei]  # Stabilisation parameter for the face

            # Assemble
            Vxe[e]  += (bc!=1) *  Γi*τi*Vxh[nodei]/α[e]
            Vye[e]  += (bc!=1) *  Γi*τi*Vyh[nodei]/α[e]
            τxxe[e] += (bc!=1) *  η/Ω*Γi*ni_x*Vxh[nodei]
            τyye[e] += (bc!=1) *  η/Ω*Γi*ni_y*Vyh[nodei]
            τxye[e] += (bc!=1) *  η*0.5*( 1.0/Ω*Γi*( ni_x*Vyh[nodei] + ni_y*Vxh[nodei] ) )
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

    @inbounds for e=1:mesh.nel 

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