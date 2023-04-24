@doc """
Solve discretised Stokes equation using either a coupled sparse-direct factorisation or Powell-Hestenes iteration.
Either LU or Cholesky factorisation (for symmetric positive definite case) can be within Powell-Hestenes iterations. 
""" StokesSolvers!

function StokesSolvers!(Vxh, Vyh, Pe, mesh, Kuu, Kup, fu, fp, M, solver; Kpu=-Kup', Kpp=spdiagm(length(fp), length(fp)), penalty=1e4, tol=1e-10)

    @printf("Start solver %s for %d DOFs\n", solver, length(fu)+length(fp))
    nVx = Int64(length(fu)/2)
    nVy = length(fu)
    
    if solver==:CoupledBackslash
        # Coupled solve
        K      = [Kuu Kup; Kpu Kpp]
        f      = [fu; fp]
        xh     = K\f
        Vxh   .= xh[1:nVx]
        Vyh   .= xh[nVx+1:nVy]
        Pe    .= xh[nVy+1:end] 
        Pe    .= Pe .- mean(Pe)
    elseif solver==:PowellHestenesCholesky || solver==:PowellHestenesLU
        # Decoupled solve
        coef  = penalty.*mesh.ke./mesh.Ω        
        Kppi  = spdiagm(coef)
        Kuusc = Kuu .- Kup*(Kppi*Kpu)
        if solver==:PowellHestenesCholesky 
            t = @elapsed Kf    = cholesky(Hermitian(Kuusc),check = false)
            @printf("Cholesky took = %02.2e s\n", t)
        elseif solver==:PowellHestenesLU
            t = @elapsed Kf    = lu(Kuusc)
            @printf("LU took = %02.2e s\n", t)
        end
        u     = zeros(length(fu), 1)
        ru    = zeros(length(fu), 1)
        fusc  = zeros(length(fu), 1)
        p     = zeros(length(fp), 1)
        rp    = zeros(length(fp), 1)
        # Iterations
        for rit=1:20
            ru        .= fu .- Kuu*u .- Kup*p
            rp        .= fp .- Kpu*u
            nrmu, nrmp = norm(ru), norm(rp)
            @printf("  --> Powell-Hestenes Iteration %02d\n  Momentum res.   = %2.2e\n  Continuity res. = %2.2e\n", rit, nrmu/sqrt(length(ru)), nrmp/sqrt(length(rp)))
            if nrmu/sqrt(length(ru)) < tol && nrmp/sqrt(length(ru)) < tol
                break
            end
            fusc .= fu  .- Kup*(Kppi*fp .+ p)
            u    .= Kf\fusc
            p   .+= Kppi*(fp .- Kpu*u)
        end
        # Post-process solve
        Vxh .= u[1:nVx]
        Vyh .= u[nVx+1:nVy]
        Pe  .= p[:]
    end
end