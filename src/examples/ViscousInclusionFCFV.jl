using FCFV_NME23, Printf

function SetUpProblem!(mesh, P, Vx, Vy, Sxx, Syy, Sxy, VxDir, VyDir, SxxNeu, SyyNeu, SxyNeu, SyxNeu, sx, sy, R, η, Formulation)
    ηm, ηc = η[1], η[2]
    # Evaluate analytic solution for boundary data
    for in=1:mesh.nf
        # Face midpoint coordinates
        x         = mesh.xf[in]
        y         = mesh.yf[in]
        # Dirichlet data
        vx, vy, p, σxx, σyy, σxy = EvalAnalDani( x, y, R, ηm, ηc )
        VxDir[in] = vx
        VyDir[in] = vy
        # Neumann data
        p, ∂vx∂x, ∂vx∂y, ∂vy∂x, ∂vy∂y = Tractions( x, y, R, ηm, ηc, 1 )
        if Formulation==:Gradient 
            SxxNeu[in] = - p + ηm*∂vx∂x 
            SyyNeu[in] = - p + ηm*∂vy∂y 
            SxyNeu[in] =       ηm*∂vx∂y
            SyxNeu[in] =       ηm*∂vy∂x
        elseif Formulation==:SymmetricGradient
            # Stress at faces - Tractions / Symmetric Cauchy stress tensor
            SxxNeu[in] = - p + 2ηm*∂vx∂x 
            SyyNeu[in] = - p + 2ηm*∂vy∂y 
            SxyNeu[in] =        ηm*(∂vx∂y + ∂vy∂x)
            SyxNeu[in] =        ηm*(∂vx∂y + ∂vy∂x)
        end
    end
    # Evaluate analytic solution on element barycentres
    for iel=1:mesh.nel
        x            = mesh.xc[iel]
        y            = mesh.yc[iel]
        vx, vy, pre, σxx, σyy, σxy = EvalAnalDani( x, y, R, ηm, ηc )
        P[iel]       = pre
        Vx[iel]      = vx
        Vy[iel]      = vy
        Sxx[iel]     = σxx
        Syy[iel]     = σyy
        Sxy[iel]     = σxy
        sx[iel]      = 0.0
        sy[iel]      = 0.0
        out          = mesh.phase[iel] == 1.0
        mesh.ke[iel] = (out==1) * 1.0*ηm + (out!=1) * 1.0*ηc         
    end
    return
end

function ViscousInclusion()

    @printf("Viscous inclusion test using first order FCFV discretisation on triangles")

    # Physics
    xmin, xmax = -3.0, 3.0    # Domain extent x
    ymin, ymax = -3.0, 3.0    # Domain extent y
    R          = 1.0          # Inclusion radius
    η          = [1.0 1e-2]   # Viscosity matrix/inclusion
    BC         = [2; 1; 1; 1] # South/East/North/West --- 1: Dirichlet / 2: Neumann

    # Numerics
    Mesher      = :AdvancingFront 
    τr          = 10
    Mesher      = :Delaunay
    τr          = 3
    solver      = 1
    Formulation = :SymmetricGradient
    # Formulation = :Gradient

    # Generate mesh 
    nx, ny = 60, 60
    if Mesher==:Delaunay       mesh = MakeTriangleMesh( nx, ny, xmin, xmax, ymin, ymax, τr, 1, R, BC, ((xmax-xmin)/nx)*((ymax-ymin)/ny), 200 ) end
    if Mesher==:AdvancingFront mesh = LoadExternalMesh(:MR, η) end
    @printf("Mesh informations:\n")
    @printf("--> %d elements\n", mesh.nel)
    @printf("--> %d faces\n",    mesh.nf)
    @printf("--> %d matrix elements\n",    sum(mesh.phase.==1))
    @printf("--> %d inclusion elements\n", sum(mesh.phase.==2))
    
    # Source term and BCs etc...
    Vxh    = zeros(mesh.nf)
    Vyh    = zeros(mesh.nf)
    Pe     = zeros(mesh.nel)
    Pa     = zeros(mesh.nel)
    Vxa    = zeros(mesh.nel)
    Vya    = zeros(mesh.nel)
    Sxxa   = zeros(mesh.nel)
    Syya   = zeros(mesh.nel)
    Sxya   = zeros(mesh.nel)
    sex    = zeros(mesh.nel)
    sey    = zeros(mesh.nel)
    VxDir  = zeros(mesh.nf)
    VyDir  = zeros(mesh.nf)
    SxxNeu = zeros(mesh.nf)
    SyyNeu = zeros(mesh.nf)
    SxyNeu = zeros(mesh.nf)
    SyxNeu = zeros(mesh.nf)

    # Set up Stokes problem
    @printf("---> Model configuration :\n")
    @time SetUpProblem!(mesh, Pa, Vxa, Vya, Sxxa, Syya, Sxya, VxDir, VyDir, SxxNeu, SyyNeu, SxyNeu, SyxNeu, sex, sey, R, η, Formulation)

    # Compute mesh properties for FCFV
    @printf("---> Compute FCFV vectors:\n")
    mesh.τ = τr.*ones(mesh.nf)  # Stabilisation per element
    @time ae, be, ze = ComputeFCFV(mesh, sex, sey, VxDir, VyDir, SxxNeu, SyyNeu, SxyNeu, SyxNeu, τr, Formulation)
    
    # Assemble element matrices and RHS
    @printf("---> Compute element matrices:\n")                      
    @time Kuu, Muu, Kup, fu, fp, tsparse = ElementAssemblyLoop(mesh, ae, be, ze, VxDir, VyDir, SxxNeu, SyyNeu, SxyNeu, SyxNeu, Formulation)
    println("---> Sparsification: ", tsparse)

    # # Solve for hybrid variables
    @printf("---> Linear solve:\n")
    Pe .= Pa
    @time StokesSolvers!(Vxh, Vyh, Pe, mesh, Kuu, Kup, fu, fp, Muu, solver; tol=5e-8, penalty=1e4)

    # Reconstruct element values
    @printf("---> Compute element values:\n")
    @time Vxe, Vye, Txxe, Tyye, Txye = ComputeElementValues(mesh, Vxh, Vyh, Pe, ae, be, ze, VxDir, VyDir, Formulation)
    
    # Evaluate residuals
    @printf("---> Evaluate residuals:\n")
    ComputeResidualsFCFV_Stokes_o1(Vxh, Vyh, Pe, mesh, ae, be, ze, sex, sey, VxDir, VyDir, SxxNeu, SyyNeu, SxyNeu, SyxNeu, Formulation)

    # Visualise
    @printf("---> Visualisation:\n")
    @time PlotMakie( mesh, Pe,  xmin, xmax, ymin, ymax; cmap=:jet1, min_v =-3, max_v =3 )

end

ViscousInclusion()