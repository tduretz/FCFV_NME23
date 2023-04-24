using FCFV_NME23, Printf

#--------------------------------------------------------------------#

@doc """
Generates model configuration (viscosity field) and set up boundary values.
""" SetUpProblem!

function SetUpProblem!(mesh, Pa, Vxa, Vya, σxxa, σyya, σxya, VxDir, VyDir, σxxNeu, σyyNeu, σxyNeu, σyxNeu, sx, sy, R, η, Formulation)
    ηm, ηc = η[1], η[2]
    # Evaluate analytic solution for boundary data
    for in=1:mesh.nf
        # Face midpoint coordinates
        x,  y        = mesh.xf[in], mesh.yf[in]
        # Dirichlet data
        vx, vy, p, σxx, σyy, σxy = EvalAnalDani( x, y, R, ηm, ηc )
        VxDir[in], VyDir[in] = vx, vy
        # Neumann data
        p, ∂vx∂x, ∂vx∂y, ∂vy∂x, ∂vy∂y = Tractions( x, y, R, ηm, ηc, 1 )
        if Formulation==:Gradient 
            σxxNeu[in] = - p + ηm*∂vx∂x 
            σyyNeu[in] = - p + ηm*∂vy∂y 
            σxyNeu[in] =       ηm*∂vx∂y
            σyxNeu[in] =       ηm*∂vy∂x
        elseif Formulation==:SymmetricGradient
            # Stress at faces - Tractions / Symmetric Cauchy stress tensor
            σxxNeu[in] = - p + 2ηm*∂vx∂x 
            σyyNeu[in] = - p + 2ηm*∂vy∂y 
            σxyNeu[in] =        ηm*(∂vx∂y + ∂vy∂x)
            σyxNeu[in] =        ηm*(∂vx∂y + ∂vy∂x)
        end
    end
    # Evaluate analytic solution on element barycentres
    for iel=1:mesh.nel
        x                               = mesh.xc[iel]
        y                               = mesh.yc[iel]
        vx, vy, pre, σxx, σyy, σxy      = EvalAnalDani( x, y, R, ηm, ηc )
        Pa[iel]                         = pre
        Vxa[iel], Vya[iel]              = vx, vy
        σxxa[iel], σyya[iel], σxya[iel] = σxx, σyy, σxy
        sx[iel],  sy[iel]               = 0.0, 0.0
        out                             = mesh.phase[iel] == 1.0
        mesh.ke[iel]                    = (out==1) * 1.0*ηm + (out!=1) * 1.0*ηc         
    end
    return
end

#--------------------------------------------------------------------#

function ViscousInclusion()

    @printf("Viscous inclusion test using first order FCFV discretisation on triangles")

    # Physics
    xmin, xmax  = -3.0, 3.0    # Domain extent x
    ymin, ymax  = -3.0, 3.0    # Domain extent y
    R           = 1.0          # Inclusion radius
    η           = [1.0 1e-2]   # Viscosity matrix/inclusion
    BC          = [2; 1; 1; 1] # South/East/North/West --- 1: Dirichlet / 2: Neumann

    # Numerics
    Mesher      = :Delaunay                # :Delaunay / :AdvancingFront (load external mesh)
    mesh_res    = :MedRes                  # :LowRes / :MedRes / :HighRes     
    solver      = :PowellHestenesCholesky  # :CoupledBackslash / :=PowellHestenesCholesky / :=PowellHestenesLU
    Formulation = :SymmetricGradient       # :Gradient / :SymmetricGradient
    τr          = 2                        # Stabilisation
    γ           = 1e5                      # Penalty factor for Powell-Hestenes solvers
    ϵ           = 1e-8                     # Tolerance of Powell-Hestenes solvers 

    # Generate mesh 
    nx, ny = 60, 60  # initial point density in x and y for triangulation 
    if Mesher==:Delaunay       mesh = MakeTriangleMesh( nx, ny, xmin, xmax, ymin, ymax, τr, 1, R, BC, ((xmax-xmin)/nx)*((ymax-ymin)/ny), 200 ) end
    if Mesher==:AdvancingFront mesh = LoadExternalMesh( mesh_res, η) end
    @printf("Mesh informations:\n")
    @printf("--> %d elements\n", mesh.nel)
    @printf("--> %d faces\n",    mesh.nf)
    @printf("--> %d matrix elements\n",    sum(mesh.phase.==1))
    @printf("--> %d inclusion elements\n", sum(mesh.phase.==2))
    if Formulation==:Gradient && solver==:PowellHestenesCholesky
        error("Gradient formulation leads to non symmetric system, can't apply Cholesky factorisation!")
    end
    
    # Source term and BCs etc...
    Vxh    = zeros(mesh.nf)   # Velocity x (hybrid variable)
    Vyh    = zeros(mesh.nf)   # Velocity y (hybrid variable)
    Pe     = zeros(mesh.nel)  # Element pressure
    σxxa   = zeros(mesh.nel)  # Total stress xx
    σyya   = zeros(mesh.nel)  # Total stress yy
    σxya   = zeros(mesh.nel)  # Total stress xy
    sex    = zeros(mesh.nel)  # Source term momentum x 
    sey    = zeros(mesh.nel)  # Source term momentum y
    VxDir  = zeros(mesh.nf)   # Dirichlet data: Face velocity x
    VyDir  = zeros(mesh.nf)   # Dirichlet data: Face velocity y
    σxxNeu = zeros(mesh.nf)   # Neuman data: Total stress xx
    σyyNeu = zeros(mesh.nf)   # Neuman data: Total stress yy
    σxyNeu = zeros(mesh.nf)   # Neuman data: Total stress xy
    σyxNeu = zeros(mesh.nf)   # Neuman data: Total stress yx

    # Analytical solution    
    Pa     = zeros(mesh.nel)  # Element pressure
    Vxa    = zeros(mesh.nel)  # Element horizontal velocity
    Vya    = zeros(mesh.nel)  # Element vertical velocity

    # Set up Stokes problem
    @printf("---> Model configuration :\n")
    @time SetUpProblem!(mesh, Pa, Vxa, Vya, σxxa, σyya, σxya, VxDir, VyDir, σxxNeu, σyyNeu, σxyNeu, σyxNeu, sex, sey, R, η, Formulation)

    # Compute mesh properties for FCFV
    @printf("---> Compute FCFV vectors:\n")
    mesh.τ = τr.*ones(mesh.nf)  # Stabilisation per element
    @time ae, be, ze = ComputeFCFV(mesh, sex, sey, VxDir, VyDir, σxxNeu, σyyNeu, σxyNeu, σyxNeu, τr, Formulation)
    
    # Assemble element matrices and RHS
    @printf("---> Compute element matrices:\n")                      
    @time Kuu, Muu, Kup, fu, fp, tsparse = ElementAssemblyLoop(mesh, ae, be, ze, VxDir, VyDir, σxxNeu, σyyNeu, σxyNeu, σyxNeu, Formulation)
    println("---> Sparsification: ", tsparse)

    # # Solve for hybrid variables
    @printf("---> Linear solve:\n")
    Pe .= Pa
    @time StokesSolvers!(Vxh, Vyh, Pe, mesh, Kuu, Kup, fu, fp, Muu, solver; tol=ϵ, penalty=γ)

    # Reconstruct element values
    @printf("---> Compute element values:\n")
    @time Vxe, Vye, Txxe, Tyye, Txye = ComputeElementValues(mesh, Vxh, Vyh, Pe, ae, be, ze, VxDir, VyDir, Formulation)
    
    # Evaluate residuals
    @printf("---> Evaluate residuals:\n")
    ComputeResidualsFCFV_Stokes_o1(Vxh, Vyh, Pe, mesh, ae, be, ze, sex, sey, VxDir, VyDir, σxxNeu, σyyNeu, σxyNeu, σyxNeu, Formulation)

    # Visualise
    @printf("---> Visualisation:\n")
    @time PlotMakie( mesh, Pe,  xmin, xmax, ymin, ymax; cmap=:turbo, min_v=-3, max_v=3, writefig=false )

end

#--------------------------------------------------------------------#

ViscousInclusion()