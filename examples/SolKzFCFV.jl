using FCFV_NME23, Printf, CairoMakie, LinearAlgebra, MathTeXEngine, Makie.GeometryBasics
using SparseArrays

#--------------------------------------------------------------------#

@doc """
Generates model configuration (viscosity field) and set up boundary values.
""" SetUpProblem!

function SetUpProblem!(mesh, P, Vx, Vy, σxx, σyy, σxy, VxDir, VyDir, σxxNeu, σyyNeu, σxyNeu, σyxNeu, sx, sy, ηip, N, nnel, nip)
    # Evaluate analytics on cell faces
    for in=1:mesh.nf
        x         = mesh.xf[in]
        y         = mesh.yf[in]
        vx, vy, p, _, η = EvalAnalSolKz(x, y)
        VxDir[in] = vx
        VyDir[in] = vy
    end
    # Evaluate analytics on barycentres
    Ni      = zeros(nnel)
    xi      = zeros(1,2)
    for e=1:mesh.nel
        x        = mesh.xc[e]
        y        = mesh.yc[e]
        vx, vy, p, ρ, η = EvalAnalSolKz(x, y)
        mesh.ke[e] = η
        P[e]   =  p
        Vx[e]  =  vx
        Vy[e]  =  vy
        σxx[e] = 0.0
        σyy[e] = 0.0
        σxy[e] = 0.0
        sx[e]  = 0.0
        sy[e]  = -ρ
        # Numbering and element nodes
        nodes    = mesh.e2v[e,:]
        x        = [mesh.xv[nodes] mesh.yv[nodes]]
        for ip=1:nip
            # Global coordinate of integration points
            Ni .= N[ip,:,:]
            xi .= Ni'*x
            _, _, _, _, η = EvalAnalSolKz(xi[1], xi[2])
            ηip[e,ip] = η
        end
    end
    return
end

#--------------------------------------------------------------------#

function SolKz(n)

    @printf("SolKz test using first order FCFV discretisation on triangles\n")

    # Physics
    xmin, xmax  = -.0, 1.0    # Domain extent x
    ymin, ymax  = -.0, 1.0    # Domain extent y
    R           = 10.0          # Inclusion radius
    η           = [1.0 1e-2]   # Viscosity matrix/inclusion
    BC          = [1; 1; 1; 1] # South/East/North/West --- 1: Dirichlet / 2: Neumann

    # Numerics
    Mesher      = :Delaunay                # :Delaunay / :AdvancingFront (load external mesh)
    # Mesher      = :AdvancingFront          # :Delaunay / :AdvancingFront (load external mesh)
    solver      = :PowellHestenesCholesky  # :CoupledBackslash / :=PowellHestenesCholesky / :=PowellHestenesLU
    # solver      = :CoupledBackslash
    # solver      = :PowellHestenesLU      # :CoupledBackslash / :=PowellHestenesCholesky / :=PowellHestenesLU
    Formulation = :SymmetricGradient       # :Gradient / :SymmetricGradient
    # Formulation = :Gradient              # :Gradient / :SymmetricGradient
    τr          = 10                       # Stabilisation
    γ           = 5e3                      # Penalty factor for Powell-Hestenes solvers
    ϵ           = 1e-8                     # Tolerance of Powell-Hestenes solvers 

    # Generate mesh 
    nx, ny = n*100, n*100  # initial point density in x and y for triangulation 
    if Mesher==:Delaunay       mesh = MakeTriangleMesh( nx, ny, xmin, xmax, ymin, ymax, τr, 0, R, BC, ((xmax-xmin)/nx)*((ymax-ymin)/ny), 200 ) end
    if Mesher==:AdvancingFront mesh = LoadExternalMesh2( nx, ny, xmin, xmax, ymin, ymax, τr, 0, R, BC, ((xmax-xmin)/nx)*((ymax-ymin)/ny), 200 ) end

    @printf("Mesh informations:\n")
    @printf("--> %d elements\n", mesh.nel)
    @printf("--> %d faces\n",    mesh.nf)
    @printf("--> %d matrix elements\n",    sum(mesh.phase.==1))
    @printf("--> %d inclusion elements\n", sum(mesh.phase.==2))
    if Formulation==:Gradient && solver==:PowellHestenesCholesky
        error("Gradient formulation leads to non symmetric system, can't apply Cholesky factorisation!")
    end
    @show sum(mesh.Ω)
    
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
    Pa     = zeros( mesh.nel)  # Element pressure
    Vxa    = zeros(mesh.nel)  # Element horizontal velocity
    Vya    = zeros(mesh.nel)  # Element vertical velocity

    # Viscosity field on integration points
    nip      = 6
    nnel     = 3
    ηip      = zeros(mesh.nel, nip) 
    ipx, ipw = IntegrationTriangle( nip )
    N, dNdX  = ShapeFunctions(ipx, nip, nnel)

    # Set up Stokes problem
    @printf("---> Model configuration :\n")
    @time SetUpProblem!(mesh, Pa, Vxa, Vya, σxxa, σyya, σxya, VxDir, VyDir, σxxNeu, σyyNeu, σxyNeu, σyxNeu, sex, sey, ηip, N, nnel, nip)

    # Viscosity integration 
    J       = zeros(2,2)
    dNdXi   = zeros(nnel, 2)
    mesh.λ .= 0.
    for e = 1:mesh.nel
        nodes   = mesh.e2v[e,:]
        x       = [mesh.xv[nodes] mesh.yv[nodes]]  
        for ip=1:nip
            # Evaluate viscosity on integration points
            η          = ηip[e,ip]
            # Jacobian
            dNdXi     .= dNdX[ip,:,:]
            mul!(J, x', dNdXi)
            # Interpolation weight
            detJ       = J[1,1]*J[2,2] - J[1,2]*J[2,1]
            w          = ipw[ip] * detJ
            # Integration
            mesh.λ[e] += w * 1.0/η
        end
    end
    @show minimum(mesh.λ)
    @show maximum(mesh.λ)
    @show minimum(mesh.δ)
    @show maximum(mesh.δ)

    # mesh.λ .= mesh.Ω./mesh.ke

    # Compute mesh properties for FCFV
    @printf("---> Compute FCFV vectors:\n")
    mesh.τe .= τr.*max.(mesh.ke, 1.0)  # Stabilisation per element
    @time ae, be, ze = ComputeFCFV(mesh, sex, sey, VxDir, VyDir, σxxNeu, σyyNeu, σxyNeu, σyxNeu, τr, Formulation)

    # Assemble element matrices and RHS
    @printf("---> Compute element matrices:\n") 
    @time Kuu, Muu, Kup, fu, fp = ElementAssemblyLoopNEW(mesh, ae, be, ze, VxDir, VyDir, σxxNeu, σyyNeu, σxyNeu, σyxNeu, Formulation);
    # @time Kuu, Muu, Kup, fu, fp = ElementAssemblyLoopToy(mesh, ae, be, ze, VxDir, VyDir, σxxNeu, σyyNeu, σxyNeu, σyxNeu, Formulation);

    ########################
    # f, ax, plt = spy(spdiagm(diag(Kuu)), markersize = 4, marker = :circle, framecolor = :lightgrey)
    
    # hidedecorations!(ax) # remove axis labeling
    # ax.title = "Visualization of a random sparse matrix"
    
    # f

    # # Assemble element matrices and RHS
    # @printf("---> Compute element matrices:\n")                      
    # @time Kuu, Muu, Kup, fu, fp, tsparse = ElementAssemblyLoopNEW(mesh, ae, be, ze, VxDir, VyDir, σxxNeu, σyyNeu, σxyNeu, σyxNeu, Formulation)
    # println("---> Sparsification: ", tsparse)

    # # Solve for hybrid variables
    @printf("---> Linear solve:\n")
    Pe .= Pa
    @time StokesSolvers!(Vxh, Vyh, Pe, mesh, Kuu, Kup, fu, fp, Muu, solver; tol=ϵ, penalty=γ)

    # Reconstruct element values
    @printf("---> Compute element values:\n")
    @time Vxe, Vye, Txxe, Tyye, Txye = ComputeElementValues(mesh, Vxh, Vyh, Pe, ae, be, ze, VxDir, VyDir, Formulation)
    
    # # Evaluate residuals
    # @printf("---> Evaluate residuals:\n")
    # ComputeResidualsFCFV_Stokes_o1(Vxh, Vyh, Pe, mesh, ae, be, ze, sex, sey, VxDir, VyDir, σxxNeu, σyyNeu, σxyNeu, σyxNeu, Formulation)

    # Visualise
    @printf("---> Visualisation:\n")

    # f = Figure(resolution = (400, 1200), fontsize=22)

    # ax1 = Axis(f[1, 1], title = L"$V_x$", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
    # p = [Polygon( Point2f0[ (mesh.xv[mesh.e2v[i,j]], mesh.yv[mesh.e2v[i,j]]) for j=1:mesh.nf_el] ) for i in 1:mesh.nel]
    # poly!(p, color = Vxe, colormap = :turbo, strokewidth = 0, strokecolor = :black, markerstrokewidth = 0, markerstrokecolor = (0, 0, 0, 0), aspect_ratio=:image, colorrange=(-2.5e-4,2.5e-4)) 
    # # scatter!(mesh.xf[mesh.bc.==-1] ,mesh.yf[mesh.bc.==-1] )
    # Colorbar(f[1, 2], colormap = :turbo, limits=(-2.5e-4, 2.5e-4), flipaxis = true, size = 25, height = Relative(2/3) )
    
    # ax2 = Axis(f[2, 1], title = L"$V_y$", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
    # p = [Polygon( Point2f0[ (mesh.xv[mesh.e2v[i,j]], mesh.yv[mesh.e2v[i,j]]) for j=1:mesh.nf_el] ) for i in 1:mesh.nel]
    # poly!(p, color = Vye, colormap = :turbo, strokewidth = 0, strokecolor = :black, markerstrokewidth = 0, markerstrokecolor = (0, 0, 0, 0), aspect_ratio=:image, colorrange=(-1.25e-4,1.25e-4)) 
    # # scatter!(mesh.xf[mesh.bc.==-1] ,mesh.yf[mesh.bc.==-1] )
    # Colorbar(f[2, 2], colormap = :turbo, limits=(-1.25e-4, 1.25e-4), flipaxis = true, size = 25, height = Relative(2/3) )
    
    # ax2 = Axis(f[3, 1], title = L"$P$", xlabel = L"$x$ [km]", ylabel = L"$y$ [km]")
    # p = [Polygon( Point2f0[ (mesh.xv[mesh.e2v[i,j]], mesh.yv[mesh.e2v[i,j]]) for j=1:mesh.nf_el] ) for i in 1:mesh.nel]
    # poly!(p, color = Pe, colormap = :turbo, strokewidth = 0, strokecolor = :black, markerstrokewidth = 0, markerstrokecolor = (0, 0, 0, 0), aspect_ratio=:image, colorrange=(-0.1,0.1)) 
    # # scatter!(mesh.xf[mesh.bc.==-1] ,mesh.yf[mesh.bc.==-1] )
    # Colorbar(f[3, 2], colormap = :turbo, limits=(-0.1, 0.1), flipaxis = true, size = 25, height = Relative(2/3) )
    
    # display(f)

    @show norm(VxDir.-Vxh)/norm(VxDir)

    # PlotMakie( mesh, Vxe,  xmin, xmax, ymin, ymax; cmap=:turbo, min_v=-2.5e-4, max_v=2.5e-4, writefig=false )
    # PlotMakie( mesh, Vya,  xmin, xmax, ymin, ymax; cmap=:turbo, min_v=-1.25e-4, max_v=1.25e-4, writefig=false )
    PlotMakie( mesh,  Pe,   xmin, xmax, ymin, ymax; cmap=:jet, min_v=-.1, max_v=.1, writefig=false )
    # # @time PlotMakie( mesh, sey,  xmin, xmax, ymin, ymax; cmap=:turbo, min_v=minimum(sey), max_v=maximum(sey), writefig=false )
    # # @time PlotMakie( mesh, log10.(mesh.ke),  xmin, xmax, ymin, ymax; cmap=:turbo, min_v=minimum(log10.(mesh.ke)), max_v=maximum(log10.(mesh.ke)), writefig=false )
end

#--------------------------------------------------------------------#

SolKz(1)