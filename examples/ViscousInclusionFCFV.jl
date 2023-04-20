using FCFV_NME23, Printf, MAT

function SetUpProblem!(mesh, P, P_f, Vx, Vy, Sxx, Syy, Sxy, VxDir, VyDir, SxxNeu, SyyNeu, SxyNeu, SyxNeu, sx, sy, R, eta)
    # Evaluate T analytic on cell faces
    etam, etai = eta[1], eta[2]
    for in=1:mesh.nf
        x        = mesh.xf[in]
        y        = mesh.yf[in]
        vx, vy, p, sxx, syy, sxy = EvalAnalDani( x, y, R, etam, etai )
        VxDir[in] = vx
        VyDir[in] = vy
        # Stress at faces - Pseudo-tractions
        p, dVxdx, dVxdy, dVydx, dVydy = Tractions( x, y, R, etam, etai, 1 )
        SxxNeu[in] = - p + etam*dVxdx 
        SyyNeu[in] = - p + etam*dVydy 
        SxyNeu[in] =       etam*dVxdy
        SyxNeu[in] =       etam*dVydx
        P_f[in]    = p
    end
    # Evaluate T analytic on barycentres
    for iel=1:mesh.nel
        x        = mesh.xc[iel]
        y        = mesh.yc[iel]
        vx, vy, pre, sxx, syy, sxy = EvalAnalDani( x, y, R, etam, etai )
        P[iel]   = pre
        Vx[iel]  = vx
        Vy[iel]  = vy
        Sxx[iel] = sxx
        Syy[iel] = syy
        Sxy[iel] = sxy
        sx[iel]  = 0.0
        sy[iel]  = 0.0
        out          = mesh.phase[iel] == 1.0
        mesh.ke[iel] = (out==1) * 1.0*etam + (out!=1) * 1.0*etai         
    end
    return
end

function LoadExternalMesh( res, η )

    mesh        = FCFV_Mesh( )
    mesh.type   = "Rubén_Mesh3"
    if res==:LR data = matread("./meshes/Mesh1Ruben.mat") end
    if res==:MR data = matread("./meshes/Mesh2Ruben.mat") end
    if res==:HR data = matread("./meshes/Mesh3Ruben.mat") end
    mesh.nel   = data["nel"]
    mesh.nf    = data["nf"]
    mesh.nn_el = data["nn_el"]
    mesh.nf_el = data["nf_el"]
    mesh.xn    = data["xn"][:]
    mesh.yn    = data["yn"][:]
    mesh.xv    = data["xn"][:]
    mesh.yv    = data["yn"][:]
    mesh.xf    = data["xf"][:]
    mesh.yf    = data["yf"][:]
    mesh.xc    = data["xc"][:]
    mesh.yc    = data["yc"][:]
    mesh.e2f   = data["e2f"]
    mesh.e2n   = data["e2n"]
    mesh.e2v   = data["e2n"]
    mesh.Γ     = data["Gamma"]
    mesh.Ω     = data["Omega"][:]
    mesh.n_x   = data["nx"]
    mesh.n_y   = data["ny"]
    mesh.phase = data["phase"][:]
    mesh.bc    = data["bc"][:]
    mesh.ke    = η[Int64.(mesh.phase)]

    return mesh
end

function ViscousInclusion()

    @printf("Viscous inclusion test using first order FCFV discretisation on triangles")

    # Physics
    xmin, xmax = -3.0, 3.0    # Domain extent x
    ymin, ymax = -3.0, 3.0    # Domain extent y
    R          = 1.0          # Inclusion radius
    η          = [1.0 1e-2]   # Viscosity matrix/inclusion
    BC         = [1; 1; 1; 1] # South/East/North/West --- 1: Dirichlet / 2: Neumann

    # Numerics
    solver     = 0

    # Load mesh 
    mesh = LoadExternalMesh(:LR, η)
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
    Pa_f   = zeros(mesh.nf)
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
    println("Model configuration :")
    @time SetUpProblem!(mesh, Pa, Pa_f, Vxa, Vya, Sxxa, Syya, Sxya, VxDir, VyDir, SxxNeu, SyyNeu, SxyNeu, SyxNeu, sex, sey, R, η)

end

ViscousInclusion()