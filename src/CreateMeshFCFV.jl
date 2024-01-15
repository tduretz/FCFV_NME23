using MAT
import Triangulate

#--------------------------------------------------------------------#

@doc """
Data structure that contains information about the mesh (number of elements, geometry, connectivity...).
""" FCFV_Mesh

Base.@kwdef mutable struct FCFV_Mesh
    type   ::Union{String, Missing}          = missing # type of mesh (UnstructTriangles)
    nel    ::Union{Int64,  Missing}          = missing # number of elements
    nf     ::Union{Int64,  Missing}          = missing # number of faces 
    nv     ::Union{Int64,  Missing}          = missing # number of vertices
    nn_el  ::Union{Int64,  Missing}          = missing # number of nodes per element
    nf_el  ::Union{Int64,  Missing}          = missing # number of faces per element
    xv     ::Union{Vector{Float64}, Missing} = missing # vertex x coordinate
    yv     ::Union{Vector{Float64}, Missing} = missing # vertex y coordinate
    xf     ::Union{Vector{Float64}, Missing} = missing # face x coordinate
    yf     ::Union{Vector{Float64}, Missing} = missing # face y coordinate
    bc     ::Union{Vector{Int64},   Missing} = missing # node tag
    xc     ::Union{Vector{Float64}, Missing} = missing # cent x coordinate
    yc     ::Union{Vector{Float64}, Missing} = missing # cent y coordinate
    e2v    ::Union{Matrix{Int64},   Missing} = missing # cell 2 node numbering
    e2f    ::Union{Matrix{Int64},   Missing} = missing # cell 2 face numbering
    Ω      ::Union{Vector{Float64}, Missing} = missing # volume of element
    n_x    ::Union{Matrix{Float64}, Missing} = missing # normal 2 face x
    n_y    ::Union{Matrix{Float64}, Missing} = missing # normal 2 face y
    Γ      ::Union{Matrix{Float64}, Missing} = missing # face length
    e2e    ::Union{Matrix{Int64},   Missing} = missing # element 2 element numbering
    f2e    ::Union{Matrix{Int64},   Missing} = missing # face 2 element numbering
    Γ_f    ::Union{Matrix{Float64}, Missing} = missing # face 2 element numbering
    Ω_f    ::Union{Matrix{Float64}, Missing} = missing # volume of element
    n_x_f  ::Union{Matrix{Float64}, Missing} = missing # normal 2 face x
    n_y_f  ::Union{Matrix{Float64}, Missing} = missing # normal 2 face y
    τ      ::Union{Vector{Float64}, Missing} = missing # face stabilisation 
    τe     ::Union{Vector{Float64}, Missing} = missing # element stabilisation 
    ke     ::Union{Vector{Float64}, Missing} = missing # element viscosity
    λ      ::Union{Vector{Float64}, Missing} = missing # viscosity integral over element
    phase  ::Union{Vector{Float64}, Missing} = missing # phase
    δ      ::Union{Vector{Float64}, Missing} = missing # anisotropy factor
    N_x    ::Union{Vector{Float64}, Missing} = missing # director vector x
    N_y    ::Union{Vector{Float64}, Missing} = missing # director vector y
end

#--------------------------------------------------------------------#

@doc """
Generate triangular mesh and populate FCFV_Mesh() structure.
""" MakeTriangleMesh

function MakeTriangleMesh( nx, ny, xmin, xmax, ymin, ymax, τr, inclusion, R, BC=[1; 1; 1; 1;], area = ((xmax-xmin)/nx)*((ymax-ymin)/ny), no_pts_incl = Int64(floor(1.0*pi*R/sqrt(((xmax-xmin)/nx)^2+((ymax-ymin)/ny)^2))); nnel=3, npel=1, xp_in=0, yp_in=0, tp_in=0  )
    # Define data required by triangle mesher
    # This function was adapted from MILAMIN (Dabrowski et al., 2008)
    regions  = Array{Float64}(undef,4,0)
    holes    = Array{Float64}(undef,2,0)
    pts_l    = 1
    pts_u    = 0
    if inclusion>=0
        # 1. Four corners of the domain
        px     = [xmin; xmax; xmax; xmin]
        py     = [ymin; ymin; ymax; ymax]
        sx     = [ 1; 2; 3; 4; ] 
        sy     = [ 2; 3; 4; 1; ]
        st     = BC  # segment markers == boundary flags
        no_pts = length(px)
        pts_l  = pts_l + no_pts
        pts_u  = pts_u + no_pts
        # Region 1
        h1     = [xmin+1e-13; ymin+1e-13; 1.0; 0.0] 
    end
    if inclusion==1
        # 2. perimeter of the inclusion
        theta0       = collect(LinRange(0.0,2.0*pi,no_pts_incl+1));
        theta        = theta0[1:end-1] # do not include last point (periodic 0 == 2pi)
        xx           = cos.(theta);
        yy           = sin.(theta);
        center_x     = (xmax + xmin)/2.0
        center_y     = (ymax + ymin)/2.0
        X            = center_x .+ R*xx;
        Y            = center_y .+ R*yy;
        no_pts       = length(X);
        st1          = -1*ones(1,no_pts);
        pts_u        = pts_u + no_pts;
        sx1          = collect(pts_l:pts_u)
        sy1          = collect(pts_l+1:pts_u+1)
        sy1[end]     = pts_l # Periodic
        h2           = [0.0; 0.0; 2.0; 0.0] # Region 2
        for i=1:no_pts_incl
            px   = push!(px, X[i])
            py   = push!(py, Y[i])
            sx   = push!(sx, sx1[i])
            sy   = push!(sy, sy1[i])
            st   = push!(st, st1[i])
        end
        regions = hcat(h1,h2)
    end
    # Concatenate arrays for better digestion
    p       = hcat(px, py)         # points
    s       = hcat(sx, sy)         # segments
    p       = p'
    s       = s'
    st      = st[:] 
    # Triangulation
    mesh                    = FCFV_Mesh()
    mesh.type               = "UnstructTriangles"
    triin                   = Triangulate.TriangulateIO()
    triin.pointlist         = Matrix{Cdouble}(vcat(px',py'))
    triin.segmentlist       = Matrix{Cint}(vcat(sx',sy'))
    triin.segmentmarkerlist = Vector{Int32}(st)
    triin.regionlist        = Matrix{Cdouble}(regions)
    astring                 = @sprintf("%0.10lf", area)
    (trimesh, vorout)       = Triangulate.triangulate("vQDpenq33o2IAa$(astring)", triin)  #
    # Read resulting mesh
    mesh.nel                = size(trimesh.trianglelist,2)
    e2v                     = trimesh.trianglelist[1:3,:]
    mesh.e2v                = e2v'
    mesh.nv                 = maximum(e2v)
    e2f                     = trimesh.trianglelist[4:6,:] .- mesh.nv
    mesh.e2f                = e2f'
    mesh.nf                 = maximum(e2f)
    mesh.xv                 = trimesh.pointlist[1,1:mesh.nv]
    mesh.yv                 = trimesh.pointlist[2,1:mesh.nv]
    mesh.xf                 = trimesh.pointlist[1,mesh.nv+1:end]
    mesh.yf                 = trimesh.pointlist[2,mesh.nv+1:end]
    mesh.bc                 = trimesh.pointmarkerlist[mesh.nv+1:end]
    mesh.phase              = trimesh.triangleattributelist[:] 
    mesh.δ                  =  ones(Float64,mesh.nel) # Default assume isotropic viscosity
    mesh.N_x                =  ones(Float64,mesh.nel) # Default assume flat anisotropy
    mesh.N_x                = zeros(Float64,mesh.nel)
    if inclusion==0 mesh.phase = ones(mesh.nel) end 
    mesh.ke                 =  ones(Float64,mesh.nel)
    mesh.τ                  = zeros(Float64,mesh.nf)   # for each face
    mesh.τe                 = zeros(Float64,mesh.nel)  # for each element
    nel                     = mesh.nel
    mesh.nn_el              = 3
    mesh.nf_el              = 3
    Ωe                      = zeros(nel)
    xc                      = zeros(nel)
    yc                      = zeros(nel)
    # Compute volumes of triangles - use vertices coordinates
    for iel=1:nel
        x1 = mesh.xv[e2v[1,iel]]
        y1 = mesh.yv[e2v[1,iel]]
        x2 = mesh.xv[e2v[2,iel]]
        y2 = mesh.yv[e2v[2,iel]]
        x3 = mesh.xv[e2v[3,iel]]
        y3 = mesh.yv[e2v[3,iel]]
        a         = sqrt((x1-x2)^2 + (y1-y2)^2)
        b         = sqrt((x2-x3)^2 + (y2-y3)^2)
        c         = sqrt((x1-x3)^2 + (y1-y3)^2)
        s         = 1/2*(a+b+c)
        Ωe[iel]   = sqrt(s*(s-a)*(s-b)*(s-c))
        xc[iel]   = 1.0/3.0*(x1+x2+x3)
        yc[iel]   = 1.0/3.0*(y1+y2+y3)
    end
    mesh.xc     = xc
    mesh.yc     = yc
    mesh.Ω      = Ωe
    # Local numbering
    nodeA = [2 3 1]
    nodeB = [3 1 2]
    nodeC = [1 2 3]
    # Compute normal to faces
    mesh.n_x = zeros(Float64, mesh.nel,mesh.nf_el)
    mesh.n_y = zeros(Float64, mesh.nel,mesh.nf_el)
    mesh.Γ   = zeros(Float64, mesh.nel,mesh.nf_el)
    mesh.e2e = zeros(  Int64, mesh.nel,mesh.nf_el)
    mesh.ke  =  ones(Float64, mesh.nel)
    mesh.λ   = zeros(Float64, mesh.nel)
    # Receive arrays from Triangulate
    vorodeges = zeros(Int64,size(vorout.edgelist))
    vorodeges[1,:] .= vorout.edgelist[1,:]
    vorodeges[2,:] .= vorout.edgelist[2,:]
    seglist   = zeros(Int64,size(trimesh.edgelist))
    seglist[1,:]   .= trimesh.edgelist[1,:]
    seglist[2,:]   .= trimesh.edgelist[2,:]

    # Assemble FCFV elements
    @inbounds for iel=1:mesh.nel 
        for ifac=1:mesh.nf_el
            # Face 
            nodei  = mesh.e2f[iel,ifac]
            # Neighbouring element for this face
            ele1 = vorodeges[1,nodei]
            ele2 = vorodeges[2,nodei]
            iel2 = (iel==ele1) * ele2 + (iel==ele2) * ele1;
            mesh.e2e[iel,ifac] = iel2;
            # Vertices
            vert1  = mesh.e2v[iel,nodeA[ifac]]
            vert2  = mesh.e2v[iel,nodeB[ifac]]
            vert3  = mesh.e2v[iel,nodeC[ifac]]
            bc     = mesh.bc[nodei]
            dx     = (mesh.xv[vert1] - mesh.xv[vert2] );
            dy     = (mesh.yv[vert1] - mesh.yv[vert2] );
            Γi     = sqrt(dx^2 + dy^2);
            # Face normal
            n_x  = -dy/Γi
            n_y  =  dx/Γi
            # Third vector
            v_x  = mesh.xf[nodei] - mesh.xc[iel]
            v_y  = mesh.yf[nodei] - mesh.yc[iel]
            # Check whether the normal points outwards
            dot                 = n_x*v_x + n_y*v_y 
            mesh.n_x[iel,ifac]  = ((dot>=0.0)*n_x - (dot<0.0)*n_x)
            mesh.n_y[iel,ifac]  = ((dot>=0.0)*n_y - (dot<0.0)*n_y)
            mesh.Γ[iel,ifac]    = Γi
            # Face stabilisation
            mesh.τ[nodei] = τr 
        end
    end

    return mesh
    
end

#--------------------------------------------------------------------#

@doc """
Reads in externally generated mesh and populate FCFV_Mesh() structure.
""" LoadExternalMesh

function LoadExternalMesh( res, η )
    mesh       = FCFV_Mesh( )
    mesh.type  = "InHouseMesher"
    if res==:LowRes  data = matread("./meshes/MeshAdvancingFrontLowRes.mat" ) end
    if res==:MedRes  data = matread("./meshes/MeshAdvancingFrontMedRes.mat" ) end
    if res==:HighRes data = matread("./meshes/MeshAdvancingFrontHighRes.mat") end
    mesh.nel   = data["nel"]
    mesh.nf    = data["nf"]
    mesh.nn_el = data["nn_el"]
    mesh.nf_el = data["nf_el"]
    mesh.xv    = data["xn"][:]
    mesh.yv    = data["yn"][:]
    mesh.xf    = data["xf"][:]
    mesh.yf    = data["yf"][:]
    mesh.xc    = data["xc"][:]
    mesh.yc    = data["yc"][:]
    mesh.e2f   = data["e2f"]
    mesh.e2v   = data["e2n"]
    mesh.Γ     = data["Gamma"]
    mesh.Ω     = data["Omega"][:]
    mesh.n_x   = data["nx"]
    mesh.n_y   = data["ny"]
    mesh.phase = data["phase"][:]
    mesh.bc    = data["bc"][:]
    mesh.ke    = η[Int64.(mesh.phase)]
    mesh.δ     =  ones(Float64,mesh.nel) # Default assume isotropic viscosity
    mesh.N_x   =  ones(Float64,mesh.nel) # Default assume flat anisotropy
    mesh.N_x   = zeros(Float64,mesh.nel)
    return mesh
end

#--------------------------------------------------------------------#

@doc """
Generate triangular mesh and reads in externally generated mesh and populate FCFV_Mesh() structure.
""" LoadExternalMesh2

function LoadExternalMesh2( nx, ny, xmin, xmax, ymin, ymax, τr, inclusion, R, BC=[1; 1; 1; 1;], area = ((xmax-xmin)/nx)*((ymax-ymin)/ny), no_pts_incl = Int64(floor(1.0*pi*R/sqrt(((xmax-xmin)/nx)^2+((ymax-ymin)/ny)^2))); nnel=3, npel=1, xp_in=0, yp_in=0, tp_in=0  )
    
    file = matopen("./meshes/square_TRI_H1.mat")
    xf  = read(file, "xf") # note that this does NOT introduce a variable ``varname`` into scope
    yf  = read(file, "yf") # note that this does NOT introduce a variable ``varname`` into scope
    e2f = read(file, "e2f") # note that this does NOT introduce a variable ``varname`` into scope
    xv  = read(file, "xv") # note that this does NOT introduce a variable ``varname`` into scope
    yv  = read(file, "yv") # note that this does NOT introduce a variable ``varname`` into scope
    e2v = read(file, "e2v") # note that this does NOT introduce a variable ``varname`` into scope
    bc  = read(file,  "bc") # note that this does NOT introduce a variable ``varname`` into scope
    close(file)

    # Read resulting mesh
    mesh                    = FCFV_Mesh()
    mesh.type               = "UnstructTriangles"
    mesh.nel                = size(e2f,1)
    mesh.e2v                = e2v
    mesh.nv                 = length(xv)
    mesh.e2f                = Int64.(e2f)
    mesh.nf                 = maximum(e2f)
    mesh.xv                 = xv
    mesh.yv                 = yv
    mesh.xf                 = xf
    mesh.yf                 = yf
    mesh.bc                 = bc
    mesh.phase              = ones(Int64,mesh.nel) 
    # if inclusion==0 mesh.phase = ones(mesh.nel) end 
    mesh.ke                 = ones(Float64,mesh.nel)
    mesh.τ                  = zeros(Float64,mesh.nf)   # for each face
    mesh.τe                 = zeros(Float64,mesh.nel)  # for each element
    nel                     = mesh.nel
    mesh.nn_el              = 3
    mesh.nf_el              = 3
    Ωe                      = zeros(nel)
    xc                      = zeros(nel)
    yc                      = zeros(nel)
    mesh.δ                  =  ones(Float64,mesh.nel) # Default assume isotropic viscosity
    mesh.N_x                =  ones(Float64,mesh.nel) # Default assume flat anisotropy
    mesh.N_x                = zeros(Float64,mesh.nel)
    # Compute volumes of triangles - use vertices coordinates
    for iel=1:nel
        xv1 = mesh.xv[mesh.e2v[iel,1]]
        yv1 = mesh.yv[mesh.e2v[iel,1]]
        xv2 = mesh.xv[mesh.e2v[iel,2]]
        yv2 = mesh.yv[mesh.e2v[iel,2]]
        xv3 = mesh.xv[mesh.e2v[iel,3]]
        yv3 = mesh.yv[mesh.e2v[iel,3]]
        a         = sqrt((xv1-xv2)^2 + (yv1-yv2)^2)
        b         = sqrt((xv2-xv3)^2 + (yv2-yv3)^2)
        c         = sqrt((xv1-xv3)^2 + (yv1-yv3)^2)
        s         = 1/2*(a+b+c)
        Ωe[iel]   = sqrt(s*(s-a)*(s-b)*(s-c))
        xc[iel]   = 1.0/3.0*(xv1 + xv2 + xv3)
        yc[iel]   = 1.0/3.0*(yv1 + yv2 + yv3)
    end
    mesh.xc     = xc
    mesh.yc     = yc
    mesh.Ω      = Ωe
    # Local numbering
    nodeA = [3 2 1]        # Modified for this specific mesh... beurk!
    nodeB = [1 3 2]
    # Compute normal to faces
    mesh.n_x = zeros(Float64, mesh.nel, mesh.nf_el)
    mesh.n_y = zeros(Float64, mesh.nel, mesh.nf_el)
    mesh.Γ   = zeros(Float64, mesh.nel, mesh.nf_el)
    mesh.e2e = zeros(  Int64, mesh.nel, mesh.nf_el)
    mesh.ke  =  ones(Float64, mesh.nel)
    mesh.λ   = zeros(Float64, mesh.nel)
    # Receive arrays from Triangulate
    # # vorodeges = zeros(Int64,size(vorout.edgelist))
    # # vorodeges[1,:] .= vorout.edgelist[1,:]
    # # vorodeges[2,:] .= vorout.edgelist[2,:]
    # # seglist   = zeros(Int64,size(trimesh.edgelist))
    # # seglist[1,:]   .= trimesh.edgelist[1,:]
    # # seglist[2,:]   .= trimesh.edgelist[2,:]

    # Assemble FCFV elements
    @inbounds for iel=1:mesh.nel 
        for ifac=1:mesh.nf_el
            # Face 
            nodei  = mesh.e2f[iel,ifac]
            # Neighbouring element for this face
            # ele1 = vorodeges[1,nodei]
            # ele2 = vorodeges[2,nodei]
            # iel2 = (iel==ele1) * ele2 + (iel==ele2) * ele1;
            # mesh.e2e[iel,ifac] = iel2;
            # Vertices
            vert1  = mesh.e2v[iel,nodeA[ifac]]
            vert2  = mesh.e2v[iel,nodeB[ifac]]
            bc     = mesh.bc[nodei]
            dx     = (mesh.xv[vert1] - mesh.xv[vert2] )
            dy     = (mesh.yv[vert1] - mesh.yv[vert2] )
            Γi     = sqrt(dx^2 + dy^2);
            # Face normal
            n_x  = -dy/Γi
            n_y  =  dx/Γi
            # Third vector
            v_x  = mesh.xf[nodei] - mesh.xc[iel]
            v_y  = mesh.yf[nodei] - mesh.yc[iel]
            # Check whether the normal points outwards
            dot                 = n_x*v_x + n_y*v_y 
            mesh.n_x[iel,ifac]  = ((dot>=0.0)*n_x - (dot<0.0)*n_x)
            mesh.n_y[iel,ifac]  = ((dot>=0.0)*n_y - (dot<0.0)*n_y)
            mesh.Γ[iel,ifac]    = Γi
            # Face stabilisation
            mesh.τ[nodei] = τr 
        end
    end    

    return mesh
    
end
