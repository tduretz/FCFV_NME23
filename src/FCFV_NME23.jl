module FCFV_NME23

using  Revise, Printf, MAT, CairoMakie, Makie.GeometryBasics
import LinearAlgebra: norm, lu, cholesky, Hermitian
import SparseArrays: spdiagm, sparse, dropzeros, spzeros, I
import Statistics: mean

include("CreateMeshFCFV.jl")
export FCFV_Mesh, LoadExternalMesh, LoadExternalMesh2, MakeTriangleMesh
include("VisuFCFV.jl")
export PlotMakie
include("DiscretisationFCFV_Stokes.jl")
export ComputeFCFV, ElementAssemblyLoop, ElementAssemblyLoopNEW, ElementAssemblyLoopToy, ComputeElementValues
include("SolversFCFV_Stokes.jl")
export StokesSolvers!
include("ComputeResidualsFCFV_Stokes.jl")
export ComputeResidualsFCFV_Stokes_o1
include("IntegrationPoints.jl")
export IntegrationTriangle, ShapeFunctions
include("EvalAnalDani.jl")
export EvalAnalDani, Tractions
include("EvalAnalSolKz.jl")
export EvalAnalSolKz

end # module FCFV_NME23