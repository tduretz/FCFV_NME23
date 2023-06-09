module FCFV_NME23

using  Revise, Printf, MAT, CairoMakie, Makie.GeometryBasics
import LinearAlgebra: norm, lu, cholesky, Hermitian
import SparseArrays: spdiagm, sparse, dropzeros
import Statistics: mean

include("CreateMeshFCFV.jl")
export FCFV_Mesh, LoadExternalMesh, MakeTriangleMesh
include("VisuFCFV.jl")
export PlotMakie
include("DiscretisationFCFV_Stokes.jl")
export ComputeFCFV, ElementAssemblyLoop, ComputeElementValues
include("SolversFCFV_Stokes.jl")
export StokesSolvers!
include("EvalAnalDani.jl")
export EvalAnalDani, Tractions
include("ComputeResidualsFCFV_Stokes.jl")
export ComputeResidualsFCFV_Stokes_o1

end # module FCFV_NME23