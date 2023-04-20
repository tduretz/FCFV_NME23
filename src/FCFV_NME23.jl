module FCFV_NME23

using  Revise, Printf, LinearAlgebra, SparseArrays, MAT, CairoMakie, Makie.GeometryBasics

include("CreateMeshFCFV.jl")
export FCFV_Mesh
# include("VisuFCFV.jl")
include("DiscretisationFCFV_Stokes.jl")
# include("SolversFCFV_Stokes.jl")
include("EvalAnalDani.jl")
export EvalAnalDani, Tractions
# include("ComputeResidualsFCFV_Stokes.jl")

end # module FCFV_NME23