using CairoMakie
using Makie.GeometryBasics

@doc """
Create patch plots, visualises constant value par element.
""" SetUpProblem!

@views function PlotMakie(mesh, v, xmin, xmax, ymin, ymax; cmap = :viridis, min_v = minimum(v), max_v = maximum(v), writefig=false)
    # min_v, max_v = minimum(v), maximum(v)
    f = Figure()
    ar = (maximum(mesh.xv) - minimum(mesh.xv)) / (maximum(mesh.xv) - minimum(mesh.yv))
    ax = Axis(f[1, 1], aspect = ar)
    p = [Polygon( Point2f0[ (mesh.xv[mesh.e2v[i,j]], mesh.yv[mesh.e2v[i,j]]) for j=1:mesh.nf_el] ) for i in 1:mesh.nel]
    poly!(p, color = v, colormap = cmap, linestyle=:none, strokewidth = 0.0, strokecolor = :black, markerstrokewidth = 0, markerstrokecolor = (0, 0, 0, 0), aspect_ratio=:image, colorrange=(min_v,max_v)) 
    # scatter!(mesh.xf[mesh.bc.==-1] ,mesh.yf[mesh.bc.==-1] )
    Colorbar(f[1, 2], colormap = cmap, limits=(min_v, max_v), flipaxis = true, size = 25, height = Relative(2/3) )
    # colgap!(f.layout, -100)
    # xlims!(ax, -0., 3.)
    # ylims!(ax, -0., 3.)
    resize_to_layout!(f)
    display(f)
    if writefig==true 
        save( string("./images/PressureField.png"), f)
    end
    return nothing
end