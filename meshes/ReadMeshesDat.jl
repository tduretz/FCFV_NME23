using DelimitedFiles, MAT, CairoMakie

function IsItNew(xf, yf, x, y)
    OK  = true
    ind = 0
    for i in eachindex(xf)
        if xf[i] == x && yf[i] == y 
            OK  = false
            ind = i
        end
    end
    return OK, ind
end

function mesh()
name = "square_TRI_H1"

# Read data
f = readdlm("./meshes/$(name).dat")
nsd       = f[1,1]
nel       = f[1,2]
nv        = f[1,3]
nextfaces = f[1,4]
i  = 2
xv =  f[i:i+nv-1,1]
yv =  f[i:i+nv-1,2]
i += nv 
e2v = f[i:i+nel-1,1:3]
i += nel
extfaces = f[i:i+nextfaces-1,1:3]

f = Figure(resolution = ( 1000,1000), fontsize=25, aspect = 2.0)

ax1 = Axis(f[1, 1], title = "h", xlabel = "x [km]", ylabel = "y [km]")

# faces
ifac = 1
xf   = []
yf   = []
e2f  = zeros(Int64, nel, 3)
for e=1:nel

    xv1 = xv[e2v[e,1]]
    yv1 = yv[e2v[e,1]]
    xv2 = xv[e2v[e,3]]
    yv2 = yv[e2v[e,3]]
    xv3 = xv[e2v[e,2]]
    yv3 = yv[e2v[e,2]]

    xf1 = 0.5*(xv1 + xv2)
    xf2 = 0.5*(xv2 + xv3)
    xf3 = 0.5*(xv1 + xv3)
    yf1 = 0.5*(yv1 + yv2)
    yf2 = 0.5*(yv2 + yv3)
    yf3 = 0.5*(yv1 + yv3)

    new, ind = IsItNew(xf, yf, xf1, yf1)
    if new
        push!(xf, xf1)
        push!(yf, yf1)
        e2f[e,1] = ifac
        ifac    +=1
        # scatter!(ax1, xf1, yf1)
    else
        e2f[e,1] = ind
    end

    new, ind = IsItNew(xf, yf, xf2, yf2)
    if new
        push!(xf, xf2)
        push!(yf, yf2)
        e2f[e,2] = ifac
        ifac    +=1
        # scatter!(ax1, xf2, yf2)
    else
        e2f[e,2] = ind
    end

    new, ind = IsItNew(xf, yf, xf3, yf3)
    if new
        push!(xf, xf3)
        push!(yf, yf3)
        e2f[e,3] = ifac
        ifac    +=1
        # scatter!(ax1, xf3, yf3)
    else
        e2f[e,3] = ind
    end
end

nf = size(xf)

# # construct face BC Array
bc = zeros(Int, nf)
numbc = 0
for i in eachindex(bc)
    if xf[i] == 0. || xf[i] == 1. || yf[i] == 0. || yf[i] == 1.
        numbc+=1
        bc[i] = 1
    end
end

@show numbc

file = matopen("./meshes/$(name).mat", "w")
write(file, "bc", bc)
write(file, "xf", xf)
write(file, "yf", yf)
write(file, "xv", xv)
write(file, "yv", yv)
write(file, "e2f", e2f)
write(file, "e2v", e2v)
close(file)

######################
nodeA = [2 3 1]
nodeB = [3 1 2]

nodeA = [3 2 1]
nodeB = [1 3 2]

    # Assemble FCFV elements
    for iel=1:1#nel 
        for ifac=1:3
            # Face 
            nodei  = e2f[iel,ifac]
            # Vertices
            vert1  = e2v[iel, nodeA[ifac]]
            vert2  = e2v[iel, nodeB[ifac]]
            @show dx     = (xv[vert1] - xv[vert2] )
            @show dy     = (yv[vert1] - yv[vert2] )
            Γi     = sqrt(dx^2 + dy^2);
            # Face normal
            n_x  = -dy/Γi
            n_y  =  dx/Γi
            # Third vector
            xc = 1//3*(sum(xv[e2v[iel,:]]))
            yc = 1//3*(sum(yv[e2v[iel,:]]))
            v_x  = xf[nodei] - xc
            v_y  = yf[nodei] - yc
            # Check whether the normal points outwards
            dot                 = n_x*v_x + n_y*v_y 
            n_x  = ((dot>=0.0)*n_x - (dot<0.0)*n_x)
            n_y  = ((dot>=0.0)*n_y - (dot<0.0)*n_y)

            lines!(ax1, [xv[vert1]; xv[vert2]], [yv[vert1]; yv[vert2]])
            scatter!(ax1, xf[ifac], yf[ifac])
            arrows!(ax1, [xf[ifac]], [yf[ifac]], [n_x], [n_y], arrowsize = 0, lengthscale=.01)

        end
    end

    display(f)

end
mesh()





