using Graphs
using GraphPlot
using Plots
#need following for GraphPlot
using Compose
import Cairo, Fontconfig
using LinearAlgebra

product = Base.Iterators.product

# some function from T² to S²
# x,y ∈ [0,1]
function myf(x, y)
  #return -sphr_euc(-x*2*π, -y*2*π)
  return sphr_euc(x*2*π, x < 0.5 ? y*2*π : -y*2*π)
end

# same function, just taking tuple as argument
myf(tup) = myf(tup[1],tup[2])

# point on S² from spherical coordinates
# θ : polar angle (z-axis = 0)
# φ : azimuthal angle (from x-z plane)
function sphr_euc(θ, φ)
  return [sin(θ)*cos(φ), sin(θ)*sin(φ), cos(θ)]
end
sphr_euc(tup) = sphr_euc(tup[1],tup[2])

dist(v,w) = sqrt(sum(((v - w) .^ 2)))
dist(v) = sqrt(sum(((v) .^ 2)))

# v,w : tuple/array (θ,φ)
function euc_dist(v,w)
  return dist(sphr_euc(v), sphr_euc(w))
end

# ============================================
# planar graph functionality

# will assume set of vertices is labeled by 1:N
# E[i] is dictionary { u => (l,r) }, where
# u,l,r are neighbors of i
# l is to the left of i->u, r is to the right
struct PlanarGraph
  N::Int
  E::Array{Dict{Int, Tuple{Int,Int}}}
end

# F = faces (u,v,w) - ordered counterclockwise
function PlanarGraph(N::Int, F::Array{Tuple{Int,Int,Int}})
  E = [Dict() for i ∈ 1:N]
  for f ∈ F
    for i ∈ 1:3
      u = f[i]
      v = f[(i+1) > 3 ? i-2 : i+1]
      w = f[(i+2) > 3 ? i-1 : i+2]
      !haskey(E[u], v) && (E[u][v] = (0,0))
      !haskey(E[u], w) && (E[u][w] = (0,0))
      E[u][v] = (w, E[u][v][2]) # tuples are not mutable
      E[u][w] = (E[u][w][1], v)
    end
  end
  return PlanarGraph(N,E)
end

neigh_left(G, u, v) = G.E[u][v][1]
neigh_right(G, u, v) = G.E[u][v][2]
# return value of 0 indicates that u->v is a boundary edge!

# test
N = 4
F = [(1,2,3),(1,3,4)]
G = PlanarGraph(N,F)

# ============================================
# set up triangulation on T²

# vertices
L, W = 100,100
N = L * W
Vcoord = [v for v ∈  product(LinRange(0,1-1/W,W),LinRange(0,1-1/L,W))]

# discretized myf
myf_dscrt(i :: Int64) = myf(Vcoord[i])

ind(i,j) = i + (j - 1)*W
# @assert Vcoord[i,j] == Vcoord[ind(i,j)]

#tri = Graph(N)
#V = 1:N
F = [(0,0,0) for i in 1:2*N]
for i ∈ 1:W
  for j ∈ 1:L
    a₁ = ind(i,j)
    a₂ = ind(i+1 > W ? 1 : i+1, j)
    a₃ = ind(i, j+1 > L ? 1 : j+1)
    a₄ = ind(i+1 > W ? 1 : i+1, j+1 > L ? 1 : j+1)
    F[a₁*2 - 1] = (a₁,a₂,a₄)
    F[a₁*2] = (a₁,a₄,a₃)
    #add_edge!(tri, a₁, a₂)
    #add_edge!(tri, a₁, a₃)
    #add_edge!(tri, a₁, a₄)
  end
end

Tri = PlanarGraph(N,F)

#draw(PNG("lattice-torus.png", 8cm, 8cm), gplot(tri, nodelabel=1:N))

# get max dist of adj vertics under myf
dd = [[(dist(myf_dscrt(u), myf_dscrt(v)), u, v) for v in keys(Tri.E[u])] for u in 1:N]
dd = collect(Iterators.flatten(dd))
max_dist_ind = argmax(dd)



# ============================================
# Choose some p₀ ∈ S², ε > 0, let U = N_ε(p₀)
# Let V₁ = { v ∈ tri | f(v) ∈ U }, V₂ = V(tri) - V₁
# ε large enough that for any triangle Δ of tri,
# if the image f(Δ) contains p₀,
# then at least one of its vertices is in N_ε(p₀)
# (the smaller |Df| is and the finer the triangulation on T² is,
# the smaller ε can be while satisfying this)
# let γ₁,... be curves at the boundary of V₂,
# i.e. they're the smallest curves in V₂ bounding the conn. comps of V₁.
# ε must also be small enough so that these γ's
# are contained in an annulus aronud U.
# Orient the γ's in T² so that they go anti-clockwise around V₁.
# The sum of the winding numbers of f(γ) around p₀
# is equal to the wrapping number of f.
#
# The (oriented) edges in such a γ is characterized by the property that
# the edges to the left are between V₁ and V₂,
# i.e. if uv ∈ γ, and uvw forms a triangle, then w ∈ V₁, u,v ∈ V₂

p₀ = sphr_euc(π/2, 0)
ε = 0.4
V₁ = Set()
V₂ = Set()
for i ∈ 1:N
  dist(myf_dscrt(i), p₀) < ε ? push!(V₁, i) : push!(V₂, i)
end


γ_list = Dict()
for u ∈ V₂
  for v ∈ keys(Tri.E[u])
    if v ∈  V₂ && Tri.E[u][v][1] ∈ V₁
      γ_list[u] = v
    end
  end
end

# split into closed curves
γs = []
visited = Set()
for u ∈ keys(γ_list)
  if u ∈ visited
    continue
  end

  new_γ = []
  v = u
  while true
    if v ∈ visited
      break
    end
    push!(visited, v)
    push!(new_γ, v)
    v = γ_list[v]
  end
  push!(γs, new_γ)
end

# compute winding numbers
# use fact that p₀ is unit vector in x-axis, so just project to y-z plane
tot_winding = 0
for γ ∈ γs
  winding_angle = 0.0
  if length(γ) ≤ 1
    continue
  end
  for i in 1:length(γ)
    p₁ = myf_dscrt(γ[i])
    p₂ = myf_dscrt(γ[i+1 > length(γ) ? 1 : i+1])
    # project in p₀ = (1,0,0) direction to y-z plane
    q₁ = cross(p₀,p₁)
    q₂ = cross(p₀,p₂)
    # get sin(θ) = prllgrm_area(p₁,p₂) / |p₁| |p₂|
    area = cross(q₁,q₂)
    sintheta = dot(area, p₀) / (dist(p₀) * dist(q₁) * dist(q₂))
    winding_angle += asin(sintheta)
  end
  num_winding = trunc(Int, round(winding_angle / (2 * π)))
  println("num turns: ", num_winding)
  tot_winding += num_winding
end


println("Wrapping number is: $tot_winding")


# end

###====================================================================
# experimental stuff
###====================================================================

# plotting projection to y-z plane
N_γ = length(γ_list)
P = zeros(N_γ, 3)
γ_ind = [0.0 for i in 1:N_γ]
count = 1
for ind in 1:length(γs)
  for i in 1:length(γs[ind])
    P[count, :] = myf_dscrt(γs[ind][i])
    γ_ind[count] = count * 1.0 / N_γ
    ## manually shift points because of overlap..
    #if ind == 1
    #  P[count, 2] += 0.005
    #end
    count += 1
  end
end

#plt = Plots.plot(P[:,2], P[:,3], col=γ_ind, seriestype=:scatter)
plt = Plots.scatter(P[:,2], P[:,3], marker_z=γ_ind)


#============================================
# primitive version of algorithm:
# mark vertices whose image under myf is near p₀
# choose ε, let U = N_ε(p₀)
# Let V₁ = { v ∈ tri | myf(v) ∈ U }, V₂ = V(tri) - V₁
# ε must be chosen such that for any triangle Δ of tri,
# if the image myf(Δ) contains p₀,
# then at least one of its vertices is in N_ε(p₀)
# also, ε must be small enough such that
# the "frontier set" of points, those of V₂ that are adjacent to V₁,
# have image in an annulus around U
# then the winding number of these frontier sets around p₀ is well-defined,
# and will agree with the wrapping number

p₀ = sphr_euc(π/2, 0)
ε = 0.4

V₁ = filter(i -> dist(myf(i), p₀) < ε, 1:N)
V₂ = filter(i -> dist(myf(i), p₀) ≥ ε, 1:N)

G₁ = induced_subgraph(G, V₁)
# tuple: graph and mapping to old vertex indices

=#


# =============================================
# ploting image of graph embedded in T² under f

N_γ = length(γ_list)
G = Graph(N_γ)
G_ind_to_Tri = [u for u in keys(γ_list)]
Tri_ind_to_G = Dict(G_ind_to_Tri[i] => i for i in 1:N_γ)
for u ∈ keys(γ_list)
  add_edge!(G, Tri_ind_to_G[u], Tri_ind_to_G[γ_list[u]])
end

P = zeros(N_γ, 3)
for i in 1:N_γ
  println(G_ind_to_Tri[i])
  P[i,:] = myf_dscrt(G_ind_to_Tri[i])
end

plt = Plots.plot(P[:,2], P[:,3], seriestype=:scatter)

=#

#=============================================
# test spherical points functions
P = zeros(121,3)
count = 1
A = LinRange(0,1,11)
product = Base.Iterators.product
for (i,j) in product(A,A)
  P[count,:] = myf(i,j)
  count += 1
end

plt = Plots.plot(P[:,1], P[:,2], P[:,3], seriestype=:scatter)

=#


#===================================================================
# some basic graph operations
# ===

# plot graph
#draw(PNG("mygraph.png", 8cm, 8cm), gplot(G₁, nodelabel=1:3))

G₁  = Graph(3)

add_edge!(G₁, 1, 2)
add_edge!(G₁, 1, 3)
add_edge!(G₁, 2, 3)

A = [
    0 1 1
    1 0 1
    1 1 0
]

G₂ = Graph(A)

@assert G₁ == G₂

=#
