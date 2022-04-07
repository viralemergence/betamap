using SimpleSDMLayers
using Statistics
using RData
using DataFrames
using DataFramesMeta
using CSV
using ProgressMeter
using GLM
using Phylo

# Get the data if they are not downloaded
if !isfile(joinpath("data", "binary.csv"))
    download("https://raw.githubusercontent.com/viralemergence/Fresnel/master/BinaryWebsite.csv", joinpath("data", "binary.csv"))
end

int = DataFrame(CSV.File(joinpath("data", "binary.csv")))
select!(int, [:Sp => :species, Symbol("Training data") => :train, Symbol("New data") => :new, :Ensemble => :model])
species = @linq int |>
                where((:train .== "Reported") .| (:new .== "New data"))
#where((:train .== "Reported").|(:new .== "New data").|(:model .== "Suspected"))
hosts = sort(species.species)

# List of species in virionette
ispath("rasters") || mkpath("rasters")

# Get the predictors
ranges = Dict{String,SimpleSDMPredictor}()
for sp in hosts
    @info sp
    fname = joinpath("rasters", replace(sp, " " => "_") * ".tif")
    if !isfile(fname)
        query = `gdal_rasterize -l "MAMMALS_TERRESTRIAL_ONLY" -a presence rangemaps/MAMMALS_TERRESTRIAL_ONLY.shp $(fname) -where "binomial LIKE '$(sp)'" -ts 1200, 600`
        run(query)
    end
    mp = geotiff(SimpleSDMResponse, fname)
    replace!(mp.grid, 0 => nothing)
    ranges[sp] = copy(mp)
    mp = nothing
    GC.gc()
end

# Richness map
richness = similar(ranges[first(hosts)])
for sp in keys(ranges)
    r = ranges[sp]
    for k in keys(r)
        if isnothing(richness[k])
            richness[k] = Float64(1.0)
        else
            richness[k] += Float64(1.0)
        end
    end
end


# Read the sharing matrix
viralsharing = DataFrame(CSV.File(joinpath("data", "PredictedNetwork.csv")))
rename!(viralsharing, "Column1" => "Ref")

# Set of species
Base.zero(::Type{Set{String}}) = Set{String}()
Base.zero(::Type{Vector{String}}) = String[]
spsets = similar(richness, Vector{String})

@showprogress for lat in latitudes(richness)
    for lon in longitudes(richness)
        if !isnothing(richness[lon, lat])
            sps = String[]
            for sp in keys(ranges)
                if !isnothing(ranges[sp][lon, lat])
                    push!(sps, sp)
                end
            end
            spsets[lon, lat] = sps
        end
    end
end

sharing = similar(richness, Float64)
msharing = similar(richness, Float64)
ssharing = similar(richness, Float64)
esharing = similar(richness, Float64)

function pielou(x)
    p = x ./ sum(x)
    return -sum(p .* log.(p)) / log(length(p))
end

@showprogress for lon in longitudes(sharing)
    for lat in latitudes(sharing)
        if !isnothing(sharing[lon, lat])
            if length(spsets[lon, lat]) == 1
                sharing[lon, lat] = 0.0
                msharing[lon, lat] = 0.0
                ssharing[lon, lat] = 0.0
                esharing[lon, lat] = 1.0
            else
                sp_codes = sort([replace(s, " " => "_") for s in spsets[lon, lat]])
                filter!(s -> s in names(viralsharing), sp_codes)
                t1 = select(viralsharing, [:Ref, Symbol.(sp_codes)...])
                sprows = findall(vec(mapslices(any, t1.Ref .== permutedims(sp_codes); dims=2)))
                S = Matrix(t1[sprows, :][1:end, 2:end])
                if size(S, 1) < 2
                    sharing[lon, lat] = 0.0
                    msharing[lon, lat] = 0.0
                    ssharing[lon, lat] = 0.0
                    esharing[lon, lat] = 1.0
                else
                    # Upper triangle
                    v = [S[i, j] for i in 1:(size(S, 1)-1) for j in (i+1):size(S, 2)]
                    sharing[lon, lat] = sum(v)
                    msharing[lon, lat] = mean(v)
                    ssharing[lon, lat] = std(v)
                    esharing[lon, lat] = pielou(v)
                end
            end
        end
    end
end

# Get the (untitled) Upham tree
tree = open(parsenexus, joinpath("data", "upham_tree.nex"))["*UNTITLED"]

# We create a layer to store the phylogenetic diversity
dphyl = similar(sharing)

@showprogress for lat in latitudes(dphyl)
    for lon in longitudes(dphyl)
        if !isnothing(spsets[lon, lat])
            pool = [replace(sp, " " => "_") for sp in spsets[lon, lat]]
            filter!(sp -> hasnode(tree, sp), pool)
            try
                dphyl[lon, lat] = sum([br.length for br in unique(vcat([branchhistory(tree, sp) for sp in pool]...))])
            catch e
                dphyl[lon, lat] = 0.0
            end
        end
    end
end

mpdist = minimum(filter(x -> x > 0.0, filter(!isnothing, unique(values(dphyl.grid)))))
for k in keys(dphyl)
    if iszero(dphyl[k])
        dphyl[k] = mpdist
    end
end

# LCBD
function LCBD(Y)
    S = (Y .- mean(Y; dims=1)) .^ 2.0
    SStotal = sum(S)
    BDtotal = SStotal / (size(Y, 1) - 1)
    SSj = sum(S; dims=1)
    SCBDj = SSj ./ SStotal
    SSi = sum(S; dims=2)
    LCBDi = SSi ./ SStotal
    return LCBDi, SCBDj, BDtotal
end

function hellinger(Y::Matrix{T}) where {T<:Number}
    yi = sum(Y; dims=2)
    return sqrt.(Y ./ yi)
end

patches = findall(!isnothing, richness.grid)

Y = zeros(Int64, (length(patches), length(ranges)))
for (i, tax) in enumerate(keys(ranges))
    sp_occ = findall(!isnothing, ranges[tax].grid)
    Y[indexin(sp_occ, patches), i] .= 1
end
Y = Y[:, findall(vec(sum(Y; dims=1) .> 0))]

lcbd = similar(richness)
lcbd.grid[patches] = LCBD(hellinger(Y))[1]

LC = rescale(lcbd, (0, 1))
PH = rescale(dphyl, (0, 1))
VC = rescale(msharing, (0, 1))

geotiff("risk_stack.tif",
    [
        LC,
        PH,
        VC
    ]
)