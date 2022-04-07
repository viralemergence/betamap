using SimpleSDMLayers
using StatsPlots
using Statistics
using Colors
using GeoJSON

# Default theme
theme(:bright)
default(; frame=:box, dpi=600)
ispath("figures") || mkpath("figures")

bv_pal = (p0=colorant"#e8e8e8", p1=colorant"#008080", p2=colorant"#fd8d1c")

# Read the risk components stacked layer
lcbd = geotiff(SimpleSDMPredictor, "risk_stack.tif", 1)
phydiv = geotiff(SimpleSDMPredictor, "risk_stack.tif", 2)
sharing = geotiff(SimpleSDMPredictor, "risk_stack.tif", 3)

# Get the landmass data from Natural Earth Data
landmass_url = "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_50m_land.geojson?raw=true"
landmass_file = joinpath("data", "land_areas.json")
if ~isfile(landmass_file)
    download(landmass_url, landmass_file)
end
io = open(landmass_file, "r")
landmass = GeoJSON.read(io)
close(io)

plotbase = plot(landmass, c=colorant"#fafafabb", lc=:grey, lw=0.6)
xaxis!(plotbase, "Longitude", (-180.0, 180.0))
yaxis!(plotbase, "Latitude", (-65.0, 90.0))

# Trivariate risk plot
triplot = deepcopy(plotbase)
plot!(triplot, lcbd, phydiv, sharing; st=:trivariate)
xaxis!(triplot, "Longitude", (-180.0, 180.0))
yaxis!(triplot, "Latitude", (-65.0, 90.0))
trileg = trivariatelegend!(
    lcbd, phydiv, sharing;
    red="High uniqueness",
    green="High diversity",
    blue="High sharing",
    inset=(1, bbox(-0.03, 0.18, 0.35, 0.35, :bottom, :left)),
    subplot=2,
    annotationfontsize=7,
    background_color=colorant"#ffffff00"
)
savefig(joinpath("figures", "risk_trivariate.png"))

# Make an overall risk layer
risk = similar(lcbd)
strength = similar(lcbd)
for k in keys(lcbd)
    r = lcbd[k]
    g = phydiv[k]
    b = sharing[k]
    c = convert(HSV, RGB(r, g, b))
    # Angular distance in the HSV space to pure yellow
    Δ = atan(cos(deg2rad(c.h)), sin(deg2rad(c.h))) - atan(cos(deg2rad(60.0)), sin(deg2rad(60.0)))
    risk[k] = 1.0 - (abs(Δ) / 2π)
    strength[k] = c.v
end

overall_risk = risk * strength

riskmap = deepcopy(plotbase)
plot!(riskmap, rescale(overall_risk, (0, 1)), clim=(0, 1), c=cgrad([bv_pal.p0, bv_pal.p1]))
xaxis!(riskmap, "Longitude", (-180.0, 180.0))
yaxis!(riskmap, "Latitude", (-65.0, 90.0))
savefig(joinpath("figures", "risk_map.png"))

# Get the land cover data
urban = SimpleSDMPredictor(EarthEnv, LandCover, 9; bottom=-56.0)

# Overkill non-zero mean function for the urban density layer
function nonzeromean(layer)
    v = collect(layer)
    filter!(!iszero, v)
    if isempty(v)
        return 0.0
    else
        return mean(v)
    end
end

# Re-project the average LC to the risk map
function coerce(template::TT, source::TS, f) where {TS<:SimpleSDMLayer,TT<:SimpleSDMLayer}
    coerced = similar(template)
    for k in keys(template)
        coerced[k] = f(clip(source, k .- stride(template), k .+ stride(template)))
    end
    return coerced
end

occupation = coerce(risk, urban, nonzeromean)

# Figures for landcover × risk
crossmap = deepcopy(plotbase)
plot!(crossmap, overall_risk, occupation; st=:bivariate, bv_pal..., classes=10, cbar=false)
xaxis!(crossmap, "Longitude", (-180.0, 180.0))
yaxis!(crossmap, "Latitude", (-65.0, 90.0))
bileg = bivariatelegend!(
    overall_risk, occupation;
    inset=(1, bbox(0.05, 0.26, 0.17, 0.17, :bottom, :left)),
    subplot=2,
    annotationfontsize=7,
    background_color=colorant"#ffffff00",
    xlab="Risk",
    ylab="Density",
    bv_pal...
)

savefig(joinpath("figures", "risk_compounded.png"))