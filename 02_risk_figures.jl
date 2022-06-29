using SimpleSDMLayers
using StatsPlots
using Statistics
using Colors, ColorSchemes, ColorBlendModes
using GeoJSON
using LinearAlgebra

# Default theme
theme(:bright)
default(; frame=:box, dpi=600, size=(700, 400), background_color=colorant"#ffffff00")
ispath("figures") || mkpath("figures")

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

# Components of risk
lcbdplot = deepcopy(plotbase)
plot!(lcbdplot, lcbd, c=:viridis)
xaxis!(lcbdplot, "Longitude", (-180.0, 180.0))
yaxis!(lcbdplot, "Latitude", (-65.0, 90.0))
savefig(joinpath("figures", "layer_uniqueness.png"))

sharingplot = deepcopy(plotbase)
plot!(sharingplot, sharing, c=:viridis)
xaxis!(sharingplot, "Longitude", (-180.0, 180.0))
yaxis!(sharingplot, "Latitude", (-65.0, 90.0))
savefig(joinpath("figures", "layer_sharing.png"))

diversityplot = deepcopy(plotbase)
plot!(diversityplot, phydiv, c=:viridis)
xaxis!(diversityplot, "Longitude", (-180.0, 180.0))
yaxis!(diversityplot, "Latitude", (-65.0, 90.0))
savefig(joinpath("figures", "layer_diversity.png"))

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

riskgrad = SimpleSDMLayers.bivariates.blue_red.grad1
plot!(riskmap, rescale(overall_risk, (0, 1)), clim=(0, 1), c=cgrad(riskgrad))
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

blendinfo = (
    SimpleSDMLayers.bivariates.blue_red,
    blendmode=ColorBlendModes.BlendColorBurn,
    classes=10
)

# Figures for landcover × risk
crossmap = deepcopy(plotbase)
plot!(crossmap, overall_risk, occupation; st=:bivariate, blendinfo..., cbar=false)
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
    blendinfo...
)
savefig(joinpath("figures", "risk_compounded.png"))

# Richness map for the bats
richness = geotiff(SimpleSDMPredictor, "richness.tif")
crossmap = deepcopy(plotbase)
plot!(crossmap, richness; cbar=false, c=cgrad(ColorSchemes.Hokusai3))
xaxis!(crossmap, "Longitude", (-180.0, 180.0))
yaxis!(crossmap, "Latitude", (-65.0, 90.0))
savefig(joinpath("figures", "bat_richness.png"))

# Distinctiveness
batdist = geotiff(SimpleSDMPredictor, "biogeo/BatDistinctiveness.tif")
virdist = geotiff(SimpleSDMPredictor, "biogeo/VirusDistinctiveness.tif")
batdist = mask(virdist, batdist)
virdist = mask(batdist, virdist)
blendinfo = (
    SimpleSDMLayers.bivariates.blue_red...,
    blendmode=ColorBlendModes.BlendMultiply,
    classes=4,
    quantiles=true
)
crossmap = deepcopy(plotbase)
plot!(crossmap, batdist, virdist; st=:bivariate, blendinfo..., cbar=false)
xaxis!(crossmap, "Longitude", (-180.0, 180.0))
yaxis!(crossmap, "Latitude", (-65.0, 90.0))
bileg = bivariatelegend!(
    batdist, virdist;
    inset=(1, bbox(0.05, 0.26, 0.17, 0.17, :bottom, :left)),
    subplot=2,
    annotationfontsize=7,
    background_color=colorant"#ffffff00",
    xlab="Bats",
    ylab="Viruses",
    blendinfo...
)
savefig(joinpath("figures", "evo_distinctiveness.png"))

# Biogeo regions
blendinfo = (
    grad1=ColorSchemes.BrBG_6,
    grad2=ColorSchemes.PiYG_6,
    blendmode=ColorBlendModes.BlendSoftLight,
    classes=6,
    quantiles=true
)

vpc1 = geotiff(SimpleSDMPredictor, "biogeo/VirusPC1.tif")
vpc2 = geotiff(SimpleSDMPredictor, "biogeo/VirusPC2.tif")
bpc1 = geotiff(SimpleSDMPredictor, "biogeo/BatPC1.tif")
bpc2 = geotiff(SimpleSDMPredictor, "biogeo/BatPC2.tif")
# This sucks but generates the same pixels everywhere
bpc1 = mask(vpc1, bpc1)
bpc2 = mask(vpc1, bpc2)
vpc2 = mask(vpc1, vpc2)
bpc1 = mask(bpc1, bpc1)
bpc2 = mask(bpc1, bpc2)
vpc2 = mask(bpc1, vpc2)
vpc1 = mask(bpc1, vpc1)

crossmap = deepcopy(plotbase)
bivariate!(crossmap, vpc1, vpc2; blendinfo..., cbar=false)
xaxis!(crossmap, "Longitude", (-180.0, 180.0))
yaxis!(crossmap, "Latitude", (-65.0, 90.0))
bileg = bivariatelegend!(
    vpc1, vpc2;
    inset=(1, bbox(0.05, 0.26, 0.17, 0.17, :bottom, :left)),
    subplot=2,
    annotationfontsize=7,
    background_color=colorant"#ffffff00",
    xlab="PCoA1",
    ylab="PCoA2",
    blendinfo...
)
savefig(joinpath("figures", "virus_biogeo.png"))

crossmap = deepcopy(plotbase)
bivariate!(crossmap, bpc1, bpc2; blendinfo..., cbar=false)
xaxis!(crossmap, "Longitude", (-180.0, 180.0))
yaxis!(crossmap, "Latitude", (-65.0, 90.0))
bileg = bivariatelegend!(
    bpc1, bpc2;
    inset=(1, bbox(0.05, 0.26, 0.17, 0.17, :bottom, :left)),
    subplot=2,
    annotationfontsize=7,
    background_color=colorant"#ffffff00",
    xlab="PCoA1",
    ylab="PCoA2",
    blendinfo...
)
savefig(joinpath("figures", "bat_biogeo.png"))

# Use convert to append the files
run(`convert figures/bat_richness.png figures/evo_distinctiveness.png -append figures/combined_richness.png`)
run(`convert figures/bat_biogeo.png figures/virus_biogeo.png -append figures/combined_biogeo.png`)