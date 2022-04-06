using SimpleSDMLayers
using StatsPlots
using Statistics
using Colors

# Default theme
theme(:bright)
default(; frame=:box, dpi=600)


# Read the risk components stacked layer
lcbd = geotiff(SimpleSDMPredictor, "risk_stack.tif", 1)
phydiv = geotiff(SimpleSDMPredictor, "risk_stack.tif", 2)
sharing = geotiff(SimpleSDMPredictor, "risk_stack.tif", 3)

# TODO get the world shapefile from github

# Trivariate risk plot
triplot = plot(lcbd, phydiv, sharing; st=:trivariate)
trileg = trivariatelegend(
    lcbd, phydiv, sharing;
    red="Unique compositions",
    green="Phylogenetic diversity",
    blue="High viral sharing"
)

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
    risk[k] = 1.0-(abs(Δ)/2π)
    strength[k] = c.v
end

plot(risk)
plot(strength)

overall_risk = risk * strength

plot(risk, strength; st=:bivariate)
plot(overall_risk)

# Get the land cover data
urban = SimpleSDMPredictor(EarthEnv, LandCover, 9; bottom=-56.0)

# Re-project the average LC to the risk map
function coerce(template::TT, source::TS, f) where {TS <: SimpleSDMLayer, TT <: SimpleSDMLayer}
    coerced = similar(template)
    for k in keys(template)
        coerced[k] = f(clip(source, k .- stride(template), k .+ stride(template)))
    end
    return coerced
end

occupation = coerce(risk, urban, maximum)

# Figures for landcover × risk
plot(overall_risk, occupation; st=:bivariate)