using YAXArrays
using Zarr
using DimensionalData
using GLMakie

tmpf = () -> begin
  era5url = "https://s3.bgc-jena.mpg.de:9000/deepextremes/v3/ERA5Cube.zarr"
  ds = open_dataset(era5url)
  sub = ds[Ti=DateTime(1998,1,1)..DateTime(2022,12,31), longitude=0..14.76,latitude=30.1..60]
  savedataset(sub,path="tmp/era5.zarr")
end

sub = open_dataset("tmp/era5.zarr");
c = Cube(sub)
  #   longitude Sampled{Float32} 0.0f0:0.25f0:14.75f0 ForwardOrdered Regular Points,
  # → latitude  Sampled{Float32} 60.0f0:-0.25f0:30.25f0 ReverseOrdered Regular Points,
  # ↗ Ti        Sampled{DateTime} [1998-01-01T00:00:00, …, 2022-12-31T00:00:00] ForwardOrdered Irregular Points,
  # ⬔ Variable  Categorical{String} ["tp", "pet", "t2mmin", "t2mmax", "t2m"] Unordered

using Statistics 
indims = InDims("Ti")
outdims = OutDims()
function apply_median(xout, xin)
    x = filter(!ismissing, xin)
    x = filter(!isnan,x)
    xout[] = isempty(x) ? missing : median(x)
end
medians = mapCube(apply_median, c[Variable=Where(contains("t2m"))];indims, outdims)
size(medians)

fig, ax, heat = heatmap(DimArray(medians[Variable=At("t2m")]))

# Lazy traditional map applied to each value
# Problem: does change the eltype to not allow missing
medians_kelvin = map(medians) do x 
  x + 273.15
end
fig, ax, heat = heatmap(DimArray(medians_kelvin[Variable=At("t2m")]))

# Moving window
function meanfilter(xout, xin)
  if ismissing(xin[2,2])
      xout .= missing
  else
  xout .= mean(skipmissing(xin))
  end
end

indims = InDims(MovingWindow("longitude", 1,1),MovingWindow("latitude", 1, 1), window_oob_value=NaN)
filteredmeans = mapCube(meanfilter, medians_kelvin, indims=indims, outdims=OutDims())

# Define new ouput dimensions
# extract a single time series near some lat/long to test
# lookup to get values of axis
using Dates
pet = c[Variable=At("pet"),
        lon=Near(11.3464),lat=Near(46.4946)]
fig,ax, pl = lines(lookup(pet, Ti).data,pet.data)
fig

using SignalDecomposition: SignalDecomposition as SD #decompose as dc
dates = lookup(pet, Ti)
stlres = SD.decompose(dates,pet.data[:], SD.TimeAnomaly()) # returns tuple

fig,ax, p = plot(dates, stlres[1]);
ax2 = Axis(fig[2,1])
plot!(ax2,stlres[2])
fig

# apply over full array
#    supply the dates as additional argument to the mapped function
import Logging
Logging.disable_logging(Logging.Info)
indims = InDims("time")
outdims = OutDims(dims(c,Ti),Dim{:Scale}(["Seasonal", "anomalies"]), 
                    path = "tmp/decomposed.zarr",backend=:zarr, overwrite=true)
function decompose_TS(xout, xin, dates)
    any(isnan,xin) && return xout .= missing
    seas, anomaly = SD.decompose(dates,xin,SD.TimeAnomaly())
    xout[:,1] = seas
    xout[:,2] = anomaly
end

Logging.disable_logging(Logging.Warn)

@time dec = mapCube(decompose_TS, 
    c[Variable=At("t2m")],
    lookup(c, Ti), # addtional variable to decompose_TS
        #lon=Near(11.3464),lat=Near(46.4946)],
    indims = indims,
    outdims = outdims)

