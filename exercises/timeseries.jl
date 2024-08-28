# read a dataset
using DataFrames, CSV
using Dates
using Statistics
using GLMakie
ds = CSV.read("examples/data/jena.csv", DataFrame)

# Mean seasonal cycle
x = ds.Variable_t2mmin;
t = ds.Ti
mean_annual_cycle = function(t,x)
  # yr = year.(t)
  # yrs = unique(yr)
  doi = dayofyear.(t) 
  dois = unique(doi)
  #doi_i = first(dois)
  mac = map(dois) do doi_i
      xi = @view x[doi .== doi_i]
      m = mean(xi)
  end
  DataFrame(doi=dois; mac)
end

mac = mean_annual_cycle(t, x)
Kelvin0 = 273.15
plot(mac.doi,mac.mac .- Kelvin0)

# Dryness Index moving average tp - pet
using RollingFunctions
ds = CSV.read("examples/data/sahara.csv", DataFrame)
# compute daily sums
dayi = year.(t) .* 1000 .+ dayofyear.(t) 
week_id = collect(zip(year.(t), dayofyear.(t) .÷ 7)) 
#week.(t) does not work because first days nad last days in 2012 get 52 assigned

ds.ai = ds.Variable_tp .- ds.Variable_pet

# Fabian's example of foldl
function streak(oldstate,newvalue)
  maxstreak, currentstreak = oldstate
  if newvalue #We extend the streak by 1
      currentstreak += 1
      maxstreak = max(currentstreak, maxstreak)
  else
      currentstreak = 0
  end
  return (maxstreak,currentstreak)
end
x = rand(Bool,1000)
foldl(streak,x,init=(0,0))

using OrderedCollections

"""
Call FUN(startindex,endindex) for each run of consecutive equal values in ids.
Returns a `NamedTuple(id, value)`.
"""
function aggregate_runs(FUN, ids)
  isempty(ids) && return((;id=[],value=[]))
  f1 = FUN(firstindex(ids),firstindex(ids))
  id = eltype(ids)[]
  value = typeof(f1)[]
  apply_fun = (oldstate,i) -> begin
    (startindex, endindex) = oldstate
    if (i == :_after_end || (ids[i] != ids[startindex])) 
        #ids[startindex] == (1998, 52) && Main.@infiltrate_main
        #DataFrame(i = [startindex,endindex,i], ids = ids[[startindex,endindex,i]])
        push!(id, ids[startindex])
        push!(value, FUN(startindex,endindex))
        (startindex, endindex) = (i,i)
    else
      endindex = i  # remember because cannot use i-1 after index changed
    end
    return (startindex, endindex)
  end
  # need to run once after the end to push last run
  ind_it = Iterators.flatten((eachindex(week_id), (:_after_end,)))
  foldl(apply_fun,ind_it,init=(firstindex(ids),firstindex(ids)))
  (;id, value)
end

tmpf = () -> begin
  # test aggregate_runs
  start_end_id = (i_start,i_end) -> (week_id[i_start], i_end - i_start + 1, week_id[i_end])
  tmp = aggregate_runs(start_end_id, week_id);
  #tmp = map(x -> first(x,5),tmp)
  all(tmp.id .== first.(tmp.value) .== last.(tmp.value))
  unique(week_id) == tmp.id
  ls = map(x -> x[2], tmp.value)
  i_single = findall(<=(1), ls) # all at least two days (might be a few at end of year)
  length(i_single) == 0
  # tmp.value[last(i_single,2)]
  # tmp.value[first(i_single,2)]
  # tmp.value[findall(ls .<= 1)]
  # unique(last.(week_id))
  # uweek_id = unique(week_id)
  # tmp2 = findall(x -> !(x ∈ uweek_id), tmp.id)
  # tmp2b = findall(x -> !(x ∈ tmp.id), uweek_id)
  # tmp3 = counter(tmp.id)
  # tmp4 = findall(collect(Base.values(tmp3)) .!= 1)
  # collect(keys(tmp3))[tmp4]
  # tmp_ds = DataFrame(Ti=ds.Ti, week_id = week_id, doi = dayofyear.(ds.Ti), i =eachindex(week_id))
  # tmp_ds2 = subset(tmp_ds, :week_id => ByRow(==((2012, 52))))
  # tmp_ds2 = tmp_ds[i_single,:]
  # tmp_ds2 = subset(tmp_ds, :week_id => ByRow(==((1998, 52))))
  #
  # empty ids
  tmp2 = aggregate_runs(start_end_id, week_id[1:0]);
  tmp2 == (id = Any[], value = Any[])
end



# compute aridity index for subset
compute_cumai = (i_start,i_end) -> begin
  l = i_end - i_start + 1
  ai = (sum(@view(ds.Variable_tp[i_start:i_end])) + sum(@view(ds.Variable_pet[i_start:i_end])))/l
end

ai_weekly = aggregate_runs(compute_cumai, week_id)
plot(last.(ai_weekly.id), ai_weekly.value)

