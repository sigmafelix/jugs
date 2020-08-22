### Initialize sf
using ArchGDAL

function get_colnames(lyr)
    feat = ArchGDAL.unsafe_getfeature(lyr, 0)
    nfield = ArchGDAL.nfield(feat)
    type_list = []
    type_coln = String[]
    for i ∈ 1:nfield
        #value = ArchGDAL.getfield(feat, i-1)
        valuedef = ArchGDAL.getname(ArchGDAL.getfielddefn(feat, i-1))
        push!(type_coln, valuedef)
        #push!(type_list, Array{typeof(value)}[])
    end
    type_coln
end

function get_datum(lyr, index)
    feat = ArchGDAL.unsafe_getfeature(lyr, index-1)
    nfield = ArchGDAL.nfield(feat)
    type_list = []
    for i ∈ 1:nfield
        value = ArchGDAL.getfield(feat, i-1)
        push!(type_list, value)
    end
    push!(type_list, get_geometry(lyr, index))
    permutedims(type_list)
end


function get_geometry(layer, index = 1)
    geom = ArchGDAL.getgeom(ArchGDAL.unsafe_getfeature(layer, index-1))
    geom
end


function jugs(gdata, num_layer = 1)
    layer = ArchGDAL.getlayer(gdata, num_layer-1)
    nfeat = ArchGDAL.nfeature(layer)
    nfield = ArchGDAL.nfield(ArchGDAL.unsafe_getfeature(layer, 0))
    colns = get_colnames(layer)
    #itype = init_type(layer)
    sf = Array{Any, nfield}
    for i ∈ 1:nfeat
        datum = get_datum(layer, i)
        sf = vcat(sf, datum)
    end
    sf = DataFrame(sf[2:size(sf)[1],:])
    names!(sf, Symbol.(push!(colns, "geom")))
    sf
end


## Compare the original and the genuine
# https://discourse.julialang.org/t/read-geopackage-layers-to-dataframe/35837
using ArchGDAL
const AG = ArchGDAL

dataset = AG.read("/c/Users/sigma/OneDrive/Data/Geo/tl_2017_us_state.shp")

function sf_init(dataset, ilayer = 0)
    layer = AG.getlayer(dataset, ilayer)

    nfeat = AG.nfeature(layer)
    nfield = AG.nfield(layer)

    # prepare Dict with empty vectors of the right type for each field
    d = Dict{String, Vector}()
    featuredefn = AG.layerdefn(layer)
    for field_no in 0:nfield-1
        field = AG.getfielddefn(featuredefn, field_no)
        name = AG.getname(field)
        typ = AG._FIELDTYPE[AG.gettype(field)]
        d[name] = typ[]
    end
    d["geometry"] = AG.IGeometry[]

    # loop over the features to fill the vectors in the Dict
    for fid in 0:nfeat-1
        AG.getfeature(layer, fid) do feature
            for (k, v) in pairs(d)
                if k == "geometry"
                    val = AG.getgeom(feature, 0)
                else
                    val = AG.getfield(feature, k)
                end
                push!(v, val)
            end
        end
    end
    df = DataFrame(d)
    df
end

###
state = AG.read("/c/Users/sigma/OneDrive/Data/Geo/tl_2017_us_state.shp")

state = sf_init(state)

county = AG.read("/c/Users/sigma/OneDrive/Data/Geo/tl_2018_us_county.shp")
county = sf_init(county)

state_sub = state[state[:GEOID] .== "12",:]

AG.intersects(county.geometry, state_sub.geometry)

using DataFrames, Plots
df = DataFrame(d)

plot(df.geometry)
###

function st_intersects(sf1, sf2)
    n1 = nrow(sf1)
    n2 = nrow(sf2)
    intersects = trues(n1, n2)

    for i ∈ 1:n1
        for j ∈ 1:n2
            ints = ArchGDAL.intersects(sf1.geometry[i], sf2.geometry[j])
            intersects[i,j] = ints
        end
    end
    intersects
end 

function contiguity_matrix(sf, Queen = true)
    nn = nrow(sf)
    contig = trues(nn, nn)
    contig0 = trues(nn, nn)

    if Queen
        for i ∈ 1:nn
            for j ∈ 1:nn
                contig[i,j] = ArchGDAL.touches(sf.geometry[i], sf.geometry[j])
            end
        end
        contig
    else # todo
        for i ∈ 1:nn
            for j ∈ 1:nn
                contig[i,j] = ArchGDAL.within(sf.geometry[i], sf.geometry[j])
                #contig0[i,j] = ArchGDAL.touches(sf.geometry[i], sf.geometry[j])
            end
        end
        #contig = contig - contig0
        contig
    end
end

check1 = contiguity_matrix(state, true)
check0 = contiguity_matrix(state, false)
sum(check0)
check2 = contiguity_matrix(county, true)



st_intersects(state, county)


function MoranI(q, w)
    denom = q'*w*q
    nom = q'*q
    nf = size(w)[1]
    expectation = -1/size(w)[1]
    wcross2 = 0.5 * sum((w+w') .^ 2)
    wcrosssum = nf * (sum(sum(w+w', dims = 2) .^ 2))
    wijsq = sum(w)^2
    qbar = sum(q)/nf
    s1 = (nf^2 - 3*nf + 3) * wcross2 - wcrosssum + (3*wijsq)
    s2 = ((1/nf)*sum((q.-qbar).^4)) / ((1/nf)*sum((q.-qbar).^2)^2) 
    s3 = wcross2 - (2*nf*wcross2) + 6*wijsq
    vI = ((nf * s1) - (s2 * s3)) / ((nf-1)*(nf-2)*(nf-3)*wijsq)
    I = denom/nom
    
    (I, expectation, vI)
end

MoranI(state.AWATER, check1)
MoranI(county.AWATER, check2)
