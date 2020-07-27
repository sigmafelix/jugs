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
