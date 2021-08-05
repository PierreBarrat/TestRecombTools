"""
    add_auspice_json_attribute(infilejson::String, outfile::String, attr_dict::Dict, attr_name::String, visible_name::String=attr_name, coloring_type::String="categorical")
"""
function add_auspice_json_attribute(infilejson::String, outfile::String, attr_dict::Dict, attr_name::String, visible_name::String=attr_name, coloring_type::String="categorical")
    json = JSON3.read(read(infilejson, String), Dict)
    add_auspice_json_attribute!(json, attr_dict, attr_name, visible_name, coloring_type)
    open(outfile, "w") do f
        JSON3.pretty(f, JSON3.write(json))
    end
end

"""
    add_auspice_json_attribute(f::Function, infilejson::String, outfile::String, tree::Tree, key, visible_name=key, coloring_type::String="categorical")

Add `f(n.data.dat[key])` to JSON file `infilejson` for all nodes `n` of `tree`. Write output to `outfile`.
"""
function add_auspice_json_attribute(f::Function, infilejson::String, outfile::String, tree::Tree, key, visible_name=string(key), coloring_type::String="categorical")
    json = JSON3.read(read(infilejson, String), Dict)
    add_auspice_json_attribute!(f, json, tree, key, visible_name, coloring_type)
    open(outfile, "w") do f
        JSON3.pretty(f, JSON3.write(json))
    end
end

"""
    add_auspice_json_attribute!(f::Function, json::Dict, tree::Tree, key, visible_name=key, coloring_type::String="categorical")

Add `f(n.data.dat[key])` to `json["tree"]` for all nodes `n` of `tree`.
"""
function add_auspice_json_attribute!(f::Function, json::Dict, tree::Tree, key, visible_name=key, coloring_type::String="categorical")
    attr_dict = Dict()
    for (label, n) in tree.lnodes
        if haskey(n.data.dat, key)
            attr_dict[label] = f(n.data.dat[key])
        end
    end
    add_auspice_json_attribute!(json, attr_dict, string(visible_name), string(visible_name), coloring_type)
end

function add_auspice_json_attribute!(json::Dict, attr_dict::Dict, attr_name::String, visible_name::String, coloring_type::String)
    _add_auspice_json_attribute!(json["tree"], attr_dict, attr_name)
    _add_auspice_json_coloring_attribute!(json, attr_name, visible_name, coloring_type)
end

function _add_auspice_json_attribute!(json::Dict, attr_dict::Dict, attr_name::String)
	if haskey(attr_dict, json["name"])
		json["node_attrs"][attr_name] = Dict("value"=>attr_dict[json["name"]])
	end
	if haskey(json, "children")
		for c in json["children"]
			_add_auspice_json_attribute!(c, attr_dict, attr_name)
		end
	end
end

function _add_auspice_json_coloring_attribute!(json::Dict, key::String, title::String, type::String)
	push!(json["meta"]["colorings"], Dict{String, Any}("key"=>key, "title"=>title, "type"=>type))
end
