module Controller
using DataFrames
using JLD2
using FileIO
export trans, load_db, save_db, add_row, getUserInput

const defaultSim = (
    radius = 10.,
    density = 0.09,
    intensity = 0.,
    parameter = "density",
    intervals = Dict(:start => 0., :step => 0.1, :stop => 0.5),
    dims = Dict(:lx => 2000, :ly => 1000),
    numberTrials = 1200,
    nCollections = 2,
)

trans(df::DataFrame) = DataFrame([[names(df)]; collect.(eachrow(df))], [:column; Symbol.(axes(df, 1))])

struct DataBase
    filename::String
    table::DataFrame
end

function load_db(filename)
    if isfile(filename)
        return DataBase(filename, load(filename, "database"))
    else
        db = DataBase(filename, DataFrame([defaultSim]))
        empty!(db.table)
        return db
    end
end

save_db(db::DataBase) = jldsave(db.filename; database = db)

function add_row(db::DataBase)
    d = Dict{Symbol,Any}()
    for (k,v) in pairs(defaultSim)

        if isa(v, Dict)
            tmp = Dict{Symbol,Any}()
            println("\tenter values for $k:")
            for (sk, sv) in v
                subinput = getUserInput(
                    typeof(sv),
                    "\t\t> enter value for $sk ($(typeof(sv))): "
                )
                tmp[sk] = subinput
            end
            # input = namedtuple(tmp)
            input = tmp
        else
            input = getUserInput(
                typeof(v),
                "\tenter value for $k ($(typeof(v))): "
            )
        end
        d[k] = input
    end
    push!(db.table, d);
end


# * Helpers
function getUserInput(T = String, msg = "")
    print("$msg ")
    if T == String
        return readline()
    else
        try
            return parse(T, readline())
        catch
            println("Sorry, I could not interpret your answer. Please try again")
            getUserInput(T, msg)
        end
    end
end

dictkeys(d::Dict) = (collect(keys(d))...,)
dictvalues(d::Dict) = (collect(values(d))...,)

function namedtuple(d::Dict{Symbol,T}) where {T}
    NamedTuple{dictkeys(d)}(dictvalues(d))
end

function namedtuple(d::Dict{String,T}) where {T}
    NamedTuple{Symbol.(dictkeys(d))}(dictvalues(d))
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    using .Controller

    ensure_load(dbpath, db) = ismissing(db) ? load_db(dbpath) : db
    

    # get path
    print("\nPlease enter database path:\n> ")
    dbpath = readline()

    let db = missing
        db = ensure_load(dbpath, db)

        while true
            print("\nwhat would you like to do: exit, reload, save, add, remove, show\n\t > ")

            input = readline()
            if contains(input, "show")
                db.table |> show
            elseif contains(input, "add")
                db |> add_row
                db.table |> show
            elseif contains(input, "remove")
                if size(db.table) |> first == 0
                    print("\t skipping, empty table.")
                    continue
                end
                idx = parse(Int, readline())
                deleteat!(db.table, idx)
                db.table |> show
            elseif contains(input, "save")
                save_db(db)
                db.table |> show
            elseif contains(input, "reload")
                println("\t ! confirm: yes / no")
                reload = contains(lowercase(readline()), "yes")
                if reload
                    db = load_db(dbpath)
                    db.table |> show
                end
            else
                println("> save? (yes / no)")
                confirm = contains(lowercase(readline()), "yes")
                if confirm
                    save_db(db)
                end
                break
            end
        end
    end
end


