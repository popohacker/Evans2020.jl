push!(LOAD_PATH,"../src/")
using Documenter, Evans2020

makedocs(modules = [Evans2020], sitename = "Evans2020.jl")

deploydocs(repo = "github.com/popohacker/Evans2020.jl.git", devbranch = "main")
