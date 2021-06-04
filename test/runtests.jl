using GridPotentialFlow
using Test
using Literate

const GROUP = get(ENV, "GROUP", "All")

ENV["GKSwstype"] = "nul"

notebookdir = "../examples"
docdir = "../docs/src/manual"
litdir = "./literate"
testdir = "./"
solverdir = "./solver"

for (root, dirs, files) in walkdir(litdir)
  for file in files
    if endswith(file,".jl")
      (GROUP == "All" || GROUP == "Notebooks") && Literate.notebook(joinpath(root, file),notebookdir)
      (GROUP == "All" || GROUP == "Documentation") && Literate.markdown(joinpath(root, file),docdir)
      (GROUP == "All" || GROUP == "Scripts") && (Literate.script(joinpath(root, file),testdir); include(file))
    end
  end
end
