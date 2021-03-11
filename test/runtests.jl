using GridPotentialFlow
using Test
using Literate

const GROUP = get(ENV, "GROUP", "Notebooks")

notebookdir = "../examples"
docdir = "../docs/src/manual"
litdir = "./literate"
solverdir = "./solver"

if GROUP == "All" || GROUP == "Solver"
  @testset "System vectors" begin
      include("solver/systemvectors.jl")
  end
  # @testset "Unregularized potential flow systems" begin
  #     include("unregularizedsystems.jl")
  # end
  @testset "Regularized potential flow systems" begin
      include("solver/regularizedsystems.jl")
  end
end

if GROUP == "All" || GROUP == "Notebooks"
  for (root, dirs, files) in walkdir(litdir)
    for file in files
      endswith(file,".jl") && Literate.notebook(joinpath(root, file),notebookdir)
    end
  end
end

if GROUP == "Documentation"
  for (root, dirs, files) in walkdir(litdir)
    for file in files
      endswith(file,".jl") && Literate.markdown(joinpath(root, file),docdir)
    end
  end
end
