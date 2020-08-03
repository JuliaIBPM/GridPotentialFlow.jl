import PotentialFlow

function compute_ẋ!(ẋ, x, t)
    plate, ambient_sys = x
    motion = ẋ[1]
    motion.ċ, motion.c̈, motion.α̇ = motion.kin(t)

    Plates.enforce_no_flow_through!(plate, motion, ambient_sys, t)
    reset_velocity!(ẋ, x)
    self_induce_velocity!(ẋ, x, t)
    display(ẋ)
end

function shed_new_vorticity!(blobs, plate, motion, t, lesp = 0.0, tesp = 0.0)
    z₊ = (blobs[end-1].z + 2plate.zs[end])/3
    z₋ = (blobs[end].z + 2plate.zs[1])/3

    blob₊ = PotentialFlow.Vortex.Blob(z₊, 1.0, δ)
    blob₋ = PotentialFlow.Vortex.Blob(z₋, 1.0, δ)
    PotentialFlow.Plates.enforce_no_flow_through!(plate, motion, blobs, t)

    Γ₊, Γ₋, _, _ = PotentialFlow.Plates.vorticity_flux!(plate, blob₊, blob₋, t, lesp, tesp);

    push!(blobs, PotentialFlow.Vortex.Blob(z₊, Γ₊, blobs[1].δ), PotentialFlow.Vortex.Blob(z₋, Γ₋, blobs[1].δ))
end

L = 1.0   # chord length
N = 128   # number of plate control points (distributed along a extrema Chebyshev grid)

ċ = -1.0L  # translation velocity
α = -π/3   # angle of attack

plate = PotentialFlow.Plate(N, L, zero(ComplexF64), α)
motion = PotentialFlow.Plates.RigidBodyMotion(ċ, 0.0);

Δt = 5e-2

δ = 0.01
lesp = 0.0
tesp = 0.0

Δz₀ = im*3Δt*exp(im*plate.α) # vectors perpendicular to the plate
z₋, z₊ = plate.zs[[1,N]] # LE and TE

blobs = PotentialFlow.Vortex.Blob.(Δz₀ .+ [z₊, z₋], 1.0, δ) # First two point vortices are placed close to the LE and TE with unit strength

PotentialFlow.Plates.enforce_no_flow_through!(plate, motion, (), 0)
Γ₊, Γ₋, _, _ = PotentialFlow.Plates.vorticity_flux!(plate, blobs[1], blobs[2], 0.0, lesp, tesp); # Determine strength of first two vortices

blobs = PotentialFlow.Vortex.Blob.(Δz₀ .+ [z₊, z₋], [Γ₊, Γ₋], δ) # Create first two point vortices now with calculated strengths

sys₀ = (plate, blobs)

sys = deepcopy(sys₀)
sys₊ = deepcopy(sys₀) # Used for storage during time-marching
ẋs = (motion, PotentialFlow.allocate_velocity(blobs))

forces = ComplexF64[];
