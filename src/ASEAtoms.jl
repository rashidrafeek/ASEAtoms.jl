module ASEAtoms

#
# Imports
#

using PyCall: PyObject, PyNULL, @py_str,
              pyisinstance
using AtomsBase: AtomsBase, AbstractSystem, FlexibleSystem, 
                 Atom, AtomView, Periodic, DirichletZero,
                 atomic_number, bounding_box, periodicity,
                 species_type
using Unitful: @u_str, ustrip

#
# Exports
# 

export ASESystem, aseread, asewrite

#
# PyCall functions
#

const aseread = PyNULL()
const asewrite = PyNULL()

function __init__()
    py"""
    from ase.io import read, write
    from ase import Atoms
    """

    copy!(aseread, py"read")
    copy!(asewrite, py"write")
end

#
# Types
#

struct ASESystem <: AbstractSystem{3}
    o::PyObject
    function ASESystem(o::PyObject; assert=true)
        if assert
            _checkisaseatoms(o)
        end
        new(o)
    end
end

function ASESystem(sys::AbstractSystem)
    pyatoms = aseatoms(sys)

    return ASESystem(pyatoms)
end

function _checkisaseatoms(o::PyObject)
    if !pyisinstance(o, py"Atoms")
        error("The given object is not an instance of ASE Atoms object")
    end
end

Base.length(sys::ASESystem) = length(sys.o)
Base.size(sys::ASESystem) = length(sys.o)

Base.position(sys::ASESystem) = collect(eachrow(sys.o.positions*u"Å"))
AtomsBase.bounding_box(sys::ASESystem) = collect(eachrow(sys.o.cell[]u"Å"))
function AtomsBase.boundary_conditions(sys::ASESystem)
    ifelse.(sys.o.pbc, Ref(Periodic()), Ref(DirichletZero()))
end

AtomsBase.atomic_symbol(sys::ASESystem) = Symbol.(sys.o.get_chemical_symbols())
AtomsBase.atomic_number(sys::ASESystem) = sys.o.get_atomic_numbers()
AtomsBase.atomic_mass(sys::ASESystem) = sys.o.get_masses()u"u"

function Base.checkbounds(sys::ASESystem, i::Integer)
    natoms = length(sys)
    if !(1 <= i <= natoms) 
        throw(BoundsError(sys, i))
    end
end

function Base.getindex(sys::ASESystem, index::Int)
    @boundscheck checkbounds(sys, index)
    AtomView(sys, index)
end

function Base.position(sys::ASESystem, i)
    @boundscheck checkbounds(sys, i)
    return sys.o[i].position*u"Å"
end
function AtomsBase.atomic_symbol(sys::ASESystem, i)
    @boundscheck checkbounds(sys, i)
    return Symbol(sys.o[i].symbol)
end
function AtomsBase.atomic_number(sys::ASESystem, i)
    @boundscheck checkbounds(sys, i)
    return convert(Int, sys.o[i].number)
end
function AtomsBase.atomic_mass(sys::ASESystem, i)
    @boundscheck checkbounds(sys, i)
    return sys.o[i].mass*u"u"
end

#
# Functions
#

"""
    ase_to_atomicsystem(atoms::PyObject; periodic=true)

Convert an ASE Atoms object, `atoms` to a `FlexibleSystem`.
"""
function ase_to_atomicsystem(atoms::PyObject; assert=true)
    if assert
        _checkisaseatoms(atoms)
    end
    box = eachrow(atoms.cell[]) .* 1u"Å"
    boundary_conditions = ifelse.(atoms.pbc, Ref(Periodic()), Ref(DirichletZero()))
    pos = atoms.positions * 1u"Å"
    sym = atoms.get_chemical_symbols()
    charges = atoms.get_initial_charges()

    if all(==(0), charges)
        atoms = Atom.(Symbol.(sym), eachrow(pos))
    else
        f = (sym, pos, charge) -> Atom(Symbol(sym), pos; charge=charge)
        atoms = f.(sym, eachrow(pos), charges)
    end

    FlexibleSystem(atoms, box, boundary_conditions)
end
function ase_to_atomicsystem(atomslist::Vector{PyObject}; assert=true)
    return ase_to_atomicsystem.(atomslist; assert)
end

"""
    ase_to_atomicsystem(path::String, args...; kwargs...)

Read structure file at `path` using `aseread` and convert it to `FlexibleSystem`. 
`args` and `kwargs` are passed on to `aseread`.
"""
function ase_to_atomicsystem(path::String, args...; kwargs...)
    atoms = aseread(path, args...; kwargs...)
    return ase_to_atomicsystem(atoms; assert=false)
end

"""
    aseatoms(sys::AbstractSystem)

Convert an `AbstractSystem` to an ASE Atoms object.
"""
function aseatoms(sys::AbstractSystem)
    numbers = atomic_number(sys)
    cell = map(bounding_box(sys)) do boxvec
        ustrip.(u"Å", boxvec)
    end
    positions = map(position(sys)) do posvec
        ustrip.(u"Å", posvec)
    end
    pbc = periodicity(sys)

    atoms = py"Atoms"(;numbers, positions, cell, pbc)

    if :data in fieldnames(typeof(sys))
        particles = collect(sys)
        datnames = collect(keys(particles[1].data))
        if :charge in datnames
            charges = getindex.(getproperty.(particles, :data), :charge)
            atoms.set_initial_charges(charges)
        end
    end
	
	return atoms
end
aseatoms(sys::ASESystem) = sys.o

end
