module ASEAtoms

#
# Imports
#

using PyCall: PyObject, PyNULL, @py_str,
              pyisinstance
using AtomsBase: AtomsBase, FastSystem, AbstractSystem, Atom,
                 AtomView, Periodic, DirichletZero,
                 atomic_number, bounding_box, periodicity
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
        _checkisaseatoms(o)
        new(o)
    end
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

Convert an ASE Atoms object, `atoms` to a `FastSystem`.
"""
function ase_to_atomicsystem(atoms::PyObject; assert=true)
    _checkisaseatoms(atoms)
    box = eachrow(atoms.cell[]) .* 1u"Å"
    boundary_conditions = ifelse.(atoms.pbc, Ref(Periodic()), Ref(DirichletZero()))
    pos = atoms.positions * 1u"Å"
    sym = atoms.get_chemical_symbols()

    atoms = Atom.(Symbol.(sym), eachrow(pos))

    FastSystem(atoms, box, boundary_conditions)
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

    return py"Atoms"(;numbers, positions, cell, pbc)
end
aseatoms(sys::ASESystem) = sys.o

end
