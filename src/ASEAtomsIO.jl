module ASEAtoms

#
# Imports
#

using PyCall: PyObject, PyNULL, @py_str,
              pyisinstance
using AtomsBase: AtomsBase, FastSystem, AbstractSystem, Atom,
                 AtomView, Periodic, DirichletZero
using Unitful: @u_str

#
# Exports
# 

export ASEAtoms, aseread, asewrite

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

struct ASEAtoms <: AbstractSystem{3}
    o::PyObject
    function ASEAtoms(o::PyObject; assert=true)
        _checkisaseatoms(o)
        new(o)
    end
end

function _checkisaseatoms(o::PyObject)
    if !pyisinstance(o, py"Atoms")
        error("The given object is not an instance of ASE Atoms object")
    end
end

Base.length(sys::ASEAtoms) = length(sys.o)
Base.size(sys::ASEAtoms) = length(sys.o)

Base.position(sys::ASEAtoms) = collect(eachrow(sys.o.positions*u"Å"))
AtomsBase.bounding_box(sys::ASEAtoms) = collect(eachrow(sys.o.cell[]u"Å"))
function AtomsBase.boundary_conditions(sys::ASEAtoms)
    ifelse.(sys.o.pbc, Ref(Periodic()), Ref(DirichletZero()))
end

AtomsBase.atomic_symbol(sys::ASEAtoms) = Symbol.(sys.o.get_chemical_symbols())
AtomsBase.atomic_number(sys::ASEAtoms) = sys.o.get_atomic_numbers()
AtomsBase.atomic_mass(sys::ASEAtoms) = sys.o.get_masses()u"u"

function Base.checkbounds(sys::ASEAtoms, i::Integer)
    natoms = length(sys)
    if !(1 <= i <= natoms) 
        error("Invalid index $i for `ASEAtoms` with indices 1:$(natoms)")
    end
end

function Base.getindex(sys::ASEAtoms, index::Int)
    @boundscheck checkbounds(sys, index)
    AtomView(sys, index)
end

function AtomsBase.atomic_symbol(sys::ASEAtoms, i)
    @boundscheck checkbounds(sys, i)
    return Symbol(sys.o[i].symbol)
end
function AtomsBase.atomic_number(sys::ASEAtoms, i)
    @boundscheck checkbounds(sys, i)
    return convert(Int, sys.o[i].number)
end
function AtomsBase.atomic_mass(sys::ASEAtoms, i)
    @boundscheck checkbounds(sys, i)
    return sys.o[i].mass*u"u"
end

#
# Functions
#

"""
    ase_to_atomicsystem(atoms; periodic=true)

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

end
