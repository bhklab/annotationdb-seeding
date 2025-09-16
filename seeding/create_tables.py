from typing import List
from typing import Optional
from sqlalchemy import String, Float, Integer, ForeignKey, LargeBinary
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column
from datetime import datetime


class Base(DeclarativeBase):
    pass


class Pubchem(Base):
    __tablename__ = "pubchem"

    cid: Mapped[int] = mapped_column(primary_key=True)
    molecular_formula: Mapped[str] = mapped_column(String(300))
    molecular_weight: Mapped[str] = mapped_column(String(50))
    smiles: Mapped[str] = mapped_column(String(2000))
    # canonical_smiles: Mapped[str] = mapped_column(String(2000))
    # isomeric_smiles: Mapped[str] = mapped_column(String(2000))
    connectivity_smiles: Mapped[str] = mapped_column(String(2000))
    inchi: Mapped[str] = mapped_column(String(1000))
    inchikey: Mapped[str] = mapped_column(String(28))
    iupac_name: Mapped[str] = mapped_column(String(1000))
    xlogp: Mapped[float] = mapped_column(Float)
    exact_mass: Mapped[str] = mapped_column(String(150))
    monoisotopic_mass: Mapped[str] = mapped_column(String(150))
    tpsa: Mapped[float] = mapped_column(Float)
    complexity: Mapped[int] = mapped_column(Integer)
    charge: Mapped[int] = mapped_column(Integer)
    h_bond_donor_count: Mapped[int] = mapped_column(Integer)
    h_bond_acceptor_count: Mapped[int] = mapped_column(Integer)
    rotatable_bond_count: Mapped[int] = mapped_column(Integer)
    heavy_atom_count: Mapped[int] = mapped_column(Integer)
    isotope_atom_count: Mapped[int] = mapped_column(Integer)
    atom_stereo_count: Mapped[int] = mapped_column(Integer)
    defined_atom_stereo_count: Mapped[int] = mapped_column(Integer)
    undefined_atom_stereo_count: Mapped[int] = mapped_column(Integer)
    bond_stereo_count: Mapped[int] = mapped_column(Integer)
    defined_bond_stereo_count: Mapped[int] = mapped_column(Integer)
    undefined_bond_stereo_count: Mapped[int] = mapped_column(Integer)
    covalent_unit_count: Mapped[int] = mapped_column(Integer)
    volume_3d: Mapped[float] = mapped_column(Float)
    x_steric_quadrupole_3d: Mapped[float] = mapped_column(Float)
    y_steric_quadrupole_3d: Mapped[float] = mapped_column(Float)
    z_steric_quadrupole_3d: Mapped[float] = mapped_column(Float)
    feature_count_3d: Mapped[int] = mapped_column(Integer)
    feature_acceptor_count_3d: Mapped[int] = mapped_column(Integer)
    feature_donor_count_3d: Mapped[int] = mapped_column(Integer)
    feature_anion_count_3d: Mapped[int] = mapped_column(Integer)
    feature_cation_count_3d: Mapped[int] = mapped_column(Integer)
    feature_ring_count_3d: Mapped[int] = mapped_column(Integer)
    feature_hydrophobe_count_3d: Mapped[int] = mapped_column(Integer)
    conformer_model_rmsd_3d: Mapped[float] = mapped_column(Float)
    effective_rotor_count_3d: Mapped[float] = mapped_column(Float)
    conformer_count_3d: Mapped[int] = mapped_column(Integer)
    fingerprint_2d: Mapped[bytes] = mapped_column(LargeBinary)
    title: Mapped[str] = mapped_column(String(100))
    patent_count: Mapped[int] = mapped_column(Integer)
    patent_family_count: Mapped[int] = mapped_column(Integer)
    literature_count: Mapped[int] = mapped_column(Integer)
    annotation_types: Mapped[str] = mapped_column(String(100))
    annotation_type_count: Mapped[int] = mapped_column(Integer)
    name: Mapped[str] = mapped_column(String(300))
    # chembl_id: Mapped[int] = mapped_column(ForeignKey("chembl.cid"))
    chembl_id: Mapped[int] = mapped_column(String(200))
    synonyms: Mapped[str] = mapped_column(String(1000))
    # date_added: Mapped[datetime] = mapped_column()


# class Chembl(Base):
#     __tablename__ = "chembl"
#     cid: Mapped[int] = mapped_column(primary_key=True)
#     pubchem_cid: Mapped[int] = mapped_column(ForeignKey("pubchem.cid"))


# class Cellosaurus(Base):
#     __tablename__ = "cellosaurus"
#     id: Mapped[int] = mapped_column(primary_key=True)
