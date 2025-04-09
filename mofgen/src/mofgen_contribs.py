from __future__ import annotations

import os
from typing import TYPE_CHECKING

from pydantic import BaseModel, Field
from tempfile import TemporaryDirectory


_installed_extra = {"mofid": True}
try:
    from mofid.run_mofid import cif2mofid
except ImportError:
    _installed_extra["mofid"] = False

from zeopp import _run_zeopp_assessment

from pymatgen.core import Structure as PmgStructure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

if TYPE_CHECKING:
    from typing_extensions import Self

class MofIdEntry(BaseModel):

    Smiles : str | None = None
    Topology : str | None = None
    SmilesLinkers : list[str] | None = None
    SmilesNodes : list[str] | None = None
    MofKey : str | None = None
    MofId : str | None = None

    @classmethod
    def from_structure(cls, structure, **kwargs) -> Self:
        old_cwd = os.getcwd()
        mof_id_output = None
        with TemporaryDirectory() as tmp_dir:
            os.chdir(tmp_dir)
            try:
                structure.to("temp.cif")
                mof_id_output = cif2mofid("temp.cif",**kwargs)
            except Exception as exc:
                print(exc)
                
        os.chdir(old_cwd)
        if mof_id_output is None:
            return cls()

        remap = {
            "Smiles": "smiles",
            "Topology": "topology",
            "SmilesLinkers": "smiles_linkers",
            "SmilesNodes": "smiles_nodes",
            "MofKey": "mofkey",
            "MofId": "mofid",
        }
        return cls(**{k: mof_id_output[v] for k, v in remap.items()})


class ZeoPPEntry(BaseModel):

    Sorbate : str | None = None
    Pld : float | None = None
    Lcd : float | None = None
    Poav : float | None = Field(None, description = "POAV in cm^3/g")
    PoavVolumetric : float | None = Field(None, description = "POAV in Å^3")
    PoavVolumeFraction : float | None = None
    Ponav : float | None = Field(None, description = "PONAV in cm^3/g")
    PonavVolumetric : float | None = Field(None, description = "PONAV in Å^3")
    PonavVolumeFraction : float | None = None
    Density : float | None = None
    UnitCellVolume : float | None = None

    @classmethod
    def from_structure(cls, structure, **zeopp_kwargs) -> list[Self]:

        output = None
        old_cwd = os.getcwd()
        with TemporaryDirectory() as tmp_dir:
            os.chdir(tmp_dir)
            try:
                output = _run_zeopp_assessment(
                    structure,
                    **zeopp_kwargs
                )
            except Exception as exc:
                print(exc)
            os.chdir(old_cwd)

        if output is None:
            return None

        mapping = {
            "Density": "Density",
            "UnitCellVolume": "Unitcell_volume",
            "Pld": "PLD",
            "Lcd": "LCD",
            "Poav": "POAV_cm^3/g",
            "PoavVolumetric": "POAV_A^3",
            "PoavVolumeFraction": "POAV_Volume_fraction",
            "Ponav": "PONAV_A^3",
            "PonavVolumetric" : "PONAV_cm^3/g",
            "PonavVolumeFraction": "PONAV_Volume_fraction",
        }
            
        return [
            ZeoPPEntry(
                Sorbate = sorbate,
                **{k: entry[v] for k, v in mapping.items()}
            )
            for sorbate, entry in output.items()
            if sorbate != "is_mof"
        ]

class MofGenEntry(BaseModel):

    Identifier : str | None = None
    Structure : PmgStructure | None = Field(None)
    Density : float | None = None
    Volume : float | None = None
    NumSites : int | None = None
    Formula : str | None = None
    FormulaReduced : str | None = None
    ChemicalSystem : str | None = None
    NumElements : int | None = None
    SpaceGroupNumber : int | None = None
    SpaceGroupSymbol : str | None = None
    Method : str | None = None
    Energy : float | None = None
    FormationEnergyPerAtom : float | None = None
    BandGap : float | None = None
    MofId : MofIdEntry | None = None
    ZeoPlusPlus : list[ZeoPPEntry] | None = None

    @classmethod
    def from_structure(cls, structure: Structure,**kwargs) -> Self:

        sg = {}
        try:
            sga = SpacegroupAnalyzer(structure)
            sg = {
                "SpaceGroupNumber": sga.get_space_group_number(),
                "SpaceGroupSymbol": sga.get_space_group_symbol(),
            }
        except Exception as exc:
            print(exc)
        
        config = {
            "Structure": structure,
            "Density": structure.density,
            "Volume": structure.volume,
            "NumSites": structure.num_sites,
            "Formula": structure.composition.formula,
            "FormulaReduced": structure.composition.reduced_formula,
            "ChemicalSystem" : "-".join(sorted([str(ele) for ele in structure.composition.elements])),
            "NumElements": len(structure.composition),
            "MofId": MofIdEntry.from_structure(structure),
            "ZeoPlusPlus": ZeoPPEntry.from_structure(structure),
            **sg
        }
        return cls(**config, **kwargs)