
from __future__ import annotations

import os
import multiprocessing
from shutil import which
import subprocess
from typing import TYPE_CHECKING

from pymatgen.core import Structure

from jobflow import job

if TYPE_CHECKING:
    from typing import Any

class ZeoPlusPlus:
    """TODO add docstr
    """

    def __init__(
        self,
        cif_path : str,
        zeopp_path : str | None = None,
        working_dir : str | None = None,
        sorbates : list[str] | str = ["N2", "CO2", "H2O"],
    ) -> None:
        
        self._cif_path = cif_path
        self.cif_name = os.path.basename(cif_path.split(".cif")[0])
        self.zeopp_path = zeopp_path or which("zeo++") or os.environ.get("ZEO_PATH")

        if isinstance(sorbates,str):
            sorbates = [sorbates]
        self.sorbates = sorbates

        self.working_dir = working_dir or os.path.dirname(cif_path)
        self._zeopp_path = zeopp_path

    @classmethod
    def from_structure(
        cls,
        structure : Structure,
        cif_path : str,
        zeopp_path : str | None = None,
        working_dir : str | None = None,
        sorbates : list[str] | str = ["N2", "CO2", "H2O"],
    ):
        structure.to(cif_path)
        return cls(
            cif_path = cif_path,
            zeopp_path = zeopp_path,
            working_dir = working_dir,
            sorbates = sorbates
        )

    def run(
        self,
        zeopp_args : list[str] | None = None,
        nproc : int = 1
    ):
        nproc = min(nproc, len(self.sorbates))

        if nproc > 1:
            sorbate_batches = [[] for _ in range(nproc)]
            iproc = 0
            for sorbate in self.sorbates:
                sorbate_batches[iproc].append(sorbate)
                iproc = (iproc + 1) % nproc
    
            manager = multiprocessing.Manager()
            output_file_path = manager.dict()
            output = manager.dict()
    
            procs = []
            for iproc in range(nproc):
                proc = multiprocessing.Process(
                    target = self._run_zeopp_many,
                    kwargs = {
                        "sorbates": sorbate_batches[iproc],
                        "file_paths_shared": output_file_path,
                        "output_shared": output,
                        "zeopp_args": zeopp_args,
                    },
                )
    
                procs.append(proc)
                proc.start()

            for proc in procs:
                proc.join()

            # convert from multiprocessing manager shared-memory object to real dict
            self.output_file_path = dict(output_file_path)
            self.output = dict(output)
            
        else:
            self.output_file_path = {}
            self.output = {}
            self._run_zeopp_many(
                sorbates = self.sorbates,
                file_paths_shared = self.output_file_path,
                output_shared = self.output,
                zeopp_args = zeopp_args,
            )

    def _run_zeopp_many(self, sorbates: list[str], file_paths_shared : dict, output_shared : dict, zeopp_args : list[str] | None = None):
        for sorbate in sorbates:
            self._run_zeopp_single(sorbate, file_paths_shared, output_shared, zeopp_args = zeopp_args)
        
    def _run_zeopp_single(self, sorbate : str, file_paths_shared : dict, output_shared : dict, zeopp_args : list[str] | None = None):

        radius_sorbate = self.get_sorbate_radius(sorbate)

        # check if zeopp_args contains -res
        output = {}
        parse_func = None
        flag_to_func = {"res": self._parse_res, "volpo": self._parse_volpo}

        zeopp_args = zeopp_args or ["-ha", "-volpo", str(radius_sorbate),  str(radius_sorbate), "50000"]

        for flag, _func in flag_to_func.items():
            if f"-{flag}" in zeopp_args:
                output_file_path = os.path.join(self.working_dir,self.cif_name) + f"_{sorbate}.{flag}"
                parse_func = _func

        zeopp_args = [self.zeopp_path] + zeopp_args + [output_file_path, self._cif_path]

        with subprocess.Popen(
            zeopp_args,
            stdout=subprocess.PIPE,
            stdin=subprocess.PIPE,
            close_fds=True,
        ) as proc:
            stdout, stderr = proc.communicate()
            if proc.returncode != 0:
                raise RuntimeError(
                    f"{self._zeopp_path} exit code: {proc.returncode}, error: {stderr!s}."
                    f"\nstdout: {stdout!s}. Please check your zeo++ installation."
                )
                    
        output = parse_func(output_file_path)

        if output == {}:
            raise ValueError(f"zeopp_args must contain either -res or -volpo, not {zeopp_args}")

        try:
            output["structure"] = Structure.from_file(self._cif_path)
        except Exception as exc:
            output["structure"] = f"Exception: {exc}"

        file_paths_shared[sorbate] = output_file_path
        output_shared[sorbate] = output
                        
    @staticmethod
    def _parse_volpo(volpo_path : str) -> dict[str,Any]:

        with open(volpo_path,"r") as f:
            data = f.read().split("\n")

        output = {}
        for line in data:
            if "PROBE_OCCUPIABLE" in line:
                # PROBE_OCCUPIABLE_VOL_CALC: not sure what the strs here mean
                # PROBE_OCCUPIABLE___RESULT: this line just copies data from previous lines
                continue

            read_value = False
            for val in line.split():

                if ":" in val:
                    key = val.split(":")[0]
                    read_value = True
                elif read_value:
                    try:
                        val = float(val)
                    except ValueError:
                        pass
                    output[key] = val
                    read_value = False

        return output
    
    @staticmethod
    def _parse_res(res_path: str) -> dict[str, Any]:
        with open(res_path, "r") as f:
            data = f.read().split()
        return {"LCD": float(data[1]), f"PLD": float(data[2])}
    
    @staticmethod
    def get_sorbate_radius(sorbate):
        # sorbate kinetic diameters in Angstrom
        # https://doi.org/10.1039/B802426J
        kinetic_diameter = {
            # noble gases
            "He": 2.551,
            "Ne": 2.82,
            "Ar": 3.542,
            "Kr": 3.655,
            "Xe": 4.047,
            # diatomic gases
            "H2": 2.8585,
            "D2": 2.8585,
            "N2": 3.72,
            "O2": 3.467,
            "Cl2": 4.217,
            "Br2": 4.296,
            # oxides
            "CO": 3.69,
            "CO2": 3.3,
            "NO": 3.492,
            "N2O": 3.838,
            "SO2": 4.112,
            "COS": 4.130,
            # others
            "H2O": 2.641,
            "CH4": 3.758,
            "NH3": 3.62,
            "H2S": 3.623,
        }

        # check sorbate is present and return radius
        try:
            return kinetic_diameter[sorbate] * 0.5
        except Exception as e:
            print("Unknown sorbate " + sorbate + ".")
        
def _run_zeopp_assessment(
    structure : Structure | str,
    zeopp_path : str | None = None,
    working_dir : str | None = None,
    sorbates : list[str] | str = ["N2", "CO2", "H2O"],
    cif_name : str | None = None,
    nproc : int = 1
) -> dict[str, Any]:
    
    if isinstance(structure,str) and os.path.isfile(structure):
        maker = ZeoPlusPlus(
            cif_path=structure,
            zeopp_path = zeopp_path,
            working_dir= working_dir,
            sorbates=sorbates
        )

    elif isinstance(structure, Structure):
        cif_name = cif_name or "structure.cif"
        maker = ZeoPlusPlus.from_structure(
            structure = structure,
            cif_path=cif_name,
            zeopp_path = zeopp_path,
            working_dir= working_dir,
            sorbates=sorbates
        )
    
    output = {sorbate: {} for sorbate in sorbates} 
    for args in [[], ["-ha", "-res"]]:
        maker.run(zeopp_args = args,nproc = nproc)
        for sorbate in sorbates:
            output[sorbate].update(maker.output[sorbate])

    output["is_mof"] = False
    if all( k in output["N2"] for k in ("PLD", "POAV_A^3", "PONAV_A^3", "POAV_Volume_fraction", "PONAV_Volume_fraction")):
        output["is_mof"] = (
            output["N2"]["PLD"] > 2.5 and
            output["N2"]["POAV_Volume_fraction"] > 0.3 and
            output["N2"]["POAV_A^3"] > output["N2"]["PONAV_A^3"] and
            output["N2"]["POAV_Volume_fraction"] > output["N2"]["PONAV_Volume_fraction"]
        )

    return output

@job
def run_zeopp_assessment(
    structure : Structure | str,
    **kwargs
) -> dict[str, Any]:
    return _run_zeopp_assessment(
        structure,**kwargs
    )