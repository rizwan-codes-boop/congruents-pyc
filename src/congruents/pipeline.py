"""Source solve -> attenuation -> reference observer component files, in Python."""
from dataclasses import dataclass
import json
from pathlib import Path
import numpy as np
from .solver import solve
from .attenuation import EBLTable
from .observer import observer_sed, distance_factor_cm2
from .spectrum import integrate_spectrum

COMPONENT_FILES = dict(zip(
    ("pion","IC_primary_disc","IC_secondary_disc","BS_primary_disc","BS_secondary_disc",
     "SY_primary_disc","SY_secondary_disc","IC_primary_halo","IC_secondary_halo",
     "SY_primary_halo","SY_secondary_halo","pion_full_calorimetry","neutrino","free_free"),
    ("spec_pi.txt","spec_IC_1_z1.txt","spec_IC_2_z1.txt","spec_BS_1_z1.txt","spec_BS_2_z1.txt",
     "spec_SY_1_z1.txt","spec_SY_2_z1.txt","spec_IC_1_z2.txt","spec_IC_2_z2.txt",
     "spec_SY_1_z2.txt","spec_SY_2_z2.txt","spec_pi_fcal1.txt","spec_nu.txt","spec_FF.txt")))


@dataclass(frozen=True)
class RunResult:
    """Complete source/observer result, with explicit optical depths and provenance."""
    source: object
    components: dict
    tau_internal: np.ndarray
    tau_ebl: np.ndarray
    distance_factor: np.ndarray
    gamma_luminosity: np.ndarray
    metadata: dict

    @property
    def total_photon_sed(self):
        """Sum physical photon components, excluding neutrinos/fcal1 diagnostic."""
        result = sum(self.components[name] for name in (*list(COMPONENT_FILES)[:11],"free_free"))
        result.flags.writeable = False
        return result

    def export(self, directory):
        """Create a NEW directory; never overwrite an existing model run."""
        folder = Path(directory)
        folder.mkdir(parents=True,exist_ok=False)
        for name,file in COMPONENT_FILES.items():
            np.savetxt(folder/file,self.components[name],fmt="%.6e",
                       header="E^2 dN/dE [GeV cmm2 sm1]",comments="")
        for file,array,label in (
                ("tau_gg.txt",self.tau_internal,"internal optical depth"),
                ("tau_EBL.txt",self.tau_ebl,"EBL optical depth"),
                ("tau_ff.txt",self.source.emission["tau_free_free"],"free-free optical depth"),
                ("distmod.txt",self.distance_factor,"dist_modulus__cmm2"),
                ("L_gamma.txt",self.gamma_luminosity,"L_gamma__GeVsm1"),
                ("E_gam.txt",self.source.photon_energy_gev,"GeV")):
            np.savetxt(folder/file,array.reshape(1,-1) if array.ndim==1 else array,fmt="%.6e",header=label,comments="")
        # NPZ retains full-precision source arrays; text files intentionally use
        # six-decimal scientific notation for compatibility with the reference writer.
        self.source.save(folder/"source.npz")
        np.savetxt(folder/"spec_total_photons.txt",self.total_photon_sed,fmt="%.6e",
                   header="Sum of photon components [GeV cm^-2 s^-1]",comments="")
        (folder/"tau_loss").mkdir()
        for name,array in self.source.diagnostics.items():
            target = folder/"tau_loss"/(name+".txt") if name.startswith("tau_loss") else folder/(name+".txt")
            np.savetxt(target,array,fmt="%.6e",header="" if name.startswith("CR_specs") else name,comments="")
        np.savetxt(folder/"T_CR.txt",self.source.kinetic_energy_gev[None,:],fmt="%.6e",header="T_CR__GeV",comments="")
        np.savetxt(folder/"fcal.txt",self.source.transport["fcal"],fmt="%.6e",header="fcal",comments="")
        (folder/"metadata.json").write_text(json.dumps(self.metadata,indent=2)+"\n")


def run(preparation, grid=None, threads=1, *, ebl_path=None, legacy_table_precision=True):
    """Same component treatment as the reference, including its neutrino writer.

    The original writer applies photon optical depths to all 14 components,
    including neutrinos. Preserve this behavior; it is not a physical endorsement.
    """
    path = Path(ebl_path) if ebl_path else Path(__file__).resolve().parents[2]/"input/tau_Eg_z_Franceschini.txt"
    ebl = EBLTable(path)
    with preparation._lock:
        preparation._ensure_open()
        source = solve(preparation,grid,threads,legacy_table_precision=legacy_table_precision,diagnostics=True)
        preparation.report("Applying EBL attenuation and observer-frame conversion")
        energy = source.photon_energy_gev
        z = np.array([g.redshift for g in preparation.galaxies])
        internal = source.tau_internal
        external = ebl.optical_depth(energy,z)
    # Optical depths retain source-grid indexing while observer_sed resamples
    # the emission. This reproduces the writer convention recorded in metadata.
    components = {name:observer_sed(energy,source.emission[name],z,
                    tau_internal=internal,tau_ebl=external) for name in COMPONENT_FILES}
    # Reference L_gamma excludes FF, neutrinos and the full-calorimetry diagnostic.
    total = sum(source.emission[name] for name in list(COMPONENT_FILES)[:11])*np.exp(-internal)
    luminosity = np.array([integrate_spectrum(energy,row,.1,100.) for row in total])
    distance = distance_factor_cm2(z)
    for a in (internal,external,luminosity,distance):
        a.flags.writeable = False
    return RunResult(source,components,internal,external,distance,luminosity,
        {**source.metadata,"ebl_path":str(path.resolve()),"ebl_sha256":ebl.sha256,
         "observer_units":"GeV cm^-2 s^-1","attenuation_indexing":"improvised legacy code writer",
         "complete_legacy_diagnostics":True})
