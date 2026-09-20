"""Python allocation/orchestration around complete OpenMP galaxy workers."""
import ctypes as ct
import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
import numpy as np
from . import constants, grids
from .provenance import python_source_provenance
from .preparation import Preparation, PROPERTY_NAMES, Table, TableData
from ._native import LOCK, IP, RunInput, TableInput, pointer, descriptor, load as _load, check, threads_value


@dataclass(frozen=True)
class SolverGrid:
    """Particle/photon sample counts and finite-volume electron cell count.

    These resolutions are independent of the shared lookup-table Grid.
    Small values are useful for smoke tests, not evidence of convergence.
    """
    cosmic_rays: int = 1000
    photons: int = 500
    cells: int = 500

    def __post_init__(self):
        """Enforce integer particle, photon and solver resolutions within supported limits."""
        for n, maximum in ((self.cosmic_rays,4096),(self.photons,4096),(self.cells,500)):
            if isinstance(n,bool) or not isinstance(n,int):
                raise TypeError("Solver grid sizes must be integers")
            if not 4 <= n <= maximum:
                raise ValueError(f"Solver grid size must be in [4, {maximum}]")


@dataclass(frozen=True)
class SolverResult:
    """Source-frame arrays and provenance returned by the galaxy workers.

    Energy axes are in GeV; electron_energy_gev is total electron energy.
    Electron populations are GeV^-1 and emission rates are GeV^-1 s^-1,
    except tau_free_free, which is dimensionless. Disc emission already
    includes free-free absorption; tau_internal is supplied separately.
    Component arrays have galaxy rows and energy columns.
    """
    kinetic_energy_gev: np.ndarray
    electron_energy_gev: np.ndarray
    photon_energy_gev: np.ndarray
    transport: dict
    electrons: dict
    emission: dict
    metadata: dict
    diagnostics: dict = field(default_factory=dict)
    tau_internal: np.ndarray = None

    def save(self,path):
        """Save component arrays and JSON metadata in a compressed NPZ archive.

        Results returned by solve contain numerical arrays loadable without pickle.
        This method does not guard against overwriting an existing destination.
        """
        arrays = {"kinetic_energy_gev":self.kinetic_energy_gev,
                  "electron_energy_gev":self.electron_energy_gev,
                  "photon_energy_gev":self.photon_energy_gev,
                  "metadata_json":np.array(json.dumps(self.metadata,sort_keys=True)),
                  "tau_internal":self.tau_internal}
        for group in ("transport","electrons","emission","diagnostics"):
            arrays.update({f"{group}__{name}":v for name,v in getattr(self,group).items()})
        np.savez_compressed(path,**arrays)


TRANSPORT = ("fcal","proton_diffusion_cm2_s","disc_diffusion_cm2_s","halo_diffusion_cm2_s",
             "primary_injection_gev_s","secondary_injection_gev_s","proton_steady_state_gev")
ELECTRONS = ("primary_disc","secondary_disc","primary_halo","secondary_halo")
EMISSION = ("IC_primary_disc","IC_secondary_disc","BS_primary_disc","BS_secondary_disc",
            "SY_primary_disc","SY_secondary_disc","IC_primary_halo","IC_secondary_halo",
            "SY_primary_halo","SY_secondary_halo","free_free","tau_free_free","pion","pion_full_calorimetry","neutrino")


def solve(preparation,grid=None,threads=1,library=None,*,legacy_table_precision=True,diagnostics=False):
    """Calculate each galaxy entirely in a native worker after serial table generation.

    Only galaxy indices are parallel: each worker performs disc then halo
    assembly/solve, emission, absorption and optional diagnostics sequentially.
    A failed worker invalidates the result; partially calculated rows never escape.
    Lookup tables remain Python-generated, cached, owned and read-only.
    """
    if not isinstance(preparation,Preparation):
        raise TypeError("Expected Preparation")
    grid = SolverGrid() if grid is None else grid
    if not isinstance(grid,SolverGrid):
        raise TypeError("Expected SolverGrid")
    if not isinstance(legacy_table_precision,bool) or not isinstance(diagnostics,bool):
        raise TypeError("Precision and diagnostics flags must be bool")
    threads_value(threads)
    with preparation._lock, LOCK:
        preparation._ensure_open()
        provenance = python_source_provenance()
        lib = _load(library if library is not None else preparation._context.library)
        n,ne,np_ = len(preparation.galaxies),grid.cosmic_rays,grid.photons
        kinetic = grids.log_grid(ne,1e-3,1e8)
        electron = grids.log_grid(ne,1e-3+constants.ME,1e8+constants.ME)
        photon = grids.log_grid(np_,1e-16,1e8)
        rows = np.ascontiguousarray([g.as_row() for g in preparation.galaxies],dtype=float)
        props = np.ascontiguousarray([[p[k] for k in PROPERTY_NAMES] for p in preparation.properties])
        transport,cp = np.zeros((n,7,ne)),np.empty(n)
        statuses = np.zeros(n,dtype=np.int32)
        preparation.report(f"Calculating transport for {n} catalogue cases with up to {min(n,threads)} OpenMP workers")
        check(lib.cg_transport(threads,n,ne,pointer(rows),pointer(props),pointer(kinetic),
              pointer(transport),pointer(cp),statuses.ctypes.data_as(IP)),statuses)
        rounded = []
        # Borrow a table descriptor, retaining any rounded copy until the native call finishes.
        def table(kind,field="3000",temperature=0.):
            original = preparation.table(kind,field,temperature)
            if legacy_table_precision and kind != "gamma":
                # Apply the six-decimal scientific-notation round trip used by table-precision
                # mode.
                def rounding(a):
                    a=np.asarray(a)
                    return np.array([float(format(float(v),".6e")) for v in a.ravel()]).reshape(a.shape)
                copy=Table(TableData(rounding(original.x),rounding(original.y),rounding(original.values)),original.metadata)
                rounded.append(copy)
                return descriptor(copy,float(format(temperature,".6e")))
            return descriptor(original,temperature)
        try:
            cmb,fir = preparation.temperature_grid("CMB"),preparation.temperature_grid("FIR")
            fields=[(f,0.) for f in ("3000","4000","7500","UV")]+[("CMB",t) for t in cmb]+[("FIR",t) for t in fir]
            preparation.report(f"Preparing {2*len(fields)+2} lookup tables (load cache or generate)")
            families=[(TableInput*len(fields))(*(table(kind,f,t) for f,t in fields)) for kind in ("emission","gamma")]
            bs,sy=table("bs"),table("sy")
            electrons,emission,internal = np.zeros((n,4,ne)),np.zeros((n,15,np_)),np.zeros((n,np_))
            # Diagnostic buffers are separate from transport and never alias inputs.
            loss,critical,budget,radio,escape=(np.zeros(shape) for shape in
                ((n,12,ne),(n,10),(n,16),(n,5),(n,2,ne)))
            args=RunInput(n,ne,np_,grid.cells,len(cmb),len(fir),
                *(pointer(a) for a in (rows,props,kinetic,electron,photon,cp,transport,electrons,emission,internal)),
                families[0],families[1],ct.pointer(bs),ct.pointer(sy),preparation.config[4],preparation.config[5],int(diagnostics),
                *(pointer(a) if diagnostics else None for a in (loss,critical,budget,radio,escape)))
            preparation.report(f"Lookup tables ready; starting galaxy spectra, matrix solves and diagnostics with up to {min(n,threads)} OpenMP workers")
            check(lib.cg_spectra(threads,ct.byref(args),statuses.ctypes.data_as(IP)),statuses)
            preparation.report("Native galaxy calculations complete")
        finally:
            for item in rounded: item.close()
        outputs={}
        if diagnostics:
            preparation.report("Preparing radiation-density report and output arrays")
            # Packing/file layout is Python-owned; all per-galaxy diagnostic
            # physics above is computed by the native spectrum workers.
            outputs={"gal_data":props[:,:9].copy(),
                "Urad_Ub":np.array([preparation.radiation(i,[])['urad_ub_ev_cm3'] for i in range(n)]),
                "E_loss_nucrit":critical,"E_loss_leptons":budget,"L_radio":radio,
                "CR_specs":np.concatenate((transport[:,6:7],electrons),axis=1).reshape(n*5,ne),
                "CR_specs_inj":np.concatenate((transport[:,4:6],escape),axis=1).reshape(n*4,ne)}
            for zone in (1,2):
                for j,process in enumerate(("SY","BS","IC","DI","IO")):
                    outputs[f"tau_loss_z{zone}_{process}"]=loss[:,(zone-1)*5+j]
            outputs["tau_loss_protons_PP"],outputs["tau_loss_protons_DI"]=loss[:,10],loss[:,11]
        for a in (kinetic,electron,photon,transport,electrons,emission,internal,*outputs.values()):
            a.flags.writeable=False
        return SolverResult(kinetic,electron,photon,
            {name:transport[:,i] for i,name in enumerate(TRANSPORT)},
            {name:electrons[:,i] for i,name in enumerate(ELECTRONS)},
            {name:emission[:,i] for i,name in enumerate(EMISSION)},
            {"solver_abi":5,"legacy_table_precision":legacy_table_precision,
             "catalogue":rows.tolist(),"catalogue_columns":["z","Mstar_Msun","Re_kpc","SFR_Msun_yr"],
             "cosmic_ray_bins":ne,"photon_bins":np_,"solver_cells":grid.cells,
             "table_grid":[preparation.grid.nx,preparation.grid.ny],"threads":threads,
             "preparation_fingerprint":preparation._preparation_hash,
             **provenance,
             "solver_python_sha256":provenance["python_source_sha256"],
             "solver_sha256":hashlib.sha256(Path(lib._name).read_bytes()).hexdigest(),
             "emission_units":"GeV^-1 s^-1 except dimensionless tau_free_free",
             "particle_units":"GeV^-1","frame":"source; disc free-free applied"},outputs,internal)
