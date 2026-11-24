# -*- coding: utf-8 -*-

# SPDX-License-Identifier: Apache-2.0
# Description: Catalogue of qubit geometries vs. charging energy (E_C), built from Palace capacitance
#              simulations of a small unit cell around a single qubit (instead of the full chip).

import copy
import hashlib
import itertools
import json
import os
import shutil
import subprocess
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import shapely
from qiskit_metal import designs
from qiskit_metal.toolbox_python.attr_dict import Dict

from SQDMetal.Comps.Qubits import TransmonTaperedInsets
from SQDMetal.Comps.Wires import WirePinStretch


def _merge(base, override):
    """Recursive dict merge: values in `override` replace those in `base` (copies, no mutation)."""
    out = Dict(copy.deepcopy(dict(base or {})))
    for k, v in (override or {}).items():
        out[k] = _merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else v
    return out


def _plain(obj):
    """Nested Dict/list -> plain JSON-serialisable structure (numbers kept, everything else str)."""
    if isinstance(obj, dict):
        return {str(k): _plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_plain(v) for v in obj]
    if isinstance(obj, (bool, int, float, str)) or obj is None:
        return obj
    if isinstance(obj, np.generic):
        return obj.item()
    return str(obj)


def _flatten(d, prefix=""):
    """{'a': {'b': 1}} -> {'a.b': 1}"""
    out = {}
    for k, v in (d or {}).items():
        path = f"{prefix}{k}"
        if isinstance(v, dict):
            out.update(_flatten(v, path + "."))
        else:
            out[path] = v
    return out


class QubitCatalogue:
    """Catalogue of qubit geometries and their charging energy E_C (plus capacitances), built from
    Palace capacitance simulations of a small unit cell around ONE qubit rather than the full chip.

    The unit cell is the qubit, a short open-ended CPW stub off its readout pin (the readout pad's
    capacitance is part of C_sigma; the rest of the resonator barely contributes) and a ground
    plane, on a small chip (cell_size x cell_size). Each geometry variant becomes one simulation
    folder; results are collected into catalogue.csv, which can then be inverted to find the
    geometry for a target E_C.

    Workflow:
        cat = QubitCatalogue("ec_catalogue", base_options=QUBIT_BASE_WITH_READOUT,
                             sim_options=Dict(palace_dir="/path/to/palace"))

        # 1. Mesh convergence on one variant (pick settings where C_sigma changes < ~0.5%)
        keys = [cat.prepare(Dict(pad_gap="100um"), mesh=Dict(min_size=m)) for m in (20e-6, 10e-6, 5e-6)]
        cat.run(keys)                    # needs Palace here; otherwise run the folders elsewhere
        cat.convergence(keys)

        # 2. Sweep (prepare + run; or prepare only, run the folders on the HPC, then collect)
        cat.sweep(Dict(pad_gap=["60um", "80um", "100um", "120um", "140um"]), run=True)

        # 3. Collect results and look up a geometry
        cat.collect(C_J_fF=2.0)
        cat.geometry_for(E_C_GHz=0.22, param="pad_gap")      # -> Dict(pad_gap="93.4um")

    Simulations made elsewhere can be added as entries with add_simulation(), e.g. single-qubit cells
    copied out of a full chip (SingleQubitCells.to_catalogue()). Those entries keep their own cell
    description (cell.builder = 'chip'), so with use_current=True geometry_for() does not mix them
    with this catalogue's own unit cells; their extra results (g, chi, kappa, ...) appear as result.*
    columns and their origin as source.* columns.

    Variants are keyed by a hash of (qubit options, cell, mesh, simulation settings), so repeated
    sweeps skip anything already prepared/simulated, and changing mesh settings makes new entries
    rather than overwriting old ones. Machine-specific settings (palace_dir, cpus, HPC JSON) are not
    part of the key.

    Inputs:
        * catalogue_dir - Folder holding sims/<cell_key>/ and catalogue.csv
        * base_options  - Qubit options shared by every variant (e.g. nqct04's QUBIT_BASE plus its
                          connection_pads); each variant's options are merged on top
        * qubit_class   - Qubit component class (default TransmonTaperedInsets)
        * readout_pin   - Qubit pin for the readout stub ('' -> no stub)
        * stub_length   - Length of the open-ended readout stub
        * cpw_width, cpw_gap - Stub CPW dimensions ('' -> taken from base_options
                          connection_pads.<readout_pin>.cpw_width/cpw_gap)
        * cell_size     - Side length of the square unit-cell chip
        * mesh          - fine_mesh_components settings for the qubit (metres): min_size,
                          max_size, taper_dist_min
        * sim_options   - Palace user_options (dielectric_material, solver_order, mesh_refinement,
                          palace_dir, ...)
        * mode          - Palace mode ('simPC', 'PC' or 'HPC')
        * cell_builder  - Optional function(design, qubit_options) -> (qubit_name, stub_name or None)
                          to build a custom unit cell (e.g. with charge/flux line stubs)
    """

    DEFAULT_MESH = Dict(min_size=10e-6, max_size=200e-6, taper_dist_min=10e-6)
    DEFAULT_SIM = Dict(dielectric_material="silicon", solver_order=2, mesh_refinement=0,
                       fillet_resolution=12, num_cpus=8)
    # Settings that do not change the physics (and so are not part of a variant's key)
    _MACHINE_KEYS = {"palace_dir", "HPC_Parameters_JSON", "num_cpus", "num_threads", "gmsh_verbosity",
                     "solns_to_save"}
    RESULT_COLUMNS = ["E_C_GHz", "C_sigma_fF", "C1_ground_fF", "C2_ground_fF", "C12_fF",
                      "C1_readout_fF", "C2_readout_fF", "C1_feed_fF", "C2_feed_fF"]
    # Files copied from an imported simulation folder (Palace field output is left behind: ~100s of MB)
    IMPORT_FILES = ("outputFiles/terminal-C.csv", "outputFiles/terminal-Cm.csv", "outputFiles/terminal-Cinv.csv",
                    "outputFiles/terminal-V.csv", "outputFiles/out.log", "outputFiles/mesh.png", "conductors.png")

    def __init__(self, catalogue_dir, base_options=None, qubit_class=TransmonTaperedInsets,
                 readout_pin="readout", stub_length="300um", cpw_width="", cpw_gap="",
                 cell_size="2mm", mesh=None, sim_options=None, mode="simPC", cell_builder=None):
        self.dir = Path(catalogue_dir).expanduser().resolve()
        self.sims_dir = self.dir / "sims"
        self.sims_dir.mkdir(parents=True, exist_ok=True)
        self.csv_path = self.dir / "catalogue.csv"
        self.base_options = Dict(base_options or {})
        self.qubit_class = qubit_class
        self.readout_pin = readout_pin
        self.stub_length = stub_length
        self.cpw_width = cpw_width
        self.cpw_gap = cpw_gap
        self.cell_size = cell_size
        self.mesh = _merge(self.DEFAULT_MESH, mesh)
        self.sim_options = _merge(self.DEFAULT_SIM, sim_options)
        self.mode = mode
        self.cell_builder = cell_builder
        self._sims = {}         # key -> prepared PALACE_Capacitance_Simulation (this session)

    # ------------------------------------------------------------------ #
    # Variants
    def variant_key(self, qubit_options, mesh=None, sim_options=None):
        """Short hash identifying a variant (geometry + cell + mesh + physics settings)."""
        spec = self._spec(qubit_options, mesh, sim_options)
        return hashlib.sha1(json.dumps(spec, sort_keys=True).encode()).hexdigest()[:12]

    def _spec(self, qubit_options, mesh=None, sim_options=None):
        sim = {k: v for k, v in _merge(self.sim_options, sim_options).items() if k not in self._MACHINE_KEYS}
        builder = getattr(self.cell_builder, "__name__", "default") if self.cell_builder else "default"
        return _plain(dict(qubit=_merge(self.base_options, qubit_options),
                           variant=Dict(qubit_options or {}),
                           cell=dict(qubit_class=self.qubit_class.__name__, readout_pin=self.readout_pin,
                                     stub_length=self.stub_length, cpw_width=self.cpw_width,
                                     cpw_gap=self.cpw_gap, cell_size=self.cell_size, builder=builder),
                           mesh=_merge(self.mesh, mesh), sim=sim, mode=self.mode))

    def make_unit_cell(self, qubit_options):
        """Builds the unit-cell design. Returns (design, qubit_name, stub_name or None)."""
        design = designs.DesignPlanar({}, overwrite_enabled=True)
        chip = design.chips.main.size
        chip.size_x = chip.size_y = self.cell_size
        chip.center_x = chip.center_y = "0mm"
        if self.cell_builder is not None:
            qname, stub = self.cell_builder(design, _merge(self.base_options, qubit_options))
            return design, qname, stub

        opts = _merge(self.base_options, qubit_options)
        opts.pos_x, opts.pos_y = "0mm", "0mm"
        self.qubit_class(design, "Q", options=opts)
        stub = None
        if self.readout_pin:
            pad = opts.get("connection_pads", Dict()).get(self.readout_pin, Dict())
            width = self.cpw_width or pad.get("cpw_width") or "10um"
            gap = self.cpw_gap or pad.get("cpw_gap") or "6um"
            WirePinStretch(design, "stub", options=Dict(
                pin_inputs=Dict(start_pin=Dict(component="Q", pin=self.readout_pin)),
                dist_extend=self.stub_length, trace_width=width, trace_gap=gap,
                end_gap="0um"))     # open end with a trace_gap cut
            stub = "stub"
        return design, "Q", stub

    # ------------------------------------------------------------------ #
    # Preparing and running simulations
    def _folder(self, key):
        return self.sims_dir / f"cell_{key}"

    def status(self, key):
        """'missing', 'prepared' or 'done' (Palace output present)."""
        folder = self._folder(key)
        if (folder / "outputFiles" / "terminal-C.csv").exists():
            return "done"
        if (folder / "variant.json").exists():
            return "prepared"
        return "missing"

    def prepare(self, qubit_options=None, mesh=None, sim_options=None, overwrite=False):
        """Builds, meshes and writes the Palace files for one variant (sims/cell_<key>/).
        Returns the variant key. Skips variants that already exist unless overwrite=True."""
        from SQDMetal.PALACE.Capacitance_Simulation import PALACE_Capacitance_Simulation

        qubit_options = Dict(qubit_options or {})
        key = self.variant_key(qubit_options, mesh, sim_options)
        if self.status(key) != "missing" and not overwrite:
            return key

        design, qname, stub = self.make_unit_cell(qubit_options)
        mesh_s = _merge(self.mesh, mesh)
        sim = PALACE_Capacitance_Simulation(
            name=f"cell_{key}", metal_design=design, sim_parent_directory=str(self.sims_dir) + "/",
            mode=self.mode, meshing="GMSH", user_options=dict(_merge(self.sim_options, sim_options)),
            create_files=True)
        sim.add_metallic(1)
        sim.add_ground_plane()
        sim.fine_mesh_components([qname], min_size=float(mesh_s.min_size), max_size=float(mesh_s.max_size),
                                 taper_dist_min=float(mesh_s.taper_dist_min), metals_only=False)
        sim.prepare_simulation()

        conductors = self._identify_conductors(design, qname, stub, sim._cur_cap_terminals)
        record = dict(key=key, created=datetime.now().isoformat(timespec="seconds"),
                      sqdmetal_git=self._git_hash(), conductor_indices=conductors,
                      spec=self._spec(qubit_options, mesh, sim_options))
        with open(self._folder(key) / "variant.json", "w") as f:
            json.dump(record, f, indent=2)
        self._sims[key] = sim
        return key

    def prepare_sweep(self, grid, fixed=None, **kwargs):
        """Prepares every combination of the grid, e.g. grid=Dict(pad_gap=["80um", "100um"]).
        Nested options use dotted names: {"connection_pads.readout.pad_gap": [...]}. `fixed`
        options are applied to every variant. Returns the list of keys of the variants prepared;
        variants that fail (e.g. pads that touch) are skipped and listed in self.failed."""
        names = list(grid.keys())
        keys = []
        self.failed = []        # (options, error) for variants that could not be prepared
        for values in itertools.product(*[grid[n] for n in names]):
            opts = Dict(copy.deepcopy(dict(fixed or {})))
            for name, value in zip(names, values):
                self._set_dotted(opts, name, value)
            try:
                keys.append(self.prepare(opts, **kwargs))
            except Exception as e:
                self.failed.append((dict(zip(names, values)), str(e)))
                print(f"QubitCatalogue: skipped {dict(zip(names, values))}: {e}")
        return keys

    def run(self, keys=None):
        """Runs Palace locally (sim_options.palace_dir) for prepared variants that have no results
        yet (all of them if keys is None). Variants prepared in an earlier session are prepared
        again first (a few seconds each)."""
        keys = keys if keys is not None else [k for k in self.keys() if self.status(k) == "prepared"]
        palace = os.path.expanduser(str(self.sim_options.get("palace_dir") or "palace"))
        if not (os.path.isfile(palace) or shutil.which(palace)):
            raise RuntimeError(f"Palace was not found ('{palace}'); set sim_options.palace_dir, or run the "
                               f"prepared folders elsewhere (see write_run_script()) and then collect().")
        for key in keys:
            if self.status(key) == "done":
                continue
            sim = self._sims.get(key)
            if sim is None:
                if self._record(key).get("source"):
                    raise ValueError(f"Variant {key} was imported with add_simulation() and has no results; "
                                     f"re-run it where it came from and add it again.")
                spec = self._record(key)["spec"]
                new_key = self.prepare(Dict(spec["variant"]), mesh=Dict(spec["mesh"]),
                                       sim_options=Dict(spec["sim"]), overwrite=True)
                if new_key != key:
                    raise ValueError(f"Variant {key} was prepared with different catalogue settings "
                                     f"(qubit class, base options, cell or simulation settings); "
                                     f"re-create the catalogue with those settings to run it.")
                sim = self._sims[key]
            sim.run()
        return keys

    def write_run_script(self, path=None, palace_cmd="palace", np_procs=4, keys=None):
        """Writes a bash script that runs Palace on every prepared variant without results (or on
        `keys`), e.g. to copy the catalogue folder to a simulation PC / HPC node and run it there.
        Each variant runs in its own folder (its config uses relative paths). Afterwards, copy the
        folder back (or point the catalogue at it) and call collect(). Returns the script path."""
        keys = keys if keys is not None else [k for k in self.keys() if self.status(k) == "prepared"]
        path = Path(path) if path else self.dir / "run_pending.sh"
        lines = ["#!/bin/bash", "# Runs Palace on pending QubitCatalogue variants", "set -e",
                 'cd "$(dirname "$0")/sims"']
        for key in keys:
            lines.append(f'(cd cell_{key} && {palace_cmd} -np {np_procs} cell_{key}.json > palace.log 2>&1) '
                         f'&& echo "done {key}"')
        path.write_text("\n".join(lines) + "\n")
        path.chmod(0o755)
        return path

    def sweep(self, grid, fixed=None, run=False, **kwargs):
        """prepare_sweep(), then (if run=True) run() and collect(). Returns the keys."""
        keys = self.prepare_sweep(grid, fixed=fixed, **kwargs)
        if run:
            self.run(keys)
            self.collect()
        return keys

    def add_simulation(self, sim_folder, qubit_options, conductor_indices, cell, mesh=None, sim_options=None,
                       qubit_class=None, source=None, results=None, components=None, mode="simPC", overwrite=False):
        """Adds a capacitance simulation that was set up and run elsewhere (e.g. a single-qubit cell
        copied out of a full chip by SingleQubitCells) as a catalogue entry.

        The small result files (terminal-*.csv, out.log, mesh.png, conductors.png) are copied from
        `sim_folder` into sims/cell_<key>/; Palace's field output is not.

        Inputs:
            * sim_folder        - Simulation folder holding outputFiles/terminal-C.csv
            * qubit_options     - The simulated qubit's full options (placement options are dropped)
            * conductor_indices - Capacitance-matrix index of ground / pad1 / pad2 (/ res / feed)
            * cell              - Description of the simulated cell; stored as cell.* columns. Include
                                  builder (e.g. 'chip') so these entries are not matched with the
                                  catalogue's own unit cells.
            * mesh, sim_options - Mesh and Palace settings used (machine-specific ones are dropped)
            * qubit_class       - Qubit class (default: the catalogue's); used for the swept/variant
                                  column, which lists options differing from the class defaults
                                  merged with the catalogue's base_options
            * source            - Where the entry came from (design, qubit, ...); stored as source.*
                                  columns. source['cell_key'], if given, distinguishes otherwise
                                  identical qubits in different cells (e.g. other resonators).
            * results           - Extra results (g, chi, kappa, ...); stored as result.* columns
            * components        - Options of the other components in the cell (kept in variant.json)
        Returns the entry's key. Existing entries are left alone unless overwrite=True."""
        src_folder = Path(sim_folder).expanduser().resolve()
        csv = src_folder / "outputFiles" / "terminal-C.csv"
        if not csv.exists():
            raise FileNotFoundError(f"No capacitance results in {sim_folder} ({csv} is missing).")
        for role in ("ground", "pad1", "pad2"):
            if role not in conductor_indices:
                raise ValueError(f"conductor_indices needs '{role}'; got {conductor_indices}.")
        qubit_class = qubit_class or self.qubit_class
        options = Dict(copy.deepcopy(dict(qubit_options)))
        for k in ("pos_x", "pos_y", "orientation"):
            options.pop(k, None)
        sim = {k: v for k, v in dict(sim_options or {}).items() if k not in self._MACHINE_KEYS}
        cell = dict(cell, qubit_class=qubit_class.__name__)
        spec = _plain(dict(qubit=options, variant=self._diff_from_defaults(options, qubit_class), cell=cell,
                           mesh=dict(mesh or {}), sim=sim, mode=mode))
        source = _plain(dict(source or {}))
        key = hashlib.sha1(json.dumps(dict(spec=spec, cell_key=source.get("cell_key", "")),
                                      sort_keys=True).encode()).hexdigest()[:12]
        folder = self._folder(key)
        if self.status(key) == "done" and not overwrite:
            return key
        for rel in self.IMPORT_FILES:
            if (src_folder / rel).exists():
                (folder / rel).parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src_folder / rel, folder / rel)
        record = dict(key=key, created=datetime.now().isoformat(timespec="seconds"),
                      sqdmetal_git=self._git_hash(), conductor_indices=dict(conductor_indices), spec=spec,
                      source=dict(source, sim_folder=str(src_folder)), results=_plain(dict(results or {})),
                      components=_plain(dict(components or {})))
        with open(folder / "variant.json", "w") as f:
            json.dump(record, f, indent=2)
        return key

    def _diff_from_defaults(self, options, qubit_class):
        """Options that differ from qubit_class's defaults merged with the catalogue's base_options."""
        try:
            defaults = qubit_class.get_template_options(designs.DesignPlanar({}, overwrite_enabled=True))
        except Exception:
            defaults = Dict()
        ref = _flatten(_plain(_merge(defaults, self.base_options)))
        out = Dict()
        for name, value in _flatten(_plain(options)).items():
            if name not in ref or str(ref[name]) != str(value):
                self._set_dotted(out, name, value)
        return out

    # ------------------------------------------------------------------ #
    # Results
    def keys(self):
        return sorted(p.name[len("cell_"):] for p in self.sims_dir.glob("cell_*") if (p / "variant.json").exists())

    def collect(self, C_J_fF=0.0):
        """Reads every variant's Palace output into catalogue.csv (pending variants are listed with
        status 'prepared'). C_J_fF is the junction capacitance added to C_sigma. Returns the table."""
        rows = []
        for key in self.keys():
            record = self._record(key)
            spec = record["spec"]
            row = dict(key=key, status=self.status(key),
                       swept=json.dumps(_flatten(spec["variant"]), sort_keys=True))
            # Every qubit option (base options merged with the variant's), so variants built with
            # different base options can be told apart
            for name, value in _flatten(spec["qubit"]).items():
                row[name] = value
                um = self._to_um(value)
                if um is not None:
                    row[f"{name}_um"] = um
            if row["status"] == "done":
                try:
                    row.update(self.ec_from_matrix(self._folder(key) / "outputFiles" / "terminal-C.csv",
                                                   record["conductor_indices"], C_J_fF))
                except Exception as e:      # e.g. an unreadable or unexpected matrix
                    row["status"] = f"error: {e}"
            row.update({f"cell.{k}": v for k, v in spec["cell"].items()})
            row.update({f"mesh.{k}": v for k, v in spec["mesh"].items()})
            row.update({f"sim.{k}": v for k, v in spec["sim"].items()})
            # Entries imported with add_simulation(): where they came from and their extra results
            row.update({f"source.{k}": v for k, v in _flatten(record.get("source", {})).items()})
            row.update({f"result.{k}": v for k, v in _flatten(record.get("results", {})).items()})
            row.update(C_J_fF=C_J_fF, sqdmetal_git=record.get("sqdmetal_git", ""), created=record.get("created", ""))
            rows.append(row)
        table = pd.DataFrame(rows)
        table.to_csv(self.csv_path, index=False)
        return table

    @staticmethod
    def ec_from_matrix(csv_path, conductor_indices, C_J_fF=0.0):
        """E_C and capacitances of a floating transmon from Palace's terminal-C.csv (Maxwell
        capacitance matrix in F), using the same expression as
        PALACE_Capacitance_Simulation.calc_params_floating_Transmon_from_files: the readout
        conductor, if present, is treated as AC ground for each pad,
            C_sigma = (C1g + C1r)(C2g + C2r) / (C1g + C1r + C2g + C2r) + C12 + C_J,
            E_C = e^2 / (2 C_sigma)."""
        e, h = 1.602176634e-19, 6.62607015e-34
        mat = pd.read_csv(csv_path).to_numpy()
        if np.allclose(mat[:, 0], np.arange(1, mat.shape[0] + 1)):
            mat = mat[:, 1:]        # first column holds the row indices
        idx = conductor_indices
        c = lambda a, b: abs(mat[idx[a], idx[b]]) if (a in idx and b in idx) else 0.0
        C1g, C2g, C12 = c("pad1", "ground"), c("pad2", "ground"), c("pad1", "pad2")
        C1r, C2r = c("pad1", "res"), c("pad2", "res")
        C_sigma = (C1g + C1r) * (C2g + C2r) / (C1g + C1r + C2g + C2r) + C12 + C_J_fF * 1e-15
        out = dict(E_C_GHz=e ** 2 / (2 * C_sigma) / h / 1e9, C_sigma_fF=C_sigma * 1e15,
                   C1_ground_fF=C1g * 1e15, C2_ground_fF=C2g * 1e15, C12_fF=C12 * 1e15,
                   C1_readout_fF=C1r * 1e15, C2_readout_fF=C2r * 1e15)
        if "feed" in idx:       # reported only (not part of C_sigma, as in calc_params_floating_Transmon)
            out.update(C1_feed_fF=c("pad1", "feed") * 1e15, C2_feed_fF=c("pad2", "feed") * 1e15)
        return out

    @property
    def table(self):
        """The last collected catalogue (catalogue.csv)."""
        return pd.read_csv(self.csv_path) if self.csv_path.exists() else pd.DataFrame()

    def convergence(self, keys):
        """Compares C_sigma and E_C across variants of the same geometry with different mesh/solver
        settings (in the order given), relative to the last (finest) one."""
        table = self.collect() if not self.csv_path.exists() else self.table
        sel = table.set_index("key").loc[keys]
        out = sel[[c for c in sel.columns if c.startswith("mesh.") or c.startswith("sim.")] + ["C_sigma_fF", "E_C_GHz"]].copy()
        ref = out["C_sigma_fF"].iloc[-1]
        out["C_sigma_change_%"] = 100 * (out["C_sigma_fF"] - ref) / ref
        return out

    _META_COLUMNS = {"key", "status", "swept", "C_J_fF", "sqdmetal_git", "created"}

    def _current_settings(self):
        """The catalogue's current base qubit options, cell, mesh and solver settings, as
        {column name: value} in the catalogue table's naming."""
        spec = self._spec(Dict())
        cur = dict(_flatten(spec["qubit"]))
        cur.update({f"cell.{k}": v for k, v in spec["cell"].items()})
        cur.update({f"mesh.{k}": v for k, v in spec["mesh"].items()})
        cur.update({f"sim.{k}": v for k, v in spec["sim"].items()})
        return cur

    def _matches(self, column, value):
        """Boolean mask of rows whose `column` equals `value` (lengths compared in um, numbers
        approximately, anything else as text; a missing column/value only matches missing)."""
        table = self._sel
        if column not in table:
            return pd.Series(value is None, index=table.index)
        if value is None or (isinstance(value, str) and value == ""):
            # empty settings are read back from the CSV as NaN
            return table[column].isna() | (table[column].astype(str) == "")
        um = self._to_um(value)
        if um is not None and f"{column}_um" in table:
            return pd.Series(np.isclose(table[f"{column}_um"].astype(float), um), index=table.index)
        if isinstance(value, (int, float)) and not isinstance(value, bool) \
                and pd.api.types.is_numeric_dtype(table[column]):
            return pd.Series(np.isclose(table[column].astype(float), value), index=table.index)
        return table[column].astype(str) == str(value)

    def geometry_for(self, E_C_GHz, param, fixed=None, use_current=True, decimals=1):
        """Geometry giving the target E_C, by interpolating catalogue results along one swept
        option (`param`, e.g. 'pad_gap').

        Only rows that share everything else are used:
            * use_current=True (default): rows built with this catalogue object's current base
              options, cell, mesh and solver settings (apart from `param`), so results from
              earlier settings in the same folder are left out;
            * fixed: explicit values that override / add to those, e.g.
              fixed={"connection_pads.readout.pad_gap": "15um", "mesh.min_size": 1e-5}.
        Anything that still varies among the selected rows raises an error naming it.

        Returns a Dict of qubit options for qubit_overrides.<Qn>.qubit: `param` plus any geometry
        given in `fixed` (the catalogue's base options are assumed to match the chip's own)."""
        table = self.table
        if table.empty:
            raise ValueError("The catalogue is empty: run collect() first.")
        table = table[table["status"] == "done"]
        col = f"{param}_um"
        if col not in table:
            raise ValueError(f"'{param}' is not a length option in the catalogue.")
        fixed = _flatten(Dict(fixed or {}))

        criteria = dict(self._current_settings()) if use_current else {}
        criteria.pop(param, None)
        criteria.update(fixed)
        self._sel = table
        excluded = {}
        for name, value in criteria.items():
            if name not in table.columns and name not in fixed and use_current:
                continue        # e.g. an option added since the catalogue was collected
            if name not in table.columns:
                raise ValueError(f"'{name}' is not a column of the catalogue.")
            mask = self._matches(name, value)
            if (~mask).any():
                excluded[name] = (int((~mask).sum()), value)
            self._sel = self._sel[mask]
        table = self._sel

        # Anything else that still varies must be fixed
        others = [c for c in table.columns
                  if c not in self._META_COLUMNS and c not in self.RESULT_COLUMNS
                  and not c.startswith(("source.", "result."))
                  and not c.endswith("_um") and c != param and table[c].nunique(dropna=True) > 1]
        if others:
            raise ValueError(f"These vary across the selected rows: {others}. Give their values in "
                             f"`fixed`.")
        table = table.dropna(subset=["E_C_GHz"]).sort_values(col)
        if len(table) < 2:
            why = "; ".join(f"{n} != {v!r} ({k} rows)" for n, (k, v) in excluded.items())
            raise ValueError(f"Need at least two simulated points along {param} with the selected "
                             f"settings, found {len(table)}." + (f" Excluded: {why}." if why else ""))
        x, ec = table[col].to_numpy(), table["E_C_GHz"].to_numpy()
        if not (np.all(np.diff(ec) > 0) or np.all(np.diff(ec) < 0)):
            pts = ", ".join(f"{xi:g}um: {ei:.4f}" for xi, ei in zip(x, ec))
            raise ValueError(f"E_C is not monotonic in {param} ({pts} GHz); narrow the range.")
        if not min(ec) <= E_C_GHz <= max(ec):
            raise ValueError(f"E_C={E_C_GHz} GHz is outside the catalogue range "
                             f"{min(ec):.4f}-{max(ec):.4f} GHz for {param} "
                             f"({x.min():g}-{x.max():g} um).")
        order = np.argsort(ec)
        value = float(np.interp(E_C_GHz, ec[order], x[order]))
        # Qubit options only (cell/mesh/solver settings are selection criteria, not options)
        out = Dict()
        for name, v in fixed.items():
            if not name.split(".")[0] in ("mesh", "sim", "cell"):
                self._set_dotted(out, name, v)
        self._set_dotted(out, param, f"{round(value, decimals)}um")
        return out

    # ------------------------------------------------------------------ #
    # Helpers
    def _record(self, key):
        with open(self._folder(key) / "variant.json") as f:
            return json.load(f)

    @staticmethod
    def _set_dotted(d, name, value):
        parts = name.split(".")
        for p in parts[:-1]:
            d = d.setdefault(p, Dict())
        d[parts[-1]] = value

    @staticmethod
    def _to_um(value):
        """'100um' / '0.1mm' / '100nm' -> 100.0 (um); anything else -> None."""
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            return None
        s = str(value).strip().lower().replace(" ", "")
        for unit, scale in (("nm", 1e-3), ("um", 1.0), ("mm", 1e3)):
            if s.endswith(unit):
                try:
                    return float(s[:-len(unit)]) * scale
                except ValueError:
                    return None
        return None

    @staticmethod
    def _identify_conductors(design, qname, stub, terminals):
        """Matches the capacitance-matrix conductors to ground / pad1 / pad2 / res by geometry:
        ground is the largest conductor, res the one containing the stub, and the pads the two
        others touching the qubit's metal (pad1 = the higher / further right one)."""
        terminals = [shapely.make_valid(t) for t in terminals]
        ground = int(np.argmax([t.area for t in terminals]))
        res = None
        if stub is not None:
            df = design.qgeometry.get_component(stub, "path")
            trace = [row["geometry"] for _, row in df.iterrows() if not bool(row["subtract"])][0]
            tip = shapely.Point(trace.coords[-1])
            res = int(np.argmin([t.distance(tip) for t in terminals]))
        q_metal = shapely.unary_union([row["geometry"] for _, row in
                                       design.qgeometry.get_component(qname, "poly").iterrows()
                                       if not bool(row["subtract"])])
        pads = [i for i, t in enumerate(terminals) if i not in (ground, res) and t.intersects(q_metal)]
        if len(pads) == 1:
            raise ValueError("The qubit's two pads touch or overlap (they form one conductor, which "
                             "would short the junction): increase pad_gap or reduce the taper sizes.")
        if len(pads) != 2:
            raise ValueError(f"Expected 2 qubit pads among the conductors, found {len(pads)}.")
        pads.sort(key=lambda i: tuple(-c for c in terminals[i].centroid.coords[0][::-1]))
        out = dict(ground=ground, pad1=pads[0], pad2=pads[1])
        if res is not None:
            out["res"] = res
        return out

    @staticmethod
    def _git_hash():
        try:
            folder = os.path.dirname(os.path.abspath(__file__))
            return subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=folder, capture_output=True,
                                  text=True, timeout=5).stdout.strip()
        except Exception:
            return ""
