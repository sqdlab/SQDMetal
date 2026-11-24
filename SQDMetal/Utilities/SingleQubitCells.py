# -*- coding: utf-8 -*-

# SPDX-License-Identifier: Apache-2.0
# Description: Copies single-qubit cells (qubit + readout resonator + coupler) out of a full chip design, with
#              the chip's own component options, and runs Palace capacitance simulations on them.

import contextlib
import copy
import hashlib
import io
import json
import os
import re
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import shapely
import qiskit_metal as metal
from qiskit_metal import Dict
from qiskit_metal.qlibrary.core import BaseQubit
from qiskit_metal.qlibrary.couplers.coupled_line_tee import CoupledLineTee

from SQDMetal.Comps.Wires import WirePinStretch
from SQDMetal.Utilities.QubitDesigner import ResonatorHalfWave


class SingleQubitCells:
    """Single-qubit capacitance simulations built from a full chip design.

    For each qubit, every component linked to it through ``pin_inputs`` (e.g. readout meander and coupling tee)
    is copied into a small design with the same class and a deep copy of the chip's options, so the cell
    cannot drift from the chip. Open-ended feedline stubs are added to coupling tees
    (``CoupledLineTee``), and the copied geometry is checked against the chip.

    Capacitance-matrix conductors are identified by geometry (ground / pad1 / pad2 / res / feed), because
    the order Palace receives depends on where the shapes sit rather than on the order components were
    created. Each simulation saves ``conductors.png`` showing the mapping.

    Typical use::

        cells = SingleQubitCells(design, "sims", user_options=dict(palace_dir="...", num_cpus=8))
        rec = cells.simulate("Q1")                       # meshes, runs Palace (cached), saves conductors.png
        p = cells.params("Q1", f_q_GHz=5.0, f_r_GHz=7.0) # calc_params_floating_Transmon_from_files output

    Parameters
    ----------
    design : QDesign
        The full chip design (already built).
    sim_dir : str or Path
        Folder holding one sub-folder per simulated cell.
    user_options : dict
        Palace ``user_options`` (palace_dir, num_cpus, solver_order, ...).
    mesh : dict, optional
        Keyword arguments for ``fine_mesh_components`` (min_size, max_size, ...), in metres.
    exclude : str
        Regex of component names never pulled into a cell (e.g. a feedline shared by all qubits).
    feedline_stub : str
        Length of the open-ended feedline kept either side of each coupling tee.
    chip_margin_mm : float
        Substrate margin around the cell.
    name : str
        Prefix for simulation folder names (``<name>_<label>_cap``).
    verbose : bool
        If False, build/meshing printouts are suppressed (Palace's own output is not).
    """

    DEFAULT_MESH = dict(mesh_sampling=160, min_size=10e-6, max_size=250e-6, taper_dist_min=12e-6,
                        taper_dist_max=250e-6)
    ROLE_LABELS = dict(ground="ground", pad1="pad 1 (readout side)", pad2="pad 2", res="readout resonator",
                       feed="feedline")
    ROLE_COLOURS = dict(ground="#c9d3e0", pad1="#d1495b", pad2="#edae49", res="#00798c", feed="#30638e")
    _MACHINE_KEYS = ("palace_dir", "num_cpus", "num_threads", "HPC_Parameters_JSON", "gmsh_verbosity")

    def __init__(self, design, sim_dir, user_options, mesh=None, exclude=r"feedline|LP_", feedline_stub="500um",
                 chip_margin_mm=0.6, name="cell", verbose=False):
        self.design = design
        self.sim_dir = Path(sim_dir).expanduser().resolve()
        self.user_options = dict(user_options)
        self.mesh = dict(self.DEFAULT_MESH, **(mesh or {}))
        self.exclude = exclude
        self.feedline_stub = feedline_stub
        self.chip_margin_mm = chip_margin_mm
        self.name = name
        self.verbose = verbose
        self._cells = {}
        self.records = {}

    # ------------------------------------------------------------------ #
    # Building cells
    def qubit_names(self):
        """Names of every qubit component (``BaseQubit`` subclasses) in the design."""
        return [n for n in self.design.components.keys() if isinstance(self.design.components[n], BaseQubit)]

    @staticmethod
    def _pin_refs(comp):
        pins = comp.options.get("pin_inputs", None) or {}
        return [v["component"] for v in pins.values() if isinstance(v, dict) and v.get("component")]

    def components(self, qubit, layer="1"):
        """Components linked to ``qubit`` through ``pin_inputs`` (either direction) on ``layer``, in the
        design's creation order. Raises if the cell reaches another qubit."""
        links = defaultdict(set)
        for n in self.design.components.keys():
            for ref in self._pin_refs(self.design.components[n]):
                links[n].add(ref)
                links[ref].add(n)
        seen, todo = {qubit}, [qubit]
        while todo:
            for nxt in links[todo.pop()]:
                if nxt not in seen and not re.search(self.exclude, nxt):
                    seen.add(nxt)
                    todo.append(nxt)
        others = (seen & set(self.qubit_names())) - {qubit}
        if others:
            raise ValueError(f"{qubit}'s cell reaches other qubits {others}: add the shared parts to `exclude`.")
        return [n for n in self.design.components.keys()
                if n in seen and str(self.design.components[n].options.get("layer", "1")) == layer]

    def cell(self, qubit):
        """The cell for ``qubit``: Dict(design, components, stubs, geometry_difference_mm2). Built once and
        checked against the full chip (raises if any copied component's geometry differs)."""
        if qubit not in self._cells:
            names = self.components(qubit)
            sub = metal.designs.DesignPlanar(metadata=Dict(design_name=f"{self.name}_{qubit}"))
            sub.overwrite_enabled = True
            for n in names:
                comp = self.design.components[n]
                type(comp)(sub, n, options=copy.deepcopy(dict(comp.options)))
            stubs = []
            for n in names:
                if isinstance(sub.components[n], CoupledLineTee):
                    p = sub.components[n].options
                    for pin in ("prime_start", "prime_end"):
                        WirePinStretch(sub, f"{n}_{pin}_stub", options=Dict(
                            pin_inputs=Dict(start_pin=Dict(component=n, pin=pin)), dist_extend=self.feedline_stub,
                            trace_width=p.prime_width, trace_gap=p.prime_gap, end_gap="0um"))
                        stubs.append(f"{n}_{pin}_stub")
            m = self.chip_margin_mm
            b = np.array([sub.components[n].qgeometry_bounds() for n in sub.components.keys()])
            x0, y0, x1, y1 = b[:, 0].min() - m, b[:, 1].min() - m, b[:, 2].max() + m, b[:, 3].max() + m
            size = sub.chips.main.size
            size.center_x, size.center_y = f"{(x0 + x1) / 2:.6f}mm", f"{(y0 + y1) / 2:.6f}mm"
            size.size_x, size.size_y = f"{x1 - x0:.6f}mm", f"{y1 - y0:.6f}mm"
            size.size_z = self.design.chips.main.size.size_z
            diff = self.check_geometry(self.design, sub, names)
            self._cells[qubit] = Dict(design=sub, components=names, stubs=stubs, geometry_difference_mm2=diff)
        return self._cells[qubit]

    # ------------------------------------------------------------------ #
    # Simulation
    def _quiet(self):
        return contextlib.nullcontext() if self.verbose else contextlib.redirect_stdout(io.StringIO())

    def _key(self, sub):
        spec = dict(components={n: dict(cls=type(sub.components[n]).__name__, options=sub.components[n].options)
                                for n in sub.components.keys()},
                    chip=dict(sub.chips.main.size), mesh=self.mesh,
                    sim={k: v for k, v in self.user_options.items() if k not in self._MACHINE_KEYS})
        return hashlib.sha1(json.dumps(spec, sort_keys=True, default=str).encode()).hexdigest()[:12]

    def _prepare(self, sim_name, sub):
        from SQDMetal.PALACE.Capacitance_Simulation import PALACE_Capacitance_Simulation
        with self._quiet():
            sim = PALACE_Capacitance_Simulation(name=sim_name, metal_design=sub,
                                                sim_parent_directory=str(self.sim_dir) + "/", mode="simPC",
                                                meshing="GMSH", user_options=dict(self.user_options),
                                                create_files=True)
            sim.add_metallic(1)
            sim.add_ground_plane()
            sim.fine_mesh_components(list(sub.components.keys()), metals_only=True, **self.mesh)
            sim.prepare_simulation()
        return sim

    def _conductors(self, qubit, sim):
        c = self.cell(qubit)
        sub, names = c.design, c.components
        resonators = [n for n in names if qubit in self._pin_refs(sub.components[n])]
        try:
            return self.identify_conductors(sub, qubit, resonators[0] if resonators else None, c.stubs,
                                            sim._cur_cap_terminals)
        except ValueError:
            for n in resonators + c.stubs:
                contact = self.ground_contact(sub, n, self.user_options.get("fillet_resolution", 12))
                if contact is not None:
                    raise ValueError(
                        f"{qubit}: {n} touches the ground plane at ({contact[0]:.4f}, {contact[1]:.4f}) mm once "
                        f"filleted, so it merges into the ground conductor. This is in the chip geometry (it would "
                        f"also short it in the GDS), e.g. a route fillet as large as the straight lead next to "
                        f"it, which leaves a sliver of ground at the pin.") from None
            raise

    @staticmethod
    def ground_contact(design, component, fillet_resolution=12):
        """Point (x, y) where ``component``'s rendered (filleted) metal touches the ground plane, or None.
        The ground plane is the chip minus every rendered cut-out, as built for simulation."""
        from SQDMetal.Utilities.QiskitShapelyRenderer import QiskitShapelyRenderer
        from SQDMetal.Utilities.ShapelyEx import ShapelyEx
        from shapely.ops import nearest_points
        g = QiskitShapelyRenderer(design=design, canvas=None, logger=None).get_net_coordinates(fillet_resolution)
        gaps = ShapelyEx.shapely_to_list(ShapelyEx.fuse_polygons_threshold(g.loc[g["subtract"]]["geometry"].buffer(0), 1e-12))
        s = design.chips.main.size
        cx, cy, sx, sy = (design.parse_value(s[k]) for k in ("center_x", "center_y", "size_x", "size_y"))
        ground = shapely.box(cx - sx / 2, cy - sy / 2, cx + sx / 2, cy + sy / 2).difference(shapely.MultiPolygon(gaps))
        rows = g.loc[(~g["subtract"]) & (g["component"] == design.components[component].id)]
        if len(rows) == 0:
            return None
        metal = shapely.unary_union(rows["geometry"].buffer(0))
        if metal.distance(ground) > 1e-9:
            return None
        p = nearest_points(metal, ground)[0]
        return (p.x, p.y)

    def simulate(self, qubit, label=None, run=True, force_rerun=False):
        """Meshes and (if ``run``) simulates the cell of ``qubit`` with Palace. Results are reused when the
        cell geometry and physics settings are unchanged (unless ``force_rerun``). Saves conductors.png
        (matrix index -> conductor) in the simulation folder.

        Returns a dict: key, qubit, components, conductor_indices, output_dir, plot, created and status
        ('simulated', 'cached' or 'prepared')."""
        label = label or qubit
        sub = self.cell(qubit).design
        sim_name = f"{self.name}_{label}_cap"
        folder = self.sim_dir / sim_name
        record_path, plot_path = folder / "cell.json", folder / "conductors.png"
        csv = folder / "outputFiles" / "terminal-C.csv"
        key = self._key(sub)
        self.sim_dir.mkdir(parents=True, exist_ok=True)

        def save_plot(terminals, idx):
            fig, ax = plt.subplots(figsize=(5, 6))
            self.plot_conductors(terminals, idx, ax, f"{label}: capacitance-matrix conductors")
            fig.savefig(plot_path, dpi=150, bbox_inches="tight")
            plt.close(fig)

        if record_path.exists() and not force_rerun:
            rec = json.loads(record_path.read_text())
            if rec["key"] == key and csv.exists():
                if not plot_path.exists():      # mesh again (no Palace run) to draw the plot
                    sim = self._prepare(sim_name, sub)
                    idx = self._conductors(qubit, sim)
                    if idx != rec["conductor_indices"]:
                        raise RuntimeError(f"{label}: conductor order changed since the cached run; "
                                           f"use force_rerun=True.")
                    save_plot(sim._cur_cap_terminals, idx)
                rec = dict(rec, label=label, plot=str(plot_path), status="cached")
                self.records[qubit] = rec
                return rec

        sim = self._prepare(sim_name, sub)
        idx = self._conductors(qubit, sim)
        save_plot(sim._cur_cap_terminals, idx)
        rec = dict(key=key, qubit=qubit, label=label, components=self.cell(qubit).components, conductor_indices=idx,
                   output_dir=str(folder / "outputFiles"), created=datetime.now().isoformat(timespec="seconds"))
        if csv.exists():
            csv.unlink()        # results of an older geometry must not be read back
        record_path.write_text(json.dumps(rec, indent=2))
        if run:
            sim.run()
            if not csv.exists():
                raise RuntimeError(f"Palace finished without {csv}; see {folder / 'outputFiles' / 'out.log'}")
        rec = dict(rec, plot=str(plot_path), status="simulated" if run else "prepared")
        self.records[qubit] = rec
        return rec

    def eigenmode(self, qubit, L_J_nH, kinetic_inductance_pH_sq=0.0, label=None, start_freq_GHz=1.0, n_modes=3,
                  C_J_fF=0.0, f_r_GHz=None, epr_interfaces=True, mesh=None, user_options=None, run=True,
                  force_rerun=False):
        """Eigenmode simulation of the cell of ``qubit`` with the junction as a lumped inductor ``L_J_nH``,
        50 Ohm ports at the ends of the feedline stubs (so the resonator and qubit are loaded as on the
        chip) and, if ``kinetic_inductance_pH_sq`` > 0, a kinetic sheet inductance (pH/square) on all metal.
        Results are reused when the cell, junction, kinetic inductance and solver settings are unchanged.

        Inputs:
            * start_freq_GHz, n_modes - the solver finds ``n_modes`` modes above ``start_freq_GHz``
            * f_r_GHz        - expected (bare) resonator frequency, used only to pick the resonator mode
            * epr_interfaces - also compute surface participations (Si-vacuum, Si-Al, Al-vacuum)
            * mesh, user_options - override this object's mesh / Palace settings for this run

        Returns a dict: modes (DataFrame: f_GHz, Q, kappa_MHz, EPR_JJ per mode), qubit_mode / resonator_mode
        (1-based), summary (dressed f_q and f_r, anharmonicity, cross-Kerr chi_qr, kappa_r, T1 limits),
        epr (the full output of calculate_hamiltonian_parameters_EPR), output_dir and status."""
        import pandas as pd
        from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation
        from SQDMetal.Utilities.Materials import MaterialInterface
        label = label or qubit
        c = self.cell(qubit)
        sub = c.design
        sim_name = f"{self.name}_{label}_eig"
        folder = self.sim_dir / sim_name
        record_path = folder / "eigenmode.json"
        eig_csv = folder / "outputFiles" / "eig.csv"
        mesh = dict(self.mesh, **(mesh or {}))
        opts = dict(self.user_options, **(user_options or {}))
        opts.update(starting_freq=start_freq_GHz * 1e9, number_of_freqs=n_modes, solns_to_save=n_modes)
        settings = dict(cell=self._key(sub), L_J_nH=L_J_nH, C_J_fF=C_J_fF, kinetic_inductance_pH_sq=kinetic_inductance_pH_sq,
                        start_freq_GHz=start_freq_GHz, n_modes=n_modes, epr_interfaces=epr_interfaces, mesh=mesh,
                        sim={k: v for k, v in opts.items() if k not in self._MACHINE_KEYS})
        key = hashlib.sha1(json.dumps(settings, sort_keys=True, default=str).encode()).hexdigest()[:12]
        self.sim_dir.mkdir(parents=True, exist_ok=True)

        cached = False
        if record_path.exists() and self._eig_rows(eig_csv) >= n_modes and not force_rerun:
            cached = json.loads(record_path.read_text()).get("key") == key
        if not cached:
            with self._quiet():
                sim = PALACE_Eigenmode_Simulation(name=sim_name, metal_design=sub,
                                                  sim_parent_directory=str(self.sim_dir) + "/", mode="simPC",
                                                  meshing="GMSH", user_options=opts, create_files=True)
                sim.add_metallic(1)
                sim.add_ground_plane()
                sim.create_port_JosephsonJunction(qubit, L_J=L_J_nH * 1e-9, C_J=C_J_fF * 1e-15)
                for stub in c.stubs:
                    sim.create_port_CPW_on_Route(stub, "end")
                sim.fine_mesh_components(list(sub.components.keys()), metals_only=True, **mesh)
                if epr_interfaces:
                    sim.setup_EPR_interfaces(metal_air=MaterialInterface("Aluminium-Vacuum"),
                                             substrate_air=MaterialInterface("Silicon-Vacuum"),
                                             substrate_metal=MaterialInterface("Silicon-Aluminium"))
                if kinetic_inductance_pH_sq:
                    sim.add_kinetic_inductance(kinetic_inductance_pH_sq * 1e-12)
                sim.prepare_simulation()
            if eig_csv.exists():
                eig_csv.unlink()        # results of other settings must not be read back
            record_path.write_text(json.dumps(dict(key=key, qubit=qubit, label=label, settings=settings,
                                                   config=str(folder / f"{sim_name}.json"),
                                                   created=datetime.now().isoformat(timespec="seconds")),
                                              indent=2, default=str))
            if not run:
                return dict(output_dir=str(folder / "outputFiles"), status="prepared", key=key)
            sim.run()
            if self._eig_rows(eig_csv) < 1:
                raise RuntimeError(f"Palace wrote no converged modes to {eig_csv}; see {folder / 'outputFiles' / 'out.log'}")
        return dict(self.eigenmode_results(folder, qubit_mode=None, f_r_GHz=f_r_GHz), key=key,
                    status="cached" if cached else "simulated")

    @staticmethod
    def _eig_rows(eig_csv):
        """Number of modes in a Palace eig.csv (0 if missing; Palace writes the header before solving,
        so an interrupted run leaves a header-only file)."""
        eig_csv = Path(eig_csv)
        if not eig_csv.exists():
            return 0
        return sum(1 for line in eig_csv.read_text().splitlines()[1:] if line.strip())

    @staticmethod
    def eigenmode_results(folder, qubit_mode=None, f_r_GHz=None):
        """Reads an eigenmode simulation folder (made by eigenmode()): mode table, EPR Hamiltonian
        parameters and a summary for the qubit mode (largest junction participation unless given) and the
        resonator mode (nearest ``f_r_GHz`` if given, else the other mode with the most junction
        participation)."""
        import pandas as pd
        from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation
        folder = Path(folder)
        out_dir = folder / "outputFiles"
        config = next(p for p in folder.glob("*.json") if p.name != "eigenmode.json")
        eig = pd.read_csv(out_dir / "eig.csv", skipinitialspace=True)
        eig.columns = [c.strip() for c in eig.columns]
        f = eig["Re{f} (GHz)"].to_numpy()
        Q = eig["Q"].to_numpy()
        with contextlib.redirect_stdout(io.StringIO()):
            epr = PALACE_Eigenmode_Simulation.calculate_hamiltonian_parameters_EPR_from_files(
                str(out_dir), str(config), print_output=False)
        p_jj = np.abs(epr["EPR"].to_numpy()[:, 0])

        def combined_Q(csv, prefix):
            """1 / sum(1/Q) over the columns named prefix[...] of a Palace postprocessing csv (inf if absent)."""
            if not (out_dir / csv).exists():
                return np.full(len(f), np.inf)
            t = pd.read_csv(out_dir / csv, skipinitialspace=True)
            t.columns = [c.strip() for c in t.columns]
            cols = [c for c in t.columns if c.startswith(prefix)]
            return 1 / (1 / t[cols].to_numpy()[: len(f)]).sum(axis=1) if cols else np.full(len(f), np.inf)

        # Palace's Q (from Im f) holds every loss in the eigenproblem: the 50 Ohm ports and the bulk dielectric
        # loss tangent. Port losses alone are in port-Q.csv; the surface (TLS) interfaces are post-processed only
        # (surface-Q.csv) and are not part of Q.
        Q_ports = combined_Q("port-Q.csv", "Q_ext")
        Q_surf = combined_Q("surface-Q.csv", "Q_surf")
        with np.errstate(divide="ignore"):
            Q_bulk = 1 / np.clip(1 / Q - 1 / Q_ports, 1e-300, None)
        modes = pd.DataFrame(dict(f_GHz=f, Q=Q, Q_ports=Q_ports, Q_bulk=Q_bulk, Q_surfaces=Q_surf,
                                  kappa_MHz=f / Q * 1e3, kappa_ports_MHz=f / Q_ports * 1e3, EPR_JJ=p_jj),
                             index=pd.Index(range(1, len(f) + 1), name="mode"))
        T1_us = lambda q, m: q[m - 1] / (2 * np.pi * f[m - 1] * 1e9) * 1e6
        qm = qubit_mode or int(np.argmax(p_jj)) + 1
        others = [m for m in modes.index if m != qm]
        if f_r_GHz is not None:
            rm = min(others, key=lambda m: abs(modes.loc[m, "f_GHz"] - f_r_GHz)) if others else None
        else:
            rm = max(others, key=lambda m: modes.loc[m, "EPR_JJ"]) if others else None
        chi = epr["Chi"].to_numpy()             # MHz; diagonal = anharmonicity, off-diagonal = cross-Kerr
        f_norm = epr["f_norms_GHz"].to_numpy().ravel()
        summary = dict(qubit_mode=qm, f_q_GHz=f[qm - 1], f_q_dressed_GHz=f_norm[qm - 1],
                       anharmonicity_MHz=-abs(chi[qm - 1, qm - 1]), Q_q=Q[qm - 1],
                       T1_q_us=T1_us(Q, qm),                 # ports + bulk dielectric
                       T1_q_ports_us=T1_us(Q_ports, qm),     # radiation into the feedline: Purcell + direct
                       T1_q_bulk_us=T1_us(Q_bulk, qm),       # substrate loss tangent
                       T1_q_surfaces_us=T1_us(Q_surf, qm),   # SA / MS / MA interfaces (post-processed)
                       EPR_JJ_q=p_jj[qm - 1])
        if rm is not None:
            summary.update(resonator_mode=rm, f_r_GHz=f[rm - 1], f_r_dressed_GHz=f_norm[rm - 1],
                           chi_qr_MHz=chi[qm - 1, rm - 1], Q_r=Q[rm - 1], Q_r_ports=Q_ports[rm - 1],
                           kappa_r_MHz=f[rm - 1] / Q_ports[rm - 1] * 1e3,     # loaded linewidth from the ports
                           kappa_r_total_MHz=f[rm - 1] / Q[rm - 1] * 1e3,
                           detuning_GHz=f[qm - 1] - f[rm - 1])
        return dict(modes=modes, qubit_mode=qm, resonator_mode=rm, summary=summary, epr=epr,
                    output_dir=str(out_dir))

    def params(self, qubit, f_q_GHz, f_r_GHz=None, res=None, C_J_fF=0.0, print_all_capacitances=False):
        """Circuit parameters from the simulated cell, via
        ``PALACE_Capacitance_Simulation.calc_params_floating_Transmon_from_files`` with the detected
        conductor indices. The resonator is ``res`` if given, else a half-wave resonator at ``f_r_GHz``."""
        from SQDMetal.PALACE.Capacitance_Simulation import PALACE_Capacitance_Simulation
        rec = self.records.get(qubit)
        if rec is None or rec["status"] == "prepared":
            raise ValueError(f"No simulation results for {qubit}: run simulate() first.")
        if res is None and f_r_GHz is not None:
            res = ResonatorHalfWave(f_r_GHz * 1e9)
        with self._quiet():
            return PALACE_Capacitance_Simulation.calc_params_floating_Transmon_from_files(
                rec["output_dir"], qubit_freq=f_q_GHz * 1e9, res=res, conductor_indices=rec["conductor_indices"],
                C_J=C_J_fF * 1e-15, print_all_capacitances=print_all_capacitances)

    # ------------------------------------------------------------------ resonator frequencies
    def resonator_length_mm(self, qubit):
        """Centre-line length (mm) of the readout resonator of ``qubit``: the routed resonator(s) in its cell,
        plus the coupled arm of each CoupledLineTee (trace ``second_cpw``, along its fillet) and the qubit's readout
        wire (``get_resonator_length_mm``) where those exist."""
        from SQDMetal.Utilities.QUtilities import QUtilities
        c = self.cell(qubit)
        sub = c.design
        unit_mm = float(sub.parse_value("1mm"))
        total = 0.0
        for n in c.components:
            comp = sub.components[n]
            if n == qubit:
                if hasattr(comp, "get_resonator_length_mm"):
                    total += comp.get_resonator_length_mm() / unit_mm
            elif type(comp).__name__ == "CoupledLineTee":
                total += QUtilities.calc_length_of_path(sub, n, trace_name="second_cpw") / unit_mm
            elif "total_length" in comp.options or type(comp).__name__.startswith("Route"):
                total += QUtilities.calc_length_of_path(sub, n) / unit_mm
        return total

    def resonator_capacitance_fF(self, qubit):
        """Total (Maxwell self-) capacitance of the resonator conductor of ``qubit`` from its capacitance
        simulation, in fF: the line's capacitance to ground plus its end loading (the qubit pad and the
        feedline coupler)."""
        import pandas as pd
        rec = self.records.get(qubit)
        if rec is None or rec["status"] == "prepared":
            raise ValueError(f"No simulation results for {qubit}: run simulate() first.")
        M = pd.read_csv(Path(rec["output_dir"]) / "terminal-C.csv").iloc[:, 1:].to_numpy()
        i = rec["conductor_indices"]["res"]
        return float(M[i, i]) * 1e15

    @staticmethod
    def qubit_pull_GHz(f_r_GHz, f_q_GHz, g_MHz):
        """Linear (normal-mode) pull of the resonator by the qubit, g^2 / (f_r - f_q), in GHz (positive when the
        resonator is above the qubit). This is the shift a linear eigenmode solve already contains."""
        return (abs(g_MHz) * 1e-3) ** 2 / (f_r_GHz - f_q_GHz)

    @staticmethod
    def resonator_calibration(f_r_eig_GHz, length_mm, C_res_fF, pull_GHz=0.0, **info):
        """Calibration of all resonators of a chip from one eigenmode simulation. A lambda/2 CPW resonator of
        inductive length l and total capacitance C (line + end loading) sits at
        f = 1 / (2 sqrt(L' l C)) = K / sqrt(l C), with K = 1 / (2 sqrt(L')). L' is the inductance per length
        (geometric + kinetic), the same for every resonator with this CPW and film. K is fixed from the resonator mode of an eigenmode run, with the qubit's
        linear pull removed. Returns a dict (JSON-serialisable) with K and the inputs; extra ``info`` is kept."""
        f_bare = f_r_eig_GHz - pull_GHz
        return dict(K=f_bare * np.sqrt(length_mm * C_res_fF), f_r_eig_GHz=f_r_eig_GHz, f_r_bare_GHz=f_bare,
                    length_mm=length_mm, C_res_fF=C_res_fF, pull_GHz=pull_GHz, **info)

    @staticmethod
    def predict_resonator(calibration, length_mm, C_res_fF, pull_GHz=0.0):
        """Resonator frequency (GHz, linear mode incl. the qubit's pull) predicted from a
        ``resonator_calibration`` for a resonator of length ``length_mm`` and capacitance ``C_res_fF``."""
        return calibration["K"] / np.sqrt(length_mm * C_res_fF) + pull_GHz

    @staticmethod
    def fit_resonator_capacitance(length_mm, C_res_fF):
        """Straight-line fit C = C' l + C_0 of simulated resonator capacitances against resonator length.

        The capacitance of one meshed resonator scatters by about 1 % from mesh to mesh (enough to move its
        frequency by tens of MHz), while C' (the line) and C_0 (end loading, nearly the same for every resonator
        of a chip) are smooth. Fitting all simulated resonators averages the mesh scatter out. Returns a dict with
        ``C_per_mm_fF``, ``C0_fF``, ``rms_pct`` (scatter of the points about the fit) and ``n``."""
        length_mm, C_res_fF = np.asarray(length_mm, float), np.asarray(C_res_fF, float)
        if length_mm.size < 2 or np.ptp(length_mm) <= 0:
            raise ValueError("Fitting C(l) needs at least two resonators of different length.")
        a, b = np.polyfit(length_mm, C_res_fF, 1)
        resid = C_res_fF / (a * length_mm + b) - 1
        return dict(C_per_mm_fF=float(a), C0_fF=float(b), rms_pct=float(np.sqrt(np.mean(resid ** 2)) * 100),
                    n=int(length_mm.size))

    @staticmethod
    def fitted_capacitance_fF(fit, length_mm):
        """Resonator capacitance (fF) at ``length_mm`` from a ``fit_resonator_capacitance`` fit."""
        return fit["C_per_mm_fF"] * np.asarray(length_mm, float) + fit["C0_fF"]

    @staticmethod
    def resonator_length_for(calibration, fit, f_GHz, pull_GHz=0.0):
        """Resonator length (mm) that puts the resonator at ``f_GHz`` (linear mode incl. the qubit's pull):
        the positive root of K / sqrt(l (C' l + C_0)) + pull = f, i.e. C' l^2 + C_0 l - (K / (f - pull))^2 = 0.
        ``calibration`` must have been made with the fitted capacitance (see ``fitted_capacitance_fF``)."""
        a, b = fit["C_per_mm_fF"], fit["C0_fF"]
        c = (calibration["K"] / (f_GHz - pull_GHz)) ** 2
        return float((-b + np.sqrt(b * b + 4 * a * c)) / (2 * a))

    @staticmethod
    def corrected_design_frequency(f_target_GHz, f_design_used_GHz, f_predicted_GHz):
        """Frequency to pass to the chip builder so the resonator lands on ``f_target_GHz``: the builder's length
        scales as 1 / f_design and f_predicted ~ f_design locally, so
        f_design_new = f_target * f_design_used / f_predicted (one step is exact to first order; re-simulate
        the capacitance once to confirm)."""
        return f_target_GHz * f_design_used_GHz / f_predicted_GHz

    def to_catalogue(self, catalogue, qubit, f_q_GHz=None, f_r_GHz=None, res=None, C_J_fF=0.0, overwrite=False):
        """Saves the simulated cell of ``qubit`` as an entry of a
        :class:`~SQDMetal.Utilities.QubitCatalogue.QubitCatalogue` (see ``QubitCatalogue.add_simulation``):
        the qubit's full options from the chip, the cell / mesh / solver settings, the conductor mapping
        and the small result files. If ``f_q_GHz`` and a resonator (``f_r_GHz`` or ``res``) are given,
        g, chi, Delta, anharmonicity, kappa and the Purcell limit are stored too (as result.* columns).
        Returns the catalogue key."""
        rec = self.records.get(qubit)
        if rec is None or rec["status"] == "prepared":
            raise ValueError(f"No simulation results for {qubit}: run simulate() first.")
        c = self.cell(qubit)
        q = c.design.components[qubit]
        results = {}
        if f_q_GHz is not None and (res is not None or f_r_GHz is not None):
            p = self.params(qubit, f_q_GHz=f_q_GHz, f_r_GHz=f_r_GHz, res=res, C_J_fF=C_J_fF)
            if res is not None:
                f_r_GHz = res.get_res_frequency() * 1e-9
            results = dict(f_q_GHz=f_q_GHz, f_r_GHz=f_r_GHz, C_J_fF=C_J_fF, g_MHz=p["g_MHz"], chi_MHz=p["chi_MHz"],
                           Delta_GHz=p["Delta_GHz"], anh_MHz=p["anh_MHz"], kappa_MHz=p["kappa_MHz"],
                           T1_purcell_ms=p["T1,p_ms"], E_J_GHz=p["E_J_GHz"], L_J_nH=p["L_J_nH"])
        others = [n for n in c.components if n != qubit]
        cell = dict(builder="chip", readout_pin="readout" if "readout" in q.pins else "",
                    components=",".join(type(c.design.components[n]).__name__ for n in others),
                    feedline_stub=self.feedline_stub if c.stubs else "", chip_margin_mm=self.chip_margin_mm)
        source = dict(type="chip", design=self.design.get_design_name(), qubit=qubit, label=rec.get("label", qubit),
                      cell_key=rec["key"])
        components = {n: dict(cls=type(c.design.components[n]).__name__, options=c.design.components[n].options)
                      for n in others}
        return catalogue.add_simulation(Path(rec["output_dir"]).parent, q.options, rec["conductor_indices"], cell,
                                        mesh=self.mesh, sim_options=self.user_options, qubit_class=type(q),
                                        source=source, results=results, components=components, overwrite=overwrite)

    # ------------------------------------------------------------------ #
    # Geometry helpers
    @staticmethod
    def metal_geom(design, name, subtract=False, layer="1"):
        """Union of a component's geometry on ``layer`` (paths buffered to their width); ``subtract=True``
        gives its ground-plane cut-outs instead."""
        shapes = []
        for kind in ("poly", "path"):
            df = design.qgeometry.get_component(name, kind)
            if df is None or len(df) == 0:
                continue
            for _, row in df.iterrows():
                if bool(row["subtract"]) != subtract or str(row["layer"]) != str(layer):
                    continue
                g = row["geometry"]
                if kind == "path":
                    g = g.buffer(row["width"] / 2, cap_style=2, join_style=2)
                shapes.append(g)
        return shapely.unary_union(shapes) if shapes else shapely.Polygon()

    @staticmethod
    def check_geometry(full, sub, names, tol_mm2=1e-9):
        """Largest symmetric-difference area (design units^2) between ``full`` and ``sub`` for each of
        ``names`` (metal and cut-outs). Raises AssertionError above ``tol_mm2``."""
        g = SingleQubitCells.metal_geom
        diff = {n: max(g(full, n, s).symmetric_difference(g(sub, n, s)).area for s in (False, True)) for n in names}
        bad = {n: a for n, a in diff.items() if a > tol_mm2}
        if bad:
            raise AssertionError(f"Cell geometry differs from the full chip: {bad}")
        return diff

    @staticmethod
    def identify_conductors(sub, qubit, resonator, stubs, terminals):
        """Maps capacitance-matrix indices (order of ``terminals``, i.e. ``sim._cur_cap_terminals``) to
        roles by geometry: ground (largest), res (nearest the ``resonator`` component), feed (nearest the
        first feedline stub) and pad1/pad2 (the remaining conductors touching the qubit; pad1 is nearer
        the qubit's 'readout' pin). ``resonator`` / ``stubs`` may be None / empty for cells without them."""
        t = [shapely.make_valid(x) for x in terminals]
        nearest = lambda geom: int(np.argmin([x.distance(geom) for x in t]))
        idx = dict(ground=int(np.argmax([x.area for x in t])))
        g = SingleQubitCells.metal_geom
        if resonator is not None:
            idx["res"] = nearest(g(sub, resonator).representative_point())
        if stubs:
            idx["feed"] = nearest(g(sub, stubs[0]).representative_point())
        taken = set(idx.values())
        q_metal = g(sub, qubit)
        pads = [i for i, x in enumerate(t) if i not in taken and x.intersects(q_metal)]
        if len(taken) != len(idx) or len(pads) != 2 or len(t) != len(idx) + 2:
            raise ValueError(f"Could not identify the conductors: {len(t)} conductors, roles {idx}, pads {pads}. "
                             f"Check sim.display_conductor_indices().")
        pins = sub.components[qubit].pins
        if "readout" in pins:
            ref = shapely.Point(np.array(pins["readout"]["middle"]))
            pads.sort(key=lambda i: t[i].distance(ref))
        idx["pad1"], idx["pad2"] = pads
        return idx

    @classmethod
    def plot_conductors(cls, terminals, idx, ax, title=""):
        """Draws each conductor coloured by role, marked with its 0-based capacitance-matrix index."""
        role_of = {i: r for r, i in idx.items()}
        for i, t in enumerate(terminals):
            role = role_of.get(i, "?")
            colour = cls.ROLE_COLOURS.get(role, "k")
            z = 1 if role == "ground" else 2
            for poly in getattr(t, "geoms", [t]):
                if poly.geom_type != "Polygon":
                    continue
                ax.fill(*poly.exterior.xy, color=colour, lw=0, zorder=z)
                for hole in poly.interiors:
                    ax.fill(*hole.xy, color="white", lw=0, zorder=z)
            ax.plot([], [], "s", color=colour, label=f"{i}: {cls.ROLE_LABELS.get(role, role)}")
            if role != "ground":
                p = t.representative_point()
                ax.annotate(str(i), (p.x, p.y), ha="center", va="center", fontsize=9, fontweight="bold",
                            color="white", bbox=dict(boxstyle="circle,pad=0.2", fc="black", ec="none"), zorder=3)
        ax.set_aspect("equal")
        ax.set_title(title)
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2, frameon=False)

    def fingerprint(self, qubit):
        """Short hash of the layer-1 geometry that sets ``qubit``'s capacitances: the class and options
        of every component in its cell (the junction, on another layer, is not included). Store it with
        simulation results to tell later whether the chip's geometry has changed since."""
        spec = {n: dict(cls=type(self.design.components[n]).__name__, options=self.design.components[n].options)
                for n in self.components(qubit)}
        return hashlib.sha1(json.dumps(spec, sort_keys=True, default=str).encode()).hexdigest()[:12]

    # ------------------------------------------------------------------ #
    # Junctions
    @staticmethod
    def junction_for(f01_GHz, E_C_GHz, ng=0.25):
        """Junction giving a transmon with charging energy ``E_C_GHz`` the transition frequency
        ``f01_GHz``, solved numerically (charge basis, see transmon_levels) rather than with the
        asymptotic sqrt(8 E_J E_C) - E_C. Returns dict(E_J_GHz, L_J_nH, E_J_over_E_C, f01_GHz, alpha_MHz)."""
        from scipy.optimize import brentq
        f01 = lambda EJ: SingleQubitCells.transmon_levels(EJ, E_C_GHz, ng)[1]
        lo, hi = 0.5 * E_C_GHz, max(10.0, (f01_GHz + E_C_GHz) ** 2 / (8 * E_C_GHz) * 4)
        if not f01(lo) < f01_GHz < f01(hi):
            raise ValueError(f"f01 = {f01_GHz} GHz is out of reach for E_C = {E_C_GHz} GHz "
                             f"(E_J/E_C between 0.5 and {hi / E_C_GHz:.0f} gives "
                             f"{f01(lo):.3g}-{f01(hi):.3g} GHz).")
        EJ = brentq(lambda x: f01(x) - f01_GHz, lo, hi, xtol=1e-12)
        E = SingleQubitCells.transmon_levels(EJ, E_C_GHz, ng)
        h, e = 6.62607015e-34, 1.602176634e-19
        L_J_nH = (h / (2 * e) / (2 * np.pi)) ** 2 / (h * EJ * 1e9) * 1e9
        return dict(E_J_GHz=EJ, L_J_nH=L_J_nH, E_J_over_E_C=EJ / E_C_GHz, f01_GHz=E[1],
                    alpha_MHz=(E[2] - 2 * E[1]) * 1e3)

    # ------------------------------------------------------------------ #
    # Misc
    @staticmethod
    def transmon_levels(EJ, EC, ng=0.25, ncut=20):
        """Lowest transmon energy levels (relative to the ground state, units of EJ/EC) from
        H = 4EC(n - ng)^2 - EJ cos(phi) in the charge basis."""
        n = np.arange(-ncut, ncut + 1)
        H = np.diag(4 * EC * (n - ng) ** 2) - EJ / 2 * (np.eye(len(n), k=1) + np.eye(len(n), k=-1))
        E = np.linalg.eigvalsh(H)
        return E - E[0]

    @staticmethod
    def design_from_notebook(path, stop_marker, verbose=False):
        """Runs a design notebook's code cells (lines starting with % or ! skipped) up to and including the
        first cell containing ``stop_marker``, in the notebook's folder, and returns its namespace (e.g.
        ``ns['design']``). Keeps the chip notebook as the single source of truth for its options."""
        import nbformat
        path = Path(path)
        nb = nbformat.read(path, as_version=4)
        ns = {"__name__": "__design_notebook__"}
        cwd = os.getcwd()
        os.chdir(path.parent)
        try:
            out = contextlib.nullcontext() if verbose else contextlib.redirect_stdout(io.StringIO())
            with out:
                for i, cell in enumerate(nb.cells):
                    if cell.cell_type != "code":
                        continue
                    src = "\n".join(l for l in cell.source.splitlines() if not l.lstrip().startswith(("%", "!")))
                    exec(compile(src, f"{path.name} [cell {i}]", "exec"), ns)
                    if stop_marker in cell.source:
                        return ns
        finally:
            os.chdir(cwd)
        raise ValueError(f"No code cell in {path.name} contains {stop_marker!r}")
