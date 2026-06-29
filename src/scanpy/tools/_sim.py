# Author: Alex Wolf (https://falexwolf.de)
"""Simulate Data.

Simulate stochastic dynamic systems to model gene expression dynamics and
cause-effect data.

Todo:
----
Beta Version. The code will be reorganized soon.

"""

from __future__ import annotations

import itertools
from pathlib import Path
from types import MappingProxyType
from typing import TYPE_CHECKING

import numpy as np
import scipy as sp

from .. import _utils, readwrite
from .. import logging as logg
from .._docs import doc_rng
from .._settings import settings
from .._utils import _doc_params
from .._utils.random import _if_legacy_apply_global, _LegacyRng

if TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import ClassVar, Literal

    from anndata import AnnData

    from .._utils.random import RNGLike, SeedLike


@_doc_params(rng=doc_rng)
def sim(  # noqa: PLR0913
    model: Literal["krumsiek11", "toggleswitch"],
    *,
    params_file: bool = True,
    tmax: int | None = None,
    branching: bool | None = None,
    nrRealizations: int | None = None,
    noiseObs: float | None = None,
    noiseDyn: float | None = None,
    step: int | None = None,
    rng: SeedLike | RNGLike | None = None,
    writedir: Path | str | None = None,
    # deprecated
    seed: int | None = None,
) -> AnnData:
    """Simulate dynamic gene expression data :cite:p:`Wittmann2009,Wolf2018`.

    Sample from a stochastic differential equation model built from
    literature-curated boolean gene regulatory networks, as suggested by
    :cite:t:`Wittmann2009`. The Scanpy implementation can be found in :cite:t:`Wolf2018`.

    Parameters
    ----------
    model
        Model file in 'sim_models' directory.
    params_file
        Read default params from file.
    tmax
        Number of time steps per realization of time series.
    branching
        Only write realizations that contain new branches.
    nrRealizations
        Number of realizations.
    noiseObs
        Observatory/Measurement noise.
    noiseDyn
        Dynamic noise.
    step
        Interval for saving state of system.
    {rng}
    writedir
        Path to directory for writing output files.

    Returns
    -------
    Annotated data matrix.

    Examples
    --------
    See this `use case <https://github.com/scverse/scanpy_usage/tree/master/170430_krumsiek11>`__

    """
    params = locals()
    if seed is not None and rng is not None:
        msg = "Cannot specify both `seed` and `rng`."
        raise TypeError(msg)
    if params_file:
        model_key = Path(model).with_suffix("").name
        from .. import sim_models

        pfile_sim = Path(sim_models.__file__).parent / f"{model_key}_params.txt"
        default_params = readwrite.read_params(pfile_sim)
        params = _utils.update_params(default_params, params)
    params["rng"] = np.random.default_rng(rng) if rng is not None else _LegacyRng(seed)
    del params["seed"]
    adata = sample_dynamic_data(**params)
    adata.uns["iroot"] = 0
    return adata


def add_args(p):
    """Update parser with tool specific arguments.

    This overwrites was is done in utils.uns_args.
    """
    # dictionary for adding arguments
    dadd_args = {
        "--opfile": {
            "default": "",
            "metavar": "f",
            "type": str,
            "help": 'Specify a parameter file (default: "sim/${exkey}_params.txt")',
        }
    }
    p = _utils.add_args(p, dadd_args)
    return p


def sample_dynamic_data(**params):  # noqa: PLR0912, PLR0915
    rng: np.random.Generator = params["rng"]
    model_key = Path(params["model"]).with_suffix("").name
    writedir = params.get("writedir")
    if writedir is None:
        writedir = settings.writedir / f"{model_key}_sim"
    else:
        writedir = Path(writedir)
    writedir.mkdir(parents=True, exist_ok=True)
    readwrite.write_params(writedir / "params.txt", params)
    # init variables
    tmax = params["tmax"]
    branching = params["branching"]
    noiseObs = params["noiseObs"]
    noiseDyn = params["noiseDyn"]
    nrRealizations = params["nrRealizations"]
    step = params["step"]  # step size for saving the figure

    nrSamples = 1  # how many files?
    maxRestarts = 1000
    maxNrSamples = 1

    # simple vector auto regressive process or
    # hill kinetics process simulation
    if "krumsiek11" not in model_key:
        # create instance, set seed
        grnsim = GRNsim(model=model_key, params=params)
        nrOffEdges_list = np.zeros(nrSamples)
        for sample in range(nrSamples):
            # random topology / for a given edge density
            if "hill" not in model_key:
                Coupl = np.array(grnsim.Coupl)
                for _sampleCoupl in range(10):
                    nrOffEdges = 0
                    for gp in range(grnsim.dim):
                        for g in range(grnsim.dim):
                            # only consider off-diagonal edges
                            if g != gp:
                                Coupl[gp, g] = 0.7 if rng.random() < 0.4 else 0
                                nrOffEdges += 1 if Coupl[gp, g] > 0 else 0
                            else:
                                Coupl[gp, g] = 0.7
                    # check that the coupling matrix does not have eigenvalues
                    # greater than 1, which would lead to an exploding var process
                    if max(sp.linalg.eig(Coupl)[0]) < 1:
                        break
                nrOffEdges_list[sample] = nrOffEdges
                grnsim.set_coupl(Coupl)
            # init type
            real = 0
            X0 = rng.random(grnsim.dim)
            Xsamples = []
            for restart in range(nrRealizations + maxRestarts):
                # slightly break symmetry in initial conditions
                if "toggleswitch" in model_key:
                    X0 = np.array([
                        0.8 for i in range(grnsim.dim)
                    ]) + 0.01 * rng.standard_normal(grnsim.dim)
                X = grnsim.sim_model(tmax=tmax, X0=X0, noiseDyn=noiseDyn)
                # check branching
                check = True
                if branching:
                    check, Xsamples = _check_branching(X, Xsamples, restart)
                if check:
                    real += 1
                    grnsim.write_data(
                        X[::step],
                        dir=writedir,
                        noiseObs=noiseObs,
                        append=restart != 0,
                        branching=branching,
                        nrRealizations=nrRealizations,
                    )
                # append some zeros
                if "zeros" in writedir.name and real == 2:
                    grnsim.write_data(
                        noiseDyn * rng.standard_normal((500, 3)),
                        dir=writedir,
                        noiseObs=noiseObs,
                        append=restart != 0,
                        branching=branching,
                        nrRealizations=nrRealizations,
                    )
                if real >= nrRealizations:
                    break
        logg.debug(
            f"mean nr of offdiagonal edges {nrOffEdges_list.mean()} "
            f"compared to total nr {grnsim.dim * (grnsim.dim - 1) / 2.0}"
        )

    # more complex models
    else:
        initType = "random"

        dim = 11
        step = 5

        grnsim = GRNsim(dim=dim, initType=initType, model=model_key, params=params)
        Xsamples = []
        for _sample in range(maxNrSamples):
            # choose initial conditions such that branchings result
            if initType == "branch":
                X0mean = grnsim.branch_init_model1(tmax)
                if X0mean is None:
                    grnsim.set_coupl()
                    continue
            real = 0
            for restart in range(nrRealizations + maxRestarts):
                if initType == "branch":
                    # vary initial conditions around mean
                    X0 = X0mean + (0.05 * rng.random(dim) - 0.025 * np.ones(dim))
                else:
                    # generate random initial conditions within [0.3,0.7]
                    X0 = 0.4 * rng.random(dim) + 0.3
                if model_key in [5, 6]:
                    X0 = np.array([0.3, 0.3, 0, 0, 0, 0])
                if model_key in [7, 8, 9, 10]:
                    X0 = 0.6 * rng.random(dim) + 0.2
                    X0[2:] = np.zeros(4)
                if "krumsiek11" in model_key:
                    X0 = np.zeros(dim)
                    X0[grnsim.varNames["Gata2"]] = 0.8
                    X0[grnsim.varNames["Pu.1"]] = 0.8
                    X0[grnsim.varNames["Cebpa"]] = 0.8
                    X0 += 0.001 * rng.standard_normal(dim)
                    if False:
                        switch_gene = restart - (nrRealizations - dim)
                        if switch_gene >= dim:
                            break
                        X0[switch_gene] = 0 if X0[switch_gene] > 0.1 else 0.8
                X = grnsim.sim_model(tmax, X0=X0, noiseDyn=noiseDyn, restart=restart)
                # check branching
                check = True
                if branching:
                    check, Xsamples = _check_branching(X, Xsamples, restart)
                if check:
                    real += 1
                    grnsim.write_data(
                        X[::step],
                        dir=writedir,
                        noiseObs=noiseObs,
                        append=restart != 0,
                        branching=branching,
                        nrRealizations=nrRealizations,
                    )
                if real >= nrRealizations:
                    break
    # load the last simulation file
    filename = max(writedir.glob("sim*.txt"))
    logg.info(f"reading simulation results {filename}")
    adata = readwrite.read(
        filename, first_column_names=True, suppress_cache_warning=True
    )
    adata.uns["tmax_write"] = tmax / step
    return adata


def write_data(  # noqa: PLR0912, PLR0913
    X,
    dir=Path("sim/test"),
    *,
    append=False,
    header="",
    varNames: Mapping[str, int] = MappingProxyType({}),
    Adj: np.ndarray | None = None,
    Coupl: np.ndarray | None = None,
    boolRules: Mapping[str, str] = MappingProxyType({}),
    model="",
    modelType="",
    invTimeStep=1,
):
    """Write simulated data.

    Accounts for saving at the same time an ID
    and a model file.
    """
    dir.mkdir(parents=True, exist_ok=True)
    # update file with sample ids
    filename = dir / "id.txt"
    if filename.is_file():
        with filename.open("r") as f:
            id = int(f.read()) + (0 if append else 1)
    else:
        id = 0
    with filename.open("w") as f:
        id = f"{id:0>6}"
        f.write(str(id))
    # dimension
    dim = X.shape[1]
    # write files with adjacancy and coupling matrices
    if not append:
        if False:
            if Adj is not None:
                # due to 'update formulation' of model, there
                # is always a diagonal dependence
                Adj = np.copy(Adj)
                if "hill" in model:
                    for i in range(Adj.shape[0]):
                        Adj[i, i] = 1
                np.savetxt(f"{dir}/adj_{id}.txt", Adj, header=header, fmt="%d")
            if Coupl is not None:
                np.savetxt(f"{dir}/coupl_{id}.txt", Coupl, header=header, fmt="%10.6f")
        # write model file
        if varNames and Coupl is not None:
            with (dir / f"model_{id}.txt").open("w") as f:
                f.write('# For each "variable = ", there must be a right hand side: \n')
                f.write(
                    "# either an empty string or a python-style logical expression \n"
                )
                f.write('# involving variable names, "or", "and", "(", ")". \n')
                f.write("# The order of equations matters! \n")
                f.write("# \n")
                f.write(f"# modelType = {modelType}\n")
                f.write(f"# invTimeStep = {invTimeStep}\n")
                f.write("# \n")
                f.write("# boolean update rules: \n")
                for k, v in boolRules.items():
                    f.write(f"{k} = {v}\n")
                # write coupling via names
                f.write("# coupling list: \n")
                names = list(varNames.keys())
                for gp in range(dim):
                    for g in range(dim):
                        if np.abs(Coupl[gp, g]) > 1e-10:
                            f.write(
                                f"{names[gp]:10} {names[g]:10} {Coupl[gp, g]:10.3} \n"
                            )
    # write simulated data
    # the binary mode option in the following line is a fix for python 3
    # variable names
    if varNames:
        header += f"{'it':>2} "
        for v in varNames:
            header += f"{v:>7} "
    with (dir / f"sim_{id}.txt").open("ab" if append else "wb") as f:
        np.savetxt(
            f,
            np.c_[np.arange(0, X.shape[0]), X],
            header=("" if append else header),
            fmt=["%4.f"] + ["%7.4f" for i in range(dim)],
        )


class GRNsim:
    """Simlulation of stochastic dynamic systems.

    Main application: simulation of gene expression dynamics.

    Also standard models are implemented.
    """

    availModels: ClassVar = dict(
        krumsiek11=(
            "myeloid progenitor network, Krumsiek et al., PLOS One 6, e22649, "
            "\n      equations from Table 1 on page 3, "
            "doi:10.1371/journal.pone.0022649 \n"
        ),
        var="vector autoregressive process \n",
        hill="process with hill kinetics \n",
    )

    writeOutputOnce = True

    def __init__(
        self,
        *,
        dim=3,
        model="ex0",
        modelType="var",
        initType="random",
        show=False,
        verbosity=0,
        Coupl=None,
        params=MappingProxyType({}),
    ):
        """Initialize.

        Params
        ------
        model
            either string for predefined model,
            or directory with a model file and a couple matrix files
        """
        # number of nodes / dimension of system
        self.dim = dim if Coupl is None else Coupl.shape[0]
        self.maxnpar = 1  # maximal number of parents
        self.p_indep = 0.4  # fraction of independent genes
        self.model = model
        self.modelType = modelType
        self.initType = initType  # string characterizing a specific initial
        self.show = show
        self.verbosity = verbosity
        # checks
        if initType not in ["branch", "random"]:
            msg = "initType must be either: branch, random"
            raise RuntimeError(msg)
        if model not in self.availModels:
            message = "model not among predefined models \n"  # noqa: F841  # TODO FIX
        # read from file
        from .. import sim_models

        model = Path(sim_models.__file__).parent / f"{model}.txt"
        if not model.is_file():
            msg = f"Model file {model} does not exist"
            raise RuntimeError(msg)
        self.model = model
        # set the coupling matrix, and with that the adjacency matrix
        self.set_coupl(Coupl=Coupl)
        # seed
        _if_legacy_apply_global(params["rng"])
        # header
        self.header = f"model = {self.model.name} \n"
        # params
        self.params = params

    def sim_model(self, tmax, X0, noiseDyn=0, restart=0):
        """Simulate the model."""
        self.noiseDyn = noiseDyn
        X = np.zeros((tmax, self.dim))
        X[0] = X0 + noiseDyn * self.params["rng"].standard_normal(self.dim)
        # run simulation
        for t in range(1, tmax):
            if self.modelType == "hill":
                Xdiff = self.Xdiff_hill(X[t - 1])
            elif self.modelType == "var":
                Xdiff = self.Xdiff_var(X[t - 1])
            else:
                msg = f"Unknown modelType {self.modelType!r}"
                raise ValueError(msg)
            X[t] = X[t - 1] + Xdiff
            # add dynamic noise
            X[t] += noiseDyn * self.params["rng"].standard_normal(self.dim)
        return X

    def Xdiff_hill(self, Xt):
        """Build Xdiff from coefficients of boolean network.

        That is, using self.boolCoeff. The employed functions
        are Hill type activation and deactivation functions.

        See Wittmann et al., BMC Syst. Biol. 3, 98 (2009),
        doi:10.1186/1752-0509-3-98 for more details.
        """
        verbosity = self.verbosity > 0 and self.writeOutputOnce
        self.writeOutputOnce = False
        Xdiff = np.zeros(self.dim)
        for ichild, child in enumerate(self.pas.keys()):
            # check whether list of parents is non-empty,
            # otherwise continue
            if self.pas[child]:
                Xdiff_syn = 0  # synthesize term
                if verbosity > 0:
                    Xdiff_syn_str = ""
            else:
                continue
            # loop over all tuples for which the boolean update
            # rule returns true, these are stored in self.boolCoeff
            for ituple, tuple in enumerate(self.boolCoeff[child]):
                Xdiff_syn_tuple = 1
                Xdiff_syn_tuple_str = ""
                for iv, v in enumerate(tuple):
                    iparent = self.varNames[self.pas[child][iv]]
                    x = Xt[iparent]
                    threshold = 0.1 / np.abs(self.Coupl[ichild, iparent])
                    Xdiff_syn_tuple *= (
                        self.hill_a(x, threshold) if v else self.hill_i(x, threshold)
                    )
                    if verbosity > 0:
                        Xdiff_syn_tuple_str += (
                            f"{'a' if v else 'i'}"
                            f"({self.pas[child][iv]}, {threshold:.2})"
                        )
                Xdiff_syn += Xdiff_syn_tuple
                if verbosity > 0:
                    Xdiff_syn_str += ("+" if ituple != 0 else "") + Xdiff_syn_tuple_str
            # multiply with degradation term
            Xdiff[ichild] = self.invTimeStep * (Xdiff_syn - Xt[ichild])
            if verbosity > 0:
                Xdiff_str = (
                    f"{child}_{child}-{child} = "
                    f"{self.invTimeStep}*({Xdiff_syn_str}-{child})"
                )
                settings.m(0, Xdiff_str)
        return Xdiff

    def Xdiff_var(self, Xt, verbosity=0):
        # subtract the current state
        Xdiff = -Xt
        # add the information from the past
        Xdiff += np.dot(self.Coupl, Xt)
        return Xdiff

    def hill_a(self, x, threshold=0.1, power=2):
        """Activating hill function."""
        x_pow = np.power(x, power)
        threshold_pow = np.power(threshold, power)
        return x_pow / (x_pow + threshold_pow)

    def hill_i(self, x, threshold=0.1, power=2):
        """Inhibiting hill function.

        Is equivalent to 1-hill_a(self,x,power,threshold).
        """
        x_pow = np.power(x, power)
        threshold_pow = np.power(threshold, power)
        return threshold_pow / (x_pow + threshold_pow)

    def nhill_a(self, x, threshold=0.1, power=2, ichild=2):
        """Normalized activating hill function."""  # noqa: D401
        x_pow = np.power(x, power)
        threshold_pow = np.power(threshold, power)
        return x_pow / (x_pow + threshold_pow) * (1 + threshold_pow)

    def nhill_i(self, x, threshold=0.1, power=2):
        """Normalized inhibiting hill function.

        Is equivalent to 1-nhill_a(self,x,power,threshold).
        """  # noqa: D401
        x_pow = np.power(x, power)
        threshold_pow = np.power(threshold, power)
        return threshold_pow / (x_pow + threshold_pow) * (1 - x_pow)

    def read_model(self):
        """Read the model and the couplings from the model file."""
        if self.verbosity > 0:
            settings.m(0, "reading model", self.model)
        # read model
        boolRules = []
        with self.model.open() as f:
            for line in f:
                if line.startswith("#") and "modelType =" in line:
                    keyval = line
                    if "|" in line:
                        keyval, _type = line.split("|")[:2]
                    self.modelType = keyval.split("=")[1].strip()
                if line.startswith("#") and "invTimeStep =" in line:
                    keyval = line
                    if "|" in line:
                        keyval, _type = line.split("|")[:2]
                    self.invTimeStep = float(keyval.split("=")[1].strip())
                if not line.startswith("#"):
                    boolRules.append([s.strip() for s in line.split("=")])
                if line.startswith("# coupling list:"):
                    break
        self.dim = len(boolRules)
        self.boolRules = dict(boolRules)
        self.varNames = {s: i for i, s in enumerate(self.boolRules.keys())}
        names = self.varNames
        # read couplings via names
        self.Coupl = np.zeros((self.dim, self.dim))
        reading = False
        with self.model.open() as f:
            for line in f:  # open(self.model.replace('/model','/couplList')):
                if line.startswith("# coupling list:"):
                    reading = True
                if not reading:
                    continue
                if not line.startswith("#"):
                    gps, gs, val = line.strip().split()
                    self.Coupl[int(names[gps]), int(names[gs])] = float(val)
        # adjancecy matrices
        self.Adj_signed = np.sign(self.Coupl)
        self.Adj = np.abs(np.array(self.Adj_signed))
        # build bool coefficients (necessary for odefy type
        # version of the discrete model)
        self.build_boolCoeff()

    def set_coupl(self, Coupl=None) -> None:
        """Construct the coupling matrix (and adjacancy matrix) from predefined models or via sampling."""
        self.varNames = {str(i): i for i in range(self.dim)}
        if self.model not in self.availModels and Coupl is None:
            self.read_model()
        elif "var" in self.model.name:
            # vector auto regressive process
            self.Coupl = Coupl
            self.boolRules = dict.fromkeys(self.varNames, "")
            names = list(self.varNames.keys())
            for gp in range(self.dim):
                pas = [
                    names[g]
                    for g in range(self.dim)
                    if np.abs(self.Coupl[gp, g] > 1e-10)
                ]
                self.boolRules[names[gp]] = "".join(
                    pas[:1] + [" or " + pa for pa in pas[1:]]
                )
                self.Adj_signed = np.sign(Coupl)
        elif self.model in ["6", "7", "8", "9", "10"]:
            self.Adj_signed = np.zeros((self.dim, self.dim))
            n_sinknodes = 2
            #             sinknodes = rng.choice(self.dim, n_sinknodes, replace=False)
            sinknodes = np.array([0, 1])
            # assume sinknodes have feeback
            self.Adj_signed[sinknodes, sinknodes] = np.ones(n_sinknodes)
            #             # allow negative feedback
            #             if self.model == 10:
            #                 plus_minus = (rng.integers(0, 2, n_sinknodes) - 0.5) * 2
            #                 self.Adj_signed[sinknodes, sinknodes] = plus_minus
            leafnodes = np.array(sinknodes)
            availnodes = np.array([i for i in range(self.dim) if i not in sinknodes])
            #             settings.m(0,leafnodes,availnodes)
            while len(availnodes) != 0:
                # parent
                parent_idx = self.params["rng"].choice(
                    np.arange(0, len(leafnodes)), size=1, replace=False
                )
                parent = leafnodes[parent_idx]
                # children
                children_ids = self.params["rng"].choice(
                    np.arange(0, len(availnodes)), size=2, replace=False
                )
                children = availnodes[children_ids]
                settings.m(0, parent, children)
                self.Adj_signed[children, parent] = np.ones(2)
                if self.model == 8:
                    self.Adj_signed[children[0], children[1]] = -1
                if self.model in [9, 10]:
                    self.Adj_signed[children[0], children[1]] = -1
                    self.Adj_signed[children[1], children[0]] = -1
                # update leafnodes
                leafnodes = np.delete(leafnodes, parent_idx)
                leafnodes = np.append(leafnodes, children)
                # update availnodes
                availnodes = np.delete(availnodes, children_ids)
        #                 settings.m(0,availnodes)
        #                 settings.m(0,leafnodes)
        #                 settings.m(0,self.Adj)
        #                 settings.m(0,'-')
        else:
            self.Adj = np.zeros((self.dim, self.dim))
            for i in range(self.dim):
                indep = self.params["rng"].binomial(1, self.p_indep)
                if indep == 0:
                    # this number includes parents (other variables)
                    # and the variable itself, therefore its
                    # self.maxnpar+2 in the following line
                    nr = self.params["rng"].integers(1, self.maxnpar + 2)
                    j_par = self.params["rng"].choice(
                        np.arange(0, self.dim), size=nr, replace=False
                    )
                    self.Adj[i, j_par] = 1
                else:
                    self.Adj[i, i] = 1
        self.Adj = np.abs(np.array(self.Adj_signed))
        # settings.m(0,self.Adj)

    def set_coupl_old(self):
        """Sample a coupling matrix using the adjacency matrix."""
        if self.model in {"krumsiek11", "var"}:
            # we already built the coupling matrix in set_coupl20()
            return
        self.Coupl = np.zeros((self.dim, self.dim))
        for i in range(self.Adj.shape[0]):
            for j, a in enumerate(self.Adj[i]):
                # if there is a 1 in Adj, specify co and antiregulation
                # and strength of regulation
                if a != 0:
                    co_anti = self.params["rng"].integers(2)
                    # set a lower bound for the coupling parameters
                    # they ought not to be smaller than 0.1
                    # and not be larger than 0.4
                    self.Coupl[i, j] = 0.0 * self.params["rng"].random() + 0.1
                    # set sign for coupling
                    if co_anti == 1:
                        self.Coupl[i, j] *= -1
        # enforce certain requirements on models
        if self.model == 1:
            self.coupl_model1()
        elif self.model == 5:
            self.coupl_model5()
        elif self.model in [6, 7]:
            self.coupl_model6()
        elif self.model in [8, 9, 10]:
            self.coupl_model8()
        # output
        if self.verbosity > 1:
            settings.m(0, self.Coupl)

    def coupl_model1(self):
        """Enforce the following signs on the couplings.

        (Model 2 has the same couplings but arbitrary signs.)
        """
        self.Coupl[0, 0] = np.abs(self.Coupl[0, 0])
        self.Coupl[0, 1] = -np.abs(self.Coupl[0, 1])
        self.Coupl[1, 1] = np.abs(self.Coupl[1, 1])

    def coupl_model5(self):
        """Toggle switch."""
        self.Coupl = -0.2 * self.Adj
        self.Coupl[2, 0] *= -1
        self.Coupl[3, 0] *= -1
        self.Coupl[4, 1] *= -1
        self.Coupl[5, 1] *= -1

    def coupl_model6(self):
        """Variant of toggle switch."""
        self.Coupl = 0.5 * self.Adj_signed

    def coupl_model8(self):
        """Variant of toggle switch."""
        self.Coupl = 0.5 * self.Adj_signed
        # reduce the value of the coupling of the repressing genes
        # otherwise completely unstable solutions are obtained
        for x in np.nditer(self.Coupl, op_flags=["readwrite"]):
            if x < -1e-6:
                x[...] = -0.2

    def coupl_model_krumsiek11(self):
        """Variant of toggle switch."""
        self.Coupl = self.Adj_signed

    def sim_model_back_help(self, Xt, Xt1):
        """Yield zero when solved for X_t given X_{t+1}."""
        return -Xt1 + Xt + self.Xdiff(Xt)

    def sim_model_backwards(self, tmax, X0):
        """Simulate the model backwards in time."""
        X = np.zeros((tmax, self.dim))
        X[tmax - 1] = X0
        for t in range(tmax - 2, -1, -1):
            sol = sp.optimize.root(
                self.sim_model_back_help, X[t + 1], args=(X[t + 1]), method="hybr"
            )
            X[t] = sol.x
        return X

    def branch_init_model1(self, tmax=100):
        # check whether we can define trajectories
        Xfix = np.array([self.Coupl[0, 1] / self.Coupl[0, 0], 1])
        if Xfix[0] > 0.97 or Xfix[0] < 0.03:
            settings.m(
                0,
                "... either no fixed point in [0,1]^2! \n"
                "    or fixed point is too close to bounds",
            )
            return None
        XbackUp = self.sim_model_backwards(
            tmax=tmax / 3, X0=Xfix + np.array([0.02, -0.02])
        )
        XbackDo = self.sim_model_backwards(
            tmax=tmax / 3, X0=Xfix + np.array([-0.02, -0.02])
        )
        Xup = self.sim_model(tmax=tmax, X0=XbackUp[0])
        Xdo = self.sim_model(tmax=tmax, X0=XbackDo[0])
        # compute mean
        X0mean = 0.5 * (Xup[0] + Xdo[0])
        if np.min(X0mean) < 0.025 or np.max(X0mean) > 0.975:
            settings.m(0, "... initial point is too close to bounds")
            return None
        if self.show and self.verbosity > 1:
            pl.figure()  # noqa: F821  TODO Fix me
            pl.plot(XbackUp[:, 0], ".b", XbackUp[:, 1], ".g")  # noqa: F821  TODO Fix me
            pl.plot(XbackDo[:, 0], ".b", XbackDo[:, 1], ".g")  # noqa: F821  TODO Fix me
            pl.plot(Xup[:, 0], "b", Xup[:, 1], "g")  # noqa: F821  TODO Fix me
            pl.plot(Xdo[:, 0], "b", Xdo[:, 1], "g")  # noqa: F821  TODO Fix me
        return X0mean

    def parents_from_boolRule(self, rule):
        """Determine parents based on boolean updaterule.

        Returns list of parents.
        """
        rule_pa = (
            rule
            .replace("(", "")
            .replace(")", "")
            .replace("or", "")
            .replace("and", "")
            .replace("not", "")
        )
        rule_pa = rule_pa.split()
        # if there are no parents, continue
        if not rule_pa:
            return []
        # check whether these are meaningful parents
        pa_old = []
        pa_delete = []
        for pa in rule_pa:
            if pa not in self.varNames:
                settings.m(0, "list of available variables:")
                settings.m(0, list(self.varNames.keys()))
                message = (
                    f"processing of rule {rule!r} yields an invalid parent: {pa} "
                    "| check whether the syntax is correct: \n"
                    'only python expressions "(",")","or","and","not" '
                    "are allowed, variable names and expressions have to be separated "
                    "by white spaces"
                )
                raise ValueError(message)
            if pa in pa_old:
                pa_delete.append(pa)
        for pa in pa_delete:
            rule_pa.remove(pa)
        return rule_pa

    def build_boolCoeff(self):
        """Compute coefficients for tuple space."""
        # coefficients for hill functions from boolean update rules
        self.boolCoeff = {s: [] for s in self.varNames}
        # parents
        self.pas = {s: [] for s in self.varNames}
        for key, rule in self.boolRules.items():
            self.pas[key] = self.parents_from_boolRule(rule)
            pasIndices = [self.varNames[pa] for pa in self.pas[key]]
            # check whether there are coupling matrix entries for each parent
            for g in range(self.dim):
                if g in pasIndices:
                    if np.abs(self.Coupl[self.varNames[key], g]) < 1e-10:
                        msg = f"specify coupling value for {key} <- {g}"
                        raise ValueError(msg)
                elif np.abs(self.Coupl[self.varNames[key], g]) > 1e-10:
                    msg = f"there should be no coupling value for {key} <- {g}"
                    raise ValueError(msg)
            if self.verbosity > 1:
                settings.m(0, "..." + key)
                settings.m(0, rule)
                settings.m(0, rule_pa)  # noqa: F821
            # now evaluate coefficients
            for tuple in list(
                itertools.product([False, True], repeat=len(self.pas[key]))
            ):
                if self.process_rule(rule, self.pas[key], tuple):
                    self.boolCoeff[key].append(tuple)
            if self.verbosity > 1:
                settings.m(0, self.boolCoeff[key])

    def process_rule(self, rule, pa, tuple):
        """Process a string that denotes a boolean rule."""
        for i, v in enumerate(tuple):
            rule = rule.replace(pa[i], str(v))
        return eval(rule)

    def write_data(
        self,
        X,
        *,
        dir=Path("sim/test"),
        noiseObs=0.0,
        append=False,
        branching=False,
        nrRealizations=1,
        seed=0,
    ):
        header = self.header
        tmax = int(X.shape[0])
        header += f"tmax = {tmax}\n"
        header += f"branching = {branching}\n"
        header += f"nrRealizations = {nrRealizations}\n"
        header += f"noiseObs = {noiseObs}\n"
        header += f"noiseDyn = {self.noiseDyn}\n"
        header += f"seed = {seed}\n"
        # add observational noise
        X += noiseObs * self.params["rng"].standard_normal((tmax, self.dim))
        # call helper function
        write_data(
            X,
            dir,
            append=append,
            header=header,
            varNames=self.varNames,
            Adj=self.Adj,
            Coupl=self.Coupl,
            model=self.model,
            modelType=self.modelType,
            boolRules=self.boolRules,
            invTimeStep=self.invTimeStep,
        )


def _check_branching(
    X: np.ndarray, Xsamples: np.ndarray, restart: int, threshold: float = 0.25
) -> tuple[bool, list[np.ndarray]]:
    """Check whether time series branches.

    Parameters
    ----------
    X
        current time series data.
    Xsamples
        list of previous branching samples.
    restart
        counts number of restart trials.
    threshold
        sets threshold for attractor identification.

    Returns
    -------
    check
        true if branching realization
    Xsamples
        updated list

    """
    check = True
    Xsamples = list(Xsamples)
    if restart == 0:
        Xsamples.append(X)
    else:
        for Xcompare in Xsamples:
            Xtmax_diff = np.absolute(X[-1, :] - Xcompare[-1, :])
            # If the second largest element is smaller than threshold
            # set check to False, i.e. at least two elements
            # need to change in order to have a branching.
            # If we observe all parameters of the system,
            # a new attractor state must involve changes in two
            # variables.
            if np.partition(Xtmax_diff, -2)[-2] < threshold:
                check = False
        if check:
            Xsamples.append(X)
    logg.debug(f"realization {restart}: {'' if check else 'no'} new branch")
    return check, Xsamples
