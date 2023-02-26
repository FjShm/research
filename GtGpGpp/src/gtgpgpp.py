import os
import argparse
import json
import yaml
import datetime
import numpy as np
import scipy.constants
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


class Stress:
    t: np.array
    pxy: np.array
    pxz: np.array
    pyz: np.array
    pxxyy: np.array
    pxxzz: np.array
    pyyzz: np.array

    def __init__(self) -> None:
        pass

    def set(self, data: np.ndarray, header: dict, unit_conversion: dict) -> None:
        ct = unit_conversion["time"]
        cp = unit_conversion["pressure"] * unit_conversion["pressure"]
        self.t = data[:, header["t"]] * ct
        self.pxy = data[:, header["pxy"]] * cp
        self.pxz = data[:, header["pxz"]] * cp
        self.pyz = data[:, header["pyz"]] * cp
        self.pxxyy = data[:, header["pxxyy"]] * cp
        self.pxxzz = data[:, header["pxxzz"]] * cp
        self.pyyzz = data[:, header["pyyzz"]] * cp
        return self


def calc_gt(stress: Stress, V: float, T: float, r: float = 0.5) -> np.ndarray:
    kB = scipy.constants.Boltzmann
    PI = stress.pxy + stress.pxz + stress.pyz
    N = stress.pxxyy + stress.pxxzz + stress.pyyzz
    gt = V / (kB * T) * ((1.0 - r) * PI / 3.0 + r * N / 12.0)
    t = stress.t
    return np.array([t, gt]).T


def calc_gpgpp(gt: np.ndarray, inputs: dict) -> np.ndarray:
    def load(gt: np.ndarray, udr: list[float], drop_minus: bool = False) -> tuple[np.array]:
        rng1 = udr[0] <= gt[:, 0]
        rng2 = gt[:, 0] <= udr[1]
        rng = rng1 * rng2
        data = gt[rng].copy()
        if drop_minus:
            rng = data[:, 1] > 0
            data = data[rng].copy()
        t = data[:, 0]
        G = data[:, 1]
        return t, G

    def windowing(G: np.array, normlized_t: np.array) -> np.array:
        return G * 0.5 * (1.0 + np.cos(normlized_t))  # hann window

    def get_freq_range(t: np.array, ndiv: int) -> tuple[float]:
        tmax, tmin = t[-1] * 1.0, t[1] - t[0]
        wmax, wmin = 1.0 * 2.0 * np.pi / tmin, 0.1 * 2.0 * np.pi / tmax

        lmax, lmin = np.ceil(np.log10(wmax)), np.floor(np.log10(wmin))
        wmax, wmin = 10**lmax, 10**lmin
        factor = 10 ** (1 / ndiv)
        return wmax, wmin, factor

    def calc_complex_modulus(G: np.array, t: np.array, w: float) -> tuple[float]:
        ta = tb = Ga = Gb = dt = 0.0

        def expi(x):
            return np.exp(complex(0.0, x))

        res = 0.0 + 0.0j
        Nd = len(t)
        for i in range(Nd):
            if i == 0:
                ta = t[i]
                Ga = G[i]
            else:
                tb = t[i]
                Gb = G[i]
                dt = tb - ta
                if w == 0:
                    res += (Ga + Gb) * 0.5 * dt
                elif w * dt < 1e-7:
                    # avoid cancellation of significant digits
                    d = complex(0.5, w * dt / 6.0) * dt
                    dc = d.conjugate()
                    res += (Ga * d + Gb * dc * complex(1.0, -w * dt)) * expi(-w * ta)
                else:
                    exwdt = expi(-w * dt)
                    d = (complex(1.0, -w * dt) - exwdt) / (w**2 * dt)
                    dc = d.conjugate()
                    res += (Ga * d + Gb * dc * exwdt) * expi(-w * ta)
                ta = tb
                Ga = Gb
        return w, -w * res.imag, +w * res.real

    ndiv = 10
    drop_minus = inputs["drop_minusGt"]
    t, G = load(gt, inputs["using_data_range"], drop_minus=drop_minus)
    G = windowing(G, (np.pi / t[-1]) * t)

    wmax, wmin, factor = get_freq_range(t, ndiv)

    omega = []
    Gp = []
    Gpp = []
    w = wmin
    while w <= wmax * 1.00001:
        r1, r2, r3 = calc_complex_modulus(G, t, w)
        omega.append(r1)
        Gp.append(r2)
        Gpp.append(r3)
        w *= factor

    gpgpp = np.array([np.array(i) for i in [omega, Gp, Gpp]]).T
    return gpgpp


def linfunc_gpp(x: np.array, x0: float, y0: float) -> np.array:
    return y0 * x / x0


def linfunc_gp(x: np.array, x0: float, y0: float) -> np.array:
    return y0 * (x / x0) ** 2


def fitting_gpgpp(gpgpp: np.ndarray, ftrng: list[float]) -> tuple:
    rng1 = ftrng[0] <= gpgpp[:, 0]
    rng2 = gpgpp[:, 0] <= ftrng[1]
    rng = rng1 * rng2
    data = gpgpp[rng].copy()

    popt_gp, pcov = curve_fit(
        linfunc_gp, data[:, 0], data[:, 1], bounds=[[ftrng[0], 0.0], [ftrng[1], np.inf]]
    )
    popt_gpp, pcov = curve_fit(linfunc_gpp, data[:, 0], data[:, 2])
    return popt_gp, popt_gpp


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="G(t), G', G'' Calculater.")
    parser.add_argument(
        "-in", "--input", help="Input file name in YAML format", default="in.yaml"
    )
    args = parser.parse_args()

    inyaml = open(args.input)
    inputs = yaml.safe_load(inyaml)
    inyaml.close()

    version = ""
    if inputs["file_version"]:
        version = "_" + datetime.datetime.now().strftime("%y%m%d-%H%M%S")
    out = open(f"log{version}.txt", mode="w")
    out.write("=== Input ===\n")
    out.write(f"{json.dumps(inputs, indent=4, sort_keys=True)}\n")

    # G(t)
    data = np.loadtxt(inputs["corr_txt_path"], comments="#")
    stress = Stress().set(data, inputs["header"], inputs["unit_conversion"])
    gt = calc_gt(stress, inputs["V"], inputs["T"], r=inputs["r"])
    np.savetxt(f"Gt{version}.txt", gt)

    # G'(w), G''(w)
    gpgpp = calc_gpgpp(gt, inputs)
    np.savetxt(f"GpGpp{version}.txt", gpgpp)

    # fitting
    popt_gp, popt_gpp = fitting_gpgpp(gpgpp, inputs["gpgpp_fitting_range"])
    cross_x = popt_gpp[1] / popt_gp[1] * (popt_gp[0] ** 2) / popt_gpp[0]
    x = np.linspace(inputs["gpgpp_fitting_range"][0] * 0.1, cross_x * 10.0)
    y_fitgp = linfunc_gp(x, *popt_gp)
    y_fitgpp = linfunc_gpp(x, *popt_gpp)
    out.write("\n=== maximum relaxation time ===\n")
    out.write(f"1/tau1: {cross_x:<.4f} [ns^-1]\n")
    out.write(f"tau1: {1./cross_x:<.4f} [ns]\n")

    # plot
    fig = plt.figure()
    ax = fig.add_subplot(
        111,
        xlabel="time [ns]",
        ylabel=r"$G(t)$ [MPa]",
        xlim=(1e-3, 1e2),
        xticks=np.logspace(-3, 2, num=6),
        xscale="log",
        ylim=(1e-3, 1e3),
        yticks=np.logspace(-3, 3, num=7),
        yscale="log",
    )
    ax.plot(gt[:, 0], gt[:, 1] * 1e-6, color="r", marker="o", mfc="none")
    fig.savefig(f"Gt{version}.eps", bbox_inches="tight")
    fig.savefig(f"Gt{version}.png", bbox_inches="tight")
    del fig, ax

    fig = plt.figure()
    ax = fig.add_subplot(
        111,
        xlabel=r"$\omega\ {\rm [ns^{-1}]}$",
        ylabel=r"$G',\ G''\ {\rm [Pa]}$",
        xscale="log",
        yscale="log",
    )
    ax.plot(gpgpp[:, 0], gpgpp[:, 1], color="r", marker="o", mfc="none", label=r"$G'$")
    ax.plot(gpgpp[:, 0], gpgpp[:, 2], color="b", marker="o", mfc="none", label=r"$G''$")
    ax.plot(x, y_fitgp, "r--", lw=1.5)
    ax.plot(x, y_fitgpp, "b--", lw=1.5)
    ax.axvline(cross_x, color="k", lw=1.0)
    ax.legend()
    fig.savefig(f"GpGpp{version}.eps", bbox_inches="tight")
    fig.savefig(f"GpGpp{version}.png", bbox_inches="tight")

    out.close()
