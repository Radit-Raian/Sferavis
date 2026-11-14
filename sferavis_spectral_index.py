import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional


def _clip01(x: float) -> float:
    return float(np.clip(x, 0.0, 1.0))

def _overlap_length(a0: float, a1: float, b0: float, b1: float) -> float:
    return max(0.0, min(a1, b1) - max(a0, b0))

def _safe_percentile(a: np.ndarray, q: float, default: float = 0.0) -> float:
    a = np.asarray(a)
    if a.size == 0:
        return default
    return float(np.percentile(a, q))

def _local_continuum_level(w: np.ndarray, f: np.ndarray, center: float, width: float = 0.03) -> float:
    
    """Estimate a local continuum level as the median flux in a window around center (excluding the center bin if possible).
    width is in microns."""
    
    mask = (w >= center - width) & (w <= center + width)
    if not np.any(mask):
        # fallback: global median
        return float(np.median(f))
    return float(np.median(f[mask]))

def _absorption_strength(w: np.ndarray, f: np.ndarray, wmin: float, wmax: float, cont_center: float) -> float:
    
    """For absorption features: strength = (cont - feature)/|cont|, clamped to [0, 1] via normalization.
    Works even if flux is not continuum-normalized, as long as cont is representative."""
    
    mask_feat = (w >= wmin) & (w <= wmax)
    if not np.any(mask_feat):
        return 0.0

    cont = _local_continuum_level(w, f, center=cont_center, width=(wmax - wmin) * 0.75 + 0.02)
    feat = float(np.median(f[mask_feat]))

    # If cont is ~0, avoid division issues
    if np.abs(cont) < 1e-12:
        return 0.0

    strength_raw = (cont - feat) / abs(cont)
    # map to [0,1] softly
    return _clip01(0.5 + strength_raw)  # heuristic mapping

def _emission_strength(w: np.ndarray, f: np.ndarray, wmin: float, wmax: float, cont_center: float) -> float:
    
    """For emission features: strength = (feature - cont)/|cont|, clamped to [0, 1]."""
    
    mask_feat = (w >= wmin) & (w <= wmax)
    if not np.any(mask_feat):
        return 0.0

    cont = _local_continuum_level(w, f, center=cont_center, width=(wmax - wmin) * 0.75 + 0.02)
    feat = float(np.median(f[mask_feat]))

    if np.abs(cont) < 1e-12:
        return 0.0

    strength_raw = (feat - cont) / abs(cont)
    return _clip01(0.5 + strength_raw)  # heuristic mapping


@dataclass(frozen=True)
class ClassSpec:
    name: str
    range_min: float
    range_max: float
    # (feature_name, feature_weight, feature_type, wmin, wmax, cont_center)
    # feature_type: 'abs' or 'em'
    features: Tuple[Tuple[str, float, str, float, float, float], ...]
    # Optional: bands that are low-transmission / noisy (used for supernova penalty)
    penalty_windows: Tuple[Tuple[float, float, float], ...]  # (wmin, wmax, penalty_weight)


def _make_class_specs() -> Dict[str, ClassSpec]:
    # Brown dwarfs (from documentation_2.docx)
    # - FeH at 1.19 μm
    # - CO at 2.295 μm
    # - CH4 initial signatures in K band at ~2.20 μm
    # - CH4 appearance in H band at 1.667 μm
    # I implement small windows around each center (configurable via feature windows).
    brown = ClassSpec(
        name="Brown dwarfs",
        range_min=1.13, range_max=2.33,
        features=(
            ("FeH_1p19", 1.0, "abs", 1.175, 1.205, 1.19),
            ("CO_2p295", 0.9, "abs", 2.270, 2.320, 2.295),
            ("CH4_K_2p20", 1.0, "abs", 2.175, 2.225, 2.20),
            ("CH4_H_1p667", 1.0, "abs", 1.63, 1.70, 1.667),
        ),
        penalty_windows=(),
    )

    # Quasars:
    # Documentation provides spectral ranges/resolutions, but no explicit line centers.
    # I use template-less heuristics: presence of broad IR coverage + generic "non-stellar"
    # via slope/curvature features would normally be used. Here we keep it range-based
    # plus a weak template-less score using "variance vs continuum".
    # To keep this code faithful to your doc, we implement quasars as:
    # - range gating (dominant)
    # - low penalty for "no strong narrow absorption at BD key bands"
    # I'll approximate that via penalizing absorption at BD centers.
    
    quasars = ClassSpec(
        name="Quasars",
        range_min=0.27, range_max=3.70,  # composite range in doc (final optical–IR composite)
        features=(
            # No explicit line/feature centers in doc; set empty.
        ),
        penalty_windows=(),
    )

    # Supernovas (explicit ion/element wavelength windows from doc)
    supernova = ClassSpec(
        name="Supernovas",
        range_min=0.8, range_max=2.5,
        features=(
            ("CaII_0p85_0p866", 0.8, "abs", 0.850, 0.866, 0.86),
            ("MgII_0p922_0p924", 0.8, "abs", 0.922, 0.924, 0.923),
            ("MgII_CoII_blend_1p091_1p095", 0.6, "abs", 1.091, 1.095, 1.093),
            ("SiII_MgII_blend_1p676_1p691", 0.6, "abs", 1.676, 1.691, 1.684),
            ("FeGroup_1p50_1p90", 0.6, "abs", 1.50, 1.90, 1.70),
        ),
        # Noise/low transmission windows (penalize features there)
        penalty_windows=(
            (1.31, 1.38, 0.6),
            (1.78, 1.88, 0.6),
            (0.93, 0.96, 0.4),
            (1.12, 1.16, 0.4),
            (2.40, 2.50, 0.4),
        )
    )

    # Asteroids: doc gives range only; no explicit bands.
    asteroids = ClassSpec(
        name="Asteroids",
        range_min=0.3, range_max=4.0,
        features=(
            # I'll treat absence of BD and SN features as "evidence" via penalties later.
        ),
        penalty_windows=(),
    )

    # Galaxies: explicit emission lines + range
    galaxies = ClassSpec(
        name="Galaxies",
        range_min=0.9, range_max=2.5,
        features=(
            ("H2_emission_2p122", 1.0, "em", 2.10, 2.14, 2.122),
            ("Brgamma_emission_2p166", 0.9, "em", 2.145, 2.185, 2.166),
        ),
        penalty_windows=(),
    )

    # Comets: emission bands at ~3.4, ~10, ~20 μm; range 3–160 μm
    comets = ClassSpec(
        name="Comets",
        range_min=3.0, range_max=160.0,
        features=(
            ("C-H_emission_3p4", 1.0, "em", 3.25, 3.55, 3.4),
            ("Si-O_broad_emission_10", 0.9, "em", 9.0, 11.5, 10.0),
            ("Si-O-Si_emission_20", 0.8, "em", 18.0, 22.5, 20.0),
        ),
        penalty_windows=(),
    )

    # Quasars range-only: store with empty features
    quasars = ClassSpec(
        name="Quasars",
        range_min=0.27, range_max=3.70,
        features=tuple(),
        penalty_windows=(),
    )

    return {
        "Brown dwarfs": brown,
        "Quasars": quasars,
        "Supernovas": supernova,
        "Asteroids": asteroids,
        "Galaxies": galaxies,
        "Comets": comets,
    }


CLASS_SPECS = _make_class_specs()


# ---------------------------
# Scoring logic
# ---------------------------

def _range_score(w: np.ndarray, c: ClassSpec) -> float:
    wmin_obs = float(np.min(w))
    wmax_obs = float(np.max(w))
    overlap = _overlap_length(wmin_obs, wmax_obs, c.range_min, c.range_max)
    denom = max(1e-12, float(wmax_obs - wmin_obs))
    # if spectrum is very narrow, overlap fraction can be small; keep it bounded
    return _clip01(overlap / denom)

def _apply_supernova_penalty(w: np.ndarray, f: np.ndarray, base: float, c: ClassSpec) -> float:
    if not c.penalty_windows:
        return base
    # penalize if observation includes penalty windows with strong "signal"
    # I'll compute average absolute deviation in those windows.
    penalty = 0.0
    for (pw0, pw1, pw) in c.penalty_windows:
        mask = (w >= pw0) & (w <= pw1)
        if np.any(mask):
            # estimate "badness" as strong variation/abs/emission
            # Here I use percentile of abs(flux) variation proxy.
            bad = float(np.std(f[mask])) / (np.std(f) + 1e-12)
            penalty += pw * bad
    return base * _clip01(1.0 - penalty)

def extract_features(wavelength_um: np.ndarray, flux: np.ndarray) -> Dict[str, float]:
    
    """Returns a dict of extracted feature strengths for all classes.
    For abs features: absorption_strength in [0,1]
    For em features: emission_strength in [0,1]"""
    
    w = np.asarray(wavelength_um, dtype=float)
    f = np.asarray(flux, dtype=float)

    feats = {}


    m = np.isfinite(w) & np.isfinite(f)
    w = w[m]
    f = f[m]

    for cls in CLASS_SPECS.values():
        for (fname, _, ftype, wmin, wmax, cont_center) in cls.features:
            if ftype == "abs":
                feats[fname] = _absorption_strength(w, f, wmin, wmax, cont_center)
            elif ftype == "em":
                feats[fname] = _emission_strength(w, f, wmin, wmax, cont_center)
            else:
                raise ValueError(f"Unknown feature type: {ftype}")
    return feats


def classify_spectrum(wavelength_um: np.ndarray, flux: np.ndarray) -> Dict:
    
    """Returns:
      - predicted_class
      - scores per class
      - extracted feature strengths (for debugging)"""
      
    w = np.asarray(wavelength_um, dtype=float)
    f = np.asarray(flux, dtype=float)

    # Basic sanitation
    m = np.isfinite(w) & np.isfinite(f)
    w = w[m]
    f = f[m]
    if w.size < 5:
        raise ValueError("Spectrum too small after removing NaNs.")

    feats = extract_features(w, f)

    scores = {}
    for cname, c in CLASS_SPECS.items():
        s = 0.0

        # Strong range prior
        rscore = _range_score(w, c)
        s += 1.5 * rscore  # range dominates

        # Add class-specific features
        if c.features:
            for (fname, wgt, ftype, wmin, wmax, cont_center) in c.features:
                s += wgt * feats.get(fname, 0.0)

        # Extra logic for classes that don't list explicit features in the doc
        # - Quasars: penalize if strong BD methane/FeH/CO absorption exists
        # - Asteroids: penalize if strong SN line windows exist, or BD methane exists
        if cname == "Quasars":
            bd_pen = (
                feats.get("FeH_1p19", 0.0) +
                feats.get("CO_2p295", 0.0) +
                feats.get("CH4_K_2p20", 0.0) +
                feats.get("CH4_H_1p667", 0.0)
            ) / 4.0
            s -= 1.0 * bd_pen

        if cname == "Asteroids":
            bd_pen = (
                feats.get("FeH_1p19", 0.0) +
                feats.get("CO_2p295", 0.0) +
                feats.get("CH4_K_2p20", 0.0) +
                feats.get("CH4_H_1p667", 0.0)
            ) / 4.0
            sn_pen = (
                feats.get("CaII_0p85_0p866", 0.0) +
                feats.get("MgII_0p922_0p924", 0.0) +
                feats.get("FeGroup_1p50_1p90", 0.0)
            ) / 3.0
            s -= 1.0 * (0.6 * bd_pen + 0.4 * sn_pen)

        if cname == "Supernovas":
            s = _apply_supernova_penalty(w, f, s, c)

        scores[cname] = float(s)

    predicted = max(scores, key=scores.get)

    return {
        "predicted_class": predicted,
        "scores": scores,
        "extracted_features": feats,
        "observed_range_um": (float(np.min(w)), float(np.max(w))),
    }


"""
# Example (to be removed)

if __name__ == "__main__":
    # Dummy test: create a synthetic BD-like absorption spectrum on 1.0–2.5 μm
    w = np.linspace(1.0, 2.5, 2000)
    f = np.ones_like(w)

    # add absorption dips
    for center, amp, width in [(1.19, 0.4, 0.015), (2.295, 0.35, 0.03), (2.20, 0.5, 0.02), (1.667, 0.45, 0.02)]:
        f -= amp * np.exp(-0.5 * ((w - center) / width) ** 2)

    res = classify_spectrum(w, f)
    print(res["predicted_class"])
    print(res["scores"])
"""
