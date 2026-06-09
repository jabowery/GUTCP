"""GUTCP cosmology audit helpers for the 2023 reconstructed Chapter 32.

The corrected GUTCP redshift counterpart implemented here is based on the
expansion extent implied by Mills's arrow-of-time section:

    E(t) = R(t) - R_min
    1 + z_GUTCP(t_emit, t_obs) = E(t_obs) / E(t_emit).

This differs from Mills's printed H(t) in Eq. (32.156), which divides the
oscillatory expansion rate by c t and therefore creates a damped sinc-like
quantity. That printed version is retained as mills_ct_path for audit
comparisons. Mills Eq. (32.165),

    lambda_infinity = lambda(r) * (1 + 2 G M / (r c^2)),

is kept as a separate endpoint gravitational-clock correction. The inverse map
from observed z to emission time is generally multi-valued over the oscillatory
history, so distance-redshift work must select a branch.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import cos, exp, isfinite, log, pi, sin, sqrt
from typing import Callable, Iterable, Literal


C = 299_792_458.0
G = 6.67430e-11
SIGMA_SB = 5.670374419e-8

SECONDS_PER_YEAR = 31_536_000.0
SECONDS_PER_BILLION_YEARS = 1.0e9 * SECONDS_PER_YEAR
LIGHT_YEAR_M = 9_460_730_472_580_800.0
MPC_KM = 3.0856775814913673e19

NOTEBOOK_CMBR_START_K = 4.4959973429079625
NOTEBOOK_T_NOW_S = 3.767087121815851e18
BOOK_M_START_KG = 2.0e54
OBSERVED_CMB_NOW_K = 2.72548
MILLS_PAGE_1579_R_SPHERE_LIGHT_YEARS = 14.02e9

PeriodConvention = Literal["numeric", "printed_32_149"]
RadiusNormalization = Literal["notebook", "page_1544"]
MassModel = Literal["constant", "book", "integrated"]
RedshiftModel = Literal[
    "expansion_extent",
    "radius_scale",
    "mills_ct_path",
    "hubble_path",
    "mills_endpoint",
    "combined_extent_endpoint",
    "combined_radius_endpoint",
    "combined",
]
LineOfSightHModel = Literal["expansion_extent", "radius_scale", "mills_ct_path", "hubble_path"]


def years_to_seconds(years: float) -> float:
    return years * SECONDS_PER_YEAR


def billion_years_to_seconds(billion_years: float) -> float:
    return billion_years * SECONDS_PER_BILLION_YEARS


def seconds_to_billion_years(seconds: float) -> float:
    return seconds / SECONDS_PER_BILLION_YEARS


def meters_to_light_years(meters: float) -> float:
    return meters / LIGHT_YEAR_M


def radial_light_travel_distance_m(t_emit_s: float, t_obs_s: float) -> float:
    """Radial light-front distance relation D = c (t_obs - t_emit)."""

    return C * (t_obs_s - t_emit_s)


def blue_shift_segment_entry_threshold_m(t_obs_s: float) -> float:
    """Distance threshold D = c t_obs for leaving the current expansion branch."""

    return C * t_obs_s


def symmetric_endpoint_neutral_threshold_m(t_obs_s: float) -> float:
    """Endpoint-only neutral threshold D = 2 c t_obs for even oscillatory branches."""

    return 2.0 * C * t_obs_s


def blueshift_threshold_status(distance_m: float, t_obs_s: float) -> str:
    """Classify a distance against the segment-entry and net endpoint thresholds."""

    entry = blue_shift_segment_entry_threshold_m(t_obs_s)
    neutral = symmetric_endpoint_neutral_threshold_m(t_obs_s)
    if distance_m < entry:
        return "current-expansion branch only"
    if distance_m < neutral:
        return "pre-current-expansion segment enters; endpoint ratio can remain red"
    if distance_m == neutral:
        return "symmetric endpoint-neutral threshold"
    return "symmetric endpoint-only net blueshift"


def blueshift_branch_constraints(
    distance_m: float,
    t_obs_s: float,
    observed_z: float | None = None,
) -> tuple[str, ...]:
    """Return branch-admissibility notes implied by distance and optional z sign."""

    entry = blue_shift_segment_entry_threshold_m(t_obs_s)
    neutral = symmetric_endpoint_neutral_threshold_m(t_obs_s)
    constraints = []
    if distance_m < entry:
        constraints.append("current-expansion branch admissible; pre-current segment not forced")
    elif distance_m == entry:
        constraints.append("on current-expansion boundary t_emit == 0")
    else:
        constraints.append("current-expansion-only branch excluded; pre-current segment required")

    if observed_z is not None:
        if observed_z >= 0.0 and distance_m > neutral:
            constraints.append("red/neutral observation excludes endpoint-symmetric net-blue branch")
        elif observed_z < 0.0 and distance_m <= neutral:
            constraints.append("blue observation excludes endpoint-symmetric red/neutral branch")
        else:
            constraints.append("observed z sign is compatible with endpoint-symmetric threshold")

    return tuple(constraints)


def cmbr_angular_view_light_years(
    r_sphere_light_years: float = MILLS_PAGE_1579_R_SPHERE_LIGHT_YEARS,
) -> float:
    """Mills page-1579 angular view A_view = 2*pi*r_sphere."""

    return 2.0 * pi * r_sphere_light_years


def cmbr_multipole_sky_fraction(ell: float) -> float:
    """Mills page-1579 sky fraction spanned by CMBR multipole ell."""

    if ell <= 0:
        raise ValueError("multipole ell must be positive")
    return 2.0 / ell


def cmbr_multipole_structure_scale_light_years(
    ell: float,
    r_sphere_light_years: float = MILLS_PAGE_1579_R_SPHERE_LIGHT_YEARS,
) -> float:
    """Mills page-1579 large-structure scale L_ell = 4*pi*r_sphere/ell."""

    return cmbr_angular_view_light_years(r_sphere_light_years) * cmbr_multipole_sky_fraction(ell)


def radius_min_from_temperature(
    temperature_k: float,
    *,
    normalization: RadiusNormalization = "notebook",
    emissivity: float = 1.0,
) -> float:
    """Return R_min from the Stefan-Boltzmann chain.

    normalization="notebook" uses the algebra present in
    GUTCPCosmologicalModel.nb:

        sqrt(c^5 / ((4*pi)^2 G e sigma T^4))

    normalization="page_1544" uses the denominator implied by reconstructed
    page 1544 Eq. (32.145):

        sqrt(c^5 / (4*pi^2 G e sigma T^4))
    """

    if normalization == "notebook":
        angular_factor = (4.0 * pi) ** 2
    elif normalization == "page_1544":
        angular_factor = 4.0 * pi**2
    else:
        raise ValueError(f"unknown radius normalization: {normalization!r}")

    return sqrt(C**5 / (angular_factor * G * emissivity * SIGMA_SB * temperature_k**4))


def mass_from_minimum_radius(radius_m: float) -> float:
    """Use Eq. (32.147), r_g = 2 G M / c^2, inverted for M."""

    return C**2 * radius_m / (2.0 * G)


def h0_km_s_mpc_to_s(h0_km_s_mpc: float) -> float:
    return h0_km_s_mpc / MPC_KM


@dataclass(frozen=True)
class CalibrationResult:
    h0_target_km_s_mpc: float
    cmb_now_target_k: float
    t_now_s: float
    cmb_start_k: float
    h0_model_km_s_mpc: float
    cmb_now_model_k: float
    iterations: int

    @property
    def t_now_billion_years(self) -> float:
        return seconds_to_billion_years(self.t_now_s)


@dataclass(frozen=True)
class FixedTimeDiagnostic:
    t_fixed_s: float
    h0_target_km_s_mpc: float
    cmb_now_target_k: float
    cmb_start_for_h0_k: float
    cmb_now_when_h0_matches_k: float
    cmb_start_for_cmb_k: float
    h0_when_cmb_matches_km_s_mpc: float

    @property
    def t_fixed_billion_years(self) -> float:
        return seconds_to_billion_years(self.t_fixed_s)


@dataclass(frozen=True)
class RedshiftBranch:
    branch_index: int
    t_emit_s: float
    t_obs_s: float
    z_model: float
    light_travel_distance_m: float

    @property
    def t_emit_billion_years(self) -> float:
        return seconds_to_billion_years(self.t_emit_s)

    @property
    def lookback_billion_years(self) -> float:
        return seconds_to_billion_years(self.t_obs_s - self.t_emit_s)

    @property
    def light_travel_distance_billion_light_years(self) -> float:
        return meters_to_light_years(self.light_travel_distance_m) / 1.0e9


@dataclass(frozen=True)
class ConsistencyReadinessRow:
    quantity: str
    paper_role: str
    gutcp_status: str
    fair_evaluation: str


@dataclass(frozen=True)
class LightBundleFrontierRow:
    source: str
    relation_fragment: str
    constraint: str
    status: str


@dataclass(frozen=True)
class ConsistencyObservableBranch:
    """Branch-specific quantities available before a GUTCP angular map exists."""

    target_z: float
    branch_index: int
    t_emit_s: float
    t_obs_s: float
    h_parallel_km_s_mpc: float
    light_travel_distance_m: float
    angular_diameter_distance_m: float | None = None
    status: str = "surrogate: d_A(z) not supplied by audited GUTCP material"

    @property
    def t_emit_billion_years(self) -> float:
        return seconds_to_billion_years(self.t_emit_s)

    @property
    def light_travel_distance_billion_light_years(self) -> float:
        return meters_to_light_years(self.light_travel_distance_m) / 1.0e9


def sine_integral(x: float, *, tolerance: float = 1.0e-16, max_terms: int = 400) -> float:
    """Return Si(x) = integral_0^x sin(u)/u du.

    The audited branch examples use |x| of order one cycle. A power series is
    accurate and dependency-free over that range; a Simpson fallback covers
    larger exploratory intervals without adding scipy/mpmath.
    """

    if x == 0.0:
        return 0.0

    sign = -1.0 if x < 0.0 else 1.0
    ax = abs(x)
    if ax <= 16.0:
        term = ax
        total = term
        for k in range(max_terms):
            ratio = -(ax * ax) * (2 * k + 1) / ((2 * k + 3) ** 2 * (2 * k + 2))
            term *= ratio
            total += term
            if abs(term) <= tolerance * max(1.0, abs(total)):
                return sign * total
        return sign * total

    n = max(4096, int(ax * 256))
    if n % 2:
        n += 1
    h = ax / n

    def f(u: float) -> float:
        if u == 0.0:
            return 1.0
        return sin(u) / u

    total = f(0.0) + f(ax)
    for i in range(1, n):
        total += (4.0 if i % 2 else 2.0) * f(i * h)
    return sign * total * h / 3.0


@dataclass(frozen=True)
class GUTCPModel:
    """Numerical GUTCP cosmology model with explicit convention switches."""

    m_start_kg: float
    period_convention: PeriodConvention = "numeric"

    @classmethod
    def from_cmbr_start(
        cls,
        temperature_k: float = NOTEBOOK_CMBR_START_K,
        *,
        normalization: RadiusNormalization = "notebook",
        period_convention: PeriodConvention = "numeric",
    ) -> "GUTCPModel":
        radius_m = radius_min_from_temperature(temperature_k, normalization=normalization)
        return cls(
            m_start_kg=mass_from_minimum_radius(radius_m),
            period_convention=period_convention,
        )

    @classmethod
    def from_book_mass(
        cls,
        m_start_kg: float = BOOK_M_START_KG,
        *,
        period_convention: PeriodConvention = "numeric",
    ) -> "GUTCPModel":
        return cls(m_start_kg=m_start_kg, period_convention=period_convention)

    @property
    def gravitational_radius_m(self) -> float:
        return 2.0 * G * self.m_start_kg / C**2

    @property
    def radius_amplitude_m(self) -> float:
        return 4.0 * pi * G * self.m_start_kg / C**2

    @property
    def average_radius_m(self) -> float:
        return self.gravitational_radius_m + self.radius_amplitude_m

    @property
    def period_s(self) -> float:
        if self.period_convention == "numeric":
            factor = 2.0
        elif self.period_convention == "printed_32_149":
            factor = 4.0
        else:
            raise ValueError(f"unknown period convention: {self.period_convention!r}")
        return factor * pi * G * self.m_start_kg / C**3

    @property
    def angular_frequency(self) -> float:
        return 2.0 * pi / self.period_s

    def theta(self, t_s: float) -> float:
        return self.angular_frequency * t_s

    def radius_m(self, t_s: float) -> float:
        """Eq. (32.153) in the notebook/numeric form."""

        return self.average_radius_m - self.radius_amplitude_m * cos(self.theta(t_s))

    def expansion_extent_m(self, t_s: float) -> float:
        """Expansion extent from the minimum-radius bounce, R(t)-R_min."""

        return self.radius_m(t_s) - self.gravitational_radius_m

    def radius_rate_m_s(self, t_s: float) -> float:
        """Eq. (32.154) generalized to the selected period convention."""

        return self.radius_amplitude_m * self.angular_frequency * sin(self.theta(t_s))

    def hubble_s(self, t_s: float) -> float:
        """Mills Eq. (32.156), H = radius_rate / (c t), in 1/s."""

        if t_s == 0.0:
            return self.radius_amplitude_m * self.angular_frequency**2 / C
        return self.radius_rate_m_s(t_s) / (C * t_s)

    def hubble_km_s_mpc(self, t_s: float) -> float:
        return self.hubble_s(t_s) * MPC_KM

    def hubble_radius_s(self, t_s: float) -> float:
        """Corrected Hubble ratio using the whole oscillatory radius."""

        return self.radius_rate_m_s(t_s) / self.radius_m(t_s)

    def hubble_radius_km_s_mpc(self, t_s: float) -> float:
        return self.hubble_radius_s(t_s) * MPC_KM

    def hubble_expansion_extent_s(self, t_s: float) -> float:
        """Corrected Hubble ratio using R(t)-R_min as expansion extent."""

        extent = self.expansion_extent_m(t_s)
        if abs(extent) <= 1.0e-15 * max(1.0, self.gravitational_radius_m):
            return float("inf")
        return self.radius_rate_m_s(t_s) / extent

    def hubble_expansion_extent_km_s_mpc(self, t_s: float) -> float:
        return self.hubble_expansion_extent_s(t_s) * MPC_KM

    def hubble_redshift_antiderivative(self, t_s: float) -> float:
        """Integral of Mills's c t denominator H(t) dt from 0 to t.

        With Eq. (32.156), H(t) = (A omega / c) sin(omega t)/t, so the
        antiderivative is (A omega / c) Si(omega t). For the notebook/numeric
        period convention, A omega / c = 4 pi.
        """

        coefficient = self.radius_amplitude_m * self.angular_frequency / C
        return coefficient * sine_integral(self.theta(t_s))

    def mass_book_kg(self, t_s: float) -> float:
        """Eq. (32.158), the book/notebook mUBook[t] matter inventory."""

        return 0.5 * self.m_start_kg * (1.0 + cos(self.theta(t_s)))

    def power_book_w(self, t_s: float) -> float:
        """Eq. (32.161), P_U(t)."""

        return C**5 * (1.0 + cos(self.theta(t_s))) / (8.0 * pi * G)

    def temperature_k(self, t_s: float, *, emissivity: float = 1.0) -> float:
        """Notebook temperature path from P_U(t)/(4 pi R(t)^2 sigma e)."""

        area = 4.0 * pi * self.radius_m(t_s) ** 2
        radiance = self.power_book_w(t_s) / area
        return (radiance / (emissivity * SIGMA_SB)) ** 0.25

    def mass_integrated_kg(self, t_s: float) -> float:
        """Notebook alternative: M0 minus integral(P_U(t)/c^2 dt).

        This is included because the notebook comments identify it as a
        correction candidate to Eq. (32.158). It is not Mills Eq. (32.158).
        """

        omega = self.angular_frequency
        radiated = C**3 * (t_s + sin(self.theta(t_s)) / omega) / (8.0 * pi * G)
        return self.m_start_kg - radiated

    def redshift_mass_kg(self, t_s: float, mass_model: MassModel) -> float:
        if mass_model == "constant":
            return self.m_start_kg
        if mass_model == "book":
            return self.mass_book_kg(t_s)
        if mass_model == "integrated":
            return self.mass_integrated_kg(t_s)
        raise ValueError(f"unknown mass model: {mass_model!r}")

    def redshift_factor_component(self, t_s: float, mass_model: MassModel = "constant") -> float:
        """B(t) = 1 + 2 G M_z(t) / (c^2 R(t)) from Eq. (32.165)."""

        radius = self.radius_m(t_s)
        mass = self.redshift_mass_kg(t_s, mass_model)
        return 1.0 + 2.0 * G * mass / (C**2 * radius)

    def one_plus_z_endpoint_candidate(
        self,
        t_emit_s: float,
        t_obs_s: float,
        *,
        mass_model: MassModel = "constant",
    ) -> float:
        """Conditional finite-r-sphere endpoint ratio from Mills Eq. (32.165).

        lambda_inf = lambda_emit B(t_emit) = lambda_obs B(t_obs), so
        1 + z = lambda_obs / lambda_emit = B(t_emit) / B(t_obs).
        The identification of this endpoint ratio with a unique cosmological
        z_GUTCP branch is an additional observational-dictionary assumption,
        not a derived result.
        """

        return self.redshift_factor_component(t_emit_s, mass_model) / self.redshift_factor_component(
            t_obs_s, mass_model
        )

    def z_endpoint_candidate(
        self,
        t_emit_s: float,
        t_obs_s: float,
        *,
        mass_model: MassModel = "constant",
    ) -> float:
        return self.one_plus_z_endpoint_candidate(t_emit_s, t_obs_s, mass_model=mass_model) - 1.0

    def one_plus_z_hubble_path(self, t_emit_s: float, t_obs_s: float) -> float:
        """Path redshift from Mills's printed c t denominator.

        For a light-front increment dD = c dt, the local Hubble redshift is
        dz ~= H(t) dD/c = H(t) dt. Multiplying adjacent increments gives
        d ln(1+z) = H(t) dt. This preserves the notebook's Eq. (32.156)
        denominator and therefore the damped sinc-like behavior under audit.
        """

        return exp(
            self.hubble_redshift_antiderivative(t_obs_s)
            - self.hubble_redshift_antiderivative(t_emit_s)
        )

    def z_hubble_path(self, t_emit_s: float, t_obs_s: float) -> float:
        return self.one_plus_z_hubble_path(t_emit_s, t_obs_s) - 1.0

    @staticmethod
    def _positive_ratio(numerator: float, denominator: float) -> float:
        scale = max(1.0, abs(numerator), abs(denominator))
        if abs(denominator) <= 1.0e-14 * scale:
            return float("inf") if numerator >= 0.0 else float("-inf")
        return numerator / denominator

    def one_plus_z_radius_scale(self, t_emit_s: float, t_obs_s: float) -> float:
        """Corrected path redshift using the whole GUTCP radius as scale."""

        return self._positive_ratio(self.radius_m(t_obs_s), self.radius_m(t_emit_s))

    def z_radius_scale(self, t_emit_s: float, t_obs_s: float) -> float:
        return self.one_plus_z_radius_scale(t_emit_s, t_obs_s) - 1.0

    def one_plus_z_expansion_extent(self, t_emit_s: float, t_obs_s: float) -> float:
        """Corrected path redshift using R(t)-R_min as expansion extent."""

        return self._positive_ratio(
            self.expansion_extent_m(t_obs_s),
            self.expansion_extent_m(t_emit_s),
        )

    def z_expansion_extent(self, t_emit_s: float, t_obs_s: float) -> float:
        return self.one_plus_z_expansion_extent(t_emit_s, t_obs_s) - 1.0

    def one_plus_z_gutcp(
        self,
        t_emit_s: float,
        t_obs_s: float,
        *,
        mass_model: MassModel = "constant",
        redshift_model: RedshiftModel = "expansion_extent",
    ) -> float:
        """GUTCP redshift counterpart for a selected path/endpoint model."""

        if redshift_model == "expansion_extent":
            return self.one_plus_z_expansion_extent(t_emit_s, t_obs_s)
        if redshift_model == "radius_scale":
            return self.one_plus_z_radius_scale(t_emit_s, t_obs_s)
        if redshift_model in ("mills_ct_path", "hubble_path"):
            return self.one_plus_z_hubble_path(t_emit_s, t_obs_s)
        if redshift_model == "mills_endpoint":
            return self.one_plus_z_endpoint_candidate(t_emit_s, t_obs_s, mass_model=mass_model)
        if redshift_model in ("combined_extent_endpoint", "combined"):
            return self.one_plus_z_expansion_extent(
                t_emit_s, t_obs_s
            ) * self.one_plus_z_endpoint_candidate(t_emit_s, t_obs_s, mass_model=mass_model)
        if redshift_model == "combined_radius_endpoint":
            return self.one_plus_z_radius_scale(t_emit_s, t_obs_s) * self.one_plus_z_endpoint_candidate(
                t_emit_s, t_obs_s, mass_model=mass_model
            )
        raise ValueError(f"unknown redshift model: {redshift_model!r}")

    def z_gutcp(
        self,
        t_emit_s: float,
        t_obs_s: float,
        *,
        mass_model: MassModel = "constant",
        redshift_model: RedshiftModel = "expansion_extent",
    ) -> float:
        return (
            self.one_plus_z_gutcp(
                t_emit_s,
                t_obs_s,
                mass_model=mass_model,
                redshift_model=redshift_model,
            )
            - 1.0
        )

    def z_radius_ratio_diagnostic(self, t_emit_s: float, t_obs_s: float) -> float:
        """Compatibility alias for the whole-radius denominator correction."""

        return self.z_radius_scale(t_emit_s, t_obs_s)

    def redshift_branch_times(
        self,
        target_z: float,
        t_obs_s: float,
        *,
        t_min_s: float | None = None,
        t_max_s: float | None = None,
        mass_model: MassModel = "constant",
        redshift_model: RedshiftModel = "expansion_extent",
        samples: int = 4096,
        tolerance: float = 1.0e-10,
    ) -> list[float]:
        """Find emission-time branches satisfying z_GUTCP(t_emit,t_obs)=target_z."""

        if target_z <= -1.0:
            raise ValueError("target_z must be greater than -1")
        if t_min_s is None:
            t_min_s = t_obs_s - self.period_s
        if t_max_s is None:
            t_max_s = t_obs_s
        if t_min_s >= t_max_s:
            raise ValueError("t_min_s must be less than t_max_s")
        if samples < 2:
            raise ValueError("samples must be at least 2")

        target = 1.0 + target_z

        def residual(t_s: float) -> float:
            return (
                self.one_plus_z_gutcp(
                    t_s,
                    t_obs_s,
                    mass_model=mass_model,
                    redshift_model=redshift_model,
                )
                - target
            )

        roots: list[float] = []

        def add_root(root: float) -> None:
            scale = max(1.0, abs(t_max_s - t_min_s))
            if not any(abs(root - old) <= tolerance * scale for old in roots):
                roots.append(root)

        step = (t_max_s - t_min_s) / samples
        prev_t = t_min_s
        prev_f = residual(prev_t)
        if abs(prev_f) <= tolerance * target:
            add_root(prev_t)

        for i in range(1, samples + 1):
            cur_t = t_min_s + i * step
            cur_f = residual(cur_t)
            if isfinite(cur_f) and abs(cur_f) <= tolerance * target:
                add_root(cur_t)
            if prev_f * cur_f < 0.0:
                lo_t, hi_t = prev_t, cur_t
                lo_f, hi_f = prev_f, cur_f
                for _ in range(90):
                    mid_t = 0.5 * (lo_t + hi_t)
                    mid_f = residual(mid_t)
                    if abs(mid_f) <= tolerance * target:
                        lo_t = hi_t = mid_t
                        break
                    if lo_f * mid_f <= 0.0:
                        hi_t, hi_f = mid_t, mid_f
                    else:
                        lo_t, lo_f = mid_t, mid_f
                add_root(0.5 * (lo_t + hi_t))
            prev_t, prev_f = cur_t, cur_f

        roots.sort()
        return roots

    def distance_redshift_branches(
        self,
        target_z: float,
        t_obs_s: float,
        *,
        t_min_s: float | None = None,
        t_max_s: float | None = None,
        mass_model: MassModel = "constant",
        redshift_model: RedshiftModel = "expansion_extent",
        samples: int = 4096,
    ) -> list[RedshiftBranch]:
        """Return light-travel-distance branches for an observed redshift."""

        roots = self.redshift_branch_times(
            target_z,
            t_obs_s,
            t_min_s=t_min_s,
            t_max_s=t_max_s,
            mass_model=mass_model,
            redshift_model=redshift_model,
            samples=samples,
        )
        return [
            RedshiftBranch(
                branch_index=index,
                t_emit_s=t_emit,
                t_obs_s=t_obs_s,
                z_model=self.z_gutcp(
                    t_emit,
                    t_obs_s,
                    mass_model=mass_model,
                    redshift_model=redshift_model,
                ),
                light_travel_distance_m=C * (t_obs_s - t_emit),
            )
            for index, t_emit in enumerate(roots)
        ]

    def line_of_sight_hubble_s(
        self,
        t_s: float,
        *,
        h_model: LineOfSightHModel = "expansion_extent",
    ) -> float:
        """Candidate line-of-sight H for the Koksbang-Heinesen test interface.

        The cited FLRW-consistency standard uses the observational line-of-sight
        expansion rate H(z), not an arbitrary time derivative. This selector
        exposes the same denominator choices audited above. Endpoint clock
        factors are not accepted here because they are finite endpoint wavelength
        ratios, not a local expansion-rate denominator.
        """

        if h_model == "expansion_extent":
            return self.hubble_expansion_extent_s(t_s)
        if h_model == "radius_scale":
            return self.hubble_radius_s(t_s)
        if h_model in ("mills_ct_path", "hubble_path"):
            return self.hubble_s(t_s)
        raise ValueError(f"unknown line-of-sight H model: {h_model!r}")

    def line_of_sight_hubble_km_s_mpc(
        self,
        t_s: float,
        *,
        h_model: LineOfSightHModel = "expansion_extent",
    ) -> float:
        return self.line_of_sight_hubble_s(t_s, h_model=h_model) * MPC_KM

    def consistency_observable_branches(
        self,
        target_z: float,
        t_obs_s: float,
        *,
        redshift_model: RedshiftModel = "expansion_extent",
        h_model: LineOfSightHModel = "expansion_extent",
        mass_model: MassModel = "constant",
        t_min_s: float | None = None,
        t_max_s: float | None = None,
        samples: int = 4096,
    ) -> list[ConsistencyObservableBranch]:
        """Return the GUTCP quantities that meet part of the C/O/M interface.

        This deliberately returns angular_diameter_distance_m=None. The
        Koksbang-Heinesen diagnostics require d_A(z), d_A'(z), d_A''(z), H(z),
        and H'(z). GUTCPz plus a branch gives t_emit(z), H_parallel(z), and a
        light-travel distance surrogate; it does not by itself give an angular
        diameter distance.
        """

        branches = self.distance_redshift_branches(
            target_z,
            t_obs_s,
            t_min_s=t_min_s,
            t_max_s=t_max_s,
            mass_model=mass_model,
            redshift_model=redshift_model,
            samples=samples,
        )
        return [
            ConsistencyObservableBranch(
                target_z=target_z,
                branch_index=branch.branch_index,
                t_emit_s=branch.t_emit_s,
                t_obs_s=branch.t_obs_s,
                h_parallel_km_s_mpc=self.line_of_sight_hubble_km_s_mpc(
                    branch.t_emit_s,
                    h_model=h_model,
                ),
                light_travel_distance_m=branch.light_travel_distance_m,
            )
            for branch in branches
        ]


def calibrate_to_h0_and_cmb(
    h0_target_km_s_mpc: float,
    *,
    cmb_now_target_k: float = OBSERVED_CMB_NOW_K,
    initial_t_now_s: float = NOTEBOOK_T_NOW_S,
    initial_cmb_start_k: float = NOTEBOOK_CMBR_START_K,
    normalization: RadiusNormalization = "notebook",
    period_convention: PeriodConvention = "numeric",
    tolerance: float = 1.0e-11,
    max_iterations: int = 50,
) -> CalibrationResult:
    """Solve the notebook calibration problem for an explicit H0 target.

    This is the Python counterpart of replacing Mathematica's Subscript[H, 0]
    entity with an explicit H0TargetKmSecMpc value:

        H0GUTCP[tNow, CMBRTempAtStart] == H0TargetKmSecMpc
        TU[tNow, CMBRTempAtStart] == observed CMB temperature
    """

    u = log(initial_t_now_s)
    v = log(initial_cmb_start_k)

    def evaluate(log_t: float, log_temp: float) -> tuple[float, float, GUTCPModel, float, float]:
        t_s = exp(log_t)
        temp_start = exp(log_temp)
        model = GUTCPModel.from_cmbr_start(
            temp_start,
            normalization=normalization,
            period_convention=period_convention,
        )
        h0 = model.hubble_km_s_mpc(t_s)
        cmb_now = model.temperature_k(t_s)
        return (
            h0 / h0_target_km_s_mpc - 1.0,
            cmb_now / cmb_now_target_k - 1.0,
            model,
            h0,
            cmb_now,
        )

    f1, f2, model, h0_model, cmb_model = evaluate(u, v)
    for iteration in range(1, max_iterations + 1):
        if max(abs(f1), abs(f2)) < tolerance:
            return CalibrationResult(
                h0_target_km_s_mpc=h0_target_km_s_mpc,
                cmb_now_target_k=cmb_now_target_k,
                t_now_s=exp(u),
                cmb_start_k=exp(v),
                h0_model_km_s_mpc=h0_model,
                cmb_now_model_k=cmb_model,
                iterations=iteration - 1,
            )

        step = 1.0e-5
        f1_u, f2_u, *_ = evaluate(u + step, v)
        f1_v, f2_v, *_ = evaluate(u, v + step)
        j11 = (f1_u - f1) / step
        j21 = (f2_u - f2) / step
        j12 = (f1_v - f1) / step
        j22 = (f2_v - f2) / step
        det = j11 * j22 - j12 * j21
        if det == 0.0:
            raise RuntimeError("singular calibration Jacobian")

        du = (j22 * f1 - j12 * f2) / det
        dv = (-j21 * f1 + j11 * f2) / det
        old_norm = abs(f1) + abs(f2)
        damping = 1.0
        for _ in range(16):
            candidate_u = u - damping * du
            candidate_v = v - damping * dv
            nf1, nf2, nmodel, nh0, ncmb = evaluate(candidate_u, candidate_v)
            new_norm = abs(nf1) + abs(nf2)
            if isfinite(new_norm) and new_norm < old_norm:
                u, v = candidate_u, candidate_v
                f1, f2, model, h0_model, cmb_model = nf1, nf2, nmodel, nh0, ncmb
                break
            damping *= 0.5
        else:
            raise RuntimeError("calibration Newton step failed to improve")

    raise RuntimeError("calibration did not converge")


def _bisect_log_temperature(
    residual: Callable[[float], float],
    *,
    low_k: float = 0.5,
    high_k: float = 20.0,
    tolerance: float = 1.0e-11,
    max_iterations: int = 100,
) -> float:
    low = log(low_k)
    high = log(high_k)
    f_low = residual(exp(low))
    f_high = residual(exp(high))
    if f_low == 0.0:
        return exp(low)
    if f_high == 0.0:
        return exp(high)
    if f_low * f_high > 0.0:
        raise RuntimeError("temperature bracket does not straddle a root")

    for _ in range(max_iterations):
        mid = 0.5 * (low + high)
        temp = exp(mid)
        f_mid = residual(temp)
        if abs(f_mid) < tolerance:
            return temp
        if f_low * f_mid <= 0.0:
            high = mid
            f_high = f_mid
        else:
            low = mid
            f_low = f_mid
    return exp(0.5 * (low + high))


def fixed_time_overconstraint_diagnostic(
    *,
    t_fixed_s: float = years_to_seconds(1.0e10),
    h0_target_km_s_mpc: float = 67.66,
    cmb_now_target_k: float = OBSERVED_CMB_NOW_K,
    normalization: RadiusNormalization = "notebook",
    period_convention: PeriodConvention = "numeric",
) -> FixedTimeDiagnostic:
    """Show why fixing Mills's 10^10 yr overconstrains the calibration.

    With t fixed, CMBRTempAtStart is the only remaining unknown. Matching H0 and
    matching the present CMB temperature generally require different values of
    that one unknown.
    """

    def model_for(temp_start_k: float) -> GUTCPModel:
        return GUTCPModel.from_cmbr_start(
            temp_start_k,
            normalization=normalization,
            period_convention=period_convention,
        )

    cmb_start_for_h0 = _bisect_log_temperature(
        lambda temp: model_for(temp).hubble_km_s_mpc(t_fixed_s) - h0_target_km_s_mpc
    )
    model_h0 = model_for(cmb_start_for_h0)

    cmb_start_for_cmb = _bisect_log_temperature(
        lambda temp: model_for(temp).temperature_k(t_fixed_s) - cmb_now_target_k
    )
    model_cmb = model_for(cmb_start_for_cmb)

    return FixedTimeDiagnostic(
        t_fixed_s=t_fixed_s,
        h0_target_km_s_mpc=h0_target_km_s_mpc,
        cmb_now_target_k=cmb_now_target_k,
        cmb_start_for_h0_k=cmb_start_for_h0,
        cmb_now_when_h0_matches_k=model_h0.temperature_k(t_fixed_s),
        cmb_start_for_cmb_k=cmb_start_for_cmb,
        h0_when_cmb_matches_km_s_mpc=model_cmb.hubble_km_s_mpc(t_fixed_s),
    )


def h0_target_sweep(
    h0_targets_km_s_mpc: Iterable[float],
    *,
    cmb_now_target_k: float = OBSERVED_CMB_NOW_K,
) -> list[CalibrationResult]:
    """Run the explicit-H0 calibration for a sequence of targets."""

    results = []
    guess_t = NOTEBOOK_T_NOW_S
    guess_temp = NOTEBOOK_CMBR_START_K
    for target in h0_targets_km_s_mpc:
        result = calibrate_to_h0_and_cmb(
            target,
            cmb_now_target_k=cmb_now_target_k,
            initial_t_now_s=guess_t,
            initial_cmb_start_k=guess_temp,
        )
        results.append(result)
        guess_t = result.t_now_s
        guess_temp = result.cmb_start_k
    return results


LIGHT_BUNDLE_FRONTIER_ROWS = (
    LightBundleFrontierRow(
        "p.1569",
        "BOOMERANG 1-degree CMBR structures imply nearly flat geometry since expansion began",
        "angular closure must separate the near-flat geometry claim from the fixed 10 Gyr approximation used in the same sentence",
        "constraining input with calibration hazard",
    ),
    LightBundleFrontierRow(
        "p.1602 footnote",
        "absolute-rest r-sphere reference and Cepheid angular-diameter/radial-velocity distance calibration",
        "optical dictionary must specify rest-frame or velocity-correction rule and how quasi-geometrical Cepheid distances anchor scale",
        "constraining input",
    ),
    LightBundleFrontierRow(
        "p.1577 Eq.32.203",
        "E-mode CMBR multipoles from Thompson scattering",
        "polarization-aware light-bundle closure must preserve the stated ell-indexed E-mode phase/amplitude convention",
        "constraining input",
    ),
    LightBundleFrontierRow(
        "p.1577 Eqs.32.204-32.205",
        "B-mode is shifted by 70 in ell with amplitude ratio DeltaAleph/(c t)",
        "candidate lensing/accelerating-spacetime closure must explain E-to-B conversion and amplitude ratio",
        "constraining input",
    ),
    LightBundleFrontierRow(
        "p.1578",
        "gravitational lensing produces rings/arcs and is linked to matter-formation patterns",
        "light-bundle relation must decide how photon bundles map into observed angular structures",
        "constraining input",
    ),
    LightBundleFrontierRow(
        "p.1579",
        "CMBR multipole scale L_ell = 4 pi r_sphere / ell",
        "large-ring/arc comparisons should use the CMBR multipole ruler relation, not BAO transverse scale by default",
        "stated relation",
    ),
)


KH_CONSISTENCY_ROWS = (
    ConsistencyReadinessRow(
        "z",
        "independent redshift coordinate for d_A(z) and H(z)",
        "available only after choosing z_GUTCP mode and branch",
        "conditional; do not import FLRW 1+z=a_o/a_e",
    ),
    ConsistencyReadinessRow(
        "t_b(z)",
        "invert the observable redshift to index model quantities",
        "implemented as a generally multi-valued branch set",
        "conditional on branch/window/endpoint convention",
    ),
    ConsistencyReadinessRow(
        "H_parallel(z)",
        "line-of-sight expansion rate; chronometers are direct, radial BAO is ruler-conditional",
        "available as H_model(t_b(z)) for selected denominator convention; BAO also needs r_d handling",
        "partial; endpoint clock factors are not local H",
    ),
    ConsistencyReadinessRow(
        "H_parallel'(z)",
        "redshift derivative used by C and M diagnostics",
        "not promoted as a stable observable until a branch rule is fixed",
        "missing for full paper-equivalent test",
    ),
    ConsistencyReadinessRow(
        "d_A(z)",
        "angular diameter distance reconstructed from supernovae",
        "not supplied by audited GUTCP equations",
        "missing; c(t_o-t_e) is only a light-travel surrogate",
    ),
    ConsistencyReadinessRow(
        "gutcpLightBundleRel constraints",
        "frontier information constraining the missing d_A optical closure",
        "pages 1569, 1577-1579, and 1602 constrain geometry/rest-frame/polarization/lensing/multipole structure but do not close d_A",
        "inventory available via light_bundle_frontier_rows()",
    ),
    ConsistencyReadinessRow(
        "d_A'(z), d_A''(z)",
        "distance derivatives required by C, O, and M",
        "not supplied because d_A(z) is not supplied",
        "missing",
    ),
    ConsistencyReadinessRow(
        "BAO ruler relation",
        "BAO constrains D_M/r_d, H_parallel*r_d, or D_V/r_d, not bare d_A",
        "not supplied by the GUTCP angular map; can be joined only as a conditional data relation",
        "conditional; do not treat transverse BAO as model-independent d_A",
    ),
    ConsistencyReadinessRow(
        "C(z), O(z), M(z)",
        "diagnostic FLRW/general-spacetime consistency relations",
        "cannot be fairly computed from current GUTCP artifact",
        "not yet evaluable at the paper's standard",
    ),
    ConsistencyReadinessRow(
        "uncertainty band",
        "bootstrap percentile envelope over symbolic reconstructions",
        "no GUTCP observational-error pipeline is present here",
        "missing; needed for sigma-level claims",
    ),
)


def light_bundle_frontier_rows() -> tuple[LightBundleFrontierRow, ...]:
    """Source constraints that any candidate GUTCP light-bundle law must face."""

    return LIGHT_BUNDLE_FRONTIER_ROWS


def flrw_consistency_readiness_rows() -> tuple[ConsistencyReadinessRow, ...]:
    """Koksbang-Heinesen standard applied to the current GUTCP artifact."""

    return KH_CONSISTENCY_ROWS


def format_consistency_readiness_table(
    rows: Iterable[ConsistencyReadinessRow] = KH_CONSISTENCY_ROWS,
) -> str:
    headers = ("Quantity", "Paper role", "GUTCP status", "Fair evaluation")
    material = [
        headers,
        *[
            (row.quantity, row.paper_role, row.gutcp_status, row.fair_evaluation)
            for row in rows
        ],
    ]
    widths = [len(header) for header in headers]
    for row in material:
        for i, cell in enumerate(row):
            widths[i] = max(widths[i], len(cell))

    def fmt(row: Iterable[str]) -> str:
        return " | ".join(cell.ljust(widths[i]) for i, cell in enumerate(row))

    lines = [fmt(headers), "-+-".join("-" * width for width in widths)]
    lines.extend(fmt(row) for row in material[1:])
    return "\n".join(lines)


def format_light_bundle_frontier_table(
    rows: Iterable[LightBundleFrontierRow] = LIGHT_BUNDLE_FRONTIER_ROWS,
) -> str:
    headers = ("Source", "Relation fragment", "Constraint", "Status")
    material = [
        headers,
        *[
            (row.source, row.relation_fragment, row.constraint, row.status)
            for row in rows
        ],
    ]
    widths = [len(header) for header in headers]
    for row in material:
        for i, cell in enumerate(row):
            widths[i] = max(widths[i], len(cell))

    def fmt(row: Iterable[str]) -> str:
        return " | ".join(cell.ljust(widths[i]) for i, cell in enumerate(row))

    lines = [fmt(headers), "-+-".join("-" * width for width in widths)]
    lines.extend(fmt(row) for row in material[1:])
    return "\n".join(lines)


AUDIT_ROWS = (
    ("32.126", "1536", "v = H R", "redshift path premise", "supporting"),
    ("32.36", "1519", "r_g = 2 G m_0 / c^2", "2 G m/c^2", "matched"),
    ("32.38", "1520", "modified Schwarzschild metric", "implicit", "supporting"),
    ("32.43", "1523", "tau = r_g/c = (v_g/c) t_i", "time metric", "supporting; OCR/dimensional risk"),
    ("32.140", "1542", "Q = c^3/(4 pi G)", "c^3/(4 Pi G)", "matched"),
    ("32.146", "1544", "R_min from Stefan-Boltzmann chain", "alephMin", "mismatch/ambiguous normalization"),
    ("32.149", "1545", "T_U = 4 pi G M/c^3, numeric implies 2 pi", "period", "mismatch"),
    ("32.153", "1546", "R(t) harmonic radius", "aleph[t_]", "matched to numeric convention"),
    ("32.154", "1546", "dR/dt = 4 pi c sin(...)", "alephRate[t_]", "matched"),
    ("32.156", "1547", "H(t) = dR/dt / (c t)", "H[t_]", "matched; denominator under audit"),
    ("32.158", "1548", "M_U(t) = M_U/2 (1 + cos(...))", "mUBook[t_]", "matched; commentary conflict"),
    ("32.165", "1551", "lambda_inf = lambda(r) (1 + 2GM/(r c^2))", "zGUTCP", "endpoint correction"),
    ("unnumbered", "1579", "L_ell = (2/ell) 2 pi r_sphere", "cmbr_multipole_structure_scale", "stated; not BAO"),
)


def format_audit_table(rows: Iterable[tuple[str, str, str, str, str]] = AUDIT_ROWS) -> str:
    headers = ("Equation", "Page", "Extracted", "Notebook", "Status")
    widths = [len(h) for h in headers]
    material = [headers, *rows]
    for row in material:
        for i, cell in enumerate(row):
            widths[i] = max(widths[i], len(cell))

    def fmt(row: Iterable[str]) -> str:
        return " | ".join(cell.ljust(widths[i]) for i, cell in enumerate(row))

    lines = [fmt(headers), "-+-".join("-" * width for width in widths)]
    lines.extend(fmt(row) for row in rows)
    return "\n".join(lines)


def demo() -> None:
    model = GUTCPModel.from_cmbr_start()
    book_model = GUTCPModel.from_book_mass()
    t_now = NOTEBOOK_T_NOW_S

    print("GUTCP 2023 redshift audit")
    print()
    print(format_audit_table())
    print()
    print("Default model: notebook CMBR-start solution")
    print(f"  M_start = {model.m_start_kg:.6e} kg")
    print(f"  R_min   = {meters_to_light_years(model.gravitational_radius_m):.6e} ly")
    print(f"  T_U     = {seconds_to_billion_years(model.period_s):.6f} billion years")
    print(f"  H(t_now notebook) = {model.hubble_km_s_mpc(t_now):.6f} km/s/Mpc")
    print()
    print("Book-mass check at t = 10^10 years")
    print(f"  H = {book_model.hubble_km_s_mpc(years_to_seconds(1.0e10)):.6f} km/s/Mpc")
    print()
    fixed_t = fixed_time_overconstraint_diagnostic()
    print("Fixed t=10^10 yr diagnostic: one unknown cannot satisfy both constraints")
    print(
        "  matching H0 requires CMBR_start_K="
        f"{fixed_t.cmb_start_for_h0_k:.6f}, then CMB_now_K="
        f"{fixed_t.cmb_now_when_h0_matches_k:.6f}"
    )
    print(
        "  matching CMB_now requires CMBR_start_K="
        f"{fixed_t.cmb_start_for_cmb_k:.6f}, then H0="
        f"{fixed_t.h0_when_cmb_matches_km_s_mpc:.6f} km/s/Mpc"
    )
    print()
    print("H0 target sweep: calibration, not independent prediction")
    print("  H0_target | t_now_Gyr | CMBR_start_K | H0_model | CMB_now_K")
    for result in h0_target_sweep((67.4, 67.66, 70.0, 73.0, 78.5)):
        print(
            f"  {result.h0_target_km_s_mpc:9.2f} | "
            f"{result.t_now_billion_years:9.3f} | "
            f"{result.cmb_start_k:12.6f} | "
            f"{result.h0_model_km_s_mpc:8.3f} | "
            f"{result.cmb_now_model_k:9.5f}"
        )
    print()
    print("Denominator audit at calibrated t_now")
    print(f"  H_ct                 = {model.hubble_km_s_mpc(t_now):.6f} km/s/Mpc")
    print(f"  H_whole_radius       = {model.hubble_radius_km_s_mpc(t_now):.6f} km/s/Mpc")
    print(
        f"  H_expansion_extent  = "
        f"{model.hubble_expansion_extent_km_s_mpc(t_now):.6f} km/s/Mpc"
    )
    print()
    print("GUTCPz denominator variants for selected emission times")
    print("  emit_Gyr | extent_z | radius_z | ct_path_z | endpoint_z | combined_extent")

    def fmt(value: float, width: int = 10) -> str:
        if value == float("inf"):
            return "inf".rjust(width)
        if value == float("-inf"):
            return "-inf".rjust(width)
        return f"{value:{width}.6f}"

    for emit_gyr in (0.0, 1.0, 10.0, 50.0, 100.0, 110.0):
        t_emit = billion_years_to_seconds(emit_gyr)
        print(
            f"  {emit_gyr:8.1f} | "
            f"{fmt(model.z_gutcp(t_emit, t_now, redshift_model='expansion_extent'), 10)} | "
            f"{model.z_gutcp(t_emit, t_now, redshift_model='radius_scale'):8.6f} | "
            f"{model.z_gutcp(t_emit, t_now, redshift_model='mills_ct_path'):9.6f} | "
            f"{model.z_gutcp(t_emit, t_now, redshift_model='mills_endpoint'):10.6f} | "
            f"{fmt(model.z_gutcp(t_emit, t_now, redshift_model='combined_extent_endpoint'), 15)}"
        )
    print()
    print("Distance-redshift branches from corrected expansion-extent GUTCPz")
    print("  target_z | branch | emit_Gyr | lookback_Gyr | light_Gly | z_model")
    for target_z in (1.0, 6.0, 10.0, 1000.0):
        branches = model.distance_redshift_branches(
            target_z,
            t_now,
            t_min_s=years_to_seconds(1.0),
            t_max_s=t_now,
            redshift_model="expansion_extent",
        )
        for branch in branches:
            print(
                f"  {target_z:8.1f} | "
                f"{branch.branch_index:6d} | "
                f"{branch.t_emit_billion_years:8.3f} | "
                f"{branch.lookback_billion_years:12.3f} | "
                f"{branch.light_travel_distance_billion_light_years:9.3f} | "
                f"{branch.z_model:8.3f}"
            )
    print()
    entry_threshold = blue_shift_segment_entry_threshold_m(t_now)
    neutral_threshold = symmetric_endpoint_neutral_threshold_m(t_now)
    print("Blueshift distance thresholds from calibrated t_now")
    print(
        "  pre-current-expansion segment enters at "
        f"D = c t_now = {meters_to_light_years(entry_threshold) / 1.0e9:.3f} Gly"
    )
    print(
        "  symmetric endpoint-neutral threshold at "
        f"D = 2 c t_now = {meters_to_light_years(neutral_threshold) / 1.0e9:.3f} Gly"
    )
    print("  distance_Gly | status")
    for distance in (0.5 * entry_threshold, 1.5 * entry_threshold, neutral_threshold, 2.5 * entry_threshold):
        print(
            f"  {meters_to_light_years(distance) / 1.0e9:12.3f} | "
            f"{blueshift_threshold_status(distance, t_now)}"
        )
    print("  redshift-sign branch pruning examples for observed z = 1")
    for distance in (1.5 * entry_threshold, 2.5 * entry_threshold):
        print(
            f"  {meters_to_light_years(distance) / 1.0e9:12.3f} | "
            f"{'; '.join(blueshift_branch_constraints(distance, t_now, observed_z=1.0))}"
        )
    print()
    print("Mills page-1579 CMBR multipole structure-scale relation")
    print(
        "  r_sphere = "
        f"{MILLS_PAGE_1579_R_SPHERE_LIGHT_YEARS / 1.0e9:.2f} Gly; "
        f"angular view = {cmbr_angular_view_light_years() / 1.0e9:.2f} Gly"
    )
    print("  ell | sky_fraction | scale_Gly")
    for ell in (15.0, 135.0, 700.0):
        print(
            f"  {ell:3.0f} | "
            f"{cmbr_multipole_sky_fraction(ell):12.5f} | "
            f"{cmbr_multipole_structure_scale_light_years(ell) / 1.0e9:9.3f}"
        )
    print()
    print("Koksbang-Heinesen FLRW-consistency standard applied to GUTCP")
    print(format_consistency_readiness_table())
    print()
    print("GUTCP light-bundle frontier constraints")
    print(format_light_bundle_frontier_table())
    print()
    print("Available branch observables over the paper's low-z domain")
    print("  target_z | branch | emit_Gyr | H_parallel | light_Gly | d_A status")
    for target_z in (0.38, 1.0, 2.0):
        observables = model.consistency_observable_branches(
            target_z,
            t_now,
            t_min_s=years_to_seconds(1.0),
            t_max_s=t_now,
            redshift_model="expansion_extent",
            h_model="expansion_extent",
        )
        for observable in observables:
            print(
                f"  {target_z:8.2f} | "
                f"{observable.branch_index:6d} | "
                f"{observable.t_emit_billion_years:8.3f} | "
                f"{observable.h_parallel_km_s_mpc:10.3f} | "
                f"{observable.light_travel_distance_billion_light_years:9.3f} | "
                f"{observable.status}"
            )


if __name__ == "__main__":
    demo()
