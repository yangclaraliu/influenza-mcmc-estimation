"""ABC-MCMC estimation of household and school influenza transmission.

The script estimates one pair of parameters for each of elementary, middle,
and high schools. Kindergarten is intentionally excluded.

For a quick smoke test:
    python "MCMC estimation(all schools).py" --quick

For the full analysis (the defaults can still take a long time):
    python "MCMC estimation(all schools).py"
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


# ------------------------------- Configuration ------------------------------

DAYS = 148
POPULATION = 1_886_100
UNDERESTIMATE = 0.354 * 0.56
VACCINE_EFFICACY = 0.60
STUDENTS_PER_CLASS = 30
CLASSES_PER_GRADE = 4
STUDENTS_PER_GRADE = STUDENTS_PER_CLASS * CLASSES_PER_GRADE
HOUSEHOLD_COUNTS = np.array([2594, 1193, 717, 222, 45, 10, 2, 6])

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR
AGE_DISTRIBUTION_FILE = DATA_DIR / "流感数据年龄分布.xlsx"
OUTPUT_DIR = SCRIPT_DIR

# The existing scripts retain the Uniform(0, 1) household prior by default.
# A school-specific launcher can switch to the truncated exponential prior.
HOUSEHOLD_PRIOR = "uniform"
HOUSEHOLD_EXPONENTIAL_SCALE = 0.32


@dataclass(frozen=True)
class SchoolConfig:
    key: str
    display_name: str
    data_file: str
    contact_probabilities: tuple[float, float, float]
    susceptibility: float
    age_distribution_column: int


SCHOOL_CONFIGS = {
    "elementary": SchoolConfig(
        key="elementary",
        display_name="小学",
        data_file="流感疫苗数据库3_小学.xlsx",
        contact_probabilities=(0.102, 0.0102, 0.00102),
        susceptibility=0.88,
        age_distribution_column=8,
    ),
    "middle": SchoolConfig(
        key="middle",
        display_name="初中",
        data_file="流感疫苗数据库4_初中.xlsx",
        contact_probabilities=(0.16, 0.016, 0.0016),
        susceptibility=0.88,
        age_distribution_column=9,
    ),
    "high": SchoolConfig(
        key="high",
        display_name="高中",
        data_file="流感疫苗数据库5_高中.xlsx",
        contact_probabilities=(0.151, 0.0151, 0.00151),
        susceptibility=0.88,
        age_distribution_column=10,
    ),
}


@dataclass
class SchoolDataset:
    school_ids: np.ndarray
    vaccination_rates: np.ndarray
    infection_rates: np.ndarray
    student_counts: np.ndarray

    @property
    def n_schools(self) -> int:
        return self.vaccination_rates.shape[0]

    @property
    def n_grades(self) -> int:
        return self.vaccination_rates.shape[1]


# ------------------------------- Data loading --------------------------------

def load_school_data(path: Path) -> SchoolDataset:
    """Aggregate student records into school-by-grade observed rates."""
    raw = pd.read_excel(path, engine="openpyxl")
    if raw.shape[1] <= 7:
        raise ValueError(f"{path.name} has fewer than eight columns.")

    frame = pd.DataFrame(
        {
            "flu": pd.to_numeric(raw.iloc[:, 1], errors="coerce"),
            "vaccinated": pd.to_numeric(raw.iloc[:, 2], errors="coerce"),
            "school": pd.to_numeric(raw.iloc[:, 5], errors="coerce"),
            "grade": pd.to_numeric(raw.iloc[:, 7], errors="coerce"),
        }
    ).dropna(subset=["school", "grade"])
    frame = frame[(frame["school"] > 0) & (frame["grade"] > 0)].copy()
    frame["school"] = frame["school"].astype(int)
    frame["grade"] = frame["grade"].astype(int)
    frame["vaccinated"] = frame["vaccinated"].fillna(0).gt(0).astype(int)
    frame["flu"] = frame["flu"].fillna(0).gt(0).astype(int)

    grouped = frame.groupby(["school", "grade"], sort=True).agg(
        students=("flu", "size"),
        vaccinated=("vaccinated", "sum"),
        infected=("flu", "sum"),
    )
    schools = np.sort(frame["school"].unique())
    n_grades = int(frame["grade"].max())
    shape = (len(schools), n_grades)
    vaccination_rates = np.full(shape, np.nan)
    infection_rates = np.full(shape, np.nan)
    student_counts = np.zeros(shape, dtype=int)
    school_to_row = {school: row for row, school in enumerate(schools)}

    for (school, grade), row in grouped.iterrows():
        i, j = school_to_row[school], grade - 1
        student_counts[i, j] = int(row["students"])
        vaccination_rates[i, j] = row["vaccinated"] / row["students"]
        infection_rates[i, j] = row["infected"] / row["students"]

    complete_schools = np.all(student_counts > 0, axis=1)
    removed_schools = schools[~complete_schools]
    if len(removed_schools):
        print(
            f"Removing schools with missing grade data from {path.name}: "
            + ", ".join(map(str, removed_schools))
        )
    schools = schools[complete_schools]
    vaccination_rates = vaccination_rates[complete_schools]
    infection_rates = infection_rates[complete_schools]
    student_counts = student_counts[complete_schools]

    if student_counts.sum() == 0:
        raise ValueError(f"No usable student records were found in {path}.")
    return SchoolDataset(
        school_ids=schools,
        vaccination_rates=vaccination_rates,
        infection_rates=infection_rates,
        student_counts=student_counts,
    )


def load_student_introduction_fraction(column: int) -> np.ndarray:
    """Load the daily fraction of community introductions affecting students."""
    raw = pd.read_excel(AGE_DISTRIBUTION_FILE, engine="openpyxl")
    if raw.shape[1] <= column:
        raise ValueError(
            f"{AGE_DISTRIBUTION_FILE.name} has no zero-based column {column}."
        )
    values = pd.to_numeric(raw.iloc[:, column], errors="coerce").fillna(0).to_numpy()
    if len(values) < DAYS:
        raise ValueError(
            f"Age-distribution column {column} has {len(values)} rows; "
            f"{DAYS} are required."
        )
    values = values[:DAYS].astype(float)
    if np.nanmax(values) > 1:
        values /= 100.0
    return np.clip(values, 0, 1)


# -------------------------- Influenza simulation -----------------------------

def community_lambda(day: float | np.ndarray) -> float | np.ndarray:
    """Piecewise community incidence curve; supports scalar and array input."""
    x = np.asarray(day, dtype=float)
    result = np.empty_like(x)
    first = x <= 71
    second = (x > 71) & (x <= 106)
    third = x > 106
    result[first] = -2.81716876e-03 * x[first] + 1.19276201e01
    result[second] = (
        -9.91888439e-03 * x[second] ** 3
        + 2.47565541 * x[second] ** 2
        - 2.00775045e02 * x[second]
        + 5.33962575e03
    )
    result[third] = (
        2.70834057e-03 * x[third] ** 3
        - 1.10831619 * x[third] ** 2
        + 1.49048335e02 * x[third]
        - 6.52465374e03
    )
    result = np.maximum(result, 0)
    return float(result) if result.ndim == 0 else result


def make_school_population(n_grades: int) -> tuple[np.ndarray, np.ndarray]:
    grades = np.repeat(np.arange(n_grades), STUDENTS_PER_GRADE)
    classes = np.concatenate(
        [
            np.repeat(
                np.arange(
                    grade * CLASSES_PER_GRADE,
                    (grade + 1) * CLASSES_PER_GRADE,
                ),
                STUDENTS_PER_CLASS,
            )
            for grade in range(n_grades)
        ]
    )
    return grades, classes


def generate_school_network(
    grades: np.ndarray,
    classes: np.ndarray,
    probabilities: tuple[float, float, float],
    rng: np.random.Generator,
) -> np.ndarray:
    p_class, p_grade, p_between = probabilities
    same_grade = grades[:, None] == grades[None, :]
    same_class = classes[:, None] == classes[None, :]
    edge_probability = np.where(
        same_class,
        p_class,
        np.where(same_grade, p_grade, p_between),
    )
    upper = np.triu(rng.random(edge_probability.shape) < edge_probability, k=1)
    return upper | upper.T


def generate_susceptibility(
    grade_vaccination_rates: np.ndarray,
    susceptibility: float,
    rng: np.random.Generator,
) -> np.ndarray:
    vaccination_rates = np.repeat(
        grade_vaccination_rates, STUDENTS_PER_GRADE
    )
    biologically_susceptible = rng.random(len(vaccination_rates)) < susceptibility
    vaccine_protected = (
        (rng.random(len(vaccination_rates)) < vaccination_rates)
        & (rng.random(len(vaccination_rates)) < VACCINE_EFFICACY)
    )
    return biologically_susceptible & ~vaccine_protected


def generate_external_introductions(
    susceptible: np.ndarray,
    student_fraction: np.ndarray,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Generate household and direct-student introductions by NHPP thinning."""
    total_students = len(susceptible)
    household_sizes = np.arange(3, 3 + len(HOUSEHOLD_COUNTS))
    size_probability = (
        household_sizes * HOUSEHOLD_COUNTS
        / np.sum(household_sizes * HOUSEHOLD_COUNTS)
    )
    sampled_sizes = rng.choice(
        household_sizes,
        size=total_students,
        replace=True,
        p=size_probability,
    )
    family_members = sampled_sizes - 1
    family_population = int(family_members.sum())

    day_grid = np.arange(1, DAYS + 1)
    scaled_rates = (
        community_lambda(day_grid)
        / UNDERESTIMATE
        * (family_population + total_students)
        / POPULATION
    )
    maximum_rate = float(np.max(scaled_rates))
    household_events = np.zeros((total_students, DAYS), dtype=np.int16)
    student_events = np.zeros((total_students, DAYS), dtype=bool)
    if maximum_rate <= 0:
        return household_events, student_events

    family_target_probability = family_members / family_members.sum()
    time = 0.0
    while True:
        time += rng.exponential(1.0 / maximum_rate)
        if time >= DAYS:
            break
        actual_rate = (
            community_lambda(time + 1)
            / UNDERESTIMATE
            * (family_population + total_students)
            / POPULATION
        )
        if rng.random() > actual_rate / maximum_rate:
            continue

        day = min(int(time), DAYS - 1)
        if rng.random() <= student_fraction[day]:
            # Target any student; run_sir infects only susceptible ones, so
            # introductions reaching immune or vaccine-protected students fail.
            target = rng.integers(total_students)
            student_events[target, day] = True
        else:
            target = rng.choice(total_students, p=family_target_probability)
            household_events[target, day] += 1

    return household_events, student_events


def run_sir(
    household_beta: float,
    school_beta: float,
    grades: np.ndarray,
    adjacency: np.ndarray,
    susceptible: np.ndarray,
    household_events: np.ndarray,
    student_events: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """Run a one-day-infectious-period SIR model and return grade attack rates."""
    n_grades = int(grades.max()) + 1
    status = np.zeros(len(grades), dtype=np.int8)  # 0=S, 1=I, 2=R
    infection_counts = np.zeros(n_grades, dtype=int)

    for day in range(DAYS):
        new_infected: set[int] = set()
        currently_susceptible = (status == 0) & susceptible

        direct_targets = np.flatnonzero(
            currently_susceptible & student_events[:, day]
        )
        new_infected.update(direct_targets.tolist())

        household_targets = np.flatnonzero(
            currently_susceptible & (household_events[:, day] > 0)
        )
        for target in household_targets:
            exposures = int(household_events[target, day])
            probability = 1 - (1 - household_beta) ** exposures
            if rng.random() < probability:
                new_infected.add(int(target))

        infected = np.flatnonzero(status == 1)
        for source in infected:
            targets = np.flatnonzero(adjacency[source] & currently_susceptible)
            for target in targets:
                if rng.random() < school_beta:
                    new_infected.add(int(target))

        status[infected] = 2
        if new_infected:
            new_array = np.fromiter(new_infected, dtype=int)
            new_array = new_array[(status[new_array] == 0) & susceptible[new_array]]
            status[new_array] = 1
            infection_counts += np.bincount(
                grades[new_array], minlength=n_grades
            )

    return infection_counts / STUDENTS_PER_GRADE


def simulate_school_type(
    parameters: np.ndarray,
    dataset: SchoolDataset,
    config: SchoolConfig,
    student_fraction: np.ndarray,
    replicates: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Simulate the school-by-grade attack-rate summary used by ABC."""
    household_beta, school_beta = parameters
    grades, classes = make_school_population(dataset.n_grades)
    simulated = np.zeros_like(dataset.vaccination_rates, dtype=float)

    for school_index in range(dataset.n_schools):
        school_total = np.zeros(dataset.n_grades)
        for _ in range(replicates):
            susceptible = generate_susceptibility(
                dataset.vaccination_rates[school_index],
                config.susceptibility,
                rng,
            )
            adjacency = generate_school_network(
                grades, classes, config.contact_probabilities, rng
            )
            household_events, student_events = generate_external_introductions(
                susceptible, student_fraction, rng
            )
            school_total += run_sir(
                household_beta,
                school_beta,
                grades,
                adjacency,
                susceptible,
                household_events,
                student_events,
                rng,
            )
        simulated[school_index] = school_total / replicates
    return simulated


def weighted_distance(simulated: np.ndarray, dataset: SchoolDataset) -> float:
    """Student-count-weighted RMSE between simulated and observed attack rates."""
    observed = dataset.infection_rates
    weights = dataset.student_counts
    valid = (weights > 0) & np.isfinite(observed)
    squared_error = (simulated[valid] - observed[valid]) ** 2
    return float(np.sqrt(np.average(squared_error, weights=weights[valid])))


def simulate_distance(
    parameters: np.ndarray,
    dataset: SchoolDataset,
    config: SchoolConfig,
    student_fraction: np.ndarray,
    replicates: int,
    rng: np.random.Generator,
) -> float:
    simulated = simulate_school_type(
        parameters, dataset, config, student_fraction, replicates, rng
    )
    return weighted_distance(simulated, dataset)


# ------------------------------- ABC-MCMC ------------------------------------

def reflect_unit_interval(value: float) -> float:
    """Reflect a random-walk proposal into the closed unit interval."""
    reflected = value % 2.0
    if reflected > 1:
        reflected = 2 - reflected
    return float(np.clip(reflected, np.finfo(float).eps, 1 - np.finfo(float).eps))


def log_uniform_prior(value: float) -> float:
    if not 0 < value < 1:
        return -np.inf
    # Uniform(0, 1) density equals one, so its log density equals zero.
    return 0.0


def log_household_prior(value: float) -> float:
    """Log density of the configured household-transmission prior."""
    if HOUSEHOLD_PRIOR == "uniform":
        return log_uniform_prior(value)
    if HOUSEHOLD_PRIOR != "truncated_exponential":
        raise ValueError(f"Unknown household prior: {HOUSEHOLD_PRIOR}")
    if not 0 < value < 1:
        return -np.inf
    scale = HOUSEHOLD_EXPONENTIAL_SCALE
    normalization = -math.expm1(-1.0 / scale)
    return -math.log(scale) - value / scale - math.log(normalization)


def sample_household_prior(rng: np.random.Generator) -> float:
    """Sample the household prior, respecting its probability support."""
    if HOUSEHOLD_PRIOR == "uniform":
        return float(rng.uniform(0, 1))
    if HOUSEHOLD_PRIOR != "truncated_exponential":
        raise ValueError(f"Unknown household prior: {HOUSEHOLD_PRIOR}")
    scale = HOUSEHOLD_EXPONENTIAL_SCALE
    truncation_mass = -math.expm1(-1.0 / scale)
    value = -scale * math.log1p(-rng.random() * truncation_mass)
    return float(np.clip(value, np.finfo(float).eps, 1 - np.finfo(float).eps))


def use_truncated_exponential_household_prior(scale: float = 0.32) -> None:
    """Use Exp(scale), conditioned on beta_household being in (0, 1)."""
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("The exponential prior scale must be positive.")
    global HOUSEHOLD_PRIOR, HOUSEHOLD_EXPONENTIAL_SCALE
    HOUSEHOLD_PRIOR = "truncated_exponential"
    HOUSEHOLD_EXPONENTIAL_SCALE = float(scale)


def log_prior(parameters: np.ndarray) -> float:
    return log_household_prior(parameters[0]) + log_uniform_prior(parameters[1])


def calibrate_epsilon(
    dataset: SchoolDataset,
    config: SchoolConfig,
    student_fraction: np.ndarray,
    draws: int,
    quantile: float,
    replicates: int,
    rng: np.random.Generator,
) -> tuple[float, pd.DataFrame]:
    """Choose epsilon from prior-predictive distance simulations."""
    records = []
    for draw in range(draws):
        parameters = np.array(
            [sample_household_prior(rng), rng.uniform(0, 1)],
            dtype=float,
        )
        distance = simulate_distance(
            parameters,
            dataset,
            config,
            student_fraction,
            replicates,
            rng,
        )
        records.append(
            {
                "draw": draw,
                "household_beta": parameters[0],
                "school_beta": parameters[1],
                "distance": distance,
            }
        )
        print(
            f"[{config.display_name}] calibration {draw + 1}/{draws}, "
            f"distance={distance:.4f}"
        )
    calibration = pd.DataFrame(records)
    epsilon = float(calibration["distance"].quantile(quantile))
    return epsilon, calibration


def abc_mcmc_chain(
    dataset: SchoolDataset,
    config: SchoolConfig,
    student_fraction: np.ndarray,
    epsilon: float,
    initial_parameters: np.ndarray,
    initial_distance: float,
    iterations: int,
    proposal_std: float,
    replicates: int,
    rng: np.random.Generator,
    chain: int,
    iteration_offset: int = 0,
    stage: str = "sampling",
) -> pd.DataFrame:
    """Run hard-threshold ABC-MCMC with a symmetric reflected random walk."""
    current = initial_parameters.astype(float).copy()
    current_distance = float(initial_distance)
    current_log_prior = log_prior(current)
    records = []
    accepted_count = 0

    for local_iteration in range(iterations):
        iteration = iteration_offset + local_iteration
        proposal = np.array(
            [
                reflect_unit_interval(
                    rng.normal(current[0], proposal_std)
                ),
                reflect_unit_interval(
                    rng.normal(current[1], proposal_std)
                ),
            ]
        )
        proposal_distance = simulate_distance(
            proposal,
            dataset,
            config,
            student_fraction,
            replicates,
            rng,
        )
        accepted = False
        if proposal_distance <= epsilon:
            proposal_log_prior = log_prior(proposal)
            log_acceptance = min(0.0, proposal_log_prior - current_log_prior)
            if math.log(rng.random()) < log_acceptance:
                current = proposal
                current_distance = proposal_distance
                current_log_prior = proposal_log_prior
                accepted = True
                accepted_count += 1

        records.append(
            {
                "chain": chain,
                "iteration": iteration,
                "proposal_household_beta": proposal[0],
                "proposal_school_beta": proposal[1],
                "proposal_distance": proposal_distance,
                "proposal_std": proposal_std,
                "household_beta": current[0],
                "school_beta": current[1],
                "distance": current_distance,
                "accepted": accepted,
                "sampling_stage": stage,
            }
        )
        if (local_iteration + 1) % max(1, iterations // 20) == 0:
            print(
                f"[{config.display_name}] {stage}, chain {chain}: "
                f"{local_iteration + 1}/{iterations}, "
                f"acceptance={accepted_count / (local_iteration + 1):.3f}"
            )
    return pd.DataFrame(records)


def gelman_rubin(samples: pd.DataFrame, parameter: str) -> float:
    chains = [
        group[parameter].to_numpy()
        for _, group in samples.groupby("chain", sort=True)
    ]
    if len(chains) < 2:
        return np.nan
    n = min(map(len, chains))
    values = np.stack([chain[:n] for chain in chains])
    within = np.mean(np.var(values, axis=1, ddof=1))
    between = n * np.var(np.mean(values, axis=1), ddof=1)
    if within == 0:
        return np.nan
    variance = ((n - 1) / n) * within + between / n
    return float(np.sqrt(variance / within))


def autocorrelation(values: np.ndarray, max_lag: int) -> np.ndarray:
    """Calculate autocorrelations from lag zero through max_lag."""
    values = np.asarray(values, dtype=float)
    max_lag = min(max_lag, len(values) - 1)
    if max_lag < 0:
        return np.array([])
    centered = values - np.mean(values)
    variance = np.dot(centered, centered)
    if variance == 0:
        return np.ones(max_lag + 1)
    result = np.ones(max_lag + 1)
    for lag in range(1, max_lag + 1):
        result[lag] = np.dot(centered[:-lag], centered[lag:]) / variance
    return result


def effective_sample_size(values: np.ndarray) -> float:
    """Estimate ESS using the initial-positive autocorrelation sequence."""
    values = np.asarray(values, dtype=float)
    if len(values) < 3 or np.var(values) == 0:
        return 1.0
    correlations = autocorrelation(values, min(1000, len(values) - 1))
    correlation_sum = 0.0
    for lag in range(1, len(correlations), 2):
        pair_sum = correlations[lag]
        if lag + 1 < len(correlations):
            pair_sum += correlations[lag + 1]
        if pair_sum <= 0:
            break
        correlation_sum += pair_sum
    integrated_time = max(1.0, 1.0 + 2.0 * correlation_sum)
    return float(min(len(values), len(values) / integrated_time))


def compute_chain_diagnostics(
    samples: pd.DataFrame,
    burn_in: int,
) -> pd.DataFrame:
    """Compute acceptance, lag-one autocorrelation, and ESS for each chain."""
    rows = []
    for chain, chain_samples in samples.groupby("chain", sort=True):
        posterior = chain_samples[chain_samples["iteration"] >= burn_in]
        acceptance_all = float(chain_samples["accepted"].mean())
        acceptance_post = float(posterior["accepted"].mean())
        for parameter in ("household_beta", "school_beta"):
            values = posterior[parameter].to_numpy()
            correlations = autocorrelation(values, 1)
            rows.append(
                {
                    "chain": chain,
                    "parameter": parameter,
                    "posterior_samples": len(values),
                    "mean": np.mean(values),
                    "std": np.std(values, ddof=1),
                    "acceptance_rate_all": acceptance_all,
                    "acceptance_rate_post_burn_in": acceptance_post,
                    "lag1_autocorrelation": (
                        correlations[1] if len(correlations) > 1 else np.nan
                    ),
                    "ess": effective_sample_size(values),
                }
            )
    return pd.DataFrame(rows)


def evaluate_pilot_acceptance(
    samples: pd.DataFrame,
    burn_in: int,
    minimum: float,
    maximum: float,
) -> tuple[bool, pd.DataFrame, str]:
    """Decide automatically whether the formal stage should be started."""
    posterior = samples[samples["iteration"] >= burn_in]
    acceptance = (
        posterior.groupby("chain", sort=True)["accepted"].mean().rename(
            "acceptance_rate"
        )
    )
    decision = acceptance.reset_index()
    decision["minimum_required"] = minimum
    decision["maximum_allowed"] = maximum
    decision["status"] = np.where(
        decision["acceptance_rate"] < minimum,
        "too_low",
        np.where(
            decision["acceptance_rate"] > maximum,
            "too_high",
            "acceptable",
        ),
    )
    decision["recommended_action"] = np.select(
        [
            decision["status"] == "too_low",
            decision["status"] == "too_high",
        ],
        [
            "decrease_proposal_std_or_loosen_epsilon",
            "increase_proposal_std_or_check_epsilon",
        ],
        default="continue_to_formal_stage",
    )
    failed = decision[decision["status"] != "acceptable"]
    if failed.empty:
        return True, decision, "all_chains_acceptable"

    reasons = [
        f"chain_{int(row.chain)}_{row.status}"
        for row in failed.itertuples(index=False)
    ]
    return False, decision, ";".join(reasons)


def summarize_samples(
    samples: pd.DataFrame,
    burn_in: int,
    config: SchoolConfig,
    epsilon: float,
    chain_diagnostics: pd.DataFrame,
) -> pd.DataFrame:
    posterior = samples[samples["iteration"] >= burn_in]
    rows = []
    for parameter in ("household_beta", "school_beta"):
        values = posterior[parameter].to_numpy()
        diagnostics = chain_diagnostics[
            chain_diagnostics["parameter"] == parameter
        ]
        r_hat = gelman_rubin(posterior, parameter)
        total_ess = float(diagnostics["ess"].sum())
        minimum_acceptance = diagnostics["acceptance_rate_post_burn_in"].min()
        maximum_acceptance = diagnostics["acceptance_rate_post_burn_in"].max()
        issues = []
        if not np.isfinite(r_hat) or r_hat >= 1.01:
            issues.append("high_r_hat")
        if total_ess < 400:
            issues.append("low_ess")
        if minimum_acceptance < 0.05:
            issues.append("low_acceptance")
        if maximum_acceptance > 0.50:
            issues.append("high_acceptance")
        rows.append(
            {
                "school_type": config.display_name,
                "parameter": parameter,
                "mean": np.mean(values),
                "std": np.std(values, ddof=1),
                "median": np.median(values),
                "lower_2.5%": np.percentile(values, 2.5),
                "upper_97.5%": np.percentile(values, 97.5),
                "r_hat": r_hat,
                "total_ess": total_ess,
                "minimum_chain_ess": diagnostics["ess"].min(),
                "mean_lag1_autocorrelation": diagnostics[
                    "lag1_autocorrelation"
                ].mean(),
                "epsilon": epsilon,
                "acceptance_rate": samples["accepted"].mean(),
                "minimum_chain_acceptance": minimum_acceptance,
                "maximum_chain_acceptance": maximum_acceptance,
                "diagnostic_status": "pass" if not issues else "check",
                "diagnostic_issues": ";".join(issues) if issues else "none",
            }
        )
    return pd.DataFrame(rows)


def estimate_school_type(
    config: SchoolConfig,
    args: argparse.Namespace,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    print(f"\nStarting {config.display_name} estimation")
    for stale_name in (
        f"{config.key}_samples.csv",
        f"{config.key}_chain_diagnostics.csv",
        f"{config.key}_diagnostics.png",
        f"{config.key}_pilot_diagnostics.png",
    ):
        (OUTPUT_DIR / stale_name).unlink(missing_ok=True)
    dataset = load_school_data(DATA_DIR / config.data_file)
    student_fraction = load_student_introduction_fraction(
        config.age_distribution_column
    )
    rng = np.random.default_rng(seed)

    epsilon, calibration = calibrate_epsilon(
        dataset,
        config,
        student_fraction,
        args.calibration_draws,
        args.epsilon_quantile,
        args.replicates,
        rng,
    )
    print(f"[{config.display_name}] calibrated epsilon={epsilon:.5f}")
    eligible = calibration[calibration["distance"] <= epsilon]
    if eligible.empty:
        eligible = calibration.nsmallest(1, "distance")

    selected_pilots = []
    all_pilot_attempts = []
    pilot_decisions = []
    chain_generators = {}
    chain_proposal_stds = {}
    for chain in range(1, args.chains + 1):
        proposal_std = args.proposal_std
        for attempt in range(1, args.max_pilot_attempts + 1):
            initial_index = (
                (chain - 1) * args.max_pilot_attempts + attempt - 1
            ) % len(eligible)
            initial = eligible.iloc[initial_index]
            attempt_seed = seed + 1_000_000 * chain + 10_000 * attempt
            chain_rng = np.random.default_rng(attempt_seed)
            attempt_samples = abc_mcmc_chain(
                    dataset=dataset,
                    config=config,
                    student_fraction=student_fraction,
                    epsilon=epsilon,
                    initial_parameters=np.array(
                        [initial["household_beta"], initial["school_beta"]]
                    ),
                    initial_distance=float(initial["distance"]),
                    iterations=args.pilot_iterations,
                    proposal_std=proposal_std,
                    replicates=args.replicates,
                    rng=chain_rng,
                    chain=chain,
                    stage=f"pilot attempt {attempt}",
                )
            passed, decision, reason = evaluate_pilot_acceptance(
                attempt_samples,
                args.pilot_burn_in,
                args.pilot_min_acceptance,
                args.pilot_max_acceptance,
            )
            attempt_samples["pilot_attempt"] = attempt
            attempt_samples["selected_for_formal"] = passed
            all_pilot_attempts.append(attempt_samples)
            decision["pilot_attempt"] = attempt
            decision["proposal_std"] = proposal_std
            decision["initial_household_beta"] = initial["household_beta"]
            decision["initial_school_beta"] = initial["school_beta"]
            decision["selected_for_formal"] = passed
            pilot_decisions.append(decision)
            acceptance = float(decision.iloc[0]["acceptance_rate"])
            print(
                f"[{config.display_name}] chain {chain}, attempt {attempt}: "
                f"acceptance={acceptance:.3f}, proposal_std={proposal_std:.5f}, "
                f"status={decision.iloc[0]['status']}"
            )
            if passed:
                selected_pilots.append(attempt_samples)
                chain_generators[chain] = chain_rng
                chain_proposal_stds[chain] = proposal_std
                break

            if reason.endswith("too_low"):
                proposal_std = max(
                    proposal_std / args.proposal_adjustment_factor, 1e-4
                )
            else:
                proposal_std = min(
                    proposal_std * args.proposal_adjustment_factor, 0.5
                )

    pilot_attempt_samples = pd.concat(all_pilot_attempts, ignore_index=True)
    pilot_decision = pd.concat(pilot_decisions, ignore_index=True)
    pilot_samples = (
        pd.concat(selected_pilots, ignore_index=True)
        if selected_pilots
        else pd.DataFrame()
    )
    accepted_chain_count = len(selected_pilots)
    pilot_target_met = accepted_chain_count == args.chains
    failed_chains = sorted(
        set(range(1, args.chains + 1)) - set(chain_generators)
    )
    pilot_reason = (
        "all_chains_acceptable"
        if pilot_target_met
        else (
            f"using_{accepted_chain_count}_of_{args.chains}_chains;"
            + "max_attempts_exceeded_for_chains_"
            + "_".join(map(str, failed_chains))
        )
    )
    pilot_attempt_samples.to_csv(
        OUTPUT_DIR / f"{config.key}_pilot_attempts.csv", index=False
    )
    pilot_samples.to_csv(
        OUTPUT_DIR / f"{config.key}_pilot_samples.csv", index=False
    )
    pilot_decision.to_csv(
        OUTPUT_DIR / f"{config.key}_pilot_decision.csv", index=False
    )
    calibration.to_csv(
        OUTPUT_DIR / f"{config.key}_epsilon_calibration.csv", index=False
    )
    if not pilot_samples.empty:
        pilot_diagnostics = compute_chain_diagnostics(
            pilot_samples, args.pilot_burn_in
        )
        pilot_diagnostics.to_csv(
            OUTPUT_DIR / f"{config.key}_pilot_chain_diagnostics.csv",
            index=False,
        )

    if accepted_chain_count == 0:
        print(
            f"[{config.display_name}] could not produce any acceptable chain "
            f"after {args.max_pilot_attempts} attempts per chain "
            f"({pilot_reason}); "
            "formal sampling was not started."
        )
        summary = pd.DataFrame(
            {
                "school_type": [config.display_name, config.display_name],
                "parameter": ["household_beta", "school_beta"],
                "epsilon": [epsilon, epsilon],
                "run_stage": ["pilot_stopped", "pilot_stopped"],
                "pilot_passed": [False, False],
                "pilot_target_met": [False, False],
                "accepted_chains": [
                    accepted_chain_count,
                    accepted_chain_count,
                ],
                "pilot_reason": [pilot_reason, pilot_reason],
            }
        )
        return pilot_samples, summary

    if pilot_target_met:
        print(
            f"[{config.display_name}] all {args.chains} chains passed; "
            f"continuing from {args.pilot_iterations} to "
            f"{args.iterations} total iterations."
        )
    else:
        print(
            f"[{config.display_name}] {accepted_chain_count} of "
            f"{args.chains} chains passed after all retries; continuing those "
            f"chains to {args.iterations} total iterations."
        )
    remaining_iterations = args.iterations - args.pilot_iterations
    continued_results = []
    for chain, pilot_chain in pilot_samples.groupby("chain", sort=True):
        last = pilot_chain.iloc[-1]
        continued_results.append(
            abc_mcmc_chain(
                dataset=dataset,
                config=config,
                student_fraction=student_fraction,
                epsilon=epsilon,
                initial_parameters=np.array(
                    [last["household_beta"], last["school_beta"]]
                ),
                initial_distance=float(last["distance"]),
                iterations=remaining_iterations,
                proposal_std=chain_proposal_stds[int(chain)],
                replicates=args.replicates,
                rng=chain_generators[int(chain)],
                chain=int(chain),
                iteration_offset=args.pilot_iterations,
                stage="formal continuation",
            )
        )

    samples = pd.concat(
        [pilot_samples, *continued_results], ignore_index=True
    ).sort_values(["chain", "iteration"], ignore_index=True)
    chain_diagnostics = compute_chain_diagnostics(samples, args.burn_in)
    summary = summarize_samples(
        samples,
        args.burn_in,
        config,
        epsilon,
        chain_diagnostics,
    )
    summary["run_stage"] = (
        "formal_completed"
        if pilot_target_met
        else "formal_completed_with_partial_chains"
    )
    summary["pilot_passed"] = True
    summary["pilot_target_met"] = pilot_target_met
    summary["accepted_chains"] = accepted_chain_count
    summary["pilot_reason"] = pilot_reason
    samples.to_csv(OUTPUT_DIR / f"{config.key}_samples.csv", index=False)
    chain_diagnostics.to_csv(
        OUTPUT_DIR / f"{config.key}_chain_diagnostics.csv", index=False
    )
    return samples, summary


# --------------------------------- Entrypoint ---------------------------------

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Estimate influenza transmission using ABC-MCMC."
    )
    parser.add_argument(
        "--schools",
        nargs="+",
        choices=tuple(SCHOOL_CONFIGS),
        default=list(SCHOOL_CONFIGS),
        help="School types to estimate (kindergarten is not included).",
    )
    parser.add_argument(
        "--iterations",
        type=int,
        default=20000,
        help="Total iterations per chain, including the pilot stage.",
    )
    parser.add_argument("--burn-in", type=int, default=4000)
    parser.add_argument("--pilot-iterations", type=int, default=2000)
    parser.add_argument("--pilot-burn-in", type=int, default=200)
    parser.add_argument("--pilot-min-acceptance", type=float, default=0.05)
    parser.add_argument("--pilot-max-acceptance", type=float, default=0.50)
    parser.add_argument("--max-pilot-attempts", type=int, default=20)
    parser.add_argument("--proposal-adjustment-factor", type=float, default=1.5)
    parser.add_argument("--chains", type=int, default=4)
    parser.add_argument("--calibration-draws", type=int, default=100)
    parser.add_argument("--epsilon-quantile", type=float, default=0.10)
    parser.add_argument("--proposal-std", type=float, default=0.05)
    parser.add_argument("--replicates", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20260825)
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Use very short settings to verify that the code runs.",
    )
    args = parser.parse_args()

    if args.quick:
        args.iterations = 20
        args.burn_in = 5
        args.pilot_iterations = 10
        args.pilot_burn_in = 2
        args.max_pilot_attempts = 3
        args.calibration_draws = 10
        args.replicates = 1
    if not 0 < args.epsilon_quantile < 1:
        parser.error("--epsilon-quantile must be between 0 and 1.")
    if not 0 <= args.burn_in < args.iterations:
        parser.error("--burn-in must be nonnegative and below --iterations.")
    if not 0 < args.pilot_iterations < args.iterations:
        parser.error("--pilot-iterations must be between 0 and --iterations.")
    if not 0 <= args.pilot_burn_in < args.pilot_iterations:
        parser.error(
            "--pilot-burn-in must be nonnegative and below --pilot-iterations."
        )
    if not (
        0
        <= args.pilot_min_acceptance
        < args.pilot_max_acceptance
        <= 1
    ):
        parser.error(
            "Pilot acceptance bounds must satisfy 0 <= minimum < maximum <= 1."
        )
    if args.chains < 1 or args.calibration_draws < 2 or args.replicates < 1:
        parser.error("chains/replicates must be >=1 and calibration draws >=2.")
    if args.max_pilot_attempts < 1:
        parser.error("--max-pilot-attempts must be at least 1.")
    if args.proposal_adjustment_factor <= 1:
        parser.error("--proposal-adjustment-factor must be greater than 1.")
    return args


def main(forced_schools: list[str] | None = None) -> None:
    args = parse_arguments()
    if forced_schools is not None:
        args.schools = forced_schools
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    summaries = []

    for index, school_key in enumerate(args.schools):
        _, summary = estimate_school_type(
            SCHOOL_CONFIGS[school_key],
            args,
            seed=args.seed + index * 1_000_000,
        )
        summaries.append(summary)

    all_summaries = pd.concat(summaries, ignore_index=True)
    summary_name = (
        f"abc_mcmc_summary_{forced_schools[0]}.xlsx"
        if forced_schools is not None and len(forced_schools) == 1
        else "abc_mcmc_summary.xlsx"
    )
    all_summaries.to_excel(OUTPUT_DIR / summary_name, index=False)
    print("\nEstimation complete")
    print(all_summaries.to_string(index=False))
    print(f"\nResults saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
