# Meta-parameters

This document records every fixed (non-estimated) parameter used in
`MCMC estimation(all schools).py`: its value, its meaning, its source, and
how it enters the code. Only `household_beta` and `school_beta` are estimated
by ABC-MCMC; everything below is held fixed, so each value needs a clear
justification.

| Parameter | Value | Code location | Status |
|---|---|---|---|
| School contact probabilities | see §1 | `SchoolConfig.contact_probabilities` | Documented |
| Susceptibility | 0.88 | `SchoolConfig.susceptibility` | **To be completed** (§2) |
| Vaccine efficacy | 0.60 | `VACCINE_EFFICACY` | **To be completed** (§3) |
| Ascertainment | 0.354 × 0.56 | `UNDERESTIMATE` | **To be completed** (§4) |
| Community incidence | piecewise curve; 1,886,100 | `community_lambda`, `POPULATION` | **To be completed** (§5) |
| Student share of community cases | age-distribution workbook | `load_student_introduction_fraction` | **To be completed** (§6) |
| Infectious period | 1 day | `run_sir` | **To be completed** (§7) |
| Static school network | same contacts all 148 days | `simulate_school_type` | **To be completed** (§8) |

---

## 1. School contact network

### Mean number of school contacts

The average numbers of school contacts for primary, middle, and high school
students (4.37669, 6.64641, and 6.10498 respectively) were taken from the
age-specific, school-setting synthetic contact matrix for China reported by
Prem et al. (2021). These values set the overall density, or expected mean
degree, of each simulated school network.

### Relative structure: within class, within grade, between grades

The relative probabilities of contact within a class, between classes in the
same grade, and between grades were derived from the SocioPatterns primary
school data (Stehlé et al., 2011). The 10-fold relationship is not stated in
the paper. It was derived from the 20-second-resolution RFID contact data, as
follows:

1. Contacts for each student pair were summed over the two observation days.
2. Pairs with at least 18 min 20 s of cumulative contact were kept as edges.
   At this threshold the empirical network has a mean degree of 4.37069,
   which closely matches the primary-school target of 4.37669 from
   Prem et al.
3. The resulting edge densities were:

   | Pair type | Edge density |
   |---|---|
   | Same class | 0.165957 |
   | Same grade, different class | 0.015985 |
   | Different grade | 0.001626 |

   The two successive ratios are about 10.38 and 9.83. Rounding both to 10
   gives the simplified structure
   **within class : within grade : between grades = 100 : 10 : 1**.

### Absolute probabilities

The absolute probabilities were then scaled so that each school type keeps
its mean number of contacts from Prem et al.:

| School type | Same class | Same grade | Different grade |
|---|---|---|---|
| Primary | 0.102 | 0.0102 | 0.00102 |
| Middle | 0.160 | 0.016 | 0.0016 |
| High | 0.151 | 0.0151 | 0.00151 |

The Stehlé et al. data directly support the relative structure for a primary
school. Applying the same 100 : 10 : 1 structure to middle and high schools
is a modelling assumption. Their overall network densities are still set by
the corresponding contact numbers from Prem et al.

### Open points

- **Grades per school type.** With 4 classes of 30 students per grade, these
  probabilities reproduce the target mean degrees for 5 primary grades,
  4 middle grades, and 3 high-school grades (the 5-4-3 system). The code
  takes the number of grades from the highest grade in each workbook. Please
  confirm that the workbooks contain 5 and 4 grades; with 6 primary grades
  the network is about 2.5% denser than intended.
- **Teachers and staff.** The Prem et al. school-setting matrix counts
  contacts with all ages, including teachers and staff. Please state whether
  the target values use only the student-age part of the matrix.
- **School layout.** The SocioPatterns school had about 2 classes of ~23
  students per grade, while the model uses 4 classes of 30. Transferring the
  density ratios assumes that pairwise contact propensity does not depend on
  school size.
- **Static network.** The same network is used on every one of the 148 days,
  whereas Prem et al. report average daily contacts.

### References

- Prem K, et al. (2021). Projecting contact matrices in 177 geographical
  regions: An update and comparison with empirical data for the COVID-19
  era. *PLOS Computational Biology*.
  <https://doi.org/10.1371/journal.pcbi.1009098>
- Stehlé J, et al. (2011). High-resolution measurements of face-to-face
  contact patterns in a primary school. *PLOS ONE*.
  <https://doi.org/10.1371/journal.pone.0023176>
- SocioPatterns primary school temporal network data:
  <https://sociopatterns.org/datasets/primary-school-temporal-network-data/>

---

## 2. Susceptibility (0.88)

**How the code uses it.** In `generate_susceptibility`, each student is
independently marked as infectable with probability 0.88. The other 12% are
fully immune for the whole season, by every route (direct community,
household, and school). The same value is used for all three school types.

**To be completed:**

- **What the value represents.** Is 0.88 the proportion of students without
  prior immunity (for example, from serology)? Or is it a relative
  susceptibility of children compared with adults? If it is a relative
  susceptibility, it should multiply transmission probabilities rather than
  remove 12% of students, and it cannot then be separated from the estimated
  betas.
- **Source.** Full reference, with the table or figure and the age group the
  value refers to.
- **Age groups.** Why is the same value used for primary, middle, and high
  school students?
- **Relationship to vaccination.** Is the source value measured in a
  population that already includes vaccinated people? If so, applying it
  independently of vaccination may count immunity twice.

---

## 3. Vaccine efficacy (0.60)

**How the code uses it.** Each student is marked as vaccinated at random,
with probability equal to their grade's observed vaccination rate. A
vaccinated student is then fully protected for all 148 days with
probability 0.6 and otherwise has no protection at all (an "all-or-nothing"
vaccine).

**To be completed:**

- **Source.** Full reference, including the season and the vaccine–strain
  match.
- **What it measures.** Protection against what outcome: any infection,
  symptomatic illness, or medically attended, lab-confirmed influenza?
- **All-or-nothing versus leaky.** Why model full protection for 60% of
  vaccinated students rather than a 60% reduction in the chance of infection
  per exposure? Because students meet the same contacts every day, the two
  choices give noticeably different results.
- **Timing.** When were students vaccinated relative to the 148-day study
  period? The model assumes protection from day 1, with no delay before
  protection develops and no waning.
- **Check against the data.** Using the student-level records, report the
  influenza rate among vaccinated and unvaccinated students, overall and by
  school type, along with the implied efficacy (1 − rate ratio). Also note
  whether the vaccination field could include students vaccinated after
  their illness.

---

## 4. Ascertainment (`UNDERESTIMATE` = 0.354 × 0.56)

**How the code uses it.** Reported community cases from `community_lambda`
are divided by 0.198 (= 0.354 × 0.56), which multiplies them by about 5, to
give true community infections. Simulated student infections are all true
infections. They are compared directly with the `flu` column in the school
workbooks.

**To be completed:**

- **Source.** Full references for 0.354 and 0.56, and what each factor
  represents (for example, the proportion seeking care, or the proportion
  tested).
- **Observed outcome.** Is the `flu` column true infection, or reported or
  diagnosed influenza? If it is reported influenza, the model compares
  simulated true infections with reported cases, and an ascertainment
  factor is also needed on the student side.

---

## 5. Community incidence (`community_lambda`, `POPULATION` = 1,886,100)

**How the code uses it.** A piecewise polynomial of daily reported cases
over 148 days, scaled to the school-plus-family population, drives
community introductions. After the ascertainment correction, it implies
about a **1.6% community attack rate** over the season. For a 720-student
school and its families, that is about 46 introductions per season, so the
household route alone can produce attack rates of only about 5%.

**To be completed:**

- **Source.** Which surveillance data, which area, and which season? How
  was the curve fitted?
- **Plausibility.** Is a 1.6% true-infection attack rate plausible for this
  season?
- **Continuity.** The fitted curve jumps at day 71 (11.7 → 14.3) and at
  day 106 (60.4 → 47.1). Can it be refitted to be continuous?

---

## 6. Student share of community cases (age-distribution workbook)

**How the code uses it.** For each community introduction, the daily value
for the school type gives the probability that it infects a student
directly. Otherwise it reaches a family member.

**To be completed:**

- **Source and definition.** How was the column derived? Is it the share of
  all city cases in that age group?
- **Population.** City-wide, students are a small share of the population.
  In the simulated school-plus-family population they are about 25%. Is
  applying the city-wide share appropriate here?
- **Double counting.** City-wide student cases already include infections
  caught at school. How much of the school transmission is therefore also
  counted as community introduction?

---

## 7. Infectious period (one day)

**How the code uses it.** An infected student is infectious on the day
after infection only, then recovers. There is no latent period, so each
school contact gets exactly one chance of being infected.

**To be completed:**

- **Meaning.** Does one day represent the biological infectious period, or
  the time a student spends at school while infectious before staying home
  sick? Influenza is typically infectious for 3–5 days.
- **Source.** If it is time at school before staying home, what supports
  it (for example, school sick-leave policy or absenteeism data)?
- **Interpretation of `school_beta`.** With a one-day period, `school_beta`
  is the probability of transmission per contact over the whole in-school
  infectious period. It is not a daily probability, so values from the
  literature are not directly comparable.

---

## 8. Static school network

**How the code uses it.** One network is generated per simulated school and
used on all 148 days, so each student has the same contacts all season.

**To be completed:**

- **Per day versus two days.** The Prem et al. targets are contacts per
  day. The SocioPatterns threshold was applied to contact summed over two
  days, so a pair that met only on one day counts as a permanent tie. Please
  comment on this mismatch.
- **Persistence of contacts.** In the SocioPatterns data, threshold day 1
  and day 2 separately. What fraction of each student's day-1 contacts
  reappear on day 2?
- **Classroom exposure.** The network includes close contacts only (about
  3 of 29 classmates). Influenza also spreads through shared classroom air.
  Is ignoring classroom-wide exposure reasonable?
