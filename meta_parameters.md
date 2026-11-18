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
