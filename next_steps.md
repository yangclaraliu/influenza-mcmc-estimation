# Next steps

The simulator and the ABC-MCMC code are in good shape mechanically. The
steps below aim to make the fit trustworthy with as few changes as possible.
They are ordered so that each step informs the next. Please work through
them in order and bring the outputs to our meetings.

## Why the trace plots look like this

The trace plots are most likely showing an **identifiability problem**,
not a sampler that needs more iterations. A parameter is identifiable when
different values of it produce noticeably different data. If many
combinations of `household_beta` and `school_beta` produce the same
simulated attack rates, the data cannot tell them apart. The posterior is
then a long, thin valley rather than a single peak, and the chains drift
slowly along it. Running the chains longer or retuning the step size does
not fix this.

### Why it arises in this model

1. **Both parameters do the same job.** Each one raises the attack rate.
   The data are school-by-grade attack rates, so what they mainly pin down
   is the total amount of transmission, not how it splits between household
   and school.
2. **The two routes cover different ranges.** The community incidence curve
   implies only about 46 introductions per season into a 720-student school
   and its families. Even with `household_beta` = 1, the household route
   gives attack rates of only about 5%. School transmission adds little
   until `school_beta` reaches about 0.2, then climbs steeply, giving
   outbreaks of 25–35% by 0.4.
   Moderate attack rates can be reached by many trade-offs between the two
   parameters.
3. **Little in the simulation differs between schools.** Every simulated
   school has the same size, network structure and introduction curve. Only
   the vaccination rates differ. The model therefore has few features of
   the data that respond differently to the two parameters.
4. **The comparison may not use all the information in the data.** Step 3
   below is designed to check this.

### An illustration

Mean simulated attack rate for a synthetic school (5 grades, 30%
vaccination coverage, 5% of community introductions reaching students
directly, 30 replicates). These are illustrative values, not results from
the study data.

| `household_beta` \ `school_beta` | 0 | 0.05 | 0.1 | 0.2 | 0.3 | 0.4 |
|---|---|---|---|---|---|---|
| 0 | 0.003 | 0.003 | 0.003 | 0.006 | 0.013 | 0.075 |
| 0.3 | 0.015 | 0.018 | 0.024 | 0.038 | 0.089 | 0.239 |
| 0.6 | 0.031 | 0.034 | 0.039 | 0.070 | 0.137 | 0.281 |
| 1.0 | 0.046 | 0.053 | 0.064 | 0.089 | 0.188 | 0.332 |

An attack rate of about 4% comes from (1.0, 0), (0.6, 0.1) and (0.3, 0.2).
An attack rate of about 9% comes from both (1.0, 0.2) and (0.3, 0.3). Along
each of these lines the data look much the same, which is the valley the
chains drift along. At higher attack rates the fit depends on `school_beta`
in the steep region, where small changes have large effects, which makes
the chains sticky.

### What this means for the next steps

The aim is not to make the chains look better, but to find out whether the
data can separate the two parameters and, if not, what extra information
could. The steps below are ordered to answer that.

## 1. Review the two new commits

- **Direct-introduction fix** (`1efdfd4`). Direct community introductions
  were drawn only from susceptible students, so vaccination coverage had no
  effect on that route. They now target any student, and only susceptible
  students become infected. Please check that you agree with the change.
- **`meta_parameters.md`** (`ae04c69`). Your explanation of the contact
  probabilities, plus sections for you to complete.

## 2. Complete `meta_parameters.md`

Complete every section marked **To be completed** (§2–§8), and confirm the
open points in §1. Each section describes what the code currently does and
lists the questions to answer. The most important for the fit are:

- **§4 Ascertainment.** Is the `flu` column true infection, or reported or
  diagnosed influenza?
- **§5 Community incidence.** Is a 1.6% community attack rate plausible?

If any value cannot be sourced, note it. It then becomes a candidate for a
sensitivity check.

## 3. Diagnostic run (before changing the fitting code)

Run the simulator on a grid of parameter values with your data:
`household_beta` from 0 to 1 in steps of 0.1, and `school_beta` from 0 to
0.3 in steps of 0.02. Produce two plots:

1. **Distance heat map** over the two parameters. Is there a single clear
   region of good fit, or a long valley?
2. **Individual simulated schools against the observed schools**, at a few
   grid points. Use single realisations, not averages over replicates.
   How does the spread of simulated attack rates across schools compare
   with the observed spread?

Questions to think about:

- What does the heat map say about whether both parameters can be
  estimated from these data?
- How do individual simulated schools compare with the 20-replicate
  average that the current distance uses?

Bring both plots to our next meeting. We will decide the changes to the
fitting code together, based on what they show.

## 4. Sensitivity check: fixed versus changing contacts

The current network keeps the same contacts for the whole season. As a
robustness check, add an option in which each infectious student's contacts
are drawn fresh on the day they are infectious, using the same
probabilities. Because each student is infectious for only one day, this is
cheap to implement in `run_sir`. Fit both versions and compare the
estimates.

If the results differ noticeably, measure how persistent contacts are in
the SocioPatterns data: threshold day 1 and day 2 separately, and record
what fraction of each student's day-1 contacts reappear on day 2. Use that
fraction to build an intermediate version, with persistent ties for that
share of contacts and fresh draws for the rest.

## Later

Once the fit works, test the sensitivity of the estimates to the other
assumptions recorded in `meta_parameters.md`. These include the infectious
period, the ascertainment factors, and susceptibility.
