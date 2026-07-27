# Continuous Ground-Truth Diagnostics Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make topomap and amplitude diagnostics use the exact saved population-level continuous effect and reject legacy truth objects that cannot support scientifically valid recovery metrics.

**Architecture:** Add two focused accessors to the shared R utilities: one validates and selects exact channel-time truth cells, and one computes peak and discrete-window truth summaries. Both diagnostic scripts consume these accessors, so neither reconstructs a waveform from hard-label onset/offset boundaries.

**Tech Stack:** R 4.4, base R assertions, dplyr/tidyr, existing brms diagnostic scripts.

## Global Constraints

- `ground_truth$continuous_effect` is the sole source of continuous diagnostic truth.
- Missing, malformed, incomplete, or duplicated truth cells cause an explicit error.
- There is no legacy Gaussian fallback.
- Existing unrelated working-tree changes must not be staged or overwritten.
- No new R package dependency is introduced.

---

### Task 1: Exact continuous-truth accessors

**Files:**
- Create: `tests/test_continuous_ground_truth.R`
- Modify: `code/utils.R`

**Interfaces:**
- Produces: `continuous_truth_cells(ground_truth, channels = NULL, times = NULL) -> data.frame(channel, time, effect)`
- Produces: `continuous_truth_summary(ground_truth, channel, peak_time, window_times) -> list(peak, window_mean)`
- Consumes: `ground_truth$continuous_effect` with unique `channel`, `time`, and `effect` columns.

- [ ] **Step 1: Write the failing accessor tests**

Create a base-R test script that sources `code/utils.R`, constructs a small truth grid, and asserts exact selection, exact discrete summaries, missing-truth failure, missing-cell failure, and duplicate-cell failure:

```r
source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

expect_error <- function(expr, pattern) {
    message_text <- tryCatch({
        force(expr)
        NA_character_
    }, error = function(e) conditionMessage(e))
    stopifnot(!is.na(message_text), grepl(pattern, message_text, fixed = TRUE))
}

truth <- list(
    true_onset = 0.20,
    true_offset = 0.30,
    continuous_effect = expand.grid(
        channel = c("CZ", "FZ"),
        time = c(0.16, 0.25, 0.34),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
)
truth$continuous_effect$effect <- c(0.121, 0.040, 1.457, 0.300, 0.121, 0.040)

selected <- continuous_truth_cells(truth, "CZ", c(0.16, 0.25, 0.34))
stopifnot(identical(selected$effect, c(0.121, 1.457, 0.121)))

summary <- continuous_truth_summary(truth, "CZ", 0.25, c(0.16, 0.25, 0.34))
stopifnot(isTRUE(all.equal(summary$peak, 1.457)))
stopifnot(isTRUE(all.equal(summary$window_mean, mean(c(0.121, 1.457, 0.121)))))

expect_error(continuous_truth_cells(list(), "CZ", 0.25),
             "Exact continuous ground truth is unavailable")

incomplete <- truth
incomplete$continuous_effect <- subset(
    incomplete$continuous_effect,
    !(channel == "CZ" & time == 0.25)
)
expect_error(continuous_truth_cells(incomplete, "CZ", 0.25),
             "Missing continuous ground-truth cells")

duplicated_truth <- truth
duplicated_truth$continuous_effect <- rbind(
    duplicated_truth$continuous_effect,
    duplicated_truth$continuous_effect[1, ]
)
expect_error(continuous_truth_cells(duplicated_truth, "CZ", 0.16),
             "Duplicated continuous ground-truth cells")

cat("continuous ground-truth tests passed\n")
```

- [ ] **Step 2: Run the test and verify RED**

Run:

```bash
PROJECT_ROOT="$PWD" Rscript --vanilla tests/test_continuous_ground_truth.R
```

Expected: FAIL because `continuous_truth_cells()` does not exist.

- [ ] **Step 3: Implement the minimal accessors**

Add to `code/utils.R`:

```r
continuous_truth_cells <- function(ground_truth, channels = NULL, times = NULL) {
    truth <- ground_truth$continuous_effect
    if (is.null(truth)) {
        stop(
            "Exact continuous ground truth is unavailable; regenerate the dataset ",
            "with the current data generator.",
            call. = FALSE
        )
    }

    required <- c("channel", "time", "effect")
    missing_columns <- setdiff(required, names(truth))
    if (length(missing_columns) > 0L) {
        stop("continuous_effect is missing required columns: ",
             paste(missing_columns, collapse = ", "), call. = FALSE)
    }

    truth <- truth[, required, drop = FALSE]
    if (!is.null(channels)) truth <- truth[truth$channel %in% channels, , drop = FALSE]
    if (!is.null(times)) truth <- truth[truth$time %in% times, , drop = FALSE]

    duplicated_cells <- duplicated(truth[c("channel", "time")])
    if (any(duplicated_cells)) {
        stop("Duplicated continuous ground-truth cells", call. = FALSE)
    }

    requested_channels <- if (is.null(channels)) unique(truth$channel) else channels
    requested_times <- if (is.null(times)) unique(truth$time) else times
    requested <- expand.grid(
        channel = requested_channels,
        time = requested_times,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    matched <- merge(requested, truth, by = c("channel", "time"), all.x = TRUE,
                     sort = FALSE)
    if (anyNA(matched$effect)) {
        stop("Missing continuous ground-truth cells", call. = FALSE)
    }

    matched[order(match(matched$channel, requested_channels),
                  match(matched$time, requested_times)), required, drop = FALSE]
}

continuous_truth_summary <- function(ground_truth, channel, peak_time,
                                     window_times) {
    cells <- continuous_truth_cells(
        ground_truth,
        channels = channel,
        times = unique(c(peak_time, window_times))
    )
    list(
        peak = cells$effect[cells$time == peak_time],
        window_mean = mean(cells$effect[cells$time %in% window_times])
    )
}
```

- [ ] **Step 4: Run the accessor test and verify GREEN**

Run:

```bash
PROJECT_ROOT="$PWD" Rscript --vanilla tests/test_continuous_ground_truth.R
```

Expected: PASS and print `continuous ground-truth tests passed`.

---

### Task 2: Make both diagnostics consume exact truth

**Files:**
- Modify: `tests/test_continuous_ground_truth.R`
- Modify: `code/08_slide_topomap_n20.R`
- Modify: `code/02_diagnostics.R`

**Interfaces:**
- Consumes: `continuous_truth_cells()` from Task 1 for topomap truth layers.
- Consumes: `continuous_truth_summary()` from Task 1 for peak and window recovery.
- Produces: unchanged diagnostic CSV schemas, but `truth_peak_uV`, `truth_mean_uV`, residuals, and recovery percentages use saved exact truth.

- [ ] **Step 1: Add failing integration guards**

Append static regression guards to `tests/test_continuous_ground_truth.R`:

```r
topomap_source <- paste(readLines("code/08_slide_topomap_n20.R"), collapse = "\n")
diagnostic_source <- paste(readLines("code/02_diagnostics.R"), collapse = "\n")

stopifnot(grepl("continuous_truth_cells", topomap_source, fixed = TRUE))
stopifnot(!grepl("temporal_weight <- function", topomap_source, fixed = TRUE))
stopifnot(grepl("continuous_truth_summary", diagnostic_source, fixed = TRUE))
stopifnot(!grepl("truth_window_mean <- function", diagnostic_source, fixed = TRUE))
```

- [ ] **Step 2: Run the test and verify RED**

Run:

```bash
PROJECT_ROOT="$PWD" Rscript --vanilla tests/test_continuous_ground_truth.R
```

Expected: FAIL because both scripts still reconstruct truth.

- [ ] **Step 3: Replace topomap truth reconstruction**

In `code/08_slide_topomap_n20.R`, remove `truth_mode`, `seed_xy`,
`temporal_weight()`, and the reconstructed spatial/temporal weights. Build the
truth layer as:

```r
gt_elec <- continuous_truth_cells(
    gt,
    channels = sensors$channel,
    times = target_times
) %>%
    rename(value = effect) %>%
    left_join(sensors, by = "channel") %>%
    mutate(layer = "0. Ground truth")
```

- [ ] **Step 4: Replace amplitude truth reconstruction**

In `code/02_diagnostics.R`:

- remove `truth_window_mean()`;
- remove the `truth_peak` argument from `amplitude_diagnostics()` and
  `run_all_diagnostics()`;
- compute `peak_time` from the same analysis grid as the spatial diagnostic;
- call `continuous_truth_summary(ground_truth, seed_ch, peak_time,
  window_times)`;
- pass `truth$peak` and `truth$window_mean` to `recovery_summary()`;
- write those exact values into `truth_peak_uV` and `truth_mean_uV`;
- remove the final assignment `truth_peak <- ground_truth$bump_amplitude`.

Use:

```r
peak_time <- pick_time((true_on + true_off) / 2, times)
window_times <- times[times >= true_on & times <= true_off]
truth <- continuous_truth_summary(
    ground_truth,
    channel = seed_ch,
    peak_time = peak_time,
    window_times = window_times
)
peak <- recovery_summary(cz_space$mean_diff, truth$peak)
window <- recovery_summary(window_draws$mean_in_window, truth$window_mean)
```

- [ ] **Step 5: Run tests and parse checks**

Run:

```bash
PROJECT_ROOT="$PWD" Rscript --vanilla tests/test_continuous_ground_truth.R
Rscript --vanilla -e 'parse(file="code/utils.R"); parse(file="code/02_diagnostics.R"); parse(file="code/08_slide_topomap_n20.R"); cat("parse checks passed\n")'
```

Expected: all tests pass and all three files parse.

- [ ] **Step 6: Verify the real high-heterogeneity truth object**

On the SSH project, after synchronising these code changes, run the accessor on
`data/montecarlo_p20_high_heterogeneity/gt_001.rds` and verify that CZ truth at
0.16, 0.25, and 0.34 seconds is approximately 0.1211, 1.4574, and 0.1211.

Expected: values come directly from `continuous_effect`; no use of
`true_onset=0.20` or `true_offset=0.30` occurs in truth construction.

- [ ] **Step 7: Review the scoped diff**

Run:

```bash
git diff --check
git diff -- code/utils.R code/02_diagnostics.R code/08_slide_topomap_n20.R tests/test_continuous_ground_truth.R
```

Expected: no whitespace errors and no unrelated edits in the scoped diff.
