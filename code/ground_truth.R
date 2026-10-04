# Population-mean effect waveform A * fbar(t) on the analysis grid (thesis Eq. 19).
population_effect_waveform <- function(bump_amplitude, onset_support_ms,
                                       downsample_factor, template_path) {
    if (!file.exists(template_path)) {
        stop("ERP template does not exist: ", template_path, call. = FALSE)
    }
    if (length(bump_amplitude) != 1L || !is.finite(bump_amplitude)) {
        stop("bump_amplitude must be one finite number", call. = FALSE)
    }
    if (!length(onset_support_ms) || any(!is.finite(onset_support_ms))) {
        stop("onset_support_ms must contain finite values", call. = FALSE)
    }
    if (length(downsample_factor) != 1L ||
        !is.finite(downsample_factor) || downsample_factor < 1L ||
        downsample_factor != as.integer(downsample_factor)) {
        stop("downsample_factor must be one positive integer", call. = FALSE)
    }

    template_environment <- new.env(parent = globalenv())
    sys.source(template_path, envir = template_environment)
    template_fields <- c("Xf", "Nf", "erp")
    missing_template_fields <- template_fields[
        !vapply(
            template_fields,
            exists,
            logical(1),
            envir = template_environment,
            inherits = FALSE
        )
    ]
    if (length(missing_template_fields)) {
        stop(
            "ERP template is missing required objects: ",
            paste(missing_template_fields, collapse = ", "),
            call. = FALSE
        )
    }
    sample_times_ms <- template_environment$Xf
    n_samples <- template_environment$Nf
    erp <- template_environment$erp
    if (length(sample_times_ms) != n_samples) {
        stop("ERP template Xf and Nf are inconsistent", call. = FALSE)
    }

    shifted_waveform <- function(onset_ms) {
        onset_index <- match(onset_ms, sample_times_ms)
        if (is.na(onset_index)) {
            stop(
                "onset_support_ms contains a value outside the ERP grid: ",
                onset_ms,
                call. = FALSE
            )
        }
        trailing_length <- n_samples - onset_index - length(erp) + 2L
        if (onset_index < 2L || trailing_length < 0L) {
            stop(
                "ERP template does not fit at onset ", onset_ms, " ms",
                call. = FALSE
            )
        }
        waveform <- c(
            rep(0, onset_index - 2L),
            erp,
            rep(0, trailing_length)
        )
        stopifnot(length(waveform) == n_samples)
        waveform
    }

    onset_waveforms <- vapply(
        onset_support_ms,
        shifted_waveform,
        numeric(n_samples)
    )
    expected_temporal <- bump_amplitude * rowMeans(onset_waveforms)
    analysis_index <- seq.int(
        from = 1L,
        to = n_samples,
        by = as.integer(downsample_factor)
    )
    analysis_times <- sample_times_ms[analysis_index] / 1000
    expected_temporal <- expected_temporal[analysis_index]

    data.frame(time = analysis_times, effect_waveform = expected_temporal)
}

paired_radius_superpopulation_truth <- function(
        ground_truth, template_path, channels = NULL, times = NULL) {
    required_fields <- c(
        "bump_amplitude", "onset_support_ms", "downsample_factor",
        "continuous_effect"
    )
    missing_fields <- setdiff(required_fields, names(ground_truth))
    if (length(missing_fields)) {
        stop(
            "ground_truth is missing required fields: ",
            paste(missing_fields, collapse = ", "),
            call. = FALSE
        )
    }
    waveform <- population_effect_waveform(
        ground_truth$bump_amplitude,
        ground_truth$onset_support_ms,
        ground_truth$downsample_factor,
        template_path
    )
    analysis_times <- waveform$time
    expected_temporal <- waveform$effect_waveform

    continuous_effect <- ground_truth$continuous_effect
    required_columns <- c("channel", "spatial_weight")
    missing_columns <- setdiff(required_columns, names(continuous_effect))
    if (length(missing_columns)) {
        stop(
            "continuous_effect is missing required columns: ",
            paste(missing_columns, collapse = ", "),
            call. = FALSE
        )
    }
    spatial <- unique(continuous_effect[required_columns])
    counts <- table(spatial$channel)
    if (any(counts != 1L)) {
        stop(
            "continuous_effect must contain one unique spatial_weight per channel",
            call. = FALSE
        )
    }
    if (any(!is.finite(spatial$spatial_weight))) {
        stop("spatial_weight must be finite", call. = FALSE)
    }

    population <- expand.grid(
        channel = spatial$channel,
        time = analysis_times,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    population$effect <-
        spatial$spatial_weight[match(population$channel, spatial$channel)] *
        expected_temporal[match(population$time, analysis_times)]

    if (!is.null(channels)) {
        unknown_channels <- setdiff(channels, spatial$channel)
        if (length(unknown_channels)) {
            stop(
                "Requested channels are absent from ground truth: ",
                paste(unknown_channels, collapse = ", "),
                call. = FALSE
            )
        }
        population <- population[population$channel %in% channels, , drop = FALSE]
    }
    if (!is.null(times)) {
        unknown_times <- setdiff(times, analysis_times)
        if (length(unknown_times)) {
            stop(
                "Requested times are absent from the analysis grid: ",
                paste(unknown_times, collapse = ", "),
                call. = FALSE
            )
        }
        population <- population[population$time %in% times, , drop = FALSE]
    }
    rownames(population) <- NULL
    population
}
