USE ProteinWSME;
GO

-- 1. All imported runs.
SELECT
    run_id,
    model_variant,
    c_value,
    min_seq_sep,
    maxfev,
    protein_count,
    source_file
FROM dbo.runs
ORDER BY model_variant, c_value, min_seq_sep;
GO

-- 2. One protein across all runs.
SELECT
    protein_id,
    model_variant,
    c_value,
    min_seq_sep,
    target_dg,
    predicted_dg,
    error,
    sq_error,
    status,
    barrier_warning
FROM dbo.protein_run_view
WHERE protein_id = N'1afk_A'
ORDER BY model_variant, c_value, min_seq_sep;
GO

-- 3. Best run for each protein by minimal squared error.
WITH ranked AS (
    SELECT
        *,
        ROW_NUMBER() OVER (
            PARTITION BY protein_id
            ORDER BY sq_error ASC, ABS(error) ASC
        ) AS rn
    FROM dbo.protein_run_view
    WHERE status = N'ok'
)
SELECT
    protein_id,
    model_variant,
    c_value,
    min_seq_sep,
    sq_error,
    predicted_dg,
    target_dg
FROM ranked
WHERE rn = 1
ORDER BY sq_error ASC;
GO

-- 4. Mean squared error grouped by run parameters.
SELECT
    model_variant,
    c_value,
    min_seq_sep,
    COUNT(*) AS protein_rows,
    AVG(sq_error) AS avg_sq_error,
    AVG(ABS(error)) AS avg_abs_error
FROM dbo.protein_run_view
WHERE status = N'ok'
GROUP BY model_variant, c_value, min_seq_sep
ORDER BY avg_sq_error ASC;
GO

-- 5. Proteins with warnings or failures.
SELECT
    protein_id,
    model_variant,
    c_value,
    min_seq_sep,
    status,
    barrier_warning,
    failure_reason
FROM dbo.protein_run_view
WHERE barrier_warning IS NOT NULL
   OR status <> N'ok'
   OR failure_reason IS NOT NULL
ORDER BY protein_id, model_variant, c_value, min_seq_sep;
GO
