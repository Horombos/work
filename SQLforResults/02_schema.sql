USE ProteinWSME;
GO

IF OBJECT_ID(N'dbo.protein_run_view', N'V') IS NOT NULL
BEGIN
    DROP VIEW dbo.protein_run_view;
END;
GO

IF OBJECT_ID(N'dbo.protein_results', N'U') IS NOT NULL
BEGIN
    DROP TABLE dbo.protein_results;
END;
GO

IF OBJECT_ID(N'dbo.runs', N'U') IS NOT NULL
BEGIN
    DROP TABLE dbo.runs;
END;
GO

CREATE TABLE dbo.runs (
    run_id INT IDENTITY(1,1) NOT NULL CONSTRAINT PK_runs PRIMARY KEY,
    source_file NVARCHAR(260) NOT NULL CONSTRAINT UQ_runs_source_file UNIQUE,
    model_variant NVARCHAR(50) NOT NULL,
    c_value INT NOT NULL,
    min_seq_sep INT NOT NULL,
    maxfev INT NULL,
    protein_count INT NOT NULL CONSTRAINT DF_runs_protein_count DEFAULT (0),
    imported_at DATETIME2(0) NOT NULL CONSTRAINT DF_runs_imported_at DEFAULT (SYSDATETIME())
);
GO

CREATE TABLE dbo.protein_results (
    result_id INT IDENTITY(1,1) NOT NULL CONSTRAINT PK_protein_results PRIMARY KEY,
    run_id INT NOT NULL,
    protein_id NVARCHAR(50) NOT NULL,
    target_dg FLOAT NULL,
    predicted_dg FLOAT NULL,
    error FLOAT NULL,
    sq_error FLOAT NULL,
    e_activation FLOAT NULL,
    q_unf INT NULL,
    q_ts INT NULL,
    q_fold INT NULL,
    barrier_warning NVARCHAR(100) NULL,
    status NVARCHAR(50) NULL,
    failure_reason NVARCHAR(400) NULL,
    CONSTRAINT FK_protein_results_runs
        FOREIGN KEY (run_id) REFERENCES dbo.runs (run_id),
    CONSTRAINT UQ_protein_results_run_protein UNIQUE (run_id, protein_id)
);
GO

CREATE INDEX IX_protein_results_protein_id
    ON dbo.protein_results (protein_id);
GO

CREATE INDEX IX_protein_results_run_id
    ON dbo.protein_results (run_id);
GO

CREATE INDEX IX_runs_parameters
    ON dbo.runs (c_value, min_seq_sep, model_variant);
GO

CREATE VIEW dbo.protein_run_view
AS
SELECT
    pr.protein_id,
    pr.target_dg,
    pr.predicted_dg,
    pr.error,
    pr.sq_error,
    pr.e_activation,
    pr.q_unf,
    pr.q_ts,
    pr.q_fold,
    pr.barrier_warning,
    pr.status,
    pr.failure_reason,
    r.run_id,
    r.model_variant,
    r.c_value,
    r.min_seq_sep,
    r.maxfev,
    r.source_file,
    r.imported_at
FROM dbo.protein_results AS pr
INNER JOIN dbo.runs AS r
    ON r.run_id = pr.run_id;
GO
