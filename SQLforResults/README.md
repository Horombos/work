SQL Server/SSMS workflow

This folder contains a SQL Server version of the protein-results database for use in SQL.

Files

01_create_database.sql creates the ProteinWSME database if it does not already exist
02_schema.sql creates tables, indexes, and the dbo.protein_run_view
03_example_queries.sql contains ready-to-run analysis queries
generate_seed_sql.py reads the CSV files and generates 04_seed_data.sql
04_seed_data.sql is the generated data-load script with plain INSERT statements

Steps in SSMS

1. run 01_create_database.sql`
2. run 02_schema.sql`
3. run 04_seed_data.sql` (make it from .py and given csv files)
4. run 03_example_queries.sql`

Main tables

dbo.runs`: one row per CSV / experiment
dbo.protein_results`: one row per protein result within a run
dbo.protein_run_view`: joined view for easier queries
