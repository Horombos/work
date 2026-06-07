import argparse
import csv
import re
from pathlib import Path


FILENAME_RE = re.compile(
    r"^global_sasa_fit_(?P<variant>.+?)_results_C(?P<c_value>\d+)_sep(?P<min_seq_sep>\d+)_maxfev(?P<maxfev>\d+)\.csv$"
)

PROJECT_DIR = Path(__file__).resolve().parent
DEFAULT_SOURCE_DIR = Path(r"C:\Users\Lenovo\Desktop\SQLdata")
DEFAULT_OUTPUT_PATH = PROJECT_DIR / "04_seed_data.sql"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a SQL Server seed script from WSME CSV result files."
    )
    parser.add_argument(
        "--source-dir",
        type=Path,
        default=DEFAULT_SOURCE_DIR,
        help=f"Directory containing the CSV files (default: {DEFAULT_SOURCE_DIR})",
    )
    parser.add_argument(
        "--output-path",
        type=Path,
        default=DEFAULT_OUTPUT_PATH,
        help=f"Output T-SQL file path (default: {DEFAULT_OUTPUT_PATH})",
    )
    return parser.parse_args()


def parse_filename(csv_path: Path) -> dict:
    match = FILENAME_RE.match(csv_path.name)
    if not match:
        raise ValueError(f"Unsupported filename format: {csv_path.name}")
    meta = match.groupdict()
    return {
        "source_file": str(csv_path),
        "model_variant": meta["variant"],
        "c_value": int(meta["c_value"]),
        "min_seq_sep": int(meta["min_seq_sep"]),
        "maxfev": int(meta["maxfev"]),
    }


def sql_string(value: str | None) -> str:
    if value is None or value == "":
        return "NULL"
    escaped = value.replace("'", "''")
    return f"N'{escaped}'"


def sql_number(value: str) -> str:
    if value == "":
        return "NULL"
    return value


def sql_int(value: str) -> str:
    if value == "":
        return "NULL"
    return str(int(float(value)))


def load_csv_rows(csv_path: Path) -> list[dict[str, str]]:
    with csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def chunked(items: list[str], size: int) -> list[list[str]]:
    return [items[index:index + size] for index in range(0, len(items), size)]


def build_insert_rows(rows: list[dict[str, str]]) -> list[str]:
    values = []
    for row in rows:
        values.append(
            "("
            + ", ".join(
                [
                    "@run_id",
                    sql_string(row["protein_id"]),
                    sql_number(row["target_dG"]),
                    sql_number(row["predicted_dG"]),
                    sql_number(row["error"]),
                    sql_number(row["sq_error"]),
                    sql_number(row["E_activation"]),
                    sql_int(row["q_unf"]),
                    sql_int(row["q_ts"]),
                    sql_int(row["q_fold"]),
                    sql_string(row["barrier_warning"]),
                    sql_string(row["status"]),
                    sql_string(row["failure_reason"]),
                ]
            )
            + ")"
        )
    return values


def generate_script(source_dir: Path) -> str:
    csv_files = sorted(source_dir.glob("*.csv"))
    if not csv_files:
        raise FileNotFoundError(f"No CSV files found in {source_dir}")

    lines: list[str] = [
        "USE ProteinWSME;",
        "GO",
        "",
        "SET NOCOUNT ON;",
        "GO",
        "",
        "BEGIN TRANSACTION;",
        "DECLARE @run_id INT;",
        "",
    ]

    for csv_path in csv_files:
        meta = parse_filename(csv_path)
        rows = load_csv_rows(csv_path)
        lines.extend(
            [
                f"-- {csv_path.name}",
                "INSERT INTO dbo.runs (",
                "    source_file,",
                "    model_variant,",
                "    c_value,",
                "    min_seq_sep,",
                "    maxfev,",
                "    protein_count",
                ")",
                "VALUES (",
                f"    {sql_string(meta['source_file'])},",
                f"    {sql_string(meta['model_variant'])},",
                f"    {meta['c_value']},",
                f"    {meta['min_seq_sep']},",
                f"    {meta['maxfev']},",
                f"    {len(rows)}",
                ");",
                "",
                "SET @run_id = SCOPE_IDENTITY();",
                "",
            ]
        )

        insert_rows = build_insert_rows(rows)
        for chunk in chunked(insert_rows, 50):
            lines.extend(
                [
                    "INSERT INTO dbo.protein_results (",
                    "    run_id,",
                    "    protein_id,",
                    "    target_dg,",
                    "    predicted_dg,",
                    "    error,",
                    "    sq_error,",
                    "    e_activation,",
                    "    q_unf,",
                    "    q_ts,",
                    "    q_fold,",
                    "    barrier_warning,",
                    "    status,",
                    "    failure_reason",
                    ")",
                    "VALUES",
                    ",\n".join(chunk) + ";",
                    "",
                ]
            )

    lines.extend(
        [
            "COMMIT TRANSACTION;",
            "GO",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    args = parse_args()
    source_dir = args.source_dir.resolve()
    output_path = args.output_path.resolve()

    if not source_dir.exists():
        raise FileNotFoundError(f"Source directory not found: {source_dir}")

    script = generate_script(source_dir)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(script, encoding="utf-8")

    print(f"Seed script created: {output_path}")


if __name__ == "__main__":
    main()
