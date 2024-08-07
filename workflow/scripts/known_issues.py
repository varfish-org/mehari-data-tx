from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.script import snakemake

from contextlib import redirect_stderr


def main(known_issues: list[str]):
    with open(snakemake.output.known_issues, "w") as file:
        print("id_type\tid\tdescription", file=file)
        for issue in known_issues:
            print(
                f"{issue['id_type']}\t{issue['id']}\t{issue['description']}",
                file=file,
            )


with open(snakemake.log[0], "w") as log, redirect_stderr(log):
    main(snakemake.params["known_issues"])
