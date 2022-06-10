import json
import sys
import subprocess
from snakemake.logging import logger

with open(snakemake.log[0], "w") as log_file:
    sys.stdout = log_file
    sys.stderr = sys.stdout

    rmd_params = {
        "sample_name": snakemake.wildcards.sample,
        "coverage_file": snakemake.input.coverage,
        "logo": snakemake.input.logo,
        "coverage_table_outfile": snakemake.output.table_outfile,
        "multiqc_ran": snakemake.params[0]["multiqc_ran"]
    }

    if snakemake.params[0]["multiqc_ran"]:
        rmd_params["multiqc_rel_path"] = snakemake.params[0]["multiqc_rel_path"]

    call = [
        * snakemake.params[0]['rscript_exec'].split(" "), snakemake.input.script,
        snakemake.input.report,
        snakemake.output.html_report,
        snakemake.input.header,
        json.dumps(rmd_params)
    ]

    logger.info(call)
    subprocess.run(
        call,
        stdout = log_file,
        stderr = log_file
    )
