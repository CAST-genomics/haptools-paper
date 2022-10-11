from pathlib import Path


out = "results/apoe4/"
# get the "sample" names from the hap filenames
config["hap_files"] = [
    x for x in Path(config["hap_files"]).glob("**/*")
    if x.is_file() and x.suffix == ".hap" and x.name[:-4] in [
        config["causal_hap"], config["snps_hap"]
    ]
]
config["samples"] = [x.name[:-4] for x in config["hap_files"]]
config["hap_files"] = dict(zip(config["samples"], config["hap_files"]))
config["genotypes"] = Path(config["genotypes"])
# create a region param that encodes a 1 Mbp region around the given site
config["region"] = tuple(map(int, config["region"].split(":")))
config["region"] = tuple(map(str, (
    config["region"][0],
    int(config["region"][1] - (1e6/2)),
    int(config["region"][1] + (1e6/2)),
)))


rule vcf2pgen:
    input:
        gts = config["genotypes"],
    params:
        prefix = lambda w, output: Path(output.pgen).with_suffix(""),
        chrom = config["region"][0],
        from_bp = config["region"][1],
        to_bp = config["region"][2],
    output:
        log = temp(out+"vcf2pgen/1000G_chr19.log"),
        pgen = out+"vcf2pgen/1000G_chr19.pgen",
        pvar = out+"vcf2pgen/1000G_chr19.pvar",
        psam = out+"vcf2pgen/1000G_chr19.psam",
    resources:
        runtime="0:04:00"
    log:
        out+"logs/vcf2pgen/1000G_chr19.log"
    benchmark:
        out+"bench/vcf2pgen/1000G_chr19.txt"
    threads: 12
    conda:
        "../envs/default.yml"
    shell:
        "plink2 --vcf {input.gts} --make-pgen --snps-only --max-alleles 2 "
        "--threads {threads} --chr {params.chrom} --from-bp {params.from_bp} "
        "--to-bp {params.to_bp} --out {params.prefix} &> {log}"

rule transform:
    input:
        pgen = rules.vcf2pgen.output.pgen,
        pvar = rules.vcf2pgen.output.pvar,
        psam = rules.vcf2pgen.output.psam,
        hap = lambda wildcards: str(config["hap_files"][wildcards.samp]),
    output:
        pgen = temp(out+"transform/{samp}.pgen"),
        pvar = temp(out+"transform/{samp}.pvar"),
        psam = temp(out+"transform/{samp}.psam"),
    resources:
        runtime="0:04:00"
    log:
        out+"logs/transform/{samp}.log"
    benchmark:
        out+"bench/transform/{samp}.txt"
    conda:
        "../envs/haptools.yml"
    shell:
        "haptools transform -v INFO -o {output.pgen} {input.pgen} {input.hap} &> {log}"

rule sim_pt:
    input:
        pgen = rules.transform.output.pgen,
        pvar = rules.transform.output.pvar,
        psam = rules.transform.output.psam,
        hap = lambda wildcards: str(config["hap_files"][wildcards.samp]),
    params:
        beta = lambda wildcards: wildcards.beta,
        h2 = lambda wildcards: wildcards.heritability,
    output:
        pts = out+"sim_pt/b{beta}/h{heritability}/{samp}.pheno",
    resources:
        runtime="0:04:00"
    log:
        out+"logs/sim_pt/b{beta}/h{heritability}/{samp}.log"
    benchmark:
        out+"bench/sim_pt/b{beta}/h{heritability}/{samp}.txt"
    conda:
        "../envs/haptools.yml"
    shell:
        "haptools simphenotype -o {output.pts} -h {params.h2} -v DEBUG {input.pgen} "
        "<( sed 's/\\t0.99$/\\t{params.beta}/' {input.hap} ) &>{log}"

rule pgen2bed:
    input:
        pgen = rules.transform.output.pgen,
        pvar = rules.transform.output.pvar,
        psam = rules.transform.output.psam,
    params:
        in_prefix = lambda w, input: Path(input.pgen).with_suffix(""),
        out_prefix = lambda w, output: Path(output.bed).with_suffix(""),
    output:
        log = temp(out+"pgen2bed/{samp}.log"),
        bed = out+"pgen2bed/{samp}.bed",
        bim = out+"pgen2bed/{samp}.bim",
        fam = out+"pgen2bed/{samp}.fam",
    resources:
        runtime="0:04:00"
    log:
        out+"logs/pgen2bed/{samp}.log"
    benchmark:
        out+"bench/pgen2bed/{samp}.txt"
    threads: 12
    conda:
        "../envs/default.yml"
    shell:
        "plink2 --pfile {params.in_prefix} --make-bed --out {params.out_prefix} &> {log}"

rule sim_gcta:
    input:
        bed = rules.pgen2bed.output.bed,
        bim = rules.pgen2bed.output.bim,
        fam = rules.pgen2bed.output.fam,
        hap = lambda wildcards: str(config["hap_files"][wildcards.samp]),
    params:
        beta = lambda wildcards: wildcards.beta,
        h2 = lambda wildcards: wildcards.heritability,
        in_prefix = lambda w, input: Path(input.bed).with_suffix(""),
        out_prefix = lambda w, output: Path(output.pts).with_suffix(""),
    output:
        par = temp(out+"sim_gcta/b{beta}/h{heritability}/{samp}.par"),
        pts = out+"sim_gcta/b{beta}/h{heritability}/{samp}.phen",
    resources:
        runtime="0:04:00"
    log:
        out+"logs/sim_gcta/b{beta}/h{heritability}/{samp}.log"
    benchmark:
        out+"bench/sim_gcta/b{beta}/h{heritability}/{samp}.txt"
    conda:
        "../envs/gcta.yml"
    shell:
        "gcta64 --bfile {params.in_prefix} --simu-qt --simu-causal-loci "
        "<(grep -E '^H' {input.hap} | cut -f5 | sed 's/$/\\t{params.beta}/') "
        "--simu-hsq {params.h2} --simu-rep 1 --out {params.out_prefix} &> {log} && "
        "sed -i 's/^0 //g' {output.pts} &>> {log}"

rule merge:
    input:
        gts = rules.vcf2pgen.output.pgen,
        gts_pvar = rules.vcf2pgen.output.pvar,
        gts_psam = rules.vcf2pgen.output.psam,
        hps = rules.transform.output.pgen,
        hps_pvar = rules.transform.output.pvar,
        hps_psam = rules.transform.output.psam,
    output:
        pgen = out+"merge/{samp}.pgen",
        pvar = out+"merge/{samp}.pvar",
        psam = out+"merge/{samp}.psam",
    resources:
        runtime="0:04:00"
    log:
        out+"logs/merge/{samp}.log"
    benchmark:
        out+"bench/merge/{samp}.txt"
    conda:
        "../envs/haptools.yml"
    shell:
        "workflow/scripts/merge.py {input.gts} {input.hps} {output.pgen} &> {log}"

def gwas_pts_input(wildcards):
    source = rules.sim_pt
    if wildcards.meth == "gcta":
        source = rules.sim_gcta
    return source.output.pts

rule gwas:
    input:
        pgen = rules.merge.output.pgen.format(samp=config["causal_hap"]),
        pvar = rules.merge.output.pvar.format(samp=config["causal_hap"]),
        psam = rules.merge.output.psam.format(samp=config["causal_hap"]),
        pts = gwas_pts_input,
    params:
        in_prefix = lambda w, input: Path(input.pgen).with_suffix(""),
        out_prefix = lambda w, output: Path(output.log).with_suffix(""),
        mv_cmd = lambda wildcards, output: (
            "mv "+str(Path(output.log).with_suffix(".PHENO1.glm.linear"))+f" {output.linear}"
            if wildcards.meth == "gcta" else "true"
        ),
    output:
        log = temp(out+"sim_{meth}/b{beta}/h{heritability}/{samp}.log"),
        linear = out+"sim_{meth}/b{beta}/h{heritability}/{samp}.{samp}.glm.linear",
    resources:
        runtime="0:04:00"
    log:
        out+"logs/gwas/{meth}/b{beta}/h{heritability}/{samp}.log"
    benchmark:
        out+"bench/gwas/{meth}/b{beta}/h{heritability}/{samp}.txt"
    threads: 1
    conda:
        "../envs/default.yml"
    shell:
        "plink2 --linear allow-no-covars --variance-standardize "
        "--pheno iid-only {input.pts} --pfile {params.in_prefix} --out {params.out_prefix} "
        "--threads {threads} &>{log} && {params.mv_cmd} &>>{log}"

rule manhattan:
    input:
        linear = expand(
            rules.gwas.output.linear,
            samp=config["samples"],
            meth="pt",
            allow_missing=True,
        ),
    params:
        linear = lambda wildcards, input: [f"-l {i}" for i in input.linear],
        red_ids = [f"-i {i}" for i in config["snps_hap"].split("-")],
        orange_ids = "-b "+config["causal_hap"],
    output:
        png = out+"sim_pt/b{beta}/h{heritability}/manhattan.pdf",
    resources:
        runtime="0:05:00"
    log:
        out+"logs/manhattan/b{beta}/h{heritability}/manhattan.log"
    benchmark:
        out+"bench/manhattan/b{beta}/h{heritability}/manhattan.txt"
    conda:
        "../envs/default.yml"
    shell:
        "workflow/scripts/manhattan.py -o {output.png} {params.linear} "
        "{params.red_ids} {params.orange_ids} &>{log}"

# rule compare_gcta:
#     input:
#         pheno = expand(
#             rules.sim_pt.output.pts,
#             samp=config["snps_hap"],
#             heritability=config["heritabilities"],
#             allow_missing=True
#         ),
#         phen = expand(
#             rules.sim_gcta.output.pts,
#             samp=config["snps_hap"],
#             heritability=config["heritabilities"],
#             allow_missing=True
#         ),
#     params:
#         pheno = lambda wildcards: expand(
#             rules.sim_pt.output.pts,
#             samp=config["snps_hap"],
#             beta=wildcards.beta,
#             allow_missing=True,
#         )[0],
#         phen = lambda wildcards: expand(
#             rules.sim_gcta.output.pts,
#             samp=config["snps_hap"],
#             beta=wildcards.beta,
#             allow_missing=True,
#         )[0],
#     output:
#         png = out+"sim_gcta/b{beta}/compare_gcta.pdf",
#     resources:
#         runtime="0:04:00"
#     log:
#         out+"logs/compare_gcta/b{beta}/compare_gcta.log"
#     benchmark:
#         out+"bench/compare_gcta/b{beta}/compare_gcta.txt"
#     conda:
#         "../envs/default.yml"
#     shell:
#         "workflow/scripts/compare_gcta.py -o {output.png} {params} &> {log}"

rule compare_gcta:
    input:
        pheno = expand(
            rules.sim_pt.output.pts,
            samp=config["snps_hap"],
            allow_missing=True
        ),
        phen = expand(
            rules.sim_gcta.output.pts,
            samp=config["snps_hap"],
            allow_missing=True
        ),
    output:
        png = out+"sim_gcta/b{beta}/h{heritability}/compare_gcta.pdf",
    resources:
        runtime="0:04:00"
    log:
        out+"logs/compare_gcta/b{beta}/h{heritability}/compare_gcta.log"
    benchmark:
        out+"bench/compare_gcta/b{beta}/h{heritability}/compare_gcta.txt"
    conda:
        "../envs/default.yml"
    shell:
        "workflow/scripts/compare_gcta.py -o {output.png} {input} &> {log}"

rule compare_gcta_manhattan:
    input:
        pheno = expand(
            rules.gwas.output.linear,
            samp=config["snps_hap"],
            meth="pt",
            allow_missing=True
        ),
        phen = expand(
            rules.gwas.output.linear,
            samp=config["snps_hap"],
            meth="gcta",
            allow_missing=True
        ),
    params:
        red_ids = [f"-i {i}" for i in config["snps_hap"].split("-")],
        orange_ids = "-b "+config["causal_hap"],
    output:
        png = out+"sim_gcta/b{beta}/h{heritability}/manhattan.pdf",
    resources:
        runtime="0:04:00"
    log:
        out+"logs/compare_gcta_manhattan/b{beta}/h{heritability}/manhattan.log"
    benchmark:
        out+"bench/compare_gcta_manhattan/b{beta}/h{heritability}/manhattan.txt"
    conda:
        "../envs/default.yml"
    shell:
        "workflow/scripts/manhattan.py -o {output.png} -l {input.pheno} "
        "-l {input.phen} {params.red_ids} {params.orange_ids} &>{log}"
