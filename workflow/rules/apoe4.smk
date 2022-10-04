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
assert config["causal_hap"] in config["samples"]
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

rule sim_pts:
    input:
        pgen = rules.transform.output.pgen,
        pvar = rules.transform.output.pvar,
        psam = rules.transform.output.psam,
        hap = lambda wildcards: str(config["hap_files"][wildcards.samp]),
    params:
        beta = lambda wildcards: wildcards.beta,
        h2 = lambda wildcards: wildcards.heritability,
    output:
        pts = out+"sim_pts/h{heritability}/b{beta}/{samp}.pheno",
    resources:
        runtime="0:04:00"
    log:
        out+"logs/sim_pts/h{heritability}/b{beta}/{samp}.log"
    benchmark:
        out+"bench/sim_pts/h{heritability}/b{beta}/{samp}.txt"
    conda:
        "../envs/haptools.yml"
    shell:
        "haptools simphenotype -o {output.pts} -h {params.h2} -v DEBUG {input.pgen} "
        "<( sed 's/EUR\\t0.99$/EUR\\t{params.beta}/' {input.hap} ) &>{log}"

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

rule gwas:
    input:
        pgen = rules.merge.output.pgen.format(samp=config["causal_hap"]),
        pvar = rules.merge.output.pvar.format(samp=config["causal_hap"]),
        psam = rules.merge.output.psam.format(samp=config["causal_hap"]),
        pts = rules.sim_pts.output.pts,
    params:
        in_prefix = lambda w, input: Path(input.pgen).with_suffix(""),
        out_prefix = lambda w, output: Path(output.log).with_suffix(""),
    output:
        log = temp(out+"sim_pts/h{heritability}/b{beta}/{samp}.log"),
        linear = out+"sim_pts/h{heritability}/b{beta}/{samp}.{samp}.glm.linear",
    resources:
        runtime="0:04:00"
    log:
        out+"logs/gwas/h{heritability}/b{beta}/{samp}.log"
    benchmark:
        out+"bench/gwas/h{heritability}/b{beta}/{samp}.txt"
    threads: 1
    conda:
        "../envs/default.yml"
    shell:
        "plink2 --linear allow-no-covars --variance-standardize "
        "--pheno {input.pts} --pfile {params.in_prefix} --out {params.out_prefix} "
        "--threads {threads} &>{log}"

rule manhattan:
    input:
        linear = expand(
            out+"sim_pts/h{heritability}/b{beta}/{samp}.{samp}.glm.linear",
            samp=config["samples"], allow_missing=True,
        ),
    params:
        linear = lambda wildcards, input: [f"-l {i}" for i in input.linear],
        red_ids = [f"-i {i}" for i in config["samples"][0].split("-")],
        orange_ids = "-b "+config["causal_hap"],
    output:
        png = out+"sim_pts/h{heritability}/b{beta}/manhattan.pdf",
    resources:
        runtime="0:05:00"
    log:
        out+"logs/manhattan/h{heritability}/b{beta}/manhattan.log"
    benchmark:
        out+"bench/manhattan/h{heritability}/b{beta}/manhattan.txt"
    conda:
        "../envs/default.yml"
    shell:
        "workflow/scripts/manhattan.py -o {output.png} {params.linear} "
        "{params.red_ids} {params.orange_ids} &>{log}"
